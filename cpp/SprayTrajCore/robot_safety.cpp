// robot_safety.cpp —— 机器人可执行性校验（VS-087）
//
// 参数与算法同源于 tools/DensoArmKine（C#）：
//   连杆表  <- Models/DensoLibrary.cs::Vs087()   （官方 vs087.urdf，reach 905）
//   FK      <- Kinematics/FkSolver.cs            T_i = T_{i-1}·Trans(offset_i)·Rot(axis_i, θ_i)
//   碰撞体  <- Cloud/ArmCollider.cs              ArmSizing.From + Build（半径唯一来源）
//   IK      <- 阻尼最小二乘（DensoArmKine 用多初值数值解，语义一致）
//
// 三条硬要求：
//   ① 非朝上法向加倾角；② 倾角改变时 TCP 到面点距离严格不变（球面旋转）；
//   ③ 枪体/连杆 vs 原始点云碰撞 + 非相邻连杆自碰撞。
//
// 优先级：**硬约束（碰撞/自碰撞/关节超限/低于绝对安全距离） > 姿态连续 > 间隙告警**。
// 实测 complex 的 idx 3163：本点最优方位下间隙 -2.35mm，被"整行统一方位"顶成 -5.81mm
// 真穿透；这种硬约束被破坏时必须回退。
// 但只要平滑后仍高于 hardClearance（默认 2mm），就不再因为"比本点最优差一点"而回退：
// 有凸点、间隙从 30mm 降到 8mm 没关系，以平滑为主；safety（默认 20mm）只负责标橙告警。
// 见 stage3 的 reverted。

#include "robot_safety.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
constexpr double kRad = 180.0 / kPi;

// ---------------------------------------------------------------- 计时/计数（RS_PROFILE=1 才输出）
// 生产不看这个；调优时靠它判断"慢在 IK 还是慢在点云查询"，否则只能瞎猜。
struct Prof {
    bool on = false;
    // 并行求解时 +1 会竞争，用 relaxed 原子（计数只在 RS_PROFILE 时看，精度要求不高）
    std::atomic<long long> ik{0}, seg{0}, evals{0};
    // eval 内部的分项耗时（纳秒，多线程累加）。只回答一个问题：
    // 「139μs/次的 eval 到底花在 IK、建碰撞体，还是点云查询」——这三条优化方向完全不同。
    std::atomic<long long> nsIk{0}, nsBuild{0}, nsSeg{0}, nsSelf{0};
    // IK 内部统计：ik() 一次调用里到底跑了几个 ikOnce、以及退化为"全种子搜"的比例。
    // 这三个数决定优化方向：若 seed 解本身已经几何收敛，只是关节超限，
    // 那么把 29 组种子全跑一遍去"换个不超限的支"就是纯浪费。
    std::atomic<long long> ikOnceN{0}, ikFast{0}, ikOver{0}, ikFail{0};
    // 「超限 → 全种子搜」这条最贵路径的投入产出：总 ikOnce、换来合法解的次数、白跑的次数
    std::atomic<long long> ikOverOnce{0}, ikOverHit{0}, ikOverMiss{0};
    // ikOnce 的迭代总次数：除以 ikOnceN 就是"平均迭代几次"，用来验证停滞退出生效。
    // **不要**在迭代循环里加 Prof::Scope 级的计时：实测即使 g_prof.on=false，
    // 三处计时也会让整个热循环慢 17%（编译器被迫跨原子写做保守假设）。
    // 内层各段占比是一次性测过的：fk 43% · 6x6 解 25% · 雅可比/正规阵 12% · 其余 20%。
    std::atomic<long long> ikIters{0};
    // 「超限多深 → 全种子搜能不能救回」直方图。阈值不靠猜：
    // 浅超限（几度）换个解支有可能滑回行程内，深超限（几十度）是"这块够不到"，
    // 再搜 28 组种子也是白跑。分桶后按实测救回率决定跳过界限。
    static constexpr int kOverBucketN = 7;
    static constexpr double kOverEdges[kOverBucketN] = {2.0, 5.0, 10.0, 20.0, 30.0, 45.0, 1e9};
    const char* const kOverNames[kOverBucketN] = {"0-2°", "2-5°", "5-10°", "10-20°",
                                                  "20-30°", "30-45°", "45°+"};
    std::atomic<long long> overCnt[kOverBucketN] = {};
    std::atomic<long long> overHitB[kOverBucketN] = {};
    // 全种子搜里有多少组种子是**重复**的（与前面某一组逐位相同）。
    // 重复种子必然给出逐位相同的结果，跑第二遍是纯粹的浪费 —— 这是唯一
    // 一种"剪掉不影响任何输出"的剪枝。实测重复率 0%，所以没有可剪的。
    std::atomic<long long> sweepSeeds{0}, sweepDup{0};
    long long lastEvals = 0;    // mark() 用：本阶段新增的 eval 次数
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point last = t0;

    void reset() {
        ik = 0; seg = 0; evals = 0;
        nsIk = 0; nsBuild = 0; nsSeg = 0; nsSelf = 0;
        ikOnceN = 0; ikFast = 0; ikOver = 0; ikFail = 0;
        ikOverOnce = 0; ikOverHit = 0; ikOverMiss = 0;
        ikIters = 0;
        sweepSeeds = 0; sweepDup = 0;
        for (int i = 0; i < kOverBucketN; ++i) { overCnt[i] = 0; overHitB[i] = 0; }
        lastEvals = 0;
    }

    // margin 为负表示超限，负数越大超得越多。桶按"超限的度数"划分。
    void bucketOver(double margin, bool hit) {
        double dep = -margin;
        int k = 0;
        while (k < kOverBucketN - 1 && dep > kOverEdges[k]) ++k;
        ++overCnt[k];
        if (hit) ++overHitB[k];
    }

    // eval 内部计时用：只在 RS_PROFILE 时读时钟（steady_clock::now 约 20ns，可接受）
    struct Scope {
        bool on;
        std::atomic<long long>& acc;
        std::chrono::steady_clock::time_point b;
        Scope(bool on_, std::atomic<long long>& a) : on(on_), acc(a) {
            if (on) b = std::chrono::steady_clock::now();
        }
        ~Scope() {
            if (!on) return;
            acc.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(
                              std::chrono::steady_clock::now() - b).count(),
                          std::memory_order_relaxed);
        }
    };
    void mark(const char* what) {
        if (!on) return;
        auto now = std::chrono::steady_clock::now();
        double dt = std::chrono::duration<double, std::milli>(now - last).count();
        double all = std::chrono::duration<double, std::milli>(now - t0).count();
        long long ev = evals.load();
        std::fprintf(stderr, "  [prof] %-14s %9.1f ms（累计 %9.1f）· eval +%lld\n",
                     what, dt, all, ev - lastEvals);
        lastEvals = ev;
        last = now;
    }
};
Prof g_prof;
constexpr double kEps = 1e-12;

inline double clampd(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }
inline double sq(double v) { return v * v; }

// ================================================================ 基础几何

struct V3 {
    double x = 0.0, y = 0.0, z = 0.0;
    V3() = default;
    V3(double a, double b, double c) : x(a), y(b), z(c) {}
    double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
    double& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
};

inline V3 operator+(const V3& a, const V3& b) { return V3(a.x + b.x, a.y + b.y, a.z + b.z); }
inline V3 operator-(const V3& a, const V3& b) { return V3(a.x - b.x, a.y - b.y, a.z - b.z); }
inline V3 operator-(const V3& a) { return V3(-a.x, -a.y, -a.z); }
inline V3 operator*(const V3& a, double s) { return V3(a.x * s, a.y * s, a.z * s); }
inline V3 operator*(double s, const V3& a) { return a * s; }
inline double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V3 cross(const V3& a, const V3& b) {
    return V3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
inline double norm(const V3& a) { return std::sqrt(dot(a, a)); }
// 零向量原样返回：调用方自己判退化，这里不静默改成 +Z（geom.unit 的同一约定）
inline V3 unit(const V3& a) {
    double n = norm(a);
    return n > kEps ? a * (1.0 / n) : a;
}

struct M3 {
    double m[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    static M3 I() {
        M3 r;
        r.m[0][0] = r.m[1][1] = r.m[2][2] = 1.0;
        return r;
    }
    V3 col(int j) const { return V3(m[0][j], m[1][j], m[2][j]); }
    void setCol(int j, const V3& v) {
        m[0][j] = v.x;
        m[1][j] = v.y;
        m[2][j] = v.z;
    }
};

inline M3 m3mul(const M3& a, const M3& b) {
    M3 r;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r.m[i][j] = a.m[i][0] * b.m[0][j] + a.m[i][1] * b.m[1][j] + a.m[i][2] * b.m[2][j];
    return r;
}
inline V3 m3mul(const M3& a, const V3& v) {
    return V3(a.m[0][0] * v.x + a.m[0][1] * v.y + a.m[0][2] * v.z,
              a.m[1][0] * v.x + a.m[1][1] * v.y + a.m[1][2] * v.z,
              a.m[2][0] * v.x + a.m[2][1] * v.y + a.m[2][2] * v.z);
}
inline M3 m3T(const M3& a) {
    M3 r;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r.m[i][j] = a.m[j][i];
    return r;
}

// 罗德里格斯（axis 内部归一化；零轴返回单位阵，与 geom.rot_axis 一致）
inline M3 rotAxis(const V3& axis, double ang) {
    V3 a = unit(axis);
    if (norm(a) < kEps) return M3::I();
    double c = std::cos(ang), s = std::sin(ang), t = 1.0 - c;
    M3 r;
    r.m[0][0] = c + a.x * a.x * t;
    r.m[0][1] = a.x * a.y * t - a.z * s;
    r.m[0][2] = a.x * a.z * t + a.y * s;
    r.m[1][0] = a.y * a.x * t + a.z * s;
    r.m[1][1] = c + a.y * a.y * t;
    r.m[1][2] = a.y * a.z * t - a.x * s;
    r.m[2][0] = a.z * a.x * t - a.y * s;
    r.m[2][1] = a.z * a.y * t + a.x * s;
    r.m[2][2] = c + a.z * a.z * t;
    return r;
}

struct Pose {
    M3 R;
    V3 t;
};
inline Pose poseMul(const Pose& a, const Pose& b) { return Pose{m3mul(a.R, b.R), m3mul(a.R, b.t) + a.t}; }
inline Pose poseTrans(double x, double y, double z) {
    Pose p;
    p.R = M3::I();
    p.t = V3(x, y, z);
    return p;
}
inline Pose poseRot(const V3& axis, double ang) {
    Pose p;
    p.R = rotAxis(axis, ang);
    return p;
}

// 旋转矩阵 → 轴角（θ≈π 用 R+I 的对角列取轴，避免除以 sinθ≈0）
inline V3 matToRotvec(const M3& R) {
    double tr = R.m[0][0] + R.m[1][1] + R.m[2][2];
    double c = clampd((tr - 1.0) * 0.5, -1.0, 1.0);
    double th = std::acos(c);
    if (th < 1e-9) return V3(0, 0, 0);
    if (std::fabs(kPi - th) > 1e-6) {
        double k = th / (2.0 * std::sin(th));
        return V3(R.m[2][1] - R.m[1][2], R.m[0][2] - R.m[2][0], R.m[1][0] - R.m[0][1]) * k;
    }
    double a00 = (R.m[0][0] + 1.0) * 0.5, a11 = (R.m[1][1] + 1.0) * 0.5, a22 = (R.m[2][2] + 1.0) * 0.5;
    V3 axis = unit(V3(std::sqrt(std::max(a00, 0.0)), std::sqrt(std::max(a11, 0.0)),
                      std::sqrt(std::max(a22, 0.0))));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (i != j && R.m[i][j] * axis[i] * axis[j] < 0) axis[j] = -axis[j];
    return axis * th;
}

inline double angleBetweenDeg(const V3& a, const V3& b) {
    double c = clampd(dot(unit(a), unit(b)), -1.0, 1.0);
    return std::acos(c) * kRad;
}

// 由枪轴（工具 +Z）构造姿态。ref 退化时换轴 —— 枪轴接近竖直是俯喷/仰喷的常态，
// 固定用 +Z 当参考会在这里除零。
inline M3 frameFromZ(const V3& zAxis) {
    V3 z = unit(zAxis);
    V3 ref(0.0, 0.0, 1.0);
    if (std::fabs(dot(z, ref)) > 0.95) ref = V3(1.0, 0.0, 0.0);
    V3 x = cross(ref, z);
    if (norm(x) < 1e-9) x = cross(V3(1.0, 0.0, 0.0), z);
    x = unit(x);
    V3 y = unit(cross(z, x));
    M3 R;
    R.setCol(0, x);
    R.setCol(1, y);
    R.setCol(2, z);
    return R;
}

// 点到线段距离
inline double pointSegDist(const V3& p, const V3& a, const V3& b) {
    V3 ab = b - a;
    double l2 = dot(ab, ab);
    if (l2 < 1e-12) return norm(p - a);
    double t = clampd(dot(p - a, ab) / l2, 0.0, 1.0);
    return norm(p - (a + ab * t));
}

// 两条线段最短距离（Ericson, Real-Time Collision Detection 5.1.9 解析解）
inline double segSegDist(const V3& p1, const V3& q1, const V3& p2, const V3& q2) {
    V3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    double a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    double s = 0.0, t = 0.0;
    if (a <= kEps && e <= kEps) return norm(r);
    if (a <= kEps) {
        s = 0.0;
        t = clampd(f / e, 0.0, 1.0);
    } else {
        double c = dot(d1, r);
        if (e <= kEps) {
            t = 0.0;
            s = clampd(-c / a, 0.0, 1.0);
        } else {
            double b = dot(d1, d2);
            double denom = a * e - b * b;
            s = denom > kEps ? clampd((b * f - c * e) / denom, 0.0, 1.0) : 0.0;
            t = (b * s + f) / e;
            if (t < 0.0) {
                t = 0.0;
                s = clampd(-c / a, 0.0, 1.0);
            } else if (t > 1.0) {
                t = 1.0;
                s = clampd((b - c) / a, 0.0, 1.0);
            }
        }
    }
    return norm((p1 + d1 * s) - (p2 + d2 * t));
}

// ================================================================ 6x6 线性代数

// 高斯-约当（部分主元）。奇异返回 false。
bool solve6(double A[6][6], const double b[6], double x[6]) {
    double a[6][7];
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) a[i][j] = A[i][j];
        a[i][6] = b[i];
    }
    for (int col = 0; col < 6; ++col) {
        int piv = col;
        double best = std::fabs(a[col][col]);
        for (int r = col + 1; r < 6; ++r) {
            if (std::fabs(a[r][col]) > best) {
                best = std::fabs(a[r][col]);
                piv = r;
            }
        }
        if (best < 1e-14) return false;
        if (piv != col)
            for (int j = col; j < 7; ++j) std::swap(a[col][j], a[piv][j]);
        double d = a[col][col];
        for (int j = col; j < 7; ++j) a[col][j] /= d;
        for (int r = 0; r < 6; ++r) {
            if (r == col) continue;
            double f = a[r][col];
            if (f == 0.0) continue;
            for (int j = col; j < 7; ++j) a[r][j] -= f * a[col][j];
        }
    }
    for (int i = 0; i < 6; ++i) x[i] = a[i][6];
    return true;
}

// 对称阵特征值（Jacobi 旋转，6x6，量级小，迭代收敛快）
void jacobiEigen(double A[6][6], double ev[6]) {
    for (int sweep = 0; sweep < 60; ++sweep) {
        double off = 0.0;
        for (int p = 0; p < 5; ++p)
            for (int q = p + 1; q < 6; ++q) off += A[p][q] * A[p][q];
        if (off < 1e-24) break;
        for (int p = 0; p < 5; ++p) {
            for (int q = p + 1; q < 6; ++q) {
                if (std::fabs(A[p][q]) < 1e-18) continue;
                double theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
                double t = (theta >= 0.0 ? 1.0 : -1.0) /
                           (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                double c = 1.0 / std::sqrt(t * t + 1.0);
                double s = t * c;
                double app = A[p][p], aqq = A[q][q], apq = A[p][q];
                A[p][p] = app - t * apq;
                A[q][q] = aqq + t * apq;
                A[p][q] = A[q][p] = 0.0;
                for (int k = 0; k < 6; ++k) {
                    if (k == p || k == q) continue;
                    double akp = A[k][p], akq = A[k][q];
                    A[k][p] = A[p][k] = c * akp - s * akq;
                    A[k][q] = A[q][k] = s * akp + c * akq;
                }
            }
        }
    }
    for (int i = 0; i < 6; ++i) ev[i] = A[i][i];
}

// ================================================================ VS-087 连杆

struct JointDef {
    const char* name;
    double ox, oy, oz;
    double ax, ay, az;
    double lo, hi;
};

// 与 DensoLibrary.Vs087 / vs087.urdf 逐字段一致
const JointDef kJoints[6] = {
    {"J1",   0.0, 0.0, 197.5,  0.0, 0.0, 1.0, -170.0, 170.0},
    {"J2",  30.0, 0.0, 197.5,  0.0, 1.0, 0.0, -100.0, 135.0},
    {"J3",   0.0, 0.0, 445.0,  0.0, 1.0, 0.0, -133.0, 150.0},
    {"J4", -20.0, 0.0, 287.0,  0.0, 0.0, 1.0, -270.0, 270.0},
    {"J5",   0.0, 0.0, 143.0,  0.0, 1.0, 0.0, -120.0, 120.0},
    {"J6",   0.0, 0.0,  80.0,  0.0, 0.0, 1.0, -360.0, 360.0},
};

constexpr double kReach = 905.0;
constexpr double kDefaultToolLen = 150.0;

struct FkOut {
    Pose frames[7];     // frames[0]=基座, frames[i]=J_i 子坐标系
    V3 origins[7];
    V3 axes[6];
    Pose flange;
    Pose tcp;
};

void fk(const double thetaDeg[6], const double toolLen, FkOut& out) {
    out.frames[0].R = M3::I();
    out.frames[0].t = V3(0, 0, 0);
    for (int i = 0; i < 6; ++i) {
        const JointDef& j = kJoints[i];
        V3 ax(j.ax, j.ay, j.az);
        out.axes[i] = unit(m3mul(out.frames[i].R, ax));
        Pose local = poseMul(poseTrans(j.ox, j.oy, j.oz), poseRot(ax, thetaDeg[i] * kDeg));
        out.frames[i + 1] = poseMul(out.frames[i], local);
    }
    out.flange = out.frames[6];
    Pose tool = poseTrans(0.0, 0.0, toolLen);   // 法兰 +Z 指向喷嘴出口
    out.tcp = poseMul(out.flange, tool);
    for (int i = 0; i < 7; ++i) out.origins[i] = out.frames[i].t;
}

// 6x6 解析雅可比（前三行线速度，后三行角速度）
void jacobian(const FkOut& f, double J[6][6]) {
    V3 p = f.tcp.t;
    for (int i = 0; i < 6; ++i) {
        V3 v = cross(f.axes[i], p - f.origins[i + 1]);
        J[0][i] = v.x;
        J[1][i] = v.y;
        J[2][i] = v.z;
        J[3][i] = f.axes[i].x;
        J[4][i] = f.axes[i].y;
        J[5][i] = f.axes[i].z;
    }
}

inline bool overLimit(const double th[6]) {
    for (int i = 0; i < 6; ++i)
        if (th[i] < kJoints[i].lo - 1e-9 || th[i] > kJoints[i].hi + 1e-9) return true;
    return false;
}
inline double limitMarginMin(const double th[6]) {
    double m = 1e18;
    for (int i = 0; i < 6; ++i) m = std::min(m, std::min(th[i] - kJoints[i].lo, kJoints[i].hi - th[i]));
    return m;
}

// 归一化雅可比条件数。**位置行必须除以特征长度**：原始 J 前三行是 mm（量级 ~500）、
// 后三行是 rad（~1），直接 SVD 得到的条件数正常姿态就有上千，与"是否接近奇异"无关。
inline double condNumber(const FkOut& f) {
    double J[6][6];
    jacobian(f, J);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 6; ++j) J[i][j] /= 500.0;
    double A[6][6];
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            double s = 0.0;
            for (int k = 0; k < 6; ++k) s += J[k][i] * J[k][j];
            A[i][j] = s;
        }
    double ev[6];
    jacobiEigen(A, ev);
    double mx = 0.0, mn = 1e18;
    for (int i = 0; i < 6; ++i) {
        mx = std::max(mx, ev[i]);
        mn = std::min(mn, ev[i]);
    }
    double svmax = std::sqrt(std::max(mx, 0.0));
    double svmin = std::sqrt(std::max(mn, 0.0));
    if (svmin <= 1e-9) return std::numeric_limits<double>::infinity();
    return svmax / svmin;
}

// ---------------------------------------------------------------- IK

struct IkSol {
    double theta[6];
    bool ok = false;
    double posErr = 0.0;
    double rotErrDeg = 0.0;
    bool over = false;
    double margin = 0.0;
    // 条件数**不在这里算**：它要一次 6x6 雅可比 + Jacobi 特征值分解（最多 60 轮），
    // 而 ik() 每次调用平均会产生 19 个候选解、其中 18 个直接丢弃。
    // 只在最终入选的那个解上算一次（eval 里已有 fk(sol.theta) 的 FK 结果，顺手就出）。
    double cond = 0.0;
};

// 固定种子：覆盖 肩前/肩后 × 肘上/肘下 × 腕翻/腕正
const double kSeeds[12][6] = {
    {0, 0, -90, 0, 0, 0},       {0, 20, -60, 0, 45, 0},
    {0, -20, -60, 0, -45, 0},   {180, 0, -90, 0, 0, 0},
    {180, 20, -60, 0, 45, 0},   {180, -20, -60, 0, -45, 0},
    {90, 10, -70, 90, 30, 90},  {-90, 10, -70, -90, 30, -90},
    {0, 60, -120, 0, 60, 0},    {180, 60, -120, 0, 60, 0},
    {0, -40, -30, 0, -30, 0},   {180, -40, -30, 0, -30, 0},
};

// 解析种子：前三个关节用两连杆定理，腕部用 ZYZ 分解。
// 为什么必须有它（实测，不是理论顾虑）：纯固定种子在此类轨迹上会收敛到 J2=154°、J3=-102°
// 的怪解 —— 残差为 0 但 J2 超限 19°，整条小臂从工件底下钻过去（穿透 24.7mm）。
// DLS 只保证"收敛到附近的解"，**种子决定解支**。
int analyticSeeds(const V3& pTcp, const M3& RTcp, double toolLen, double out[][6]) {
    V3 d = RTcp.col(2);
    V3 W = pTcp - d * (toolLen + 80.0);   // 腕心近似（VS-087 腕非理想球腕，只当种子）
    double r = std::hypot(W.x, W.y);
    double wz = W.z;
    const double lat = 30.0, a1 = 445.0, a2 = 430.0, baseH = 197.5;
    int cnt = 0;
    for (int si = 0; si < 2; ++si) {
        double s1 = (si == 0) ? 1.0 : -1.0;
        double th1 = std::atan2(W.y, W.x) * kRad + (s1 > 0 ? 0.0 : 180.0);
        double dx = r * s1 - lat;
        double dz = wz - baseH;
        double D2 = dx * dx + dz * dz;
        double c3 = (D2 - a1 * a1 - a2 * a2) / (2.0 * a1 * a2);
        if (std::fabs(c3) > 1.0) continue;
        for (int sj = 0; sj < 2; ++sj) {
            double sgn = (sj == 0) ? 1.0 : -1.0;
            double th3 = sgn * std::acos(clampd(c3, -1.0, 1.0)) * kRad;
            double th2 = std::atan2(dx, dz) * kRad -
                         std::atan2(a2 * std::sin(th3 * kDeg), a1 + a2 * std::cos(th3 * kDeg)) * kRad;
            M3 R03 = m3mul(rotAxis(V3(0, 0, 1), th1 * kDeg), rotAxis(V3(0, 1, 0), (th2 + th3) * kDeg));
            M3 R36 = m3mul(m3T(R03), RTcp);
            for (int sb = 0; sb < 2; ++sb) {
                double bsign = (sb == 0) ? 1.0 : -1.0;
                double sbv = std::hypot(R36.m[0][2], R36.m[1][2]) * bsign;
                double bb = std::atan2(sbv, R36.m[2][2]);
                double aa, cc;
                if (std::fabs(sbv) < 1e-9) {
                    aa = 0.0;
                    cc = std::atan2(-R36.m[1][0], R36.m[0][0]);
                } else {
                    aa = std::atan2(R36.m[1][2], R36.m[0][2]);
                    cc = std::atan2(R36.m[2][1], -R36.m[2][0]);
                }
                if (cnt >= 16) return cnt;
                double* s = out[cnt++];
                s[0] = th1;
                s[1] = th2;
                s[2] = th3;
                s[3] = aa * kRad;
                s[4] = bb * kRad;
                s[5] = cc * kRad;
            }
        }
    }
    return cnt;
}

// 阻尼最小二乘的单次迭代
IkSol ikOnce(const V3& pTarget, const M3& RTarget, double toolLen, const double seed[6]) {
    double th[6];
    for (int i = 0; i < 6; ++i) th[i] = seed[i];
    const double wRot = 100.0;   // 1rad 姿态误差 ≈ 100mm 位置误差，把两行拉到同一量纲
    double lam = 30.0;
    double bestErr = 1e300;
    int stall = 0;               // 连续"没把误差压下去"的迭代次数
    int iters = 0;
    // FkOut 约 1.2KB（7 个 Pose + origins + axes），声明在循环里等于每轮先把
    // 这些 double 全部置零一遍，而 fk() 每轮会把每个字段都写满 —— 这层清零是纯浪费。
    // 提到循环外不改变任何数值：fk 只读自己上一轮写过的 frames[i]，不读初值。
    FkOut f;

    for (int it = 0; it < 90; ++it) {
        ++iters;
        fk(th, toolLen, f);
        V3 ep = pTarget - f.tcp.t;
        V3 er = matToRotvec(m3mul(RTarget, m3T(f.tcp.R)));
        double err[6] = {ep.x, ep.y, ep.z, wRot * er.x, wRot * er.y, wRot * er.z};
        double en = 0.0;
        for (int i = 0; i < 6; ++i) en += err[i] * err[i];
        en = std::sqrt(en);
        const double prevBest = bestErr;
        if (en < bestErr) {
            bestErr = en;
            lam = std::max(lam * 0.7, 0.5);
        } else {
            lam = std::min(lam * 0.5 + lam, 5e3);
        }
        // 停滞按**相对**进步判定：一次只压下去 0.5% 也算停住。
        // 只看"绝对没进步"抓不到真实情况 —— 实测界面口径下 ikOnce 平均要跑 64.8 次，
        // 就是因为每十来次挤出一点点进步，把绝对版 stall 反复清零。
        if (en < prevBest * 0.995) stall = 0;
        else ++stall;
        if (norm(ep) < 0.02 && norm(er) * kRad < 0.02) break;

        // ---- 停滞提前退出（不改变任何返回值的语义）
        // 这个 DLS 的阻尼是"不进步就 ×1.5"单调涨的，连涨十几次就接近上限 5e3，
        // 步长已经小到 1e-3 量级，后面几十次迭代只是在原地磨 —— 这类种子本来就收敛不了。
        // 提前退出对判据**无副作用**：
        //   · 未收敛的解 ok=false，进不了 valid()，也永远排不到 rankCmp 的已收敛解之前；
        //   · bestErr 只用于调阻尼，不参与返回，th 仍是当前迭代值（和跑满 90 次一样是它）。
        // 门槛 bestErr > 0.5（加权单位，≈ 位置 0.5mm）保证"已经快收敛的不动"。
        if (stall >= 8 && bestErr > 0.5) break;

        double J[6][6];
        double JJt[6][6];
        jacobian(f, J);
        for (int i = 3; i < 6; ++i)
            for (int j = 0; j < 6; ++j) J[i][j] *= wRot;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                double s = 0.0;
                for (int k = 0; k < 6; ++k) s += J[i][k] * J[j][k];
                JJt[i][j] = s + (i == j ? lam * lam : 0.0);
            }
        double y[6];
        if (!solve6(JJt, err, y)) break;
        double dth[6];
        for (int i = 0; i < 6; ++i) {
            double s = 0.0;
            for (int k = 0; k < 6; ++k) s += J[k][i] * y[k];
            dth[i] = s;
        }
        // 单步限幅 30°：线性化在限位附近会失效，不夹住会来回震荡
        double mx = 0.0;
        for (int i = 0; i < 6; ++i) mx = std::max(mx, std::fabs(dth[i]));
        if (mx > 30.0 * kDeg) {
            double s = 30.0 * kDeg / mx;
            for (int i = 0; i < 6; ++i) dth[i] *= s;
        }
        for (int i = 0; i < 6; ++i) th[i] += dth[i] * kRad;
    }
    if (g_prof.on) g_prof.ikIters += iters;   // 常开会在 8 线程间抢同一条缓存行

    // 退出循环时 th 已经走了一步（最后一轮的 dth 已加），所以 f 是过期的，必须重算。
    fk(th, toolLen, f);
    IkSol sol;
    for (int i = 0; i < 6; ++i) sol.theta[i] = th[i];
    sol.posErr = norm(pTarget - f.tcp.t);
    sol.rotErrDeg = norm(matToRotvec(m3mul(RTarget, m3T(f.tcp.R)))) * kRad;
    sol.ok = (sol.posErr <= std::max(0.02 * 20, 0.5) && sol.rotErrDeg <= std::max(0.02 * 20, 0.5));
    sol.over = overLimit(th);
    sol.margin = limitMarginMin(th);
    // cond 不在这里算：见 IkSol::cond 的注释（19 个候选里 18 个会被丢弃）。
    return sol;
}

inline int rankCmp(const IkSol& a, const IkSol& b) {
    // 收敛 > 不限位 > 单轴超出总量少 > 残差小 > 余量大
    auto key = [](const IkSol& s, double& k3, double& k4, double& k5) {
        double overSum = 0.0;
        if (s.margin < 0.0) overSum = -s.margin;
        k3 = overSum;
        k4 = s.posErr + 100.0 * s.rotErrDeg;
        k5 = -s.margin;
    };
    int ka = (a.ok ? 0 : 1), kb = (b.ok ? 0 : 1);
    if (ka != kb) return ka < kb ? -1 : 1;
    int oa = (a.over ? 1 : 0), ob = (b.over ? 1 : 0);
    if (oa != ob) return oa < ob ? -1 : 1;
    double a3, a4, a5, b3, b4, b5;
    key(a, a3, a4, a5);
    key(b, b3, b4, b5);
    if (a3 != b3) return a3 < b3 ? -1 : 1;
    if (a4 != b4) return a4 < b4 ? -1 : 1;
    if (a5 != b5) return a5 < b5 ? -1 : 1;
    return 0;
}

// 关节空间距离（平方，度²）：用来在多个合法 IK 解里挑“离上一点最近”的那支。
// 只用于比较，不开方。
inline double jointDist2(const double a[6], const double b[6]) {
    double s = 0.0;
    for (int i = 0; i < 6; ++i) {
        double d = a[i] - b[i];
        s += d * d;
    }
    return s;
}

inline double maxJointDelta(const double a[6], const double b[6]) {
    double m = 0.0;
    for (int i = 0; i < 6; ++i) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
}

// 上一点解如果单关节跳超过这个角度，就怀疑换了 IK 解支；再多试解析种子，
// 选离上一点最近的合法解。低于它直接返回，绝大多数点不额外花时间。
constexpr double kBranchContinuityTrigDeg = 10.0;
// 阶段 3 复用阶段 1 缓存的连续性阈值。比 IK 解支阈值更严：
// 缓存是“另一条轨迹”的解，只要和最终上一点差得多就重解，避免 5~10° 的隐形换支。
constexpr double kCacheContinuityTrigDeg = 5.0;

// 种子顺序：上一点解（连续性）→ 解析种子（解支正确）→ 固定种子表（兜底）
IkSol ik(const V3& pTarget, const M3& RTarget, double toolLen, const double* seed) {
    if (g_prof.on) ++g_prof.ik;   // 仅诊断计数：常开时 8 线程抢同一条缓存行，实测有 2~3% 开销
    double seeds[1 + 16 + 12][6];
    int ns = 0;
    if (seed) {
        for (int i = 0; i < 6; ++i) seeds[ns][i] = seed[i];
        ++ns;
    }
    ns += analyticSeeds(pTarget, RTarget, toolLen, seeds + ns);
    for (int i = 0; i < 12; ++i) {
        for (int k = 0; k < 6; ++k) seeds[ns][k] = kSeeds[i][k];
        ++ns;
    }

    auto valid = [](const IkSol& s) {
        return s.ok && s.posErr < 0.05 && s.rotErrDeg < 0.05 && !s.over;
    };
    IkSol best;
    bool hasBest = false;
    IkSol bestValid;
    bool hasValid = false;
    double bestValidDist = 1e300;
    // 本点是否走的是"种子收敛但超限"这条最贵的路，以及最终有没有换来合法分支。
    // 只用于决定要不要继续压缩这条路径，不参与任何判据。
    bool seedOver = false;
    double seedMargin = 0.0;
    int nOnce = 0;

    // 先跑上一点种子：连续且合法就直接返回，保持主路径速度。
    if (seed) {
        IkSol sol = ikOnce(pTarget, RTarget, toolLen, seed);
        if (g_prof.on) ++g_prof.ikOnceN;
        ++nOnce;
        best = sol;
        hasBest = true;
        if (valid(sol)) {
            if (maxJointDelta(sol.theta, seed) <= kBranchContinuityTrigDeg) {
                if (g_prof.on) ++g_prof.ikFast;
                return sol;
            }
            // 跳得大：先把它作为“最近的合法解”的兜底，再找有没有更近的解支。
            bestValid = sol;
            hasValid = true;
            bestValidDist = jointDist2(sol.theta, seed);
        } else if (sol.ok) {
            // 几何收敛但关节超限：后面的 28 组种子就是为了找"换个支不超限"。
            if (g_prof.on) ++g_prof.ikOver;
            seedOver = true;
            seedMargin = sol.margin;
        } else {
            if (g_prof.on) ++g_prof.ikFail;
        }
    }

    // 解析种子 / 固定种子。有 seed 且上一点解跳得大时，选关节距离最近的合法解；
    // 没有 seed 时保持旧行为（第一个合法解直接返回）。
    for (int i = (seed ? 1 : 0); i < ns; ++i) {
        // 诊断：这组种子是否与前面某组逐位相同（相同则结果必然相同，跑它也白跑）
        if (g_prof.on) {
            ++g_prof.sweepSeeds;
            for (int q = 0; q < i; ++q) {
                bool same = true;
                for (int k = 0; k < 6; ++k)
                    if (seeds[q][k] != seeds[i][k]) { same = false; break; }
                if (same) { ++g_prof.sweepDup; break; }
            }
        }
        IkSol sol = ikOnce(pTarget, RTarget, toolLen, seeds[i]);
        if (g_prof.on) ++g_prof.ikOnceN;
        ++nOnce;
        if (valid(sol)) {
            if (!seed) { if (g_prof.on) ++g_prof.ikFast; return sol; }
            double dist = jointDist2(sol.theta, seed);
            if (!hasValid || dist < bestValidDist) {
                bestValid = sol;
                bestValidDist = dist;
                hasValid = true;
            }
        }
        if (!hasBest || rankCmp(sol, best) < 0) {
            best = sol;
            hasBest = true;
        }
    }
    if (seedOver) {
        // 「超限 → 全种子搜」这条路的投入产出：跑了多少 ikOnce、有没有换来合法解。
        if (g_prof.on) {
            g_prof.ikOverOnce += nOnce;
            if (hasValid) ++g_prof.ikOverHit; else ++g_prof.ikOverMiss;
            g_prof.bucketOver(seedMargin, hasValid);
        }
    }
    if (hasValid) return bestValid;
    return best;
}

// ---------------------------------------------------------------- 碰撞体（胶囊）

struct Sizing {
    double jointR, joint1R, linkR, wristLinkR, baseColumnR, baseNeckR;
    double gripR, toolPipeR;
};

Sizing sizingFrom(double reach) {
    Sizing s;
    s.jointR = clampd(reach * 0.038, 12.0, 46.0);
    s.joint1R = s.jointR * 1.25;
    s.linkR = s.jointR * 0.62;
    s.wristLinkR = s.jointR * 0.72;
    s.baseColumnR = s.jointR * 2.6;
    s.baseNeckR = s.jointR * 1.15;
    s.gripR = s.jointR * 0.62;
    s.toolPipeR = s.gripR * 0.45;
    return s;
}

constexpr double kGripFraction = 0.72;    // 握把占「法兰→TCP」全长比例

enum SegGroup { SG_TOOL = 0, SG_LINK = 1, SG_JOINT = 2, SG_BASE = 3 };

struct Seg {
    V3 a, b;
    double r;
    int group;
    int j0, j1;    // 关联关节号（自碰撞排除用；-1 表示无）
    const char* part;
    int kind;      // RS_SEG_*，用于把「最近部位」报成人能看懂的名字
};

// 只造工具段（握把 + 细管）。
//
// **枪长 toolLen 已经把喷嘴包含在内**：toolLen 的定义是「法兰面 → 喷嘴出口」，喷嘴出口
// 就是 TCP。所以枪体几何只有「法兰 → TCP」这 450mm 一段实体，**TCP 之外没有任何东西**。
//
// 这里踩过一次严重的坑：曾在 TCP 之外再加一段 "喷嘴" 胶囊（tip = tcp + u*(L*0.18)），
// 等于沿枪轴朝工件又伸了 0.18×450 = 81mm，半径 8.56mm。喷距 standoff=100 是从 TCP 量的，
// 于是 100mm 里有 81mm 被这段凭空的胶囊吃掉，剩下 19mm 再减半径——
// 槽型盘悬挂 1467 点报出「枪体撞 398 点、最小表面距离 −5.33mm」，而 3D 视图里枪离工件
// 明明有 100mm。改成不虚构后，同一份数据最差间隙 +69.39mm、0 穿透，完全印证那 398 个
// P0 全是假阳性。不要再把"喷嘴"当成枪长之外的一段加回来。
//
// 为什么单独抽出来：**枪体的位置由「目标 TCP + 枪轴」解析确定，与机器人型号/基座/IK 全都无关**
// —— 法兰就在 TCP 沿枪轴退 toolLen 的地方。所以"枪会不会戳进工件"这一条在
// 「要求可达」和「不要求可达」两种模式下是同一个判据、同一份几何，不依赖 IK 解。
// buildSegments 也从这里取工具段，避免两处各写一份（改一处漏一处的老坑）。
// 用「法兰点 + TCP 点」而不是「TCP + 轴向 + 长度」：FK 那条路必须原样用 f.flange.t，
// 自己按轴向反推法兰会引入 1ulp 差异，等距时最近的段名会串（工具细管 ↔ 工具握把），
// 让"改动前后逐字节一致"这条回归红线失效。
void buildToolSegs(const V3& flange, const V3& tcp, const Sizing& sz, std::vector<Seg>& out) {
    V3 d = tcp - flange;
    double L = norm(d);
    if (L <= 0.5) return;
    V3 u = d * (1.0 / L);
    V3 neck = flange + u * (L * kGripFraction);   // 握把末端：枪身由粗变细的地方
    out.push_back({flange, neck, sz.gripR, SG_TOOL, -1, -1, "工具握把", RS_SEG_GRIP});
    out.push_back({neck, tcp, sz.toolPipeR, SG_TOOL, -1, -1, "工具细管", RS_SEG_PIPE});
}

void buildSegments(const FkOut& f, const Sizing& sz, bool includeBase, std::vector<Seg>& out) {
    out.clear();
    if (includeBase) {
        double baseH = std::max(4.0, f.origins[1].z * 0.55);
        out.push_back({V3(0, 0, 0), V3(0, 0, baseH), sz.baseColumnR, SG_BASE, -1, -1, "底座立柱", RS_SEG_BASE_COL});
        out.push_back({V3(0, 0, baseH), f.origins[1], sz.baseNeckR, SG_BASE, -1, -1, "底座颈", RS_SEG_BASE_NECK});
    }
    for (int i = 1; i <= 6; ++i) {
        double jr = (i == 1) ? sz.joint1R : sz.jointR;
        out.push_back({f.origins[i], f.origins[i], jr, SG_JOINT, i, i, "关节", RS_SEG_JOINT});
        if (i + 1 < 7) {
            double rr = (i >= 4) ? sz.wristLinkR : sz.linkR;
            out.push_back({f.origins[i], f.origins[i + 1], rr, SG_LINK, i, i + 1, "连杆", RS_SEG_LINK});
        }
    }
    // 工具段：法兰→TCP 的方向就是工具 +Z，展开成握把/细管/喷嘴三段
    buildToolSegs(f.flange.t, f.tcp.t, sz, out);
}

// 部位名（含关节号），给控制台/CSV/界面用
std::string segName(int kind, int j0, int j1) {
    char buf[48];
    switch (kind) {
        case RS_SEG_GRIP:      return "工具握把";
        case RS_SEG_PIPE:      return "工具细管";
        case RS_SEG_BASE_COL:  return "底座立柱";
        case RS_SEG_BASE_NECK: return "底座颈";
        case RS_SEG_JOINT:
            std::snprintf(buf, sizeof(buf), "关节J%d", j0);
            return buf;
        case RS_SEG_LINK:
            std::snprintf(buf, sizeof(buf), "连杆J%d-J%d", j0, j1);
            return buf;
        default: return "—";
    }
}

// 自碰撞段对：非相邻、且不共享关节号（同一关节球与它的两根连杆本来就贴在一起）
void buildSelfPairs(const std::vector<Seg>& segs, std::vector<std::pair<int, int>>& out) {
    out.clear();
    std::vector<int> idx;
    for (int i = 0; i < (int)segs.size(); ++i)
        if (segs[i].group == SG_LINK || segs[i].group == SG_JOINT) idx.push_back(i);
    for (size_t a = 0; a < idx.size(); ++a) {
        for (size_t b = a + 1; b < idx.size(); ++b) {
            const Seg& sa = segs[idx[a]];
            const Seg& sb = segs[idx[b]];
            bool share = false;
            if (sa.j0 >= 0 && (sa.j0 == sb.j0 || sa.j0 == sb.j1)) share = true;
            if (sa.j1 >= 0 && (sa.j1 == sb.j0 || sa.j1 == sb.j1)) share = true;
            if (share) continue;
            out.push_back({idx[a], idx[b]});
        }
    }
}

// ================================================================ 点云（KD-tree）

struct KdNode {
    float mn[3], mx[3];
    int left = -1, right = -1;
    int start = 0, count = 0;   // count>0 → 叶
};

// 体素降采样：每格保留一个代表点。
// 碰撞判定只关心"表面在哪"，5mm 体素相对 20mm 安全间隙是二阶小量；而 60 万点降到 2 万点，
// KD 构建从 500ms 降到 20ms（big 样本实测：KD 占了整条流程的 70%）。
// 想保精度就不开；产线预判开 3~5mm。
void voxelDecimate(const double* xyz, int n, double voxel, std::vector<double>& out) {
    out.clear();
    if (n <= 0 || voxel <= 1e-6) return;
    const double inv = 1.0 / voxel;
    struct Entry { long long key; int idx; };
    std::vector<Entry> ents;
    ents.reserve((size_t)n);
    for (int i = 0; i < n; ++i) {
        long long ix = (long long)std::floor(xyz[i * 3 + 0] * inv);
        long long iy = (long long)std::floor(xyz[i * 3 + 1] * inv);
        long long iz = (long long)std::floor(xyz[i * 3 + 2] * inv);
        // 各 21 位打包（±100 万格，5mm 体素下覆盖 ±5km，足够）
        long long key = ((ix & 0x1FFFFF) << 42) | ((iy & 0x1FFFFF) << 21) | (iz & 0x1FFFFF);
        ents.push_back({key, i});
    }
    // 排序比哈希快：60 万条 12 字节的记录排序 ~60ms，哈希插入反而更慢
    std::sort(ents.begin(), ents.end(), [](const Entry& a, const Entry& b) {
        return a.key != b.key ? a.key < b.key : a.idx < b.idx;
    });
    out.reserve((size_t)n * 3 / 4);
    long long prev = 0;
    bool first = true;
    for (const auto& e : ents) {
        if (!first && e.key == prev) continue;
        prev = e.key;
        first = false;
        out.push_back(xyz[e.idx * 3 + 0]);
        out.push_back(xyz[e.idx * 3 + 1]);
        out.push_back(xyz[e.idx * 3 + 2]);
    }
}

class CloudIndex {
public:
    void build(const double* xyz, int n) {
        px_.resize(n);
        py_.resize(n);
        pz_.resize(n);
        for (int i = 0; i < n; ++i) {
            px_[i] = (float)xyz[i * 3 + 0];
            py_[i] = (float)xyz[i * 3 + 1];
            pz_[i] = (float)xyz[i * 3 + 2];
        }
        raw_ = n;
        idx_.resize(n);
        for (int i = 0; i < n; ++i) idx_[i] = i;
        nodes_.clear();
        nodes_.reserve(2 * (size_t)n / 8 + 16);
        if (n > 0) root_ = buildNode(0, n);
        // 质心 / AABB
        double sx = 0, sy = 0, sz = 0;
        for (int i = 0; i < n; ++i) {
            sx += px_[i];
            sy += py_[i];
            sz += pz_[i];
        }
        if (n > 0) {
            centroid = V3(sx / n, sy / n, sz / n);
            bmin = V3(*std::min_element(px_.begin(), px_.end()),
                      *std::min_element(py_.begin(), py_.end()),
                      *std::min_element(pz_.begin(), pz_.end()));
            bmax = V3(*std::max_element(px_.begin(), px_.end()),
                      *std::max_element(py_.begin(), py_.end()),
                      *std::max_element(pz_.begin(), pz_.end()));
        }
    }

    int size() const { return (int)px_.size(); }
    V3 point(int i) const { return V3(px_[i], py_[i], pz_[i]); }

    // 线段到点云的最短距离；超过 radius 直接返回 radius（= 只报下界）
    double segmentDist(const V3& a, const V3& b, double radius) const {
        if (nodes_.empty()) return radius;
        if (g_prof.on) ++g_prof.seg;
        double best = radius;
        queryRec(root_, a, b, best);
        return best;
    }

    // 平均点间距（抽样 200 个点的最近邻距离中位数）—— 给"沿枪轴投影"的管道半径用
    double estimateSpacing() const {
        int n = size();
        if (n < 8) return 1.0;
        int k = std::min(200, n);
        std::vector<double> ds;
        ds.reserve(k);
        for (int s = 0; s < k; ++s) {
            int i = (int)((long long)s * (n - 1) / (k - 1));
            V3 p = point(i);
            double best = 1e300;
            // 小规模暴力 + 邻域剪枝：只需近似值，不值得另建一棵树
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                double dx = px_[j] - p.x, dy = py_[j] - p.y, dz = pz_[j] - p.z;
                double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < best) best = d2;
            }
            ds.push_back(std::sqrt(best));
        }
        std::sort(ds.begin(), ds.end());
        return ds[ds.size() / 2];
    }

    V3 centroid;
    V3 bmin, bmax;
    int raw_ = 0;

    std::vector<float> px_, py_, pz_;

private:
    int buildNode(int start, int count) {
        int ni = (int)nodes_.size();
        nodes_.push_back(KdNode());
        float mn[3] = {1e30f, 1e30f, 1e30f}, mx[3] = {-1e30f, -1e30f, -1e30f};
        for (int i = start; i < start + count; ++i) {
            int p = idx_[i];
            float c[3] = {px_[p], py_[p], pz_[p]};
            for (int k = 0; k < 3; ++k) {
                mn[k] = std::min(mn[k], c[k]);
                mx[k] = std::max(mx[k], c[k]);
            }
        }
        for (int k = 0; k < 3; ++k) {
            nodes_[ni].mn[k] = mn[k];
            nodes_[ni].mx[k] = mx[k];
        }
        if (count <= 12) {
            nodes_[ni].start = start;
            nodes_[ni].count = count;
            return ni;
        }
        int axis = 0;
        double ext = -1;
        for (int k = 0; k < 3; ++k) {
            double e = (double)mx[k] - mn[k];
            if (e > ext) {
                ext = e;
                axis = k;
            }
        }
        int mid = count / 2;
        auto cmp = [&](int a, int b) {
            float ca = axis == 0 ? px_[a] : (axis == 1 ? py_[a] : pz_[a]);
            float cb = axis == 0 ? px_[b] : (axis == 1 ? py_[b] : pz_[b]);
            return ca < cb;
        };
        std::nth_element(idx_.begin() + start, idx_.begin() + start + mid,
                         idx_.begin() + start + count, cmp);
        int l = buildNode(start, mid);
        int r = buildNode(start + mid, count - mid);
        nodes_[ni].left = l;
        nodes_[ni].right = r;
        return ni;
    }

    // 节点包围盒到线段的距离下界。两个下界取大者：
    //   (a) AABB(线段) vs 节点 AABB —— 线段在自己的 AABB 里，所以它 ≤ 真值；
    //   (b) dist(节点中心, 线段) − 半对角 —— 盒内任一点距中心不超过半对角。
    // 长线段时 (a) 很松，(b) 补上；两个都成立，取 max 合法。
    static double lowerBound(const KdNode& n, const V3& A, const V3& B) {
        double gap2 = 0.0;
        for (int k = 0; k < 3; ++k) {
            double smin = std::min(A[k], B[k]), smax = std::max(A[k], B[k]);
            if (smax < n.mn[k]) gap2 += sq((double)n.mn[k] - smax);
            else if (smin > n.mx[k]) gap2 += sq(smin - (double)n.mx[k]);
        }
        double lb1 = std::sqrt(gap2);
        V3 c((n.mn[0] + n.mx[0]) * 0.5, (n.mn[1] + n.mx[1]) * 0.5, (n.mn[2] + n.mx[2]) * 0.5);
        double dc = pointSegDist(c, A, B);
        double half = 0.5 * std::sqrt(sq((double)n.mx[0] - n.mn[0]) + sq((double)n.mx[1] - n.mn[1]) +
                                      sq((double)n.mx[2] - n.mn[2]));
        double lb2 = dc - half;
        return std::max(lb1, std::max(lb2, 0.0));
    }

    void queryRec(int node, const V3& A, const V3& B, double& best) const {
        if (best <= 0.0) return;
        const KdNode& nd = nodes_[node];
        if (lowerBound(nd, A, B) >= best) return;
        if (nd.count > 0) {
            for (int i = nd.start; i < nd.start + nd.count; ++i) {
                int p = idx_[i];
                double d = pointSegDist(V3(px_[p], py_[p], pz_[p]), A, B);
                if (d < best) best = d;
            }
            return;
        }
        double lbL = nd.left >= 0 ? lowerBound(nodes_[nd.left], A, B) : 1e300;
        double lbR = nd.right >= 0 ? lowerBound(nodes_[nd.right], A, B) : 1e300;
        if (lbL <= lbR) {
            if (nd.left >= 0) queryRec(nd.left, A, B, best);
            if (nd.right >= 0) queryRec(nd.right, A, B, best);
        } else {
            if (nd.right >= 0) queryRec(nd.right, A, B, best);
            if (nd.left >= 0) queryRec(nd.left, A, B, best);
        }
    }

    std::vector<int> idx_;
    std::vector<KdNode> nodes_;
    int root_ = -1;
};

// ================================================================ 轨迹输入

struct TrajIn {
    std::vector<V3> pos;    // TCP
    std::vector<V3> gun;    // 枪轴朝向（指向工件 = -n）
    std::vector<unsigned char> spray;
    std::vector<double> speed;
    int n() const { return (int)pos.size(); }
    V3 normal(int i) const { return -gun[i]; }
};

// 「行」的定义：**一笔**（走完一趟、掉头之前的那段连续运动），不是"一段连续 spray"。
//
// 这个区别是致命的：母版生成器跨折返**不关喷**，所以 longstripe 的 45 笔在 CSV 里只有
// 2 个 spray 连续段（11044 + 511 点）。按 spray 段当行，等于把 22 笔当成 1 行 ——
// 方位角只能取一个（22 笔朝向不同的壁共用一个进刀侧）；姿态又是面组级常量，跨折返
// 合并会把两笔朝向不同的壁强行统一成一个姿态（上一笔的高倾角把下一笔整段抬起来）。
//
// 折返判据用**水平投影**方向掉头：与"本笔起始方向"夹角 > 120°（dot < −0.5）即断行。
// 不用三维方向是因为曲面/法向变化会让三维步长抖动；而折返在水平面上一定是掉头。
// 断点落在折返圆弧里（≈120° 处），正是两笔之间该断开的位置。
// 阈值取 120° 而不是 90°：曲面上的笔画本身会连续转向，90° 会把一笔正常弯道切成两行。
void computeRows(const TrajIn& t, std::vector<std::pair<int, int>>& rows) {
    rows.clear();
    const double kRevCos = -0.5;      // cos120°
    const int kMinRowPts = 4;         // 短于此不再断，避免折返圆弧里连着断出碎片行
    int n = t.n(), i = 0;
    while (i < n) {
        if (!t.spray[i]) {
            ++i;
            continue;
        }
        int a = i, j = i;
        double hx = 0.0, hy = 0.0;
        bool hasH = false;
        while (j + 1 < n && t.spray[j + 1]) {
            V3 st = t.pos[j + 1] - t.pos[j];
            double L = std::sqrt(st.x * st.x + st.y * st.y);
            if (L > 1e-9) {
                double ux = st.x / L, uy = st.y / L;
                if (!hasH) {
                    hx = ux; hy = uy; hasH = true;
                } else if (ux * hx + uy * hy < kRevCos && (j - a + 1) >= kMinRowPts) {
                    rows.push_back({a, j});
                    a = j + 1;
                    hx = ux; hy = uy;      // 新一笔的起始方向
                    ++j;
                    continue;
                }
            }
            ++j;
        }
        rows.push_back({a, j});
        i = j + 1;
    }
}

// ================================================================ 评价器

struct Params {
    double standoff = 100.0;
    double safety = 20.0;          // 间隙告警阈值（软约束，只影响 P1 报告）
    double hardClearance = 2.0;    // 绝对安全距离（硬约束，低于它才回退）
    double search = 120.0;
    double toolLen = kDefaultToolLen;
    double tiltNonUpDeg = 60.0;      // δ=90°（竖直侧壁）时的倾角，实际按 δ/90 线性缩放（2026-10-26 由 30 改 60）
    double upThresholdDeg = 45.0;    // 保留字段（倾角已改线性，不再依赖阈值判"非朝上"）
    double tiltMaxDeg = 60.0;
    double tiltStepDeg = 5.0;
    double tiltRateDeg = 3.0;
    double maxSprayAngleDeg = 30.0;
    double jointMarginMin = 1.0;
    double condMax = 30.0;
    double transitLiftMax = 400.0;
    double voxelMm = 0.0;
    int linkStride = 1;
    int threads = 1;
    bool azStance = true;    // true = 按机器人站位逐行择优；false = 只用「向上」（旧口径）
    bool includeBase = false;
    bool fixTransit = true;
    bool azFallback = true;
    // false = 不要求可达：跳过 IK / 连杆 / 自碰撞 / 关节限位，只留几何层。
    // 判据分家见 robot_safety.h 里 requireReach 的注释。
    bool requireReach = true;
    // true = 单次逆解（fast）：姿态纯解析、每点只解 1 次 IK、不回退不救援。
    // 与 requireReach 正交：fast + requireReach 就是"查安全、但只解一次"。
    bool fastPose = false;
    // true = 倾角逐点连续：倾角只由该点法向决定（tilt = tiltNonUpDeg × δ/90），
    // 不做面级常量。见 robot_safety.h 里 RS_POSE_TILT_PERP 的说明。
    bool tiltPerPoint = false;
    // true = 逐点解析救 J5：倾角与方位都由「前臂方向 ⊕ 法向」解析反推，
    // 目标是让 |θ5| 留够 jointMarginMin 余量。见 robot_safety.h 里 RS_POSE_TILT_J5 的说明。
    // 与 tiltPerPoint 同时置位时优先本模式。
    bool tiltJ5 = false;
    // 点云系里的「重力向上」单位向量。**默认 (0,0,1) = 旧口径，逐字节不变**。
    //
    // 为什么需要显式给：倾角口径里"正位置/朝上"是相对**重力**定义的，而点云是相机系，
    // 相机 +Z 是光轴不是重力。举本次实测的煤气灶背板：主面法向与点云 +Z 夹角只有 5°，
    // 旧口径判它"基本朝上、只需 3.3° 倾角"；但工件是立着挂的，真实 α=86°，
    // 需求 57.5°。差的不是精度，是整条口径的前提。
    V3 upW = V3(0, 0, 1);
    // 倾角基准方向（点云系单位向量）。**不是**重力，见上面 tiltAxisAndDir 前的说明。
    // 默认与 upW 同值，基座位姿已知时在 impl 里改成基座 ±X。**每点倾角只读这一个**。
    V3 tiltRefW = V3(0, 0, 1);
    // 平放面倾角方向（基座系）：0=不指定（按法向自动选 ±X）；+1=往基座 +X 倾；−1=往基座 −X 倾。
    // 只在基座位姿已知时生效，见 impl 里对 tiltRefW 的赋值。
    double tiltSideX = 0.0;
};

struct PR {
    int idx = -1;
    double tiltDeg = 0.0;
    double azDeg = 0.0;
    double clearance = 0.0;   // min(clrTool, clrLink)
    double clrTool = 1e300;   // 仅枪体
    double clrLink = 1e300;   // 仅连杆/关节/底座
    double selfGap = 1e300;
    double margin = 0.0;
    double cond = 0.0;
    double standoffErr = 0.0;
    double lifted = 0.0;
    double theta[6] = {0, 0, 0, 0, 0, 0};
    V3 pos, gun;
    double sprayAngle = 0.0;
    bool ikOk = false;
    bool over = false;
    bool reverted = false;
    unsigned int viol = 0;
    const char* part = "";
    int partKind = -1, partJ0 = -1, partJ1 = -1;   // 全局最近部位（RS_SEG_*）
    int linkKind = -1, linkJ0 = -1, linkJ1 = -1;   // 连杆/关节里最近的那一段
};

// 位姿优劣（越大越好）：自碰撞通过 > 间隙大 > 关节余量大。
// 碰撞(P0) 永远压过姿态连续(P2)，所以 selfGap<0 直接判负。
// requireReach=false 时没有自碰撞可言，直接比间隙。
inline int poseCmp(const PR& a, const PR& b, bool requireReach) {
    if (requireReach) {
        int oa = (a.selfGap >= 0.0) ? 1 : 0, ob = (b.selfGap >= 0.0) ? 1 : 0;
        if (oa != ob) return oa < ob ? -1 : 1;
    }
    if (a.clearance != b.clearance) return a.clearance < b.clearance ? -1 : 1;
    if (requireReach && a.margin != b.margin) return a.margin < b.margin ? -1 : 1;
    return 0;
}

// 硬约束：IK 成功、不自碰撞、间隙不低于绝对安全距离、关节不超限。
// 注意 margin < jointMarginMin 只是 P1 软告警，不在这里判死；over（真正超限）才是硬约束。
// requireReach=false：只剩「间隙不低于绝对安全距离」——这里的 clearance 此时就是
// min(枪体)，因为连杆/关节根本没建碰撞体（见 eval 的 requireReach 分支）。
inline bool hardOk(const PR& q, const Params& pp) {
    if (!(q.clearance >= pp.hardClearance)) return false;
    if (!pp.requireReach) return true;
    return q.ikOk && q.selfGap >= 0.0 && !q.over;
}

// 两个“都不满足硬约束”的位姿二选一：先比硬约束通过项，再比自间隙/间隙/关节余量。
// 与 poseCmp 的区别：selfGap<0 时还要比谁的负值更小；关节真正超限（over）
// 也单独判负，不能因为另一个位姿间隙大就选一个关节超限的位姿。
inline int hardCmp(const PR& a, const PR& b, bool requireReach) {
    if (requireReach) {
        int ia = a.ikOk ? 1 : 0, ib = b.ikOk ? 1 : 0;
        if (ia != ib) return ia < ib ? -1 : 1;
        int sa = (a.selfGap >= 0.0) ? 1 : 0, sb = (b.selfGap >= 0.0) ? 1 : 0;
        if (sa != sb) return sa < sb ? -1 : 1;
        int oa = a.over ? 0 : 1, ob = b.over ? 0 : 1;
        if (oa != ob) return oa < ob ? -1 : 1;
        if (a.selfGap != b.selfGap) return a.selfGap < b.selfGap ? -1 : 1;
    }
    if (a.clearance != b.clearance) return a.clearance < b.clearance ? -1 : 1;
    if (requireReach && a.margin != b.margin) return a.margin < b.margin ? -1 : 1;
    return 0;
}

class Evaluator {
public:
    Evaluator(const CloudIndex& cloud, const Pose& base, const Params& p)
        : cloud_(cloud), p_(p) {
        baseR_ = base.R;
        baseRT_ = m3T(base.R);
        baseT_ = base.t;
        sz_ = sizingFrom(kReach);
        double z0[6] = {0, 0, 0, 0, 0, 0};
        FkOut f;
        fk(z0, p.toolLen, f);
        std::vector<Seg> segs;
        buildSegments(f, sz_, p.includeBase, segs);
        buildSelfPairs(segs, pairs_);
    }

    V3 toBaseP(const V3& p) const { return m3mul(baseRT_, p - baseT_); }
    // 方向向量：点云系 → 基座系（只转不平移）。toCloudR 的反向。
    V3 toBaseV(const V3& v) const { return m3mul(baseRT_, v); }
    M3 toBaseR(const M3& R) const { return m3mul(baseRT_, R); }
    V3 toCloudP(const V3& p) const { return m3mul(baseR_, p) + baseT_; }
    // 方向向量：基座系 → 点云系（只有旋转，不平移）。refX 这类「来自 FK 的轴」必须过这一步。
    V3 toCloudR(const V3& v) const { return m3mul(baseR_, v); }
    // 基座原点（点云系）。方位决策要用它判"机器人从哪一侧够过来"。
    V3 baseOrigin() const { return baseT_; }
    // 基座→点云系的旋转矩阵。逐点救 J5 要把基座系里的前臂方向转到点云系。
    const M3& baseR() const { return baseR_; }

    // 给定枪身朝外方向 d（倾斜后的法向 n'），求 IK + 碰撞 + 自碰撞。
    // **工具 +Z 是 −d**：d 是从工件指向枪的外法向，工具 +Z 是喷射方向（法兰→喷嘴）= −d。
    // 这条搞反过一整轮：位置全对、喷距守恒，但 IK 只能选"反腕 + 大臂下探"的解
    // （θ2=154°、J2 超限 19°、小臂从工件底下穿过去 24.7mm）。
    bool eval(const V3& face, const V3& d, const V3* refX, const double* seed,
              double tiltDeg, double azDeg, bool checkLinks, PR& out) const {
        if (g_prof.on) ++g_prof.evals;
        V3 zTool = -unit(d);
        V3 pC = face + d * p_.standoff;

        // ---- 不要求可达：不解 IK，只查枪体与点云的间隙
        //
        // 这一步是"点云不适配这台机器人"时的唯一正确口径：枪体的三个碰撞体由
        // 「TCP + 枪轴」解析确定（法兰 = TCP − toolLen·枪轴），**与机器人无关**，
        // 所以它不是"省略了检查"，而是把与型号无关的那部分单独算准。
        // 不解 IK 的副作用全是正向的：没有解支噪声 → 面级姿态不会被"这个点冷启解不出"
        // 带偏；省掉 8 组种子 × 90 次迭代 → 19k 点从半分钟降到几秒。
        if (!p_.requireReach) {
            std::vector<Seg> segs;
            // 碰撞段一律建在**基座系**：后面的 toCloudP 只负责 基座系→点云系。
            // 这里踩过一次：直接把点云系的 TCP/枪轴喂进 buildToolSegs，再被 toCloudP
            // 变换一次，等于整根枪被平移了一个基座偏移（自动基座是 −252,150,−197），
            // 于是"枪体撞了 1512 点、最深 −21.2mm"，全是假的。
            buildToolSegs(toBaseP(pC) - toBaseV(zTool) * p_.toolLen, toBaseP(pC), sz_, segs);
            double worstTool = 1e300;
            const Seg* segTool = nullptr;
            Prof::Scope _p(g_prof.on, g_prof.nsSeg);
            for (const Seg& s : segs) {
                double limit = std::min(p_.search, (worstTool < 1e299) ? worstTool + s.r : p_.search);
                double surface = cloud_.segmentDist(toCloudP(s.a), toCloudP(s.b), limit) - s.r;
                if (surface < worstTool) { worstTool = surface; segTool = &s; }
            }
            out.idx = -1;
            out.tiltDeg = tiltDeg;
            out.azDeg = azDeg;
            out.clearance = worstTool;
            out.clrTool = worstTool;
            out.clrLink = 1e300;              // 没建连杆碰撞体，别报 0 让人以为量到了
            out.part = segTool ? segTool->part : "";
            out.partKind = segTool ? segTool->kind : -1;
            out.partJ0 = segTool ? segTool->j0 : -1;
            out.partJ1 = segTool ? segTool->j1 : -1;
            out.linkKind = -1; out.linkJ0 = -1; out.linkJ1 = -1;   // -1 = 没查连杆
            out.ikOk = false;                 // 如实：没有解过 IK
            out.over = false;
            out.margin = 0.0;
            out.cond = 0.0;
            for (int i = 0; i < 6; ++i) out.theta[i] = 0.0;
            out.pos = pC;
            out.gun = -d;
            out.sprayAngle = 0.0;
            out.standoffErr = std::fabs(norm(pC - face) - p_.standoff);
            out.selfGap = 1e300;
            return true;
        }

        M3 Rc = frameFromZ(zTool);
        if (refX) {
            // refX 来自上一点 fk() 的工具 x 轴，是**基座系**的；zTool 是**点云系**的。
            // 不转坐标系直接投影 = 拿两个不同坐标系的向量做正交分解，工具 x 轴每步都会被
            // 额外偏转一份「基座偏航」的角量，并且逐点累加 → 手腕绕枪轴持续自旋。
            // 实测（complex，自动基座偏航 18.08°）：J6 每点 -17.6°、整行拧到 J6 限位；
            // 而 tray（自动基座偏航 0）同一段代码完全平滑（每点 0.5°）——正是坐标系混用的铁证。
            V3 rx = toCloudR(*refX);
            V3 x0 = rx - zTool * dot(rx, zTool);
            if (norm(x0) > 1e-6) {
                x0 = unit(x0);
                V3 y0 = unit(cross(zTool, x0));
                Rc.setCol(0, x0);
                Rc.setCol(1, y0);
                Rc.setCol(2, zTool);
            }
        }
        IkSol sol;
        {
            Prof::Scope _p(g_prof.on, g_prof.nsIk);
            sol = ik(toBaseP(pC), toBaseR(Rc), p_.toolLen, seed);
        }
        if (!sol.ok) return false;

        FkOut f;
        fk(sol.theta, p_.toolLen, f);
        std::vector<Seg> segs;
        {
            Prof::Scope _p(g_prof.on, g_prof.nsBuild);
            buildSegments(f, sz_, p_.includeBase, segs);
        }

        double worstTool = 1e300, worstLink = 1e300;
        const Seg* segTool = nullptr;
        const Seg* segLink = nullptr;
        {
        Prof::Scope _p(g_prof.on, g_prof.nsSeg);
        for (const Seg& s : segs) {
            if (!checkLinks && s.group != SG_TOOL) continue;
            bool isTool = (s.group == SG_TOOL);
            double cur = isTool ? worstTool : worstLink;
            // 剪枝上限：若这一段连「当前最差表面距离」都够不着，它不可能更新最差值。
            // 不这么做就只能靠抽样（原来每 3 点才判一次连杆），那等于漏判 2/3 的点。
            double limit = std::min(p_.search, (cur < 1e299) ? cur + s.r : p_.search);
            double dist = cloud_.segmentDist(toCloudP(s.a), toCloudP(s.b), limit);
            double surface = dist - s.r;
            if (surface < cur) {
                if (isTool) { worstTool = surface; segTool = &s; }
                else { worstLink = surface; segLink = &s; }
            }
        }
        }
        const Seg* worstSeg = (worstTool <= worstLink) ? segTool : segLink;
        double worst = std::min(worstTool, worstLink);
        double selfGap = 1e300;
        {
        Prof::Scope _p(g_prof.on, g_prof.nsSelf);
        for (const auto& pr : pairs_) {
            const Seg& sa = segs[pr.first];
            const Seg& sb = segs[pr.second];
            double dax = segSegDist(sa.a, sa.b, sb.a, sb.b);
            selfGap = std::min(selfGap, dax - sa.r - sb.r);
        }
        }

        out.idx = -1;
        out.tiltDeg = tiltDeg;
        out.azDeg = azDeg;
        out.clearance = worst;
        out.clrTool = worstTool;
        out.clrLink = worstLink;
        out.part = worstSeg ? worstSeg->part : "";
        out.partKind = worstSeg ? worstSeg->kind : -1;
        out.partJ0 = worstSeg ? worstSeg->j0 : -1;
        out.partJ1 = worstSeg ? worstSeg->j1 : -1;
        out.linkKind = segLink ? segLink->kind : -1;
        out.linkJ0 = segLink ? segLink->j0 : -1;
        out.linkJ1 = segLink ? segLink->j1 : -1;
        out.ikOk = true;
        out.over = sol.over;
        out.margin = sol.margin;
        out.cond = condNumber(f);      // 只对最终入选的这个解算一次
        for (int i = 0; i < 6; ++i) out.theta[i] = sol.theta[i];
        out.pos = pC;
        out.gun = -d;
        out.sprayAngle = 0.0;
        out.standoffErr = std::fabs(norm(pC - face) - p_.standoff);
        out.selfGap = selfGap;
        return true;
    }

private:
    const CloudIndex& cloud_;
    Params p_;
    M3 baseR_, baseRT_;
    V3 baseT_;
    Sizing sz_;
    std::vector<std::pair<int, int>> pairs_;
};

// ---------------------------------------------------------------- 倾角方向

// 把法向 n 绕倾斜轴转 tiltDeg，方位可绕 n 偏 azDeg。tilt=0 严格返回 n（不倾斜的点一字不变）。
inline V3 dirOf(const V3& n, const V3& uUp, const V3& axis, double tiltDeg, double azDeg) {
    if (tiltDeg <= 1e-9 || norm(axis) < 1e-6) return unit(n);
    V3 u = uUp;
    if (std::fabs(azDeg) > 1e-9) {
        u = m3mul(rotAxis(n, azDeg * kDeg), uUp);
        u = unit(u - n * dot(u, n));
    }
    return unit(n * std::cos(tiltDeg * kDeg) + u * std::sin(tiltDeg * kDeg));
}

// 「把法向往上扳」的旋转轴 + 切平面内朝上的方向。
// ⚠ 退化判据必须用**归一化之前**的长度 |u| = sin(n 与 up 的夹角)。
//   先 unit 再比阈值是错的（unit 之后恒为 1，判据永不成立）。这个 bug 的表现很有欺骗性：
//   平整朝上的面被判成"非朝上"，于是往切平面里的噪声方向扳，倾角加到 40°、枪位横甩，
//   而喷距仍然守恒 —— 数字上一切"正常"。
//
// 倾角基准方向（"正对工件"的那个正方向），与「重力向上」分开。
//
// ⚠ 这是本文件唯一一处会把两者搞混的地方，务必分清：
//   · Params::upW      = 重力向上。抬枪沿它、自动猜基座靠它定"下方"，**不是**倾角基准。
//   · Params::tiltRefW = 倾角基准 = 「基座 → 工件」那条轴线，基座系里就是 ±X。
// 为什么要拆：槽型盘悬挂实测，枪轴 n=(1.00,0.01,-0.03) 本来正对工件（与 +X 差 1.97°），
// 但旧代码拿重力 +Z 当基准，量出 acos(n·Z)=91.88°，判定"这是竖直侧壁"，
// 于是按 tiltNonUpDeg=60 给出 61.3° 倾角，把枪扳到与 +X 差 59.4°（IJK 变成
// [0.509,0.005,-0.861]，整条修正轨迹变成斜的、整体下坠）。基准换成 +X 后
// 同一个点算出 60×1.97/90 = 1.3°，实测倾斜 0 点、倾角均值 0.0°，枪保持正对。
//
// 取值：基座位姿已知 → 基座 ±X（符号由主法向自动选，tiltSideX 可强制）；否则退回重力向上
// （自动假设的基座是猜的，拿它当基准等于把猜测当标定；此时重力向上是唯一可信的"上"）。

inline void tiltAxisAndDir(const V3& nIn, const V3& ref, V3& axis, V3& uUp) {
    V3 n = unit(nIn);
    V3 u = ref - n * dot(ref, n);
    // ⚠ 退化判据必须用归一化之前的长度 |u| = sin(n 与 ref 的夹角)。
    if (norm(u) < std::sin(5.0 * kDeg)) {
        axis = V3(0, 0, 0);
        uUp = V3(0, 0, 0);
        return;
    }
    uUp = unit(u);
    axis = unit(cross(n, uUp));
}

// 永不退化的切平面基。与 tiltAxisAndDir 的唯一区别：ref 与法向几乎平行时**不返回零向量**，
// 而是换一个世界轴兜底，保证「倾斜」在这个面上仍然有定义。
//
// 为什么必须有这个：tiltAxisAndDir 退化时返回零向量，调用方一律把 tilt 归 0、az 只有 0
// —— 于是朝上的大平面（法向 ∥ 重力）上"倾角"这个自由度被彻底关掉。实测 big（煤气灶背板
// 平放，枪轴与 −Z 夹角中位仅 3.9°）退化率 **100%**，tilt 唯一值只有 0°、
// 9761 个喷涂点里非零倾角的点 0 个；而它恰恰是 J5 贴限位最严重的一件（−5.23°）。
// 也就是说：唯一能救它的杠杆，被退化判据关掉了。
inline void tangentBasisRobust(const V3& nIn, const V3& ref, V3& axis, V3& uUp) {
    V3 n = unit(nIn);
    V3 u = ref - n * dot(ref, n);
    if (norm(u) < std::sin(5.0 * kDeg)) {
        // ref 与 n 近平行 → 取三个世界轴里「投影最长」的那个当参考，必定可用
        const V3 cands[3] = {V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1)};
        double bestL = -1.0;
        V3 bestU(0, 0, 0);
        for (int k = 0; k < 3; ++k) {
            V3 w = cands[k] - n * dot(cands[k], n);
            double L = norm(w);
            if (L > bestL) { bestL = L; bestU = w; }
        }
        if (bestL < 1e-9) {           // 理论上不可能（n 是单位向量）
            axis = V3(0, 0, 0);
            uUp = V3(0, 0, 0);
            return;
        }
        u = bestU;
    }
    uUp = unit(u);
    axis = unit(cross(n, uUp));
}

// ---------------------------------------------------------------- J5 的可解方向（解析）
//
// ⚠ 这是本文件唯一一处「用运动学结构反推姿态」的地方，别的都是几何/试解。务必读懂再改。
//
// 结论（推导与实测见 outputs/diag_关节限位_J5_2026-11-10.md）：
//   |θ5| = arccos(u_z)，u = R03ᵀ·g，R03 = Rz(θ1)·Ry(θ2+θ3)
//   R03 的第三列就是**前臂方向** a = (sin(θ2+θ3)cosθ1, sin(θ2+θ3)sinθ1, cos(θ2+θ3))
//   于是 u_z = g·a  ⟹  **|θ5| = ∠(a, g)**
// 也就是说：θ5 只取决于「前臂方向」与「枪轴方向」的夹角，**与绕枪轴滚转完全无关**
// （滚转只改 θ4/θ6 的分配）。这就是为什么"多解"救不了 J5：解支只换臂构型，
// θ5 由枪轴⊕前臂共同决定，枚举 28 个种子后最优支与当前支差中位 0.00°。
//
// 反过来用：**想改 θ5，唯一有效的杠杆是改枪轴方向 g**（= 改倾斜方向/倾角）。
// 固定 TCP 与喷距时，枪轴只能在「绕面点、半径 standoff」的球面上动：
//   g = −(cosT·n + sinT·U)
// 最大化 |a·g| 关于 U 的一阶条件给出
//   U* = −a_T/|a_T|，  a_T = a − n(a·n)   （前臂在切平面上的投影，**取反**）
// 代回得到闭式：
//   θ5(T) = arccos(−a_n·cosT + |a_T|·sinT)
// 注意 **不取绝对值**：u_z 可为负，θ5∈[0,180]，取绝对值会把 −0.577 的 125.2° 误算成 54.8°。
// 物理直觉：|θ5| 小 = 腕子伸直 = 工具轴与前臂同向；所以要让枪轴"顺着前臂伸出去"。
//
// 退化：|a_T| ≈ 0（前臂垂直于面）时倾斜对 θ5 完全无效 —— 这是结构死结，
// 只能换站位/换喷距，不是角度能救的。此时返回 norm(axis)=0，调用方按老路径走。

// 由关节角（度）取基座系前臂方向 a。θ4/θ5/θ6 不参与 —— 这正是"滚转无效"的体现。
inline V3 forearmDirOf(const double thetaDeg[6]) {
    const double t1 = thetaDeg[0] * kDeg;
    const double t23 = (thetaDeg[1] + thetaDeg[2]) * kDeg;
    const double s23 = std::sin(t23);
    return V3(s23 * std::cos(t1), s23 * std::sin(t1), std::cos(t23));
}

// 给定法向与目标关节余量，解析求「最省倾角的救 J5 倾斜方向」。
//   nC  : 点云系外法向（= −gun）
//   nB  : 同一法向在**基座系**下的表示
//   aB  : 基座系前臂方向（forearmDirOf 的结果）
//   marginTargetDeg : 期望把 |θ5| 压到多少度以内（限位 120°，通常给 120−余量）
//   maxTiltDeg      : 工艺/安全允许的最大倾角
// 输出：outTiltNeed（解除 θ5 危机所需最小倾角；-1 表示无解/退化）
// 返回 true 表示算出了有效倾角需求。**最优倾斜方向由调用方按 nC 现算**
// （切平面内 a_T 方向），因为把基座系方向转回点云系需要旋转矩阵，绕一圈只是引入误差。
//
// 方向约定：dirOf 里 uUp 是"绕 n 的起始方向"、az 是相对它的右手转角。调用方把
// uUp 直接设成 U*、az 设 0，就等价于「方位指向前臂反方向」，不必再走 stanceAzDeg 的启发式。
inline bool j5RescueTilt(const V3& nB, const V3& aB,
                         double marginTargetDeg, double maxTiltDeg,
                         double& outTiltNeed) {
    // 前臂在**基座系**切平面上的投影（推导全在基座系里成立）
    const double an = dot(aB, nB);
    const double at = norm(aB - nB * an);
    // 退化判据：切向分量 < cos(85°) ⟹ 前臂几乎垂直于面，倾斜改不动 θ5（结构死结）
    if (at < 0.0872) return false;

    const double R = std::sqrt(an * an + at * at);   // = |aB|
    if (R < 1e-9) return false;
    double C = std::cos(marginTargetDeg * kDeg);
    C = std::max(-1.0, std::min(1.0, C / R));
    // θ5(T) = arccos(R·cos(T − φ))，φ = atan2(at, −an)
    const double phiDeg = std::atan2(at, -an) * kRad;
    const double dPhiDeg = std::acos(C) * kRad;
    // 两个根里取最小的非负根 = 最省倾角的解
    double best = 1e300;
    for (int s = -1; s <= 1; s += 2) {
        double T = phiDeg + s * dPhiDeg;
        if (T > 180.0) T -= 360.0;
        if (T >= -1e-9 && T < best) best = T;
        if (T + 360.0 <= maxTiltDeg + 1e-9 && T + 360.0 < best) best = T + 360.0;
    }
    if (best > maxTiltDeg) {
        // 目标余量在倾角预算内够不着：退而求其次，取预算上限（能救多少救多少）
        best = maxTiltDeg;
    }
    outTiltNeed = best;
    return true;
}

// 该法向**需要多大倾角**（0 = 正对基准方向，α=90° 侧壁时给满 tiltNonUpDeg）。
// 抽成独立函数是因为行级方位打分也要用它 —— 两处若各写一份，
// 打分用的倾角和真正执行的倾角会悄悄分叉，A/B 结论就不可信了。
//
// α 与"悬垂面"都用 dot(n, ref)，不再写死 n.z。ref=(0,0,1) 时两者逐位等价。
// ref 传的是**倾角基准**（Params::tiltRefW），不是重力向上 —— 见上面 tiltRefW 的说明。
inline double tiltNeedDeg(const V3& nIn, const V3& ref, const Params& p) {
    V3 n = unit(nIn);
    double c = dot(n, ref);
    c = std::max(-1.0, std::min(1.0, c));
    double alphaDeg = std::acos(c) * kRad;
    double t = p.tiltNonUpDeg * (alphaDeg / 90.0);
    if (t > p.tiltMaxDeg) t = p.tiltMaxDeg;
    // 背离基准方向的面：扳过去会把枪顶进工件本体，这条路物理上不成立，不做线性外推
    if (c < -0.7) t = 0.0;
    return t;
}

// 「站位侧」方位角：机器人基座相对该加工点在哪一侧，枪就往那一侧倾。
//   u  = (base - face) 在切平面上的投影（指向基座）
//   方位角按 dirOf 的约定给：az 是 uUp 绕 n 右手旋转的角度。
// 为什么不用"向上"当默认：向上只是"重力上方通常空旷"的经验代理，
// 内腔里上方往往正对着对面的壁；而基座必然立在件外，从基座侧压进来
// 天然不需要把腕子伸进腔体里。两种都可能是错的，所以下面按**实测**择优。
// 退化（基座几乎在法向延长线上，如平放件的正上方）→ 返回 nan，表示无倾向。
inline double stanceAzDeg(const V3& nIn, const V3& face, const V3& base,
                          const V3& uUp, const V3& axis) {
    if (norm(axis) < 1e-6 || norm(uUp) < 1e-6) return std::numeric_limits<double>::quiet_NaN();
    V3 n = unit(nIn);
    V3 toBase = base - face;
    double L = norm(toBase);
    if (L < 1e-6) return std::numeric_limits<double>::quiet_NaN();
    V3 u = toBase - n * dot(toBase, n);
    // 切向分量 < sin(20°)·|toBase| ⇔ 基座基本在法向上，没有"哪一侧"可言
    if (norm(u) < 0.34 * L) return std::numeric_limits<double>::quiet_NaN();
    u = unit(u);
    double az = std::atan2(dot(u, cross(n, uUp)), dot(u, uUp)) * kRad;
    // 归一化到 (-180, 180]
    while (az <= -180.0) az += 360.0;
    while (az > 180.0) az -= 360.0;
    return az;
}

// ================================================================ 主规划

// ---------------------------------------------------------------- 逐点解析救 J5
//
// 本函数是整个 RS_POSE_TILT_J5 模式的决策核心：给定一个喷涂点，直接算出
// 「该用哪个倾斜方向、倾多大角」，使得 |θ5| 留够余量、且倾角尽量小。
//
// 三个输入决定一切：
//   · 外法向 n          → 枪轴只能绕它所在的切平面偏
//   · 前臂方向 a        → |θ5| = ∠(a, g)，a 由 θ1、θ2+θ3 给出（与腕部无关）
//   · 关节余量目标      → 把 |θ5| 压到 (120° − jointMarginMin) 以内
//
// 输出 (tilt, az)：az 是相对**稳健切平面基**（tangentBasisRobust 的 uUp）的方位角，
// 由 U* = −a_T 反算得到。因为 a 逐点连续，U* 与 tilt 也逐点连续 —— 姿态可以随点变，
// 但不会来回甩腕。这是它相对"面级常量"的核心优势，也是用户明确要求的口径。
//
// 注意：本函数**不解 IK**。它是纯几何+运动学的解析解，所以可以放心在阶段 1 逐点调用
// （10 万点也只是一次乘加）。能不能真解出来，仍由阶段 3 的 IK 说了算。
struct J5Plan {
    bool valid = false;      // false = 退化（前臂垂直于面）或已在目标内
    double tiltDeg = 0.0;
    double azDeg = 0.0;
    double j5PredDeg = 180.0;   // 预测的 |θ5|（便于诊断输出）
};

// ⚠ 三个方向参数的坐标系归属必须分清（这里踩过一次）：
//   nCloud = 点云系外法向（t.normal(i)）；aBase/nBase = 基座系。
//   把基座系的 n 当点云系用，会让 aT 的切向分解在错误平面上做，方位角整体偏掉。
inline J5Plan planTiltForJ5(const V3& nCloud, const V3& aBase, const V3& nBase,
                            const M3& baseR,   // 基座系→点云系的旋转
                            const Params& p) {
    J5Plan out;
    if (p.jointMarginMin <= 0.0 && p.tiltMaxDeg <= 0.0) return out;
    // J5 限位 ±120°（上下对称，取小的那侧）
    const double j5Lim = std::min(std::fabs(kJoints[4].hi), std::fabs(kJoints[4].lo));
    // 目标：留出 jointMarginMin + 1° 的缓冲。多给的这 1° 是留给「IK 收敛后 θ5
    // 与解析值有零点几度偏差」的，不然算出来刚好卡线、实解又掉出去。
    const double tgt = j5Lim - std::max(0.0, p.jointMarginMin) - 1.0;
    if (tgt <= 1.0) return out;

    const V3 aU = unit(aBase), nU = unit(nBase);
    const V3 nC = unit(nCloud);
    const double an = dot(aU, nU);
    const double at = std::sqrt(std::max(0.0, 1.0 - an * an));

    // ---- 0) 先看**当前会执行的姿态**下 θ5 是多少。不紧张就别动 —— 这是本次改动的关键。
    //
    // 判据用 tilt=0 的枪轴（g = −n），因为「非危机点就 tilt=0」是本模式对非危机点的
    // 唯一承诺（见阶段 2 的说明）。用别的东西当"现状"会和执行口径分叉。
    const double c0 = dot(aU, -nU);
    const double j5Now = std::acos(std::max(-1.0, std::min(1.0, c0))) * kRad;
    if (j5Now <= tgt) return out;                // 已在目标内 → 不倾斜（valid=false）

    // ---- 1) 解析求「把 |θ5| 压到 tgt 所需的最小倾角」（全在基座系里算）
    double tiltNeed = 0.0;
    if (!j5RescueTilt(nU, aU, tgt, p.tiltMaxDeg, tiltNeed)) return out;

    // ---- 2) 最优倾斜方向 U* = −a_T（切平面内、指向前臂反方向）
    // baseR 的列就是基座轴在点云系下的表示，所以点云系 a = baseR · aBase。
    const V3 aC = unit(m3mul(baseR, aU));
    V3 aT = aC - nC * dot(aC, nC);
    const double atLen = norm(aT);
    if (atLen < 0.0872) return out;          // 与 j5RescueTilt 同口径的退化判据
    aT = aT * (-1.0 / atLen);                // U*

    // ---- 3) 用稳健切平面基表达 U* 的方位角
    V3 axis, uUp;
    tangentBasisRobust(nC, p.upW, axis, uUp);
    if (norm(uUp) < 1e-6) return out;
    double az = std::atan2(dot(aT, axis), dot(aT, uUp)) * kRad;
    while (az <= -180.0) az += 360.0;
    while (az > 180.0) az -= 360.0;

    // ---- 4) 倾角：解析值优先，但受工艺喷角上限约束
    // 工艺上喷角 ≡ 倾角（d 是 n 绕切向轴转 tilt），不能为了救 J5 把喷角无限放大。
    // 解析需求超过上限时取上限：能救多少救多少，救不回来的由阶段 3 / 报告如实暴露。
    double tilt = tiltNeed;
    if (tilt > p.maxSprayAngleDeg) tilt = p.maxSprayAngleDeg;
    if (tilt > p.tiltMaxDeg) tilt = p.tiltMaxDeg;

    // ---- 5) 预测该姿态下的 |θ5|，并做「不许更差」的单调性保护
    // U* 与 a_T 反向，故 a·g = −a_n·cosT + |a_T|·sinT
    double c = -an * std::cos(tilt * kDeg) + at * std::sin(tilt * kDeg);
    c = std::max(-1.0, std::min(1.0, c));
    out.j5PredDeg = std::acos(c) * kRad;

    // ⚠ 保护：解析解只在**真的把 θ5 拉回来**时才采纳。
    // 目标值 tgt 是"上限"，但倾角受 maxSprayAngleDeg 钳制后可能落在
    // 「把 θ5 推得更远」的那一侧（闭式解有正负两支，取小正值不一定单调改善）。
    // 出现这种情况说明该点靠倾斜救不了（前臂方向不对），应当保持基线姿态、
    // 让报告如实暴露，而不是执行一个更糟的姿态。实测 big 就有点落在这一侧。
    if (out.j5PredDeg >= j5Now - 1e-6) return out;

    out.valid = true;
    out.tiltDeg = tilt;
    out.azDeg = az;
    return out;
}

struct VioCount {
    int byRule[RS_VIOL_COUNT] = {0};
    int total = 0;
    void add(unsigned int bit) {
        for (int i = 0; i < RS_VIOL_COUNT; ++i)
            if (bit & (1u << i)) ++byRule[i];
        ++total;
    }
};

// 第 i 点「至少需要多大倾角 / 要不要换方位」
//   azPref：本行的进刀侧方位（行级决策的结果；0 = 向上）。逐点不再自己决定方位，
//           因为方位角是"枪从哪一侧压过来"，逐点变等于沿轨迹反复拧腕，实机不可执行。
void needTilt(const Evaluator& ev, const Params& p, const TrajIn& t, int i,
              const double* seed, const V3* refX, bool checkLinks, double azPref,
              double& outTilt, double& outAz, PR& outPR, bool& hasPR) {
    V3 n = t.normal(i);
    V3 face = t.pos[i] + t.gun[i] * p.standoff;

    // 倾角需求按法向**线性**给，而不是"非朝上就一律 ±20°"的阶梯：
    //   δ = 法向与 +Z（"正位置"，即采集/相机拍过去的方向）的夹角
    //       （0° = 正朝上，90° = 竖直侧壁，180° = 完全倒挂）
    //   tilt0 = tiltNonUpDeg · (δ / 90°)，封顶 tiltMaxDeg
    // 系数 2026-10-26 由 30 提到 60（用户口径"60°×δ/90"）：侧壁（δ=90）直接给满 60°，
    // 即"把枪尽量扳到正位置再压下去"；δ<45 的面仍只给 <30° 的轻微补偿。
    // 取 90° 作满量程锚点是因为它物理上就是"面立起来了"的界限：此时给满基准倾角，
    // 比侧壁更倒挂的面（δ>90°）继续线性增长，直到上限。
    // 阶梯做法在 45° 边界上产生 20° 的阶跃，反映到轨迹上就是相邻点倾角突跳，
    // 实机过不了；线性化之后法向连续 → δ 连续 → 倾角需求连续。
    double tilt0 = tiltNeedDeg(n, p.tiltRefW, p);

    V3 axis, uUp;
    tiltAxisAndDir(n, p.tiltRefW, axis, uUp);
    bool degenerate = norm(uUp) < 1e-6;
    if (degenerate) {
        axis = V3(0, 0, 0);
        uUp = V3(0, 0, 0);
        tilt0 = 0.0;
    }

    std::vector<std::pair<double, double>> cand;
    cand.push_back({tilt0, azPref});       // 线性基准值 + 本行进刀侧，先试这一发（贴合需求，不无谓加大）
    if (!degenerate) {
        // 之后的搜索走**固定 5° 网格**，不能锚在 tilt0 上等距偏移。
        // tilt0 随法向连续变化，若候选 = tilt0 + k·step，候选值会随 α 一起漂移：
        // 同一个物理点可能从"刚好过关的 40°"滑到"刚好撞上的 36.48°"，纯数值抖动。
        // curved idx 2990 实测：网格版 40° 自碰撞间隙 +1.18mm，漂移版 36.48° 为 −1.16mm。
        double g = std::ceil(tilt0 / p.tiltStepDeg - 1e-9) * p.tiltStepDeg;
        while (g <= tilt0 + 1e-9) g += p.tiltStepDeg;   // 网格首档必须严格大于 tilt0
        for (double tt = g; tt <= p.tiltMaxDeg + 1e-9; tt += p.tiltStepDeg)
            cand.push_back({tt, azPref});
    }

    bool hasBest = false;
    PR bestAny;
    double bestClrKey = -1e300, bestMarginKey = -1e300;
    int flatRun = 0;
    bool hasPrev = false;
    double prevClr = 0.0, prevSelf = 0.0, prevMargin = 0.0;

    auto tryCand = [&](double tiltDeg, double az, PR& got) -> bool {
        V3 d = dirOf(n, uUp, axis, tiltDeg, az);
        PR pr;
        if (!ev.eval(face, d, refX, seed, tiltDeg, az, checkLinks, pr)) return false;
        pr.sprayAngle = angleBetweenDeg(d, n);
        if (pr.clearance >= p.safety && pr.selfGap >= 0.0 && !pr.over && pr.margin >= p.jointMarginMin) {
            got = pr;
            return true;
        }
        double clrKey = (pr.selfGap >= 0.0) ? pr.clearance : -1e300;
        if (clrKey > bestClrKey || (clrKey == bestClrKey && pr.margin > bestMarginKey)) {
            bestClrKey = clrKey;
            bestMarginKey = pr.margin;
            bestAny = pr;
            hasBest = true;
        }
        // 早停计数：只有**三个指标都不动**才算"倾角对这条方向的碰撞毫无影响"。
        // 只比 clearance 会误停：间隙常常早早饱和到 search 上限（120mm），
        // 而自碰撞/关节余量仍在随倾角变化，此时停掉就漏掉了后面能过关的档位。
        // curved idx 2990 实测：饱和区里 36.48° 自碰撞 −1.16mm、40° 反而 +1.18mm。
        if (hasPrev && std::fabs(pr.clearance - prevClr) < 0.5
                    && std::fabs(pr.selfGap - prevSelf) < 0.5
                    && std::fabs(pr.margin - prevMargin) < 0.5) ++flatRun;
        else flatRun = 0;
        prevClr = pr.clearance;
        prevSelf = pr.selfGap;
        prevMargin = pr.margin;
        hasPrev = true;
        return false;
    };

    for (const auto& c : cand) {
        PR got;
        if (tryCand(c.first, c.second, got)) {
            outTilt = c.first;
            outAz = c.second;
            outPR = got;
            hasPR = true;
            return;
        }
        // 连续 4 个候选间隙变化 < 0.5mm → 倾角对这个方向的碰撞毫无影响，
        // 再加也是白加（典型是碰撞来自基座/连杆布局，那是布局问题，倾角救不了）
        if (flatRun >= 4) break;
    }

    // 兜底：换倾斜方位（有些内壁往上扳会顶到对面的壁）
    // 方位清单围绕**本行进刀侧**展开；进刀侧本身不是"向上"时，把向上也留作最后一档
    // —— 站位侧可能正对着对面的壁，那时候"往上"反而是唯一空的方向。
    if (p.azFallback && !degenerate && flatRun < 4) {
        std::vector<double> azs;
        azs.push_back(azPref + 30.0);
        azs.push_back(azPref - 30.0);
        azs.push_back(azPref + 60.0);
        azs.push_back(azPref - 60.0);
        if (std::fabs(azPref) > 1e-9) azs.push_back(0.0);
        const double ttBase[4] = {10.0, 25.0, 40.0, 55.0};
        for (double az : azs) {
            // 归一到 (-180,180]，并夹掉因 ±180 撞车的重复候选（如 azPref=180 时 +30/-30 同一条线）
            while (az <= -180.0) az += 360.0;
            while (az > 180.0) az -= 360.0;
            for (int m = 0; m < 4; ++m) {
                double tt = std::max(std::ceil(tilt0 / p.tiltStepDeg - 1e-9) * p.tiltStepDeg, ttBase[m]);
                if (tt > p.tiltMaxDeg + 1e-9) continue;
                PR got;
                if (tryCand(tt, az, got)) {
                    outTilt = tt;
                    outAz = az;
                    outPR = got;
                    hasPR = true;
                    return;
                }
            }
        }
    }

    if (hasBest) {
        outTilt = bestAny.tiltDeg;
        outAz = bestAny.azDeg;
        outPR = bestAny;
        hasPR = true;
        return;
    }
    outTilt = tilt0;
    outAz = azPref;
    hasPR = false;
}

double percentile(std::vector<double> v, double q) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    double pos = (double)(v.size() - 1) * q;
    int lo = (int)std::floor(pos), hi = (int)std::ceil(pos);
    double f = pos - lo;
    return v[lo] * (1.0 - f) + v[hi] * f;
}

// 沿枪轴从轨迹点打到点云的距离（管道内点的投影中位数）。
// 不能用"最近点距离"：斜面/曲面上的最近点根本不在法向上，实测 complex 的最近点中位数
// 只有 72mm，会把正确的 standoff=100 误判成参数不一致。换成沿轴投影后全部样本落在 98.5~100.5。
bool rayStandoff(const CloudIndex& cloud, const V3& p, const V3& gun, double tube, double& out) {
    std::vector<double> tvec;
    double tube2 = tube * tube;
    for (int i = 0; i < cloud.size(); ++i) {
        V3 rel = cloud.point(i) - p;
        double tt = dot(rel, gun);
        if (tt <= 0.0 || tt >= 300.0) continue;
        double perp2 = dot(rel, rel) - tt * tt;
        if (perp2 < tube2) tvec.push_back(tt);
    }
    if (tvec.size() < 3) return false;
    std::sort(tvec.begin(), tvec.end());
    out = tvec[tvec.size() / 2];
    return true;
}

}  // namespace

// ================================================================ C ABI

const char* robot_safety_build_stamp(void) { return __DATE__ " " __TIME__; }

RobotSafetyParams robot_safety_default_params(void) {
    RobotSafetyParams p;
    std::memset(&p, 0, sizeof(p));
    p.standoff = 100.0;
    p.safety = 20.0;
    p.hardClearance = 2.0;   // 绝对安全距离：只在这个硬下限被破坏时才回退
    p.search = 120.0;
    p.toolLen = kDefaultToolLen;
    p.tiltNonUpDeg = 60.0;   // δ=90°（竖直侧壁）时的倾角。30→60 见 needTilt 的说明
    p.upThresholdDeg = 45.0;
    p.tiltMaxDeg = 60.0;
    p.tiltStepDeg = 5.0;
    p.tiltRateDeg = 3.0;
    p.maxSprayAngleDeg = 30.0;
    p.jointMarginMin = 1.0;
    p.condMax = 30.0;
    p.transitLiftMax = 400.0;
    p.voxelMm = 0.0;
    // 默认逐点全判连杆。原来的 3 点抽 1 点是为省时间加的，代价是 2/3 的点根本没查连杆，
    // 而"连杆会不会撞工件"恰恰是现场最关心的一条。eval 里已用剪枝上限把开销压回来。
    p.linkStride = 1;
    p.includeBase = 0;
    p.fixTransit = 1;
    p.azFallback = 1;
    // 默认要求可达（= 完整口径）。关掉它只影响"谁来对结果负责"，不影响几何判据。
    p.requireReach = 1;
    p.threads = 1;           // 默认单线程：结果可复现，便于对拍/回归
    p.azMode = 0;            // 默认「按站位逐行择优」；1 = 只用向上（旧口径）
    p.fastPose = 0;          // 默认完整搜索；1 = 单次逆解（见 fastPose 注释）
    return p;
}

const char* robot_safety_rule_name(unsigned int bit) {
    switch (bit) {
        case RS_VIOL_IK_FAIL: return "IK_FAIL";
        case RS_VIOL_COLLISION: return "COLLISION";
        case RS_VIOL_CLEARANCE: return "CLEARANCE";
        case RS_VIOL_SELF_COLLISION: return "SELF_COLLISION";
        case RS_VIOL_JOINT_LIMIT: return "JOINT_LIMIT";
        case RS_VIOL_SPRAY_ANGLE: return "SPRAY_ANGLE";
        case RS_VIOL_NEAR_SINGULAR: return "NEAR_SINGULAR";
        case RS_VIOL_TILT_JUMP: return "TILT_JUMP";
        case RS_VIOL_LINK_HIT: return "LINK_HIT";
        case RS_VIOL_TOOL_HIT: return "TOOL_HIT";
        case RS_VIOL_STEP_JUMP: return "STEP_JUMP";
        default: return nullptr;
    }
}

void robot_safety_free(RobotSafetyResult* res) {
    if (!res) return;
    delete[] res->points;
    res->points = nullptr;
    if (res->error) {
        delete[] res->error;
        res->error = nullptr;
    }
    res->count = 0;
}

// 真实实现在这里：比公开入口多一个 upCloud。公开的 robot_safety_check 转手调它并传 NULL,
// 这样 RobotSafetyParams 的布局一个字节都不用动（sizeof 仍 152），旧调用方二进制兼容。
static RobotSafetyResult robot_safety_check_impl(const SprayTrajPoint* pts, int n,
                                                 const double* cloudXyz, int cloudCount,
                                                 const double* basePose16,
                                                 const RobotSafetyParams* paramsIn,
                                                 const double* upCloud);

RobotSafetyResult robot_safety_check(const SprayTrajPoint* pts, int n,
                                     const double* cloudXyz, int cloudCount,
                                     const double* basePose16,
                                     const RobotSafetyParams* paramsIn) {
    return robot_safety_check_impl(pts, n, cloudXyz, cloudCount, basePose16, paramsIn, nullptr);
}

RobotSafetyResult robot_safety_check_up(const SprayTrajPoint* pts, int n,
                                        const double* cloudXyz, int cloudCount,
                                        const double* basePose16,
                                        const RobotSafetyParams* paramsIn,
                                        const double* upCloud) {
    return robot_safety_check_impl(pts, n, cloudXyz, cloudCount, basePose16, paramsIn, upCloud);
}

static RobotSafetyResult robot_safety_check_impl(const SprayTrajPoint* pts, int n,
                                                 const double* cloudXyz, int cloudCount,
                                                 const double* basePose16,
                                                 const RobotSafetyParams* paramsIn,
                                                 const double* upCloud) {
    RobotSafetyResult R;
    std::memset(&R, 0, sizeof(R));
    R.minClearanceMm = std::numeric_limits<double>::infinity();
    R.minClearanceToolMm = std::numeric_limits<double>::infinity();
    R.minClearanceLinkMm = std::numeric_limits<double>::infinity();
    R.count = 0;

    auto fail = [&](const char* msg) -> RobotSafetyResult {
        size_t len = std::strlen(msg) + 1;
        char* e = new char[len];
        std::memcpy(e, msg, len);
        R.error = e;
        return R;
    };

    if (!pts || n <= 0) return fail("轨迹为空");
    if (!cloudXyz || cloudCount <= 0) return fail("点云为空");

    Params p;
    {
        const RobotSafetyParams d = robot_safety_default_params();
        const RobotSafetyParams& s = paramsIn ? *paramsIn : d;
        p.standoff = s.standoff > 0 ? s.standoff : d.standoff;
        p.safety = s.safety > 0 ? s.safety : d.safety;
        // 0 是合法值（只要不穿透）；负数才回落到默认 2mm。
        p.hardClearance = s.hardClearance >= 0 ? s.hardClearance : d.hardClearance;
        p.search = s.search > 0 ? s.search : d.search;
        p.toolLen = s.toolLen > 0 ? s.toolLen : d.toolLen;
        p.tiltNonUpDeg = s.tiltNonUpDeg >= 0 ? s.tiltNonUpDeg : d.tiltNonUpDeg;
        p.upThresholdDeg = s.upThresholdDeg > 0 ? s.upThresholdDeg : d.upThresholdDeg;
        p.tiltMaxDeg = s.tiltMaxDeg > 0 ? s.tiltMaxDeg : d.tiltMaxDeg;
        p.tiltStepDeg = s.tiltStepDeg > 0 ? s.tiltStepDeg : d.tiltStepDeg;
        p.tiltRateDeg = s.tiltRateDeg >= 0 ? s.tiltRateDeg : d.tiltRateDeg;
        p.maxSprayAngleDeg = s.maxSprayAngleDeg > 0 ? s.maxSprayAngleDeg : d.maxSprayAngleDeg;
        p.jointMarginMin = s.jointMarginMin;
        p.condMax = s.condMax > 0 ? s.condMax : d.condMax;
        p.transitLiftMax = s.transitLiftMax > 0 ? s.transitLiftMax : d.transitLiftMax;
        p.voxelMm = s.voxelMm > 0 ? s.voxelMm : d.voxelMm;
        p.linkStride = s.linkStride > 0 ? s.linkStride : d.linkStride;
        p.includeBase = s.includeBase != 0;
        p.fixTransit = s.fixTransit != 0;
        p.azFallback = s.azFallback != 0;
        p.threads = s.threads > 0 ? s.threads : d.threads;
        // 语义：0 才是「不要求可达」；非 0（含 default_params 给的 1）= 完整口径。
        p.requireReach = (s.requireReach != 0);
        // 语义：0（含零初始化结构）= 默认＝按站位择优；只有显式填 1 才退回旧口径
        p.azStance = (s.azMode != 1);
        p.fastPose = (s.fastPose & RS_POSE_FAST) != 0;
        p.tiltPerPoint = (s.fastPose & RS_POSE_TILT_PERP) != 0;
        p.tiltJ5 = (s.fastPose & RS_POSE_TILT_J5) != 0;
        p.tiltSideX = s.tiltSideX;
    }

    // ---- 重力向上方向（点云系）
    // NULL / 零向量 / 长度过小 ⇒ 保持旧的 (0,0,1)，旧调用方结果逐字节不变。
    // 给了就归一化：允许调用方传未归一化的方向（例如从矩阵行直接取），不要求它先算模长。
    if (upCloud) {
        V3 u(upCloud[0], upCloud[1], upCloud[2]);
        double L = norm(u);
        // 阈值取 1e-6 而不是 0：界面上"没填"常是 0,0,0，而 1e-12 这种残值更可能是
        // 单位换算留下的垃圾。两者都该退回默认，不该被当成"朝某个方向"。
        if (L > 1e-6) p.upW = u * (1.0 / L);
    }

    // ---- 点云索引
    g_prof.on = std::getenv("RS_PROFILE") != nullptr;
    g_prof.reset();
    g_prof.t0 = g_prof.last = std::chrono::steady_clock::now();
    CloudIndex cloud;
    {
        // 可选体素降采样：只影响"点云有多细"，不影响判据口径
        std::vector<double> dec;
        const double* src = cloudXyz;
        int srcN = cloudCount;
        if (p.voxelMm > 0.5) {
            voxelDecimate(cloudXyz, cloudCount, p.voxelMm, dec);
            if ((int)(dec.size() / 3) >= 64) {
                src = dec.data();
                srcN = (int)(dec.size() / 3);
            }
        }
        cloud.build(src, srcN);
    }
    double spacing = cloud.estimateSpacing();
    g_prof.mark("点云KD构建");

    // ---- 轨迹
    TrajIn t;
    t.pos.resize(n);
    t.gun.resize(n);
    t.spray.resize(n);
    t.speed.resize(n);
    for (int i = 0; i < n; ++i) {
        t.pos[i] = V3(pts[i].x, pts[i].y, pts[i].z);
        V3 g(pts[i].i, pts[i].j, pts[i].k);
        double gl = norm(g);
        if (gl < 1e-9) return fail("轨迹里有零长度的 IJK，无法定姿态");
        t.gun[i] = g * (1.0 / gl);
        t.spray[i] = pts[i].spray ? 1 : 0;
        t.speed[i] = pts[i].speed;
    }
    std::vector<std::pair<int, int>> rows;
    computeRows(t, rows);

    // ---- 基座位姿：T_base_in_cloud
    // 喷涂点法向的均值（归一化）。两处要用：自动猜基座时定"外侧"，定倾角基准时选 ±X。
    V3 nMean = V3(0, 0, 1);
    {
        V3 nm(0, 0, 0);
        int cnt = 0;
        for (int i = 0; i < n; ++i)
            if (t.spray[i]) {
                nm = nm + t.normal(i);
                ++cnt;
            }
        if (cnt && norm(nm) > 1e-9) nMean = unit(nm);
    }

    Pose base;
    bool baseAuto = false;
    if (basePose16) {
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) base.R.m[r][c] = basePose16[r * 4 + c];
            base.t[r] = basePose16[r * 4 + 3];
        }
    } else {
        baseAuto = true;
        // 猜法：基座沿工件外侧退 0.50·reach、下沉 0.35·reach，让工件质心落在基座系
        // (0.5·reach, 0, 0.35·reach) 的舒适区。⚠ 这是假设，不是标定值。
        //
        // 「外侧」与「下方」都相对**重力**定义：基座立在地面上，工件在它前方上方。
        // 所以水平投影必须去掉 up 方向的分量，下沉沿 -up —— up=(0,0,1) 时与旧版逐位一致。
        V3 h = nMean - p.upW * dot(nMean, p.upW);   // 水平投影（去掉重力分量）
        if (norm(h) < 0.2) h = V3(1, 0, 0) - p.upW * dot(V3(1, 0, 0), p.upW);
        if (norm(h) < 1e-6) h = V3(0, 1, 0) - p.upW * dot(V3(0, 1, 0), p.upW);
        h = unit(h);
        V3 origin = cloud.centroid - h * (0.50 * kReach) - p.upW * (0.35 * kReach);
        V3 x = unit(cloud.centroid - origin);
        V3 z = p.upW;                                // 基座 +Z = 重力向上
        V3 y = unit(cross(z, x));
        x = unit(cross(y, z));
        base.R.setCol(0, x);
        base.R.setCol(1, y);
        base.R.setCol(2, z);
        base.t = origin;
    }

    // ---- 竖直向上 & 倾角方向：安全层是在**基座系**里算的
    // 基座位姿已知（手眼 / 恒等 / 手动）时，基座 +Z 就是竖直向上，倾角大小（δ）必须以它
    // 为准。此前这里沿用界面在**点云系**里填的方向，而那个方向是相机姿态的函数：
    // 相机侧装/仰视时它和真实竖直差 90° 以上，表现为"平放件被当成倒挂天花板"——
    // δ 判成 174°，tiltNeedDeg 命中 c<−0.7 直接给 0 倾角，整条倾角口径反了。
    if (basePose16) {
        p.upW = unit(base.R.col(2));
        if (std::getenv("RS_DBG_UP"))
            std::fprintf(stderr, "[safety] 竖直向上取基座 +Z = (%.3f,%.3f,%.3f)（点云系）\n",
                         p.upW.x, p.upW.y, p.upW.z);
    }
    // ---- 倾角基准：基座 ±X（"基座 → 工件"那条轴线），**不是**重力 +Z。
    // 上文 p.upW 仍是重力向上（抬枪、判"下方"要用），两者从此分家。
    //
    // 为什么默认就必须是 ±X：安全层是在基座系里算的，"枪正对工件"这个正方向在基座系里
    // 就是 X 轴（工件在基座 +X 方向、挂在前面）。拿重力 +Z 当基准等于问"枪是否朝上"，
    // 而工件是立着挂的、枪是横着打的，两者天然差 90° —— 每个点都会被判成"竖直侧壁"，
    // 一律给满 tiltNonUpDeg，整条轨迹被扳歪（实测枪轴从 [0.999,0.010,-0.033] 扳到
    // [0.509,0.005,-0.861]，与 +X 差 59.4°，轨迹整体下坠）。
    //
    // 符号（+X 还是 −X）：取与主法向更"同向"的那一侧，使"正对"落在 0° 而不是 180°。
    //   工件挂在基座 +X 侧时主法向≈+X，选 +X → α≈0 → 倾角 0，正对，符合预期；
    //   若某件主面朝 −X（背面朝基座），选 −X → 同样 α≈0，不会被误判成"背离"（c<−0.7 归零）。
    // tiltSideX 非 0 时按用户强制，不自动判。
    p.tiltRefW = p.upW;
    if (basePose16) {
        V3 ax = unit(base.R.col(0));
        double side = 1.0;
        if (std::fabs(p.tiltSideX) > 1e-9) {
            side = (p.tiltSideX > 0 ? 1.0 : -1.0);
        } else if (dot(nMean, ax) < 0.0) {
            // 主法向朝 −X：基准跟着翻，否则 α 会算成 180°−ε 而命中 c<−0.7 归零
            side = -1.0;
        }
        p.tiltRefW = ax * side;
        if (std::getenv("RS_DBG_UP"))
            std::fprintf(stderr, "[safety] 倾角基准 = 基座 %cX = (%.3f,%.3f,%.3f)（点云系）"
                                 " %s；重力向上 = (%.3f,%.3f,%.3f)\n",
                         side > 0 ? '+' : '-', p.tiltRefW.x, p.tiltRefW.y, p.tiltRefW.z,
                         std::fabs(p.tiltSideX) > 1e-9 ? "用户强制" : "按主法向自动选侧",
                         p.upW.x, p.upW.y, p.upW.z);
    } else {
        // 没有基座位姿时退回重力向上：自动假设的基座是猜的，不能拿它的 X 轴当标定。
        p.tiltRefW = p.upW;
        if (std::getenv("RS_DBG_UP"))
            std::fprintf(stderr, "[safety] 无基座位姿，倾角基准退回重力向上 = (%.3f,%.3f,%.3f)\n",
                         p.tiltRefW.x, p.tiltRefW.y, p.tiltRefW.z);
    }

    Evaluator ev(cloud, base, p);

    int nThreads = std::max(1, std::min(p.threads, 64));

    // 点 → 行 的归属在建行之后就定死，阶段 0/1/3 都按它取行级方位
    std::vector<int> rowOf(n, -1);
    for (size_t ri = 0; ri < rows.size(); ++ri)
        for (int i = rows[ri].first; i <= std::min(rows[ri].second, n - 1); ++i) rowOf[i] = (int)ri;

    // ---- 姿态（倾角 + 进刀方位）的决策在「阶段 2」做，且作用域是「面组」而不是点：
    // 同一个面上姿态恒定（工艺要求），只有换面时才允许改变。见阶段 2 的说明。
    int azStanceRows = 0;

    // ---- 阶段 1：每点的倾角需求（带 IK 种子连续性）
    std::vector<double> need(n, 0.0), azOf(n, 0.0);
    // 逐点救 J5 用的「基线 |θ5|」记录（诊断用）。
    // ⚠ 这里**只留诊断数组**：J5 的倾角/方位不再由阶段 1 预算并存下来，
    //   而是由阶段 3.6 在拿到真实执行位姿之后现场算（见阶段 3.6 的说明）。
    //   历史教训：第一版在这里预算 (tiltJ5, azJ5) 并在阶段 2 覆盖基线，
    //   等于让 J5 取代基线 —— complex 关节步长 47→93、longstripe P0 2→22。
    std::vector<float> j5Base(n, 0.0f);
    std::vector<char> hasCache(n, 0);
    std::vector<PR> cache(n);
    std::vector<char> hasPrevSol(n, 0);

    // 逐点独立求解（只借用"上一点的解"当 IK 种子加速收敛），所以可以按块并行：
    // 每块内部仍串行传递种子，块首用解析种子自启。块大小取 n/(线程数·4)，
    // 让尾块不至于拖尾（长件上点耗时差异很大：臂展中心的点快、边缘的点慢）。
    int chunk = std::max(1, n / (nThreads * 4));

    auto runRange = [&](int b0, int e0) {
        double prevTheta[6] = {0, 0, 0, 0, 0, 0};
        bool hasSeed = false;
        V3 prevX;
        bool hasRefX = false;
        for (int i = b0; i < e0; ++i) {
            bool checkLinks = (i % std::max(1, p.linkStride) == 0);
            double tilt = 0.0, az = 0.0;
            PR pr;
            bool has = false;
            // ---- 单次逆解模式 / 逐点连续倾角：倾角取**纯几何需求**（法向线性），
            // 不解 IK、不试方位、不留缓存。代价是"这个倾角到底解不解得出来"要等阶段 3 才知道；
            // 收益是这里一次 IK 都不花。这一句就是 fast 模式省掉 60% 时间的全部原因
            // （阶段 1 原本占 58.4s/83s）。tiltPerPoint 走同一条：它的倾角同样与 IK 无关。
            if (p.fastPose || p.tiltPerPoint) {
                need[i] = tiltNeedDeg(t.normal(i), p.tiltRefW, p);
                azOf[i] = 0.0;
                continue;
            }
            // 阶段 1 一律按「向上」找倾角：方位要等倾角定下来才谈得上（见阶段 2 末尾）。
            // 这也保证了"旧口径安全的样本"在阶段 1 得到与旧版完全相同的倾角与缓存。
            needTilt(ev, p, t, i, hasSeed ? prevTheta : nullptr, hasRefX ? &prevX : nullptr,
                     checkLinks, 0.0, tilt, az, pr, has);
            need[i] = tilt;
            azOf[i] = az;
            if (has) {
                pr.idx = i;
                cache[i] = pr;
                hasCache[i] = 1;
                for (int k = 0; k < 6; ++k) prevTheta[k] = pr.theta[k];
                hasSeed = true;
                FkOut f;
                fk(pr.theta, p.toolLen, f);
                prevX = f.tcp.R.col(0);
                hasRefX = true;
            }
        }
    };

    if (nThreads <= 1) {
        runRange(0, n);
    } else {
        std::atomic<int> next{0};
        std::vector<std::thread> pool;
        pool.reserve(nThreads);
        for (int k = 0; k < nThreads; ++k) {
            pool.emplace_back([&] {
                for (;;) {
                    int b = next.fetch_add(chunk);
                    if (b >= n) break;
                    runRange(b, std::min(n, b + chunk));
                }
            });
        }
        for (auto& th : pool) th.join();
    }

    // ---- 阶段 2：面级姿态常量（倾角 + 进刀方位）
    //
    // 工艺要求（用户口径）：**同一个面上喷涂的过程中，喷角即姿态尽量不要改变**。
    // 旧口径是「逐点求最小可行倾角 + 只升不降的限幅平滑」，倾角是逐点自由量：
    // 只要它逐点自由，就一定会随局部位姿抖动。开口盒子实测一笔 77 点里倾角从 30°
    // 走到 55° 再落回 30°，相邻点 TCP 沿球面来回摆（位移可达 40mm），画出来是锯齿，
    // 落盘执行是一路拧腕。
    //
    // 现在把「倾角 + 方位」一起升级成**面级常量**：相邻、法向一致、共面的行归成一个
    // 面组，每个面组只允许一个 (倾角, 方位)，姿态因此只可能在换面时改变。
    // 面内个别点在这个常量下过不了硬约束时，仍由阶段 3 退回本点最优并计入 reverted
    // —— 那是物理上必须偏的点（否则就撞了），不是决策抖动。
    //
    // 候选倾角直接取该面内阶段 1 求出的 need 值（本来就在 tiltStepDeg 网格上）+ tilt0，
    // 不另造一套网格：两套网格会让「打分用的倾角」和「执行用的倾角」悄悄分叉。
    // 若这批候选全都有 P0（典型是"整面只有某个不直观的方位才够得着"），再做一轮
    // 广搜（倾角 × 每 30° 方位）—— 面级常量只有一份，选错了整面都废，值得多花这点时间。
    const double kSameFaceCos = std::cos(2.0 * kDeg);   // 行法向夹角 ≤2° 视为同一个面
    const double kSameFacePlaneMm = 2.0;                // 共面判据：沿法向偏移 ≤2mm

    std::vector<int> groupOfRow(rows.size(), -1);
    std::vector<std::pair<int, int>> groups;            // 每个面组的行区间 [r0, r1]
    {
        std::vector<V3> rowNrm(rows.size(), V3(0, 0, 0));
        for (size_t r = 0; r < rows.size(); ++r) {
            V3 s(0, 0, 0);
            int c = 0;
            for (int i = rows[r].first; i <= rows[r].second; ++i)
                if (t.spray[i]) { s = s + t.normal(i); ++c; }
            if (c) rowNrm[r] = unit(s);
        }
        int g = -1;
        for (size_t r = 0; r < rows.size(); ++r) {
            bool same = false;
            if (r > 0) {
                int m0 = (rows[r - 1].first + rows[r - 1].second) / 2;
                int m1 = (rows[r].first + rows[r].second) / 2;
                V3 f0 = t.pos[m0] + t.gun[m0] * p.standoff;
                V3 f1 = t.pos[m1] + t.gun[m1] * p.standoff;
                same = dot(rowNrm[r], rowNrm[r - 1]) >= kSameFaceCos &&
                       std::fabs(dot(rowNrm[r], f1 - f0)) <= kSameFacePlaneMm;
            }
            if (!same) { ++g; groups.push_back({(int)r, (int)r}); }
            else groups.back().second = (int)r;
            groupOfRow[r] = g;
        }
    }

    std::vector<double> finalTilt = need;
    // 空移点（spray=0）：倾角恒为 0，姿态就是轨迹自带的 IJK，一列都不许改。
    // 理由：空移时枪是关的，"为了更好解而倾转"没有工艺意义，反而会把算法自己设计好的
    // 避障路径（抬到安全高度、沿法向退出）扭成另一条 —— 实测开口盒子 1958 个空移点里
    // 1669 个被倾转（最大 30°），IJK 变化到 0.87，等于换了一条空移路径。
    // 置 0 后 dirOf(nrm,uUp,axis,0,az) == nrm == -gun，逐位回到原轨迹姿态。
    for (int i = 0; i < n; ++i)
        if (!t.spray[i]) finalTilt[i] = 0.0;
    std::vector<double> azRow(rows.size(), 0.0);
    std::vector<double> faceTilt(groups.size(), 0.0), faceAz(groups.size(), 0.0);
    if (p.tiltPerPoint) {
        // ---- 逐点连续口径（RS_POSE_TILT_PERP）
        //
        // 倾角只由**该点法向**决定：δ = 法向与「正位置」(+Z，即采集/相机拍过去的方向)的夹角，
        // tilt = tiltNonUpDeg × δ/90，封顶 tiltMaxDeg。法向连续 ⟹ 倾角连续，
        // 一圈桶壁就是一条平滑的倾角曲线，不会再出现 0°/25° 两块台阶。
        //
        // 为什么不走面级择优：面级常量的唯一好处是"面内姿态绝对不动"，代价是**换面即跳变**。
        // 垃圾桶实测（面级口径）：tilt 只有 0°(306 点)/25°(95)/19.7°(69)/20.5°(50) 几个离散值，
        // 相邻点枪轴夹角最大 48.7°、>5° 的有 17 段；槽型盘 348 点全 0°、其余侧 23~35°。
        // 用户口径明确：「法向偏离正位置越大就倾得越大，线性给差值」—— 要的是连续，
        // 不是面内恒定。安全性交给阶段 3 逐点兜底（±1 步倾角救援 + 过不了照实报 P0），
        // 不再靠"整面一起歪"去救个别点。
        for (int i = 0; i < n; ++i)
            if (t.spray[i]) finalTilt[i] = tiltNeedDeg(t.normal(i), p.tiltRefW, p);
        azRow.assign(rows.size(), 0.0);   // 方位固定「向上」，与倾角同一条连续规则
    } else {
        // 面组 → 喷涂点索引（抽样打分用；过渡点位置可改，不参与姿态决策）
        std::vector<std::vector<int>> gp(groups.size());
        for (size_t g = 0; g < groups.size(); ++g)
            for (int r = groups[g].first; r <= groups[g].second; ++r)
                for (int i = rows[r].first; i <= rows[r].second; ++i)
                    if (t.spray[i]) gp[g].push_back(i);

        struct FaceScore {
            int hit = 0;        // 枪体/连杆 穿透点云（P0）
            int unreach = 0;    // IK 解不出（P0）
            int selfHit = 0;    // 自碰撞（P0）
            // 关节超出限位。**只统计、不参与择优**，原因写在下面（这行是刻意的，别"修"）。
            //
            // 看起来它该进 P0：hardOk 把 over 当硬约束。但实测把 over 计入 P0 之后，
            // 开口盒子 face3 会从 (60°, 0°) 换成 (60°, −60°) —— 关节超限点数确实从
            // 1111 降到 835，代价是整面 2080 点的枪轴方位歪 60°：与相邻面组差 60°、
            // 交界处 TCP 一次甩 100mm（选 −120° 时是 152mm），而且喷涂膜厚会沿面倾斜。
            // 也就是说，"减少逐点回退"的唯一手段是让整面斜着进枪，这在工艺上是硬伤。
            // 关节超限是个点级的腕部问题，换 5~10° 倾角就能救（位移只有 9mm），
            // 交给阶段 3 逐点处理，比让整面歪掉便宜得多。
            int over = 0;
            int shortClr = 0;   // 间隙 < safety 的点数（工艺软目标）
            double minClr = 1e300;
            int p0() const { return hit + unreach + selfHit; }
        };

        auto scoreFace = [&](size_t g, double tilt, double az, int S) -> FaceScore {
            FaceScore s;
            const std::vector<int>& P = gp[g];
            int len = (int)P.size();
            if (len <= 0) return s;
            if (S > len) S = len;
            // 热种子链：阶段 3 是「接着上一点继续解」的，冷启动（全零种子）会把本来就够得着
            // 的点判成不可达。开口盒子 +X 内壁整面 2104 点就是这么被判死的：所有候选都是
            // 冷启动 400/400 不可达，选优只能退回「倾角最小」，整面 30° 全不可达。
            double seed[6] = {0, 0, 0, 0, 0, 0};
            bool hasSeed = false;
            for (int k = 0; k < S; ++k) {
                int i = P[(S <= 1) ? 0 : (size_t)k * (len - 1) / (S - 1)];
                V3 n = t.normal(i);
                V3 face = t.pos[i] + t.gun[i] * p.standoff;
                V3 axis, uUp;
                tiltAxisAndDir(n, p.tiltRefW, axis, uUp);
                if (norm(uUp) < 1e-6) continue;
                V3 d = dirOf(n, uUp, axis, tilt, az);
                PR pr;
                // 抽样本来就少，连杆一起查：抽错一次就是整面选错，省这点时间不划算
                if (!ev.eval(face, d, nullptr, hasSeed ? seed : nullptr, tilt, az, true, pr)) {
                    ++s.unreach;
                    hasSeed = false;      // 断链后重新冷启，别拿远处的解硬顶
                    continue;
                }
                for (int q = 0; q < 6; ++q) seed[q] = pr.theta[q];
                hasSeed = true;
                if (pr.clrTool < 0.0 || pr.clrLink < 0.0) ++s.hit;
                if (pr.selfGap < 0.0) ++s.selfHit;
                if (pr.over) ++s.over;
                if (std::isfinite(pr.clearance)) {
                    if (pr.clearance < p.safety) ++s.shortClr;
                    if (pr.clearance < s.minClr) s.minClr = pr.clearance;
                }
            }
            return s;
        };

        // 选优次序：致命点数少 > 不破工艺喷角 > 间隙不足点数少 > 倾角小 > 方位取 0（向上）
        //
        // 第一档必须比"撞/不可达/自碰撞"更粗：只按撞数排，会把"整面都不可达但一点不撞"
        // 排在"大部分够得着、少数点蹭边"前面（开口盒子 +X 面就这么选错过）。
        //
        // 第二档（工艺喷角）单独拿出来是必须的：倾角就是喷角（= 枪轴与法向的夹角），
        // 而间隙不足只是 P1 告警。若不设这一档，"整面把间隙从 8mm 拉到 20mm"会赢过
        // "整面守住 30° 喷角"，代价是**整个面**常年超工艺喷角（实测 longstripe 从
        // 30~35° 被抬到 59°，P0 反而多了 4 个）。面级常量只有一份，这种全局违规
        // 不能为了软告警而接受：P0 该救就救，救完就退回工艺基准。
        auto betterThan = [&](const FaceScore& a, double ta, double aza,
                              const FaceScore& b, double tb, double azb) {
            // 不要求可达时，"不可达/自碰撞"这两个计数恒为 0，P0 退化成纯粹的
            // 「枪体撞了几点」。这正是"只要不撞、姿态越规整越好"的那把尺子。
            int pa = p.requireReach ? a.p0() : a.hit;
            int pb = p.requireReach ? b.p0() : b.hit;
            if (pa != pb) return pa < pb;
            bool oka = ta <= p.maxSprayAngleDeg + 1e-9;
            bool okb = tb <= p.maxSprayAngleDeg + 1e-9;
            if (oka != okb) return oka;
            if (a.hit != b.hit) return a.hit < b.hit;
            if (p.requireReach) {
                if (a.unreach != b.unreach) return a.unreach < b.unreach;
                if (a.selfHit != b.selfHit) return a.selfHit < b.selfHit;
            }
            if (a.shortClr != b.shortClr) return a.shortClr < b.shortClr;
            if (std::fabs(ta - tb) > 1e-6) return ta < tb;
            return std::fabs(aza) < std::fabs(azb) - 1e-9;
        };

        struct Cand { double tilt, az; FaceScore s; };
        auto cmpCand = [&](const Cand& a, const Cand& b) {
            return betterThan(a.s, a.tilt, a.az, b.s, b.tilt, b.az);
        };

        auto decideGroup = [&](size_t g) {
            const std::vector<int>& P = gp[g];
            if (P.empty()) return;
            int ra = groups[g].first, rb = groups[g].second;

            // ---- 单次逆解模式：面级姿态用「轻量择优」定，不做抽样打分池、不做广搜、不做密集复核。
            //
            // 为什么不能纯解析（倾角 = 法向需求）：实测开口盒子，纯解析时整件最大倾角 30°，
            // 但 +X 内壁在 30° 下 2104 点一个都解不出（这面需要 60° 才能把腕部翻进去），
            // 可达率从 99.3% 掉到 89.0%。倾角需求是**几何**需求，能不能解是**运动学**问题，
            // 两者在深腔件上差 30°。所以这里保留"试几个倾角、看哪个解得出"，只是把
            // 采样从 S=32/96/400 三档（约 3.5 万次 eval）压到每组 ≤7 倾角 × ≤8 点（约 500 次）。
            //
            // 打分是词典型的：**先看 ikOk 数**（解不出 = 这点根本执行不了，权重最高），
            // 再看 hardOk 数（够得着但撞了 = 站位/尺寸问题，不是倾角能救的），
            // 平手取小倾角（工艺喷角小、膜厚更均匀）。方位仍只用「向上」——
            // 翻方位一次 TCP 要绕面点转 2·standoff·sin(Δaz/2)（60° 时 100mm），fast 模式不开这个口子。
            if (p.fastPose) {
                double need0 = 0.0;
                for (int i : P) need0 = std::max(need0, need[i]);
                std::vector<double> cand;
                cand.push_back(need0);
                double g0 = std::ceil(need0 / p.tiltStepDeg - 1e-9) * p.tiltStepDeg;
                if (g0 <= need0 + 1e-9) g0 += p.tiltStepDeg;   // 网格首档严格大于需求（同 needTilt 口径）
                for (double tt = g0; tt <= p.tiltMaxDeg + 1e-9; tt += p.tiltStepDeg)
                    cand.push_back(tt);

                int K = std::min<int>(24, (int)P.size());
                int bestC = 0, bestIk = -1, bestHard = -1;
                double bestClr = -1e300;
                for (int c = 0; c < (int)cand.size(); ++c) {
                    int ikN = 0, hardN = 0;
                    double minClr = 1e300;
                    double seed[6] = {0, 0, 0, 0, 0, 0};
                    bool hasSeed = false;
                    V3 prevX;
                    bool hasRefX = false;
                    for (int k = 0; k < K; ++k) {
                        int i = P[K > 1 ? (int)((long long)(P.size() - 1) * k / (K - 1)) : 0];
                        V3 nrm = t.normal(i);
                        V3 face = t.pos[i] + t.gun[i] * p.standoff;
                        V3 axis, uUp;
                        tiltAxisAndDir(nrm, p.tiltRefW, axis, uUp);
                        if (norm(uUp) < 1e-6) { axis = V3(0, 0, 0); uUp = V3(0, 0, 0); }
                        V3 d = dirOf(nrm, uUp, axis, cand[c], 0.0);
                        PR pr;
                        if (!ev.eval(face, d, hasRefX ? &prevX : nullptr,
                                     hasSeed ? seed : nullptr, cand[c], 0.0, true, pr)) {
                            hasSeed = false;
                            continue;
                        }
                        ++ikN;
                        if (hardOk(pr, p)) ++hardN;
                        minClr = std::min(minClr, pr.clearance);
                        for (int q = 0; q < 6; ++q) seed[q] = pr.theta[q];
                        hasSeed = true;
                        FkOut f;
                        fk(pr.theta, p.toolLen, f);
                        prevX = f.tcp.R.col(0);
                        hasRefX = true;
                    }
                    // 词序比较：解得出 > 过得了硬约束 > 间隙大 > 倾角小（倾角小由 bestC 先到先得保证）
                    if (ikN > bestIk ||
                        (ikN == bestIk && (hardN > bestHard ||
                                           (hardN == bestHard && minClr > bestClr + 1e-9)))) {
                        bestIk = ikN;
                        bestHard = hardN;
                        bestClr = minClr;
                        bestC = c;
                    }
                }
                double faceT = cand[bestC];
                faceTilt[g] = faceT;
                faceAz[g] = 0.0;
                for (int r = ra; r <= rb; ++r) {
                    azRow[r] = 0.0;
                    for (int i = rows[r].first; i <= rows[r].second; ++i)
                        if (t.spray[i]) finalTilt[i] = faceT;
                }
                if (g_prof.on) {
                    std::fprintf(stderr, "  [face%2zu] fast 行 %3d..%3d 点%5d 倾角%5.1f 候选%zu 采样%d ikOk%d hardOk%d 最小%7.1f\n",
                                 g, ra, rb, (int)P.size(), faceT, cand.size(), K, bestIk, bestHard, bestClr);
                }
                return;
            }

            V3 nsum(0, 0, 0);
            for (int i : P) nsum = nsum + t.normal(i);
            V3 nFace = unit(nsum);
            V3 axis, uUp;
            tiltAxisAndDir(nFace, p.tiltRefW, axis, uUp);
            bool degenerate = norm(uUp) < 1e-6;       // 朝上的面：倾角无意义，恒 0
            double tilt0 = degenerate ? 0.0 : tiltNeedDeg(nFace, p.tiltRefW, p);

            // 候选倾角：tilt0 + 该面内阶段 1 实际需要过的倾角（都在网格上）
            std::vector<double> tilts;
            tilts.push_back(tilt0);
            for (int i : P) {
                double v = std::round(need[i] * 1000.0) / 1000.0;
                if (std::fabs(v - tilt0) > 1e-6) tilts.push_back(v);
            }
            std::sort(tilts.begin(), tilts.end());
            tilts.erase(std::unique(tilts.begin(), tilts.end()), tilts.end());
            if (degenerate) tilts.assign(1, 0.0);

            // 候选方位：向上(0) + 站位侧 + 阶段 1 兜底投出来的方位（保持旧口径能救的行仍能救）
            std::vector<double> azs;
            azs.push_back(0.0);
            // 方位候选的白名单边界：**工艺上「上翻方向」的左右极限是 ±60°**。
            // 为什么必须卡：方位决定枪从哪一侧进，超过 90° 就是"从侧面/下面斜插"，
            // 开口盒子 face3 实测会选出方位 −120° 的姿态 —— 它的关节超限点数确实最少
            // （S=400 下 107，方位 0 是 162），但整面 2080 点的 TCP 要绕面点转 120°，
            // 与相邻面组（方位 0）交界处一次甩 152mm，比它省下的那点关节余量贵得多。
            // 可涂性可以用倾角换，不该用"整面歪着进"换。
            const double kAzLimitDeg = 60.0;
            auto pushAz = [&](double v) {
                if (!std::isfinite(v)) return;
                while (v <= -180.0) v += 360.0;
                while (v > 180.0) v -= 360.0;
                if (std::fabs(v) > kAzLimitDeg + 1e-9) return;
                for (double e : azs) if (std::fabs(e - v) < 1e-6) return;
                azs.push_back(v);
            };
            if (!degenerate && p.azStance) {
                int im = P[P.size() / 2];
                V3 fm = t.pos[im] + t.gun[im] * p.standoff;
                pushAz(stanceAzDeg(nFace, fm, ev.baseOrigin(), uUp, axis));
            }
            if (!degenerate && p.azFallback) {
                for (int i : P) {
                    if (std::fabs(azOf[i]) > 1e-9) {
                        pushAz(azOf[i]);
                        if (azs.size() >= 4) break;
                    }
                }
            }
            std::sort(azs.begin(), azs.end());
            if (degenerate) azs.assign(1, 0.0);

            // ---- 候选评估：每一档都只在**同一采样数**内做比较 ----
            // 混着比会挑出"只在小样本上好看"的姿态：开口盒子 +X 面实测把 (55°,60°) 的 P0=10
            // 让位给了 (30°,0°) 的"8 点全不可达"，最后整面 2103 点全部不可达。
            const int S1 = 32;
            std::vector<Cand> pool;
            pool.reserve(tilts.size() * azs.size());
            for (double tt : tilts)
                for (double aa : azs) pool.push_back({tt, aa, scoreFace(g, tt, aa, S1)});
            std::sort(pool.begin(), pool.end(), cmpCand);

            // 常规候选全都过不了 P0 → 广搜：倾角整档 × 每 30° 方位。
            // 实测开口盒子 +X 面就是这种情况：常规候选里要么整面够不着，要么整片穿透，
            // 而"够得着且只蹭到几个点"的那个姿态 (55°, 方位 60°) 根本不在候选集里
            // （阶段 1 的兜底只在需用倾角上试方位 ±30/±60，且从没试过 55° 与方位 60° 的组合）。
            if (!degenerate && pool.front().s.p0() > 0) {
                std::vector<Cand> wide;
                for (double tt = tilt0; tt <= p.tiltMaxDeg + 1e-9; tt += p.tiltStepDeg) {
                    for (int a = -180; a < 180; a += 30) {
                        double aa = (double)a;
                        // 与 pushAz 同一条工艺边界：只允许"上翻方向左右 60° 内"
                        if (std::fabs(aa) > kAzLimitDeg + 1e-9) continue;
                        bool dup = false;
                        for (const Cand& c : pool) {
                            if (std::fabs(c.tilt - tt) < 1e-6 && std::fabs(c.az - aa) < 1e-6) { dup = true; break; }
                        }
                        if (dup) continue;
                        wide.push_back({tt, aa, scoreFace(g, tt, aa, S1)});
                    }
                }
                for (const Cand& c : wide) pool.push_back(c);   // 同为 S1，可直接同比
                std::sort(pool.begin(), pool.end(), cmpCand);
            }

            // P0 差在"2% 抽样 + 2 点"以内的候选视为同一档。抽样噪声不该决定整面的姿态：
            // 先把这一档捞到前面，再在档内按（工艺喷角 → 间隙 → 倾角）排；
            // 否则一个大倾角的"碰巧少撞两点"会把小倾角的好候选挤出复核名单
            // （longstripe 实测：整面 11290 点，59° P0=9、49° P0=4、30° 更差，
            //  但 59° 让整面喷角超工艺上限 —— 这种整面违规不能由 1~2 个点的噪声决定）。
            {
                int minP0c = pool.front().s.p0();
                for (const Cand& c : pool) minP0c = std::min(minP0c, c.s.p0());
                int bandC = std::max(2, (int)(0.02 * S1 + 0.5));
                std::stable_partition(pool.begin(), pool.end(), [&](const Cand& c) {
                    return c.s.p0() <= minP0c + bandC;
                });
                std::sort(pool.begin(), pool.end(), cmpCand);
            }

            // 密集复核：粗采样分不出的差距，密集采样常常分得出来。
            //
            // keep 必须足够宽。开口盒子 face3 的候选池里，真正的第一名是
            // `(55°, 30°)`（S=400 下 P0=6），但它在 S=16 粗筛时是 6/16 失败、
            // 排到第 9 —— 而终选只在 keep 名前挑，于是选走了 (60°, −120°)
            // （P0=107，整面斜着进、和相邻面组差 120°，最大突跳 152mm）。
            // 采样越少，"失败 6 个"这种绝对数越不可比；把粗筛抬到 32、复核扩到 16
            // 名，就是为了让这类候选有机会被同一把尺子重新量一次。
            int keep = std::min<int>(16, (int)pool.size());
            int Sfinal = std::min(96, (int)P.size());
            if (Sfinal > S1) {
                for (int k = 0; k < keep; ++k) pool[k].s = scoreFace(g, pool[k].tilt, pool[k].az, Sfinal);
                std::sort(pool.begin(), pool.begin() + keep, cmpCand);
            }
            // 还有 P0 → 再上一档采样复核（只在同一档内部比较，量纲才一致）
            if (pool[0].s.p0() > 0 && Sfinal < (int)P.size()) {
                int Sbig = std::min(400, (int)P.size());
                if (Sbig > Sfinal) {
                    for (int k = 0; k < keep; ++k) pool[k].s = scoreFace(g, pool[k].tilt, pool[k].az, Sbig);
                    std::sort(pool.begin(), pool.begin() + keep, cmpCand);
                    Sfinal = Sbig;
                }
            }
            // ---- 最终定档：同上，P0 在容忍带内视为同档，同档内先守工艺喷角
            int minP0 = pool[0].s.p0();
            for (int k = 1; k < keep; ++k) minP0 = std::min(minP0, pool[k].s.p0());
            int band = std::max(2, (int)(0.02 * Sfinal + 0.5));
            auto pickBetter = [&](const Cand& x, const Cand& y) {
                int px = x.s.p0(), py = y.s.p0();
                bool inx = px <= minP0 + band, iny = py <= minP0 + band;
                if (inx != iny) return inx;
                bool ox = x.tilt <= p.maxSprayAngleDeg + 1e-9;
                bool oy = y.tilt <= p.maxSprayAngleDeg + 1e-9;
                if (ox != oy) return ox;
                if (x.s.shortClr != y.s.shortClr) return x.s.shortClr < y.s.shortClr;
                if (px != py) return px < py;
                if (std::fabs(x.tilt - y.tilt) > 1e-6) return x.tilt < y.tilt;
                return std::fabs(x.az) < std::fabs(y.az) - 1e-9;
            };
            // ---- 抽样打分毕竟不是真解：对前 3 名再做一次「连续段热种子复核」
            //
            // 抽样点彼此间隔上百毫米，"接着上一个抽样点当种子"和阶段 3 里"接着上一个真实点
            // 当种子"根本不是一回事。这一步的代价是实测出来的：`longwide` 单层 12476 点，
            // 抽样口径下 10° 优于 20°，但真解里 10° 有 15 点不可达（HEAD 是 PASS，它变
            // REJECT）。所以前 3 名各取 3 段 64 个连续点，按阶段 3 的方式解一遍，硬约束
            // 全过的优先；都不全过就仍按抽样排序取第一（不因为复核失败而放弃整面）。
            auto verifyFace = [&](double vt, double vz) -> int {
                const std::vector<int>& Q = gp[g];
                const int len2 = (int)Q.size();
                if (len2 <= 0) return 0;
                const int kW = 64;
                int fails = 0;
                for (int w = 0; w < 3; ++w) {
                    int c = (len2 <= kW) ? 0 : (int)((long long)len2 * (w * 2 + 1) / 6);
                    int a2 = std::max(0, std::min(len2 - kW, c - kW / 2));
                    if (len2 <= kW) a2 = 0;
                    int b2 = std::min(len2 - 1, a2 + kW - 1);
                    double seed2[6] = {0, 0, 0, 0, 0, 0};
                    bool has2 = false;
                    for (int k = a2; k <= b2; ++k) {
                        int i = Q[k];
                        V3 n = t.normal(i);
                        V3 face = t.pos[i] + t.gun[i] * p.standoff;
                        V3 axis, uUp;
                        tiltAxisAndDir(n, p.tiltRefW, axis, uUp);
                        if (norm(uUp) < 1e-6) continue;
                        V3 d = dirOf(n, uUp, axis, vt, vz);
                        PR pr;
                        if (!ev.eval(face, d, nullptr, has2 ? seed2 : nullptr, vt, vz, true, pr)) {
                            ++fails;
                            has2 = false;
                            continue;
                        }
                        if (!hardOk(pr, p)) ++fails;
                        for (int q = 0; q < 6; ++q) seed2[q] = pr.theta[q];
                        has2 = true;
                    }
                }
                return fails;
            };

            int bestK = 0;
            for (int k = 1; k < keep; ++k)
                if (pickBetter(pool[k], pool[bestK])) bestK = k;

            std::vector<int> order(keep);
            for (int k = 0; k < keep; ++k) order[k] = k;
            std::stable_sort(order.begin(), order.end(),
                             [&](int x, int y) { return pickBetter(pool[x], pool[y]); });
            int verifyFail = -1, verified = 0;
            for (int t2 = 0; t2 < std::min(3, keep); ++t2) {
                int vf = verifyFace(pool[order[t2]].tilt, pool[order[t2]].az);
                ++verified;
                if (t2 == 0) verifyFail = vf;
                if (vf == 0) { bestK = order[t2]; verifyFail = 0; break; }
            }
            if (bestK != 0) std::swap(pool[0], pool[bestK]);
            faceTilt[g] = pool[0].tilt;
            faceAz[g] = pool[0].az;
            if (const char* dbgG = std::getenv("RS_DBG_FACE")) {
                if (atoi(dbgG) == (int)g) {
                    printf("[pool%zu] 候选池前 10（前 %d 名按 S=%d 复核，其余仍是 S=%d 粗筛，两者不可比）:\n",
                           g, keep, Sfinal, S1);
                    for (int k = 0; k < (int)pool.size() && k < 10; ++k)
                        printf("   [S=%3d] tilt=%5.1f az=%7.1f P0=%4d 撞%4d 不可达%4d 自撞%3d 关节%4d 间隙不足%3d 最小%7.1f\n",
                               (k < keep ? Sfinal : S1),
                               pool[k].tilt, pool[k].az, pool[k].s.p0(), pool[k].s.hit,
                               pool[k].s.unreach, pool[k].s.selfHit, pool[k].s.over,
                               pool[k].s.shortClr, pool[k].s.minClr);
                    fflush(stdout);
                }
            }
            for (int r = ra; r <= rb; ++r) {
                azRow[r] = faceAz[g];
                for (int i = rows[r].first; i <= rows[r].second; ++i)
                    if (t.spray[i]) finalTilt[i] = faceTilt[g];
            }
            if (g_prof.on) {
                std::fprintf(stderr, "  [face%2zu] 行 %3d..%3d 点%5d 倾角%5.1f 方位%7.1f | 抽%d: P0=%d(撞%d/不可达%d/自撞%d) 间隙不足%d 最小%6.1f | 复核%d名 首选失败%d\n",
                             g, ra, rb, (int)P.size(), faceTilt[g], faceAz[g], Sfinal,
                             pool[0].s.p0(), pool[0].s.hit, pool[0].s.unreach, pool[0].s.selfHit,
                             pool[0].s.shortClr, pool[0].s.minClr, verified, verifyFail);
            }
        };

        if (nThreads <= 1 || groups.size() <= 2) {
            for (size_t g = 0; g < groups.size(); ++g) decideGroup(g);
        } else {
            std::atomic<size_t> nextGroup{0};
            std::vector<std::thread> pool;
            int nt = std::max(1, std::min<int>(nThreads, (int)groups.size()));
            pool.reserve(nt);
            for (int k = 0; k < nt; ++k) {
                pool.emplace_back([&] {
                    for (;;) {
                        size_t g = nextGroup.fetch_add(1);
                        if (g >= groups.size()) break;
                        decideGroup(g);
                    }
                });
            }
            for (auto& th : pool) th.join();
        }
    }
    g_prof.mark("阶段2 面级姿态");

    int azChangedRows = 0;
    for (size_t r = 0; r < azRow.size(); ++r)
        if (std::fabs(azRow[r]) > 1e-9) ++azChangedRows;

    // 统一取「第 i 点的方位角」。
    // ⚠ 这里**不区分 J5 模式**：J5 现在是叠加层（见阶段 3.6），基线阶段一律走行级 azRow。
    //   第一版在这里按 p.tiltJ5 分流到逐点数组，等于让 J5 模式绕过基线方位决策 ——
    //   那是"取代"而不是"叠加"，正是复杂件退化的根源。
    auto azAt = [&](int i) -> double {
        return (rowOf[i] >= 0) ? azRow[rowOf[i]] : 0.0;
    };
    // 统一取「切平面基」。基线一律用原判据（退化→倾角无意义，恒 0），
    // 与 needTilt/tiltNeedDeg 的退化口径一致。J5 叠加层自己用 tangentBasisRobust。
    auto basisAt = [&](const V3& nrm, V3& axis, V3& uUp) {
        tiltAxisAndDir(nrm, p.tiltRefW, axis, uUp);
    };
    if (g_prof.on) {
        int nonZeroPts = 0;
        for (int i = 0; i < n; ++i)
            if (rowOf[i] >= 0 && std::fabs(azRow[rowOf[i]]) > 1e-9) ++nonZeroPts;
        // 行内倾角最大步进：面级常量的直接验收指标，正常应恒为 0
        double maxStep = 0.0;
        for (size_t r = 0; r < rows.size(); ++r)
            for (int i = rows[r].first + 1; i <= rows[r].second; ++i)
                maxStep = std::max(maxStep, std::fabs(finalTilt[i] - finalTilt[i - 1]));
        std::fprintf(stderr, "  [az] rows=%zu 面组=%zu 非零方位行=%d 受影响点=%d · 行内倾角最大步进=%.1f°\n",
                     rows.size(), groups.size(), azChangedRows, nonZeroPts, maxStep);
    }
    // ---- 阶段 2.5：倾角的连续性收尾
    //
    // 这一步存在的唯一理由：**倾斜会移动位置**。倾斜 ≡ 保持面点、保持喷距、把枪轴向
    // 「正位置」扳 tilt 度，于是 TCP = 面点 + 枪轴·standoff 必然被推出去
    // standoff·2sin(tilt/2)（standoff=100、tilt=43° → 73.5mm）。原始轨迹的 TCP 本来就
    // 落在面点上，所以"手腕原地翻"在 raw 里几乎不动位置；一旦倾斜，它就变成绕面点走
    // 半径 standoff·sin(tilt) 的球面，**raw 的姿态噪声被 1:1 放大成位置噪声**。
    //
    // 垃圾桶实测（734 点）：raw 相邻步长中位 4.00、最大 5.08，一条跳点都没有；
    // 修正后出现 30 段 >10mm、最大 74.90mm。两个来源：
    //   ① 换段边界：喷涂末点被推出 73.5mm，紧邻的空移点按设计 tilt 恒 0、位置一字不动
    //      （空移点不倾斜是刻意的，见上文"置 0 后逐位回到原轨迹姿态"），一步 74.9mm。
    //   ② 桶底 Ø9mm 锥形小坑：raw 沿坑壁绕一圈，位置只走 3mm 但枪轴逐点翻 20~58°、
    //      方位绕满 360°（枪轴 z 分量全程恒定 -0.450，是锥面无疑）。倾斜后 TCP 落到
    //      Ø110mm 的圈上，步长 22~55mm。
    //
    // 两处都不撞件（步长 >10mm 的 30 段，段到点云最小距离 73.2mm 起），是节拍/关节
    // 平滑问题，不是碰撞问题。所以下面只做"让轨迹连续"这一件事，不改任何安全结论。
    //
    // RS_NO_CONT=1 可关掉整段（**只用于 A/B 回归**：同一份二进制跑出"收尾前/后"两组，
    // 否则对比只能跨版本比，机型/参数一变动结论就不可信）。生产不要设这个变量。

    // 枪体专用评估器：requireReach=0 → 只解析计算「TCP + 枪体三段碰撞体」与点云的间隙，
    // 不解 IK（闸门用）。与机器人构型无关，天然无解支噪声。仅供阶段 2.5 的闸门使用。
    Params pTool = p;
    pTool.requireReach = 0;
    Evaluator evTool(cloud, base, pTool);
    if (std::getenv("RS_NO_CONT") == nullptr) {
        // 喷涂连续段：**全局**按 spray 序列切，不能按行切。
        // onepass 模式下整条轨迹是一条连续喷涂路径（行与行之间没有空移点），
        // 按行切会把一条连续线切成 17 段，于是每段首尾各做一次 ramp ——
        // 实测垃圾桶 #88~#102 上直接出现 V 形低谷（tilt 40→6.6→21），
        // 本来连续的一条线被人为压成锯齿。空移点才是真正的分段依据。
        std::vector<std::pair<int, int>> runs;
        for (int i = 0; i < n;) {
            while (i < n && !t.spray[i]) ++i;
            if (i >= n) break;
            const int a = i;
            while (i < n && t.spray[i]) ++i;
            runs.push_back({a, i - 1});
        }
        const double Sw = p.standoff > 1e-6 ? p.standoff : 100.0;
        const double kStepBudget = 12.0;   // 相邻点 TCP 位移预算 mm
        // 把"位移预算"折成"允许的倾角步长"：倾斜矢量模 = 2·standoff·sin(tilt/2)
        const double dTiltCap = (kStepBudget < 2.0 * Sw)
                                    ? 2.0 * std::asin(0.5 * kStepBudget / Sw) * kRad
                                    : 90.0;

        const std::vector<double> origTilt = finalTilt;   // 闸门要用

        // ① 步长预算（默认关，RS_SHRINK=1 才启用）：按「法向变化率」反解每点的倾角上限，
        //    再对倾角序列做双向限速。
        //
        // 位移的准确表达式只有一个：TCP = pos + standoff·(d − n)，所以
        //     ΔTCP = Δpos + standoff·Δ(d − n)
        // 前一项是 raw 的步长（本来就平滑，垃圾桶中位 4.00mm）；后一项是**倾斜矢量**
        // v = standoff·(d − n) 的变化，它有两部分：
        //   模：|v| = standoff·2sin(tilt/2)  —— 随倾斜本身长出来
        //   方向：跟着法向方位转，转角 Δφ 与法向完全相同
        // 所以 |Δv| ≈ standoff·sin(tilt)·Δφ + standoff·Δtilt。**两个方向都是 tilt 越小越好**：
        // tilt→0 时 v→0，TCP 逐位回到 raw 位置（pos + 0 = 原始轨迹点，而原始轨迹是平滑的）。
        //
        // **默认关闭**的原因：倾角在这台臂上很多时候不只是"扳正"，它就是抬枪避让的手段。
        // 18 件 A/B 实测，无条件收缩会把间隙整片搞坏（垃圾桶 32.5→2.4mm、
        // 开口盒子 −14.2→−32.2mm、槽型盘 20.6→6.8mm），而且逐点收缩会在序列上
        // 打出台阶，反而造出新的步长跳。边界斜坡（②）没有这个问题 —— 它只动喷涂段
        // 首尾，那些点间隙通常很空，而且下面有整段闸门兜底。
        //
        // ⚠ 也不要改成"抬高倾角"。抬 tilt 会让 d 更贴近竖直、d 的夹角变小，看着像在
        //   变平滑，但倾斜矢量反而更大、TCP 被推得更远 —— 实测把 #96 抬到 60° 后它
        //   比 raw 偏了 100mm，段内最大位移从 66mm 涨到 71mm，越修越糟。
        int shrunkPts = 0;
        if (std::getenv("RS_SHRINK") != nullptr) {
            for (const auto& rn : runs) {
                const int a = rn.first, b = rn.second;
                // a) 按法向转角反解倾角上限：standoff·sin(t)·Δφ ≤ budget
                for (int i = a + 1; i <= b; ++i) {
                    const double dphi = angleBetweenDeg(t.normal(i), t.normal(i - 1)) * kDeg;
                    if (dphi < 1e-6) continue;
                    const double s = kStepBudget / (Sw * dphi);
                    const double capT = (s >= 1.0) ? 90.0 : std::asin(s) * kRad;
                    if (finalTilt[i] > capT) {
                        if (finalTilt[i] - capT > 0.5) ++shrunkPts;
                        finalTilt[i] = capT;
                    }
                }
                // b) 双向限速（先前向再后向，等价于对"V 形低谷"做两侧斜坡）
                for (int i = a + 1; i <= b; ++i)
                    finalTilt[i] = std::min(finalTilt[i], finalTilt[i - 1] + dTiltCap);
                for (int i = b - 1; i >= a; --i)
                    finalTilt[i] = std::min(finalTilt[i], finalTilt[i + 1] + dTiltCap);
            }
        }

        // ② 边界斜坡：每个喷涂段首尾强制 tilt=0，再用 ±dTiltCap 双向限速生成斜坡。
        //
        // 为什么不是"削掉首尾各 N 点、中间用 smoothstep"：N 是按 tmax/dTiltCap 算的，
        // 短段根本装不下 N 个点，于是出现「首点削到 1.7°、第二点还是 49.67°」——
        // 实测复杂点云 #1540→#1541 反而从 74.10mm 涨到 81.40mm，比不修还差。
        // 「首尾置 0 + 双向限速」自动适应段长：
        //   长段 → 斜坡-平台-斜坡，中间倾角一字不动；
        //   短段 → 整体被压成三角形，峰值 = dTiltCap × 段长/2，绝不会出台阶。
        //
        // tilt=0 的点位置 = 面点 + n·standoff = raw 位置，与紧邻的空移点（tilt 恒 0、
        // 位置一字不动）天然接上 —— 换段那一跳就是这么消掉的。实测垃圾桶
        // #309→#310 从 74.91mm 消掉（#303~#309 七点 43→0 渐变，间隙 63.5→66.0 反而变好）。
        //
        // ⚠ 轨迹首点/末点**不做斜坡**：它们前面/后面没有空移点要衔接，削了纯属白削。
        // 实测过反例：垃圾桶 #0 被削到 21°，连杆间隙 46.6→6.9mm —— 桶腔里枪更竖直
        // 意味着手腕要探得更深。斜坡的收益只存在于「喷涂段 ↔ 空移段」的接缝上。
        int rampRuns = 0, rampPts = 0, raisedPts = 0;
        for (const auto& rn : runs) {
            const int a = rn.first, b = rn.second;
            if (b - a < 1) continue;
            const std::vector<double> keep(finalTilt.begin() + a, finalTilt.begin() + b + 1);
            const bool rampA = (a > 0);         // a>0 时 t.spray[a-1] 必为 0（run 是极大喷涂段）
            const bool rampB = (b < n - 1);
            if (rampA) finalTilt[a] = 0.0;
            if (rampB) finalTilt[b] = 0.0;
            for (int i = a + 1; i <= b; ++i)
                finalTilt[i] = std::min(finalTilt[i], finalTilt[i - 1] + dTiltCap);
            for (int i = b - 1; i >= a; --i)
                finalTilt[i] = std::min(finalTilt[i], finalTilt[i + 1] + dTiltCap);

            // 闸门（逐点，口径 = **枪体间隙**，用 requireReach=0 的轻量评估器）：
            //
            // ⚠ 为什么要闸门：tilt=0 的点位置 = raw 位置、枪轴 = raw 枪轴，但**枪体不一定安全**
            // ——垃圾桶就是现成的反例：空腔件里枪沿法向伸进去，枪身会蹭到对面桶壁，
            // #0 的 raw 姿态枪体间隙 −20.7mm，阶段 2 给的 41° 倾角恰恰是把枪从桶里拉出来。
            // 所以「斜坡把端点 tilt 归 0」在这个件上天然行不通，必须逐点验枪体。
            //
            // ⚠ 为什么只查枪体、不解 IK：第一版用全位姿评估（冷启 IK + frameFromZ 正则手腕），
            // 同一个 #0 被报 −20.7mm「穿件」之外还多了一层假阳性——冷启选了另一支肘。
            // 手臂构型是「连续种子」的领地（阶段 3），这里插手只会用阶段 3 根本不会
            // 采用的构型否决正确的斜坡。枪体三段碰撞体由「TCP+枪轴」解析确定，
            // 与机器人无关，是这里唯一可信的口径。
            //
            // 处理方式不是整段回滚（实测过：回滚会在序列上留台阶，反而造新跳），
            // 而是**就地二分上搜**：在 [目标值, 原值] 里找能过安全线的最小倾角。
            // 间隙随倾角基本单调（倾角越大枪越出腔），原值是阶段 2 选的、必然过线，
            // 所以二分必有解；找到的就是「安全前提下最贴近斜坡目标」的值。
            int changed = 0, raisedPts = 0;
            for (int i = a; i <= b; ++i) {
                const double o = keep[i - a];
                if (std::fabs(finalTilt[i] - o) < 1e-6) continue;
                const V3 nrm = t.normal(i);
                const V3 face = t.pos[i] + t.gun[i] * Sw;
                V3 axis, uUp;
                basisAt(nrm, axis, uUp);
                if (norm(uUp) < 1e-6) continue;   // 法向退化时 dirOf 恒返回 n，tilt 无效
                ++changed;
                const double az = azAt(i);
                auto clrAt = [&](double tt) {
                    PR pr;
                    evTool.eval(face, dirOf(nrm, uUp, axis, tt, az),
                                nullptr, nullptr, tt, az, false, pr);
                    return pr.clearance;
                };
                // 地板：最多花 2mm 余量，且不许低于 min(原值, safety)。
                // 第一版只用了 min(原值, safety)，18 件实测暴露出它太松：
                // 12-34-54 #1146 原本 71.7mm 的枪体间隙被削到 22.4mm（贴着 safety=20 的地板）——
                // 为节拍/观感花掉 49mm 安全余量是烂交易，平滑是锦上添花，余量是保命钱。
                // 收紧成 clrA−2 后，二分会自动把倾角顶回「间隙刚好还剩 clrA−2」的位置：
                // 平滑改善保留能白拿的部分（垃圾桶 #303~#309 削倾角间隙反而 63.5→66.0，全保留），
                // 要花钱的部分自动放弃（那一跳留给 STEP_JUMP 如实报出来）。
                const double clrA = clrAt(o);
                const double floorClr = std::max(std::min(clrA, p.safety), clrA - 2.0);
                if (clrAt(finalTilt[i]) < floorClr) {
                    double lo = finalTilt[i], hi = o;   // 不变量：lo 过不了线，hi 过得了
                    for (int it = 0; it < 20; ++it) {
                        const double mid = 0.5 * (lo + hi);
                        if (clrAt(mid) < floorClr) lo = mid;
                        else hi = mid;
                    }
                    finalTilt[i] = hi;
                    ++raisedPts;
                }
            }
            if (changed > 0) {
                ++rampRuns;
                rampPts += changed;
            }
        }
        if (g_prof.on) {
            double maxStep = 0.0;
            for (const auto& rn : runs) {
                bool hasPrev = false;
                V3 prevT(0, 0, 0);
                for (int i = rn.first; i <= rn.second; ++i) {
                    const V3 nrm = t.normal(i);
                    const double az = azAt(i);
                    V3 axis, uUp;
                    basisAt(nrm, axis, uUp);
                    if (norm(uUp) < 1e-6) { hasPrev = false; continue; }
                    const V3 dd = dirOf(nrm, uUp, axis, finalTilt[i], az);
                    const V3 T = t.pos[i] + (dd - nrm) * Sw;
                    if (hasPrev) maxStep = std::max(maxStep, norm(T - prevT));
                    prevT = T;
                    hasPrev = true;
                }
            }
            std::fprintf(stderr,
                         "  [cont] 喷涂段=%zu 生效斜坡=%d 段/%d 点 安全回抬=%d 点 收缩=%d 点 · "
                         "段内相邻位移 max %.1fmm（预算 %.1f）\n",
                         runs.size(), rampRuns, rampPts, raisedPts, shrunkPts, maxStep, kStepBudget);
        }
    }
    g_prof.mark("阶段2.5 连续性收尾");

    // ---- 阶段 3：按 final 重新求解 + 归集违规
    std::vector<PR> results(n);
    std::vector<char> hasResult(n, 0);
    VioCount vc;
    double standoffMaxErr = 0.0;
    int recomputed = 0, reverted = 0;
    int reachable = 0;

    // 逐行独立（行内仍需种子连续性，行间互不依赖）→ 可以按行并行。
    // 每行行首没有"上一点解"，用解析种子自启；行首本来就是换行处，姿态本来就会跳。
    struct Acc3 {
        int recomputed = 0, reverted = 0, reachable = 0;
        double standoffMaxErr = 0.0;
        // 面级姿态在这些点上「为什么过不了硬约束」的分类。
        // 用来回答一个反复被问的问题：逐点回退到底是「真撞了不得不救」，
        // 还是「只差一点余量就被当成撞」。分类互斥，按严重度优先计入一条。
        int ffPierce = 0;   // 间隙 < 0：真的穿进去了
        int ffSelf = 0;     // 自碰撞（胶囊近似间的重叠）
        int ffClr = 0;      // 0 <= 间隙 < hardClearance：只是余量不足，没撞
        int ffJoint = 0;    // 关节超限
        int ffIk = 0;       // 解不出逆解
        // 「过不了硬约束 → 同方位 ±倾角步长 救援」的产出比。
        // 只回答一个问题：这些救援 eval 有多少真的把点救回来了（其余的纯烧时间）。
        int rescueEvals = 0, rescueHit = 0;
        // 按**救援前的失败类型**分桶（互斥，口径与 hardOk 一致，同 ff* 那套）。
        // 超限深度不预示救援成败（实测各深度档都是 0.2% 左右），
        // 但失败类型可能预示得很准：IK 解不出来 / 真穿透 / 自碰撞 / 关节超限 / 只差余量，
        // 五种情况里"改倾角能不能救"完全不同。分桶后按救回率决定跳过哪一类。
        int rkEvals[5] = {0, 0, 0, 0, 0};
        int rkHit[5] = {0, 0, 0, 0, 0};
        // 被救回的点里，救援前超限最深的那个（margin 最负）。用来定"超限跳过"的阈值：
        // 阈值只要不比它更负，就一个救得回来的点都不会漏。
        double hitMarginWorst = 0.0;
    };

    // 违规位单点判定。阶段 3 和后面的「回退段统一姿态」都会改位姿，两处必须同口径，
    // 所以集中在这里算一次，别各写一份。
    //
    // isSpray=false（空移点）：只判「撞不撞、走不走得到」——IK / 工具撞 / 连杆撞 /
    // 自碰撞 / 关节限位 / 奇异。**不判** CLEARANCE（喷距）与 SPRAY_ANGLE（喷角）：
    // 这两个是"喷出来的膜厚对不对"，而空移段枪是关的，报出来全是假告警
    // （实测开口盒子空移点被误报 SPRAY_ANGLE 636 点、CLEARANCE 一整片）。
    auto setViolBits = [&](PR& pr, bool isSpray) {
        const unsigned int kMask = RS_VIOL_IK_FAIL | RS_VIOL_TOOL_HIT | RS_VIOL_LINK_HIT |
                                   RS_VIOL_CLEARANCE | RS_VIOL_SELF_COLLISION |
                                   RS_VIOL_JOINT_LIMIT | RS_VIOL_SPRAY_ANGLE |
                                   RS_VIOL_NEAR_SINGULAR;
        pr.viol &= ~kMask;
        // 不要求可达：IK_FAIL / SELF_COLLISION / JOINT_LIMIT / LINK_HIT 这四位一律不置。
        // 它们全都是「这台机器人 + 这条站位」的结论，点云本身没有对错。
        // 几何层照旧：枪体碰撞、间隙不足、工艺喷角都还照实报（这三条换任何机器人都成立）。
        if (!p.requireReach) {
            if (pr.clrTool < 0.0) pr.viol |= RS_VIOL_TOOL_HIT;
            if (isSpray) {
                if (pr.clearance >= 0.0 && pr.clearance < p.safety) pr.viol |= RS_VIOL_CLEARANCE;
                if (pr.sprayAngle > p.maxSprayAngleDeg + 1e-9) pr.viol |= RS_VIOL_SPRAY_ANGLE;
            }
            return;
        }
        if (!pr.ikOk) { pr.viol |= RS_VIOL_IK_FAIL; return; }
        // 穿透按部位细分：连杆/关节撞 与 枪体/喷嘴撞 是两类完全不同的现场问题
        // （前者要改站位或姿态，后者只要抬高喷距/换角度），混在一个 COLLISION 里没法决策。
        // RS_VIOL_COLLISION 保留定义但不再置位，避免 byRule 重复计数。
        if (pr.clrTool < 0.0) pr.viol |= RS_VIOL_TOOL_HIT;
        if (pr.clrLink < 0.0) pr.viol |= RS_VIOL_LINK_HIT;
        if (isSpray && pr.clearance >= 0.0 && pr.clearance < p.safety) pr.viol |= RS_VIOL_CLEARANCE;
        if (pr.selfGap < 0.0) pr.viol |= RS_VIOL_SELF_COLLISION;
        if (pr.over || pr.margin < p.jointMarginMin) pr.viol |= RS_VIOL_JOINT_LIMIT;
        if (isSpray && pr.sprayAngle > p.maxSprayAngleDeg + 1e-9) pr.viol |= RS_VIOL_SPRAY_ANGLE;
        if (pr.cond > p.condMax) pr.viol |= RS_VIOL_NEAR_SINGULAR;
    };

    auto runRange3 = [&](int a, int b, Acc3& acc) {
        double prevTheta[6] = {0, 0, 0, 0, 0, 0};
        bool hasSeed = false;
        V3 prevX;
        bool hasRefX = false;
        for (int i = a; i <= b; ++i) {
            V3 nrm = t.normal(i);
            V3 face = t.pos[i] + t.gun[i] * p.standoff;
            double tilt = finalTilt[i];
            double az = azAt(i);

            V3 axis, uUp;
            basisAt(nrm, axis, uUp);
            if (norm(uUp) < 1e-6) {
                axis = V3(0, 0, 0);
                uUp = V3(0, 0, 0);
            }

            bool reuse = hasCache[i] && std::fabs(cache[i].tiltDeg - tilt) < 1e-6 &&
                         std::fabs(cache[i].azDeg - az) < 1e-6;
            // 阶段 1 的缓存是在“阶段 1 的倾角 + 阶段 1 的上一点解”下算出来的；
            // 阶段 3 的 final 倾角可能已经变了，上一点也可能被重解。若直接复用缓存，
            // 会出现“上一点已换到 A 支、本点缓存仍在 B 支”的关节跳变（complex idx 3161：
            // J6 从 32° 跳到 8°，但目标只转了 3°）。缓存与当前上一点不连续时重解一次，
            // 再用硬约束决定是否退回缓存。
            if (reuse && hasSeed &&
                maxJointDelta(cache[i].theta, prevTheta) > kCacheContinuityTrigDeg)
                reuse = false;
            PR pr;
            if (reuse) {
                pr = cache[i];
                pr.idx = i;
            } else {
                V3 d = dirOf(nrm, uUp, axis, tilt, az);
                bool checkLinks = (i % std::max(1, p.linkStride) == 0);
                bool ok = ev.eval(face, d, hasRefX ? &prevX : nullptr, hasSeed ? prevTheta : nullptr,
                                  tilt, az, checkLinks, pr);
                ++acc.recomputed;
                if (!ok) {
                    // ---- 面级姿态在这个点解不出逆解：先按「同方位、倾角最小偏离」救援，
                    // 再判不可达。面级常量是工艺要求，但不该让个别点直接失去可涂性。
                    // 实测 longwide（单层 12476 点）：面级 10° 下有 15 点解不出，20° 下全好；
                    // 抽样打分看不见 0.12% 的失败（400 抽样期望命中 0.5 个），只能在逐点解时救。
                    bool rescued = false;
                    // 单次逆解模式（fastPose）：救援只试**一步**（±1 个倾角步长），
                    // 保证每点最多 2 次 IK。完整的 ±1/±2/±3 三级搜救是质量优化，
                    // 不属于安全判定 —— 关掉不影响任何一条违规的判据。
                    // 空移点（spray=0）不救援：mMax=0 让循环直接不跑。
                    // 它的姿态是轨迹生成器为避障算好的，靠倾转"救"成可解等于把避障路径改掉
                    // （可能直接撞件）。解不出就照实报 IK_FAIL。
                    // 逐点连续口径下救援给到**上限**：线性值只是"法向要求的姿态"，
                    // 个别点过不了时允许单独往上抬（最多到 tiltMax），否则长板这类件
                    // 实测 P0 会从 93 涨到 825 —— 倾角在这些件上不只是工艺量，还是避障手段。
                    int mMax = 0;
                    if (t.spray[i]) {
                        if (p.fastPose) mMax = 1;
                        else if (p.tiltPerPoint)
                            mMax = (int)std::floor((p.tiltMaxDeg - tilt) / p.tiltStepDeg + 1e-9);
                        else mMax = 3;
                    }
                    for (int m = 1; m <= mMax && !rescued; ++m) {
                        for (int sgn = 1; sgn >= -1 && !rescued; sgn -= 2) {
                            double t2 = tilt + sgn * m * p.tiltStepDeg;
                            if (t2 < -1e-9 || t2 > p.tiltMaxDeg + 1e-9) continue;
                            V3 d2 = dirOf(nrm, uUp, axis, t2, az);
                            PR pr2;
                            if (!ev.eval(face, d2, hasRefX ? &prevX : nullptr,
                                         hasSeed ? prevTheta : nullptr, t2, az, checkLinks, pr2))
                                continue;
                            if (!hardOk(pr2, p)) continue;
                            pr2.idx = i;
                            pr2.sprayAngle = angleBetweenDeg(-pr2.gun, nrm);
                            pr2.reverted = true;
                            pr = pr2;
                            ++acc.reverted;
                            rescued = true;
                        }
                    }
                    if (!rescued) {
                        // 所有候选倾角都解不出逆解 → 不可达（P0）
                        PR bad;
                        bad.idx = i;
                        bad.tiltDeg = tilt;
                        // 不可达点没有真实位姿，但"指令倾角"是确定值：写进去，别让 tilt_deg 列
                        // 落成 0 —— 那会让"整行 60° 倾斜、其中 45 点不可达"看起来像"倾角在 0/60 之间跳"。
                        bad.sprayAngle = tilt;
                        bad.clearance = -std::numeric_limits<double>::infinity();
                        // 自碰撞标 nan 而不是 inf：inf 经 %f 落到 CSV 会变成一个天文数字，
                        // 看着像"量到了"，实际是没测（没有位姿就没有碰撞体）。
                        bad.selfGap = std::numeric_limits<double>::quiet_NaN();
                        bad.ikOk = false;
                        // 位置写「本该倾斜到的 TCP」= face + d·standoff，不写原始未倾斜点。
                        // 写成 t.pos[i] 会与相邻的倾斜点差 2·standoff·sin(tilt/2)（60° 时
                        // 整整 100mm），一条轨迹被拆成「倾斜壳」和「原始壳」两套，画出来
                        // 全是长刺，落盘执行时也是一串不可插补的突跳。
                        bad.pos = face + d * p.standoff;
                        bad.gun = -d;
                        bad.viol = RS_VIOL_IK_FAIL;
                        results[i] = bad;
                        hasResult[i] = 1;
                        continue;   // 违规位在输出阶段按 pr.viol 统一统计，这里不必重复累加
                    }
                }
                pr.idx = i;
                // ---- 平滑回退：只回退硬约束，不再“差一点就回退”
                // 用户现场口径：有凸点没关系，只要不碰撞就以平滑为主。
                // 所以这里把“平滑后间隙比本点最优小”从回退条件里拿掉：
                //   · 平滑位姿满足硬约束（IK / 自碰撞 / 间隙≥hardClearance / 关节不超限）
                //     → 保留平滑，哪怕间隙比本点最优小几毫米；
                //   · 只有硬约束被破坏，且本点最优更安全时，才退回本点最优。
                // 注意：margin < jointMarginMin 仍只是 P1 软告警，不触发回退；
                //       safety（默认 20mm）也只影响报告标橙，不再是回退阈值。
                if (p.fastPose || p.tiltPerPoint) {
                    // ---- 单次逆解 / 逐点倾角逐点救援：姿态过不了硬约束时，**只**在同方位上
                    // 试 ±1 个倾角步长（最多 2 次额外 eval），不做 ±2/±3、不退回本点最优、不翻方位。
                    //
                    // 逐点倾角模式的救援幅度同样是 ±1 步：倾角本来就是法向的连续函数，
                    // 大范围搜救会把"线性给的差值"改回台阶，正是这次要消灭的东西。
                    //
                    // 为什么保留这一步：开口盒子实测，完全不救援时 P0 从 868 涨到 4013
                    // （这些点的冲突是"换 5° 就能避开"的，不是几何死结）；加上 ±1 步回到可接受范围，
                    // 代价仍然每点 ≤2 次 IK。±2/±3 的收益递减、甩腕风险递增，fast 模式明确不做。
                    // 空移点（spray=0）**不救援**：它的姿态是轨迹生成器的避障结果，
                    // 靠倾转让它"更好解"等于换一条空移路径 —— 实测开口盒子 819 个空移点
                    // 就是这样被退回 35°/55° 的倾斜位姿（TCP 偏 57~105mm，还导致喷角假告警）。
                    // 真的撞就交给下面的「过渡点抬升」按 +Z 抬枪解决，姿态一列不改。
                    // ---- 关节已经超出限位过多、倾角救援够不着的点：不再试解
                    // 倾角改的是腕部朝向，改不了臂展。关节超限是"这块够不到"，
                    // 同方位的小幅倾角只能让最差关节再挪一点，救不回来。
                    // 实测背板件（9598 点）：不跳过时救援 20608 次 eval、救回 142 点；
                    // 跳过深超限的点后只剩 6623 次 eval、仍救回 139 点 ——
                    // 三分之二的救援开销换 3 个点，而整件本来就是 REJECT。
                    // 代价是实测的：背板 5/9678 点由"已救回"退回违规（其中 1 点 P0），
                    // cx 1/5560 点。判定结论不变（两边都是 REJECT）。
                    //
                    // 阈值不写死，按**本次救援的实际幅度**算：
                    //   救援最多把倾角改 ±rescueSpan（面级常量模式 = 2×步长 = 10°；
                    //   逐点倾角模式 = 一直到 tiltMax），倾角到关节的放大实测不到 2 倍，
                    //   所以超限超过 2×rescueSpan 就不可能靠这次救援回到限位内。
                    // 面级常量模式 → 阈值 20°，正是"滑回行程内"能成立的上限；
                    // 逐点倾角模式救援幅度大，阈值自动放宽（多数情况等于不跳过，行为不变）。
                    // 超限的点照实报 JOINT_LIMIT，判据一条不少，只是不再反复试解。
                    const double rescueSpan = p.tiltPerPoint
                        ? std::max(p.tiltStepDeg, p.tiltMaxDeg - tilt)
                        : 2.0 * p.tiltStepDeg;
                    const bool overFar = pr.over && pr.margin < -2.0 * rescueSpan;
                    if (!hardOk(pr, p) && t.spray[i] && !overFar) {
                        // fast 模式的救援：同方位、倾角 ±1/±2 步（最多 4 次 eval，只在失败点上花）。
                        //
                        // 为什么给到 ±2：面级姿态是**抽样**选的（每组 ≤8 点），单点级别的冲突
                        // （典型是 1 个点的自碰撞）抽样根本看不见。开口盒子实测 ±1 步 P0≈3969、
                        // ±1..±3 步 3869（多花 0.8s 只换 100 点，不值）；但槽型盘只有 1 个点自碰撞，
                        // ±1 步救不到（差 10°），±2 步正好够 —— 这是"非 REJECT 变 REJECT"的差别，
                        // 99% 的时间里用不上，用上了就值。
                        //
                        // 仍然不做：退回阶段 1 的逐点最优姿态（那必须先把每个点都解一遍）。
                        // 逐点连续口径：允许一路抬到 tiltMax —— 线性值只决定"起手姿态"，
                        // 这个点真过不了（自碰撞/连杆撞件）就单独往上抬，别为了一行连续
                        // 让整条轨迹不可用（长板实测：只 ±2 步时 P0=825，面级常量时 93）。
                        const int mMax2 = p.tiltPerPoint
                            ? (int)std::floor((p.tiltMaxDeg - tilt) / p.tiltStepDeg + 1e-9)
                            : 2;
                        const double marginBefore = pr.margin;   // 救援前超限有多深（诊断用）
                        // 救援前失败类型（与 ff* 同口径，互斥按严重度）
                        // 实测（背板件）：IK失败 0/0 · 穿透 163 点/救回 0 · 自碰撞 8/2
                        //               · 关节超限 2017/137 · 仅余量不足 13/0
                        // 也就是说**只有关节超限这一类救得回来**。但不要因此就
                        // "只在关节超限时救援"：穿透那 163 点只占 489 次 eval（7%），
                        // 剪掉它省不了什么，却要新增一条会漏救的判据。
                        const int rk = !pr.ikOk ? 0
                                     : (pr.clearance < 0.0 ? 1
                                     : (pr.selfGap < 0.0 ? 2
                                     : (pr.over ? 3 : 4)));
                        ++acc.rkEvals[rk];
                        for (int m = 1; m <= mMax2; ++m) {
                            bool fixed = false;
                            for (int sgn = -1; sgn <= 1 && !fixed; sgn += 2) {
                                double t2 = tilt + sgn * m * p.tiltStepDeg;
                                if (t2 < -1e-9 || t2 > p.tiltMaxDeg + 1e-9) continue;
                                V3 d2 = dirOf(nrm, uUp, axis, t2, az);
                                PR pr2;
                                ++acc.rescueEvals;
                                if (!ev.eval(face, d2, hasRefX ? &prevX : nullptr,
                                             hasSeed ? prevTheta : nullptr, t2, az, checkLinks, pr2))
                                    continue;
                                if (!hardOk(pr2, p)) continue;
                                pr2.idx = i;
                                pr2.sprayAngle = angleBetweenDeg(-pr2.gun, nrm);
                                pr2.reverted = true;
                                pr = pr2;
                                ++acc.reverted;
                                ++acc.rescueHit;
                                ++acc.rkHit[rk];
                                if (marginBefore < acc.hitMarginWorst) acc.hitMarginWorst = marginBefore;
                                fixed = true;
                            }
                            if (fixed) break;
                        }
                    }
                }
                else if (hasCache[i] && t.spray[i]) {
                    // 空移点不参与这条"退回本点最优"的回退链：cache[i] 是阶段 1 用
                    // 面级需求倾角算的解，退回去等于把避障姿态换成喷涂姿态（见上一条注释）。
                    // 空移点的 saf 修正只走「过渡点抬升」，位置改变、姿态不变。
                    bool smoothHard = hardOk(pr, p);
                    if (!smoothHard) {
                        // 按严重度归类「面级姿态为什么过不了」，互斥计数。
                        // 口径与 hardOk 严格一致：margin < jointMarginMin 只是 P1 软告警，
                        // 不触发回退，所以不能算进来（否则会把软告警说成"必须救"）。
                        if (!pr.ikOk) ++acc.ffIk;
                        else if (pr.clearance < 0.0) ++acc.ffPierce;
                        else if (pr.selfGap < 0.0) ++acc.ffSelf;
                        else if (pr.over) ++acc.ffJoint;
                        else if (pr.clearance < p.hardClearance) ++acc.ffClr;
                        // 先在「面级常量」附近找**偏离最小**的合法姿态（同方位，倾角 ±step/±2step/±3step）。
                        // 直接退回本点最优的话，倾角可能从 60° 跳回 35°，TCP 一次甩出上百毫米
                        // —— 那正是"修正后轨迹一团乱刺"的来源。物理上必须偏的点也要尽量少偏。
                        bool fixed = false;
                        for (int m = 1; m <= 3 && !fixed; ++m) {
                            for (int sgn = -1; sgn <= 1 && !fixed; sgn += 2) {
                                double t2 = tilt + sgn * m * p.tiltStepDeg;
                                if (t2 < -1e-9 || t2 > p.tiltMaxDeg + 1e-9) continue;
                                V3 d2 = dirOf(nrm, uUp, axis, t2, az);
                                PR pr2;
                                if (!ev.eval(face, d2, hasRefX ? &prevX : nullptr, hasSeed ? prevTheta : nullptr,
                                             t2, az, checkLinks, pr2)) continue;
                                if (!hardOk(pr2, p)) continue;
                                pr2.idx = i;
                                pr2.sprayAngle = angleBetweenDeg(-pr2.gun, nrm);
                                pr2.reverted = true;
                                pr = pr2;
                                ++acc.reverted;
                                fixed = true;
                            }
                        }
                        if (!fixed) {
                            // 退回「逐点最优」只在两种情况下做，且都不许白退：
                            //
                            // 1) 先试「面级方位 + 该点最优倾角」。方位不变，位移只有倾角那一点
                            //    （2·standoff·sin(Δ/2)，5° → 8.7mm）。
                            // 2) 缓存原样退回，但**只认方位相同的**，而且必须真过得了硬约束。
                            //
                            // 为什么卡这么死：方位是"从哪一侧进枪"，翻一次 TCP 就绕面点转
                            // 2·standoff·sin(Δaz/2)（60° → 100mm），画出来是长刺，落盘是一串
                            // 不可插补的突跳。而实测开口盒子 7070..7075 这 6 个点，面级姿态的
                            // self_gap 是 −0.38~−0.55、7072 的连杆还穿进去 −20.3mm —— 缓存原样
                            // 退回也只是把"撞"换成"另一种撞"（间隙 0.25~1.15mm，自碰撞照样不过），
                            // 却为此付出了 71mm 的甩腕。既然两个都过不了硬约束，就保持面级姿态：
                            // 姿态规整、违规照实报，由站位/基座去解决，而不是在这里甩姿态。
                            double cTilt = cache[i].tiltDeg;
                            if (std::fabs(cache[i].azDeg - az) > 1e-6) {
                                V3 d2 = dirOf(nrm, uUp, axis, cTilt, az);
                                PR pr2;
                                if (ev.eval(face, d2, hasRefX ? &prevX : nullptr,
                                            hasSeed ? prevTheta : nullptr, cTilt, az, checkLinks, pr2) &&
                                    hardOk(pr2, p)) {
                                    pr2.idx = i;
                                    pr2.sprayAngle = angleBetweenDeg(-pr2.gun, nrm);
                                    pr2.reverted = true;
                                    pr = pr2;
                                    ++acc.reverted;
                                    fixed = true;
                                }
                            }
                            if (!fixed && hardOk(cache[i], p)) {
                                pr = cache[i];
                                pr.idx = i;
                                pr.reverted = true;
                                ++acc.reverted;
                            }
                        }
                    }
                }
            }
            pr.sprayAngle = angleBetweenDeg(-pr.gun, nrm);
            acc.standoffMaxErr = std::max(acc.standoffMaxErr, pr.standoffErr);

            setViolBits(pr, t.spray[i] != 0);

            if (pr.ikOk) ++acc.reachable;
            results[i] = pr;
            hasResult[i] = 1;
            for (int k = 0; k < 6; ++k) prevTheta[k] = pr.theta[k];
            hasSeed = true;
            FkOut f;
            fk(pr.theta, p.toolLen, f);
            prevX = f.tcp.R.col(0);
            hasRefX = true;
        }
    };

    if (nThreads <= 1) {
        Acc3 acc;
        runRange3(0, n - 1, acc);
        recomputed = acc.recomputed;
        reverted = acc.reverted;
        reachable = acc.reachable;
        standoffMaxErr = acc.standoffMaxErr;
    } else {
        // 任务 = 每行一个区间；再补上没被任何行覆盖的孤立点（正常轨迹没有，防御性处理）
        std::vector<std::pair<int, int>> tasks;
        tasks.reserve(rows.size() + 8);
        std::vector<char> covered(n, 0);
        for (const auto& r : rows) {
            if (r.second < r.first) continue;
            tasks.push_back(r);
            for (int i = r.first; i <= std::min(r.second, n - 1); ++i) covered[i] = 1;
        }
        for (int i = 0; i < n; ++i)
            if (!covered[i]) tasks.push_back({i, i});

        std::vector<Acc3> accs(nThreads);
        std::atomic<int> nextTask{0};
        std::vector<std::thread> pool;
        pool.reserve(nThreads);
        for (int k = 0; k < nThreads; ++k) {
            pool.emplace_back([&, k] {
                for (;;) {
                    int t = nextTask.fetch_add(1);
                    if (t >= (int)tasks.size()) break;
                    runRange3(tasks[t].first, tasks[t].second, accs[k]);
                }
            });
        }
        for (auto& th : pool) th.join();
        int ffPierce = 0, ffSelf = 0, ffClr = 0, ffJoint = 0, ffIk = 0;
        int rescueEvals = 0, rescueHit = 0;
        int rkEvals[5] = {0, 0, 0, 0, 0};
        int rkHit[5] = {0, 0, 0, 0, 0};
        double hitMarginWorst = 0.0;
        for (const auto& a : accs) {
            recomputed += a.recomputed;
            reverted += a.reverted;
            reachable += a.reachable;
            standoffMaxErr = std::max(standoffMaxErr, a.standoffMaxErr);
            ffPierce += a.ffPierce; ffSelf += a.ffSelf; ffClr += a.ffClr;
            ffJoint += a.ffJoint;   ffIk += a.ffIk;
            rescueEvals += a.rescueEvals; rescueHit += a.rescueHit;
            for (int k = 0; k < 5; ++k) { rkEvals[k] += a.rkEvals[k]; rkHit[k] += a.rkHit[k]; }
            if (a.hitMarginWorst < hitMarginWorst) hitMarginWorst = a.hitMarginWorst;
        }
        if (g_prof.on && rescueEvals > 0) {
            std::fprintf(stderr, "  [rs] 救援 eval %d 次 · 救回 %d 点（%.1f%%）· 白跑 %d 次"
                                 " · 被救回点救援前最深超限 %.1f°\n",
                         rescueEvals, rescueHit, 100.0 * rescueHit / rescueEvals,
                         rescueEvals - rescueHit, hitMarginWorst);
            // 失败类型 → 救回率。顺序与 rk 编码一致。
            static const char* const kNm[5] = {"IK失败", "穿透", "自碰撞", "关节超限", "仅余量不足"};
            std::fprintf(stderr, "  [rs] 失败类型 → 救回：");
            for (int k = 0; k < 5; ++k)
                std::fprintf(stderr, "%s%s %d 点/救回 %d", k ? " · " : "",
                             kNm[k], rkEvals[k], rkHit[k]);
            std::fprintf(stderr, "\n");
        }
        if (g_prof.on)
            std::fprintf(stderr, "  [fb] 面级姿态过不了硬约束 %d 点：穿透 %d · 自碰撞 %d · 间隙不足 %d · 关节 %d · IK %d\n",
                         ffPierce + ffSelf + ffClr + ffJoint + ffIk,
                         ffPierce, ffSelf, ffClr, ffJoint, ffIk);
    }

    // ---- 阶段 3.5：回退段统一姿态（把「阶梯」压成「平台」）
    //
    // 逐点回退是必要的（硬约束不通过就得偏），但它会在同一个面上造出阶梯：开口盒子实测
    // 1417 个回退点落在 81 段里，其中 18 段在十几个点内连着出现 30°→40°→55°→30° 三种姿态。
    // 视觉上仍是刺，节拍上仍是连续甩腕 —— 而要求是「同一个面喷涂过程中姿态尽量不要变」。
    //
    // 所以把相邻回退点（允许夹 ≤2 个非回退点、且不跨行）归成一段，段内只用一个姿态：
    // 候选 = 段内逐点用过的姿态 + 面级姿态；按「段内硬约束失败点数 → 偏离面级姿态的大小」
    // 择优。只在硬约束层面取舍：仍有失败点的那些点保留各自的最佳位姿，绝不为了整齐而降标准。
    int segUnified = 0;
    {
        const int kGap = 2;
        auto isRev = [&](int k) { return hasResult[k] && results[k].reverted && t.spray[k]; };
        int i = 0;
        while (i < n) {
            if (!isRev(i)) { ++i; continue; }
            const int a = i;
            const int ri = rowOf[a];
            int b = i;
            while (b + 1 < n) {
                int k = b + 1;
                while (k < n && !isRev(k)) ++k;
                if (k < n && k - b - 1 <= kGap && rowOf[k] == ri) b = k;
                else break;
            }
            i = b + 1;
            if (b - a < 1) continue;   // 单点段没有阶梯可压

            const double faceTilt = finalTilt[a];
            const double faceAz = azAt(a);
            std::vector<std::pair<double, double>> cands;
            auto addCand = [&](double t2, double z2) {
                for (const auto& c : cands)
                    if (std::fabs(c.first - t2) < 1e-9 && std::fabs(c.second - z2) < 1e-9) return;
                cands.emplace_back(t2, z2);
            };
            for (int k = a; k <= b; ++k)
                if (isRev(k)) addCand(results[k].tiltDeg, results[k].azDeg);
            addCand(faceTilt, faceAz);

            int bestCi = -1, bestFail = n + 1;
            double bestDev = 1e300;
            for (int ci = 0; ci < (int)cands.size(); ++ci) {
                int fail = 0;
                // 不能提前 break：bestFail 降到 0 之后再"只在 fail<bestFail 时评估"等于
                // 把后面所有候选都当成 0 失败（没评估），于是按"偏离最小"挑回了面级姿态 ——
                // 那恰恰是过不了硬约束才对退的那个。必须每个候选都完整评一遍。
                for (int k = a; k <= b; ++k) {
                    V3 nrm = t.normal(k);
                    V3 fpt = t.pos[k] + t.gun[k] * p.standoff;
                    V3 ax, up;
                    basisAt(nrm, ax, up);
                    if (norm(up) < 1e-6) { ax = V3(0, 0, 0); up = V3(0, 0, 0); }
                    V3 d2 = dirOf(nrm, up, ax, cands[ci].first, cands[ci].second);
                    bool checkLinks = (k % std::max(1, p.linkStride) == 0);
                    const double* seed = (k > a && hasResult[k - 1]) ? results[k - 1].theta : nullptr;
                    PR tmp;
                    if (!ev.eval(fpt, d2, nullptr, seed, cands[ci].first, cands[ci].second,
                                 checkLinks, tmp) ||
                        !hardOk(tmp, p))
                        ++fail;
                }
                double dev = std::fabs(cands[ci].first - faceTilt) +
                             0.5 * std::fabs(cands[ci].second - faceAz);
                if (fail < bestFail || (fail == bestFail && dev < bestDev)) {
                    bestCi = ci; bestFail = fail; bestDev = dev;
                }
            }
            // 段内允许极少数点仍过不了硬约束（那几点保留自己的位姿，见下面的逐点兜底）；
            // 容忍带与面组择优同口径：2% 段长 + 1 点。全段都必须过才算统一，会一段都统一不了 ——
            // 因为"阶梯"本身就是逐点需求单调变化造成的，最长那段 187 点里必然有几个点要更高倾角。
            const int band = std::max(1, (int)(0.02 * (b - a + 1) + 0.5));
            if (std::getenv("RS_DBG_SEG"))
                std::fprintf(stderr, "  [seg] %d..%d 长%d 候选%zu 最优失败%d 容忍%d\n",
                             a, b, b - a + 1, cands.size(), bestFail, band);
            if (bestCi < 0 || bestFail > band) continue;

            const double bTilt = cands[bestCi].first;
            const double bAz = cands[bestCi].second;
            int applied = 0;
            for (int k = a; k <= b; ++k) {
                V3 nrm = t.normal(k);
                V3 fpt = t.pos[k] + t.gun[k] * p.standoff;
                V3 ax, up;
                basisAt(nrm, ax, up);
                if (norm(up) < 1e-6) { ax = V3(0, 0, 0); up = V3(0, 0, 0); }
                V3 d2 = dirOf(nrm, up, ax, bTilt, bAz);
                bool checkLinks = (k % std::max(1, p.linkStride) == 0);
                const double* seed = (k > a && hasResult[k - 1]) ? results[k - 1].theta : nullptr;
                PR tmp;
                if (!ev.eval(fpt, d2, nullptr, seed, bTilt, bAz, checkLinks, tmp) ||
                    !hardOk(tmp, p)) {
                    continue;   // 这一点在本段姿态下过不了硬约束：原样保留，不动它
                }
                const bool moved = norm(tmp.pos - results[k].pos) > 1e-6;
                tmp.idx = k;
                tmp.sprayAngle = angleBetweenDeg(-tmp.gun, nrm);
                tmp.reverted = (std::fabs(bTilt - faceTilt) > 1e-6 ||
                                std::fabs(bAz - faceAz) > 1e-6);
                setViolBits(tmp, true);   // isRev 已过滤掉空移点，这里必然是喷涂点
                results[k] = tmp;
                if (moved) ++applied;
            }
            // 只要有一段真的改了位姿就计一段；段内个别点没跟上（见上面的 continue）不算失败，
            // 那些点仍各自满足硬约束 —— 这里绝不为整齐牺牲硬约束。
            if (applied > 0) ++segUnified;
        }

        // 段统一之后，回退点数与可达点数都得按最终位姿重数一遍，否则报的是"过程账"。
        reverted = 0;
        reachable = 0;
        for (int k = 0; k < n; ++k) {
            if (!hasResult[k]) continue;
            if (results[k].reverted) ++reverted;
            if (results[k].ikOk) ++reachable;
        }
        if (g_prof.on)
            std::fprintf(stderr, "  [seg] 回退段统一=%d 段 · 回退点=%d 可达=%d\n",
                         segUnified, reverted, reachable);
    }

    // 生效行数：只有"确实改变了个点位姿"的行才计入 azStanceRows。
    // 阶段 3 的硬约束回退会把低于 hardClearance / 自碰撞 / 关节超限的位姿逐点退回，
    // 于是会出现"改判 N 行、实际一点没变"（complex 上就是：改判 1 行、输出与旧口径逐字节相同）。
    // 报数字要报生效的。
    {
        std::vector<char> applied(rows.size(), 0);
        for (int i = 0; i < n; ++i) {
            if (rowOf[i] < 0) continue;
            int ri = rowOf[i];
            if (std::fabs(azRow[ri]) < 1e-9) continue;
            if (!results[i].reverted && std::fabs(results[i].azDeg - azRow[ri]) < 1e-6) applied[ri] = 1;
        }
        azStanceRows = 0;
        for (char c : applied) if (c) ++azStanceRows;
    }

    // ---- 阶段 3.6：逐点解析救 J5（RS_POSE_TILT_J5）—— **叠加层**
    //
    // 这一段的定位必须先说清楚，否则很容易又写回"另一套姿态方案"：
    //
    //   **基线（面级常量 / 逐点连续）负责"够得着、不撞、喷角合规"。J5 只在
    //   「θ5 贴限位」这一件事上做局部修补。修补不了的保持基线姿态、如实报 P0。**
    //
    // 为什么必须是叠加而不是并列：第一版把 J5 做成阶段 2 的头号分支，非危机点一律
    // tilt=0，等于用一次局部修补推翻整套基线。代价实测（9 算例 A/B）：
    //   · complex：基线全点 6.21° 倾角被抹成 0，阶段 3 靠回退硬撑，关节步长>5° 47→93；
    //   · longstripe：基线 40° 倾角是把枪伸进窄槽的避障刚需，抹成 0 后连杆撞件，
    //     P0 从 2 涨到 22（LINK_HIT 20 点）；
    //   · longwide：基线 2.81° 是抬出浅槽用的，抹 0 后 7 点新增超限。
    // 三个算例的退化源头完全一样：**把非危机点的基线姿态丢掉**。
    //
    // 所以现在的顺序是：基线跑完 → 阶段 3 已拿到每点真实执行位姿（results[]）→
    // 这里只在「θ5 超出目标」的点上，用解析解**替换该点的 (tilt, az)**，
    // 其余点一个字节都不动。
    //
    // 判据用 results[i].theta —— 这是**真正会下发给机器人**的关节角，不是种子、
    // 不是近似。第一版图省时间用 analyticSeeds 估前臂方向，实测与真值差中位 28.6°、
    // 最大 59°，直接导致 big 上 8005/9761 点被误判成危机点、关节步长>5° 从 3 涨到 407。
    //
    // 三个已知边界（踩过，别再踩）：
    //   ① 退化面（前臂近乎垂直于面，|a_T|<cos85°）倾斜对 θ5 无效，是结构死结 → 跳过；
    //   ② 危机点方位逐点独立算会跨 ±180 翻转、实机拧腕 → 段内取同一空间方向（圆均值）；
    //   ③ 替换后的姿态必须重新过 hardOk，过不了就退回基线（不许为了救 J5 制造撞件）。
    int j5Rescued = 0, j5Rejected = 0, j5Segs = 0;
    if (p.tiltJ5) {
        std::vector<double> j5T(n, 0.0), j5A(n, 0.0);
        std::vector<char> j5Hit(n, 0);
        // |θ5| 的目标：限位 120° 减去要求余量，再留 1° 给 IK 收敛误差
        const double j5Lim = std::min(std::fabs(kJoints[4].hi), std::fabs(kJoints[4].lo));
        const double j5Tgt = j5Lim - std::max(0.0, p.jointMarginMin) - 1.0;

        if (j5Tgt > 1.0) {
            // ---- 步骤 1：逐点判定 + 解析规划（只读 results[]，不写任何东西）
            for (int i = 0; i < n; ++i) {
                if (!t.spray[i] || !hasResult[i] || !results[i].ikOk) continue;
                // 基线下的真实 |θ5|
                const double j5Now = std::fabs(results[i].theta[4]);
                j5Base[i] = (float)j5Now;
                if (j5Now <= j5Tgt) continue;          // 基线已经够安全 → 不动

                const V3 aBase = forearmDirOf(results[i].theta);
                const V3 nBase = ev.toBaseV(t.normal(i));
                // 第 1 参 = 点云系法向（算方位要用点云系切平面基），第 3 参 = 基座系（解析倾角）
                const J5Plan pl = planTiltForJ5(t.normal(i), aBase, nBase, ev.baseR(), p);
                if (!pl.valid) { ++j5Rejected; continue; }
                j5T[i] = pl.tiltDeg;
                j5A[i] = pl.azDeg;
                j5Hit[i] = 1;
            }

            // ---- 步骤 2：危机点按连续段分组，段内统一倾斜方向（枪轴必须连续）
            //
            // 逐点独立用 U* 会让方位在相邻点间跳几十度（实测同一面上 74.5°→69.0°→111.2°→59.1°，
            // 还有跨 ±180 的翻转），落到实机上就是沿轨迹反复拧腕。但**倾角不能统一**：
            // 每个点需要的倾角是它自己解出来的，强行统一会把"刚好够"的点推到撞件。
            const int kGap = 2;                        // 允许夹 ≤2 个非危机点仍算同段
            int i = 0;
            while (i < n) {
                if (!j5Hit[i]) { ++i; continue; }
                const int a = i;
                int b = i;
                while (b + 1 < n) {
                    int k = b + 1;
                    while (k < n && !j5Hit[k]) ++k;
                    if (k < n && k - b - 1 <= kGap && rowOf[k] == rowOf[a]) b = k;
                    else break;
                }
                i = b + 1;

                // 段内方位用**空间方向**求圆均值，而不是角度平均
                // （179° 与 −179° 直接平均得 0°，方向完全反了）
                V3 usum(0, 0, 0);
                int cnt = 0;
                for (int k = a; k <= b; ++k) {
                    if (!j5Hit[k]) continue;
                    const V3 nrm = t.normal(k);
                    V3 axis, uUp;
                    tangentBasisRobust(nrm, p.upW, axis, uUp);
                    if (norm(uUp) < 1e-6) continue;
                    const double r = j5A[k] * kDeg;
                    usum = usum + uUp * std::cos(r) + axis * std::sin(r);
                    ++cnt;
                }
                if (cnt == 0 || norm(usum) < 1e-6) continue;
                usum = unit(usum);
                // 各点法向不同 → 同一个空间方向的 az 也不同，逐点反算回来
                for (int k = a; k <= b; ++k) {
                    if (!j5Hit[k]) continue;
                    const V3 nrm = t.normal(k);
                    V3 axis, uUp;
                    tangentBasisRobust(nrm, p.upW, axis, uUp);
                    if (norm(uUp) < 1e-6) continue;
                    V3 uu = usum - nrm * dot(usum, nrm);
                    if (norm(uu) < 1e-6) continue;
                    uu = unit(uu);
                    j5A[k] = std::atan2(dot(uu, axis), dot(uu, uUp)) * kRad;
                }
                ++j5Segs;
            }

            // ---- 步骤 3：逐点试执行，过 hardOk 才采纳
            //
            // 这是「不许把 J5 危机换成撞件」的唯一保障。注意必须用**本点前一点的
            // 已定结果**当 IK 种子（种子链），冷启会让本来连续的解跳到另一支。
            //
            // 采纳失败分因记录（RS_DBG_J5=1 打印），因为三类失败的处置完全不同：
            //   IK 失败   → 解析方向本身不可达，只能退回基线；
            //   硬约束    → 解析方向撞件/推别的关节超限，值得在附近搜一次方位（见下）；
            //   θ5 未改善 → 该点前臂方向不对，倾斜救不了，老实退回。
            int rejIk = 0, rejHard = 0, rejWorse = 0, rejSing = 0, azFb = 0;
            const bool dbgJ5 = std::getenv("RS_DBG_J5") != nullptr;

            // 尝试一个 (tilt, az) 并判定是否可采纳。
            // 返回 0=可采纳 1=IK失败 2=硬约束 3=θ5未改善 4=腕部奇异。可采纳时写入 out。
            auto tryOne = [&](int idx, double tilt, double az, PR& out) -> int {
                const V3 nrm = t.normal(idx);
                const V3 face = t.pos[idx] + t.gun[idx] * p.standoff;
                V3 axis, uUp;
                tangentBasisRobust(nrm, p.upW, axis, uUp);
                if (norm(uUp) < 1e-6) return 1;              // 退化基
                const double* seed = (idx > 0 && hasResult[idx - 1]) ? results[idx - 1].theta : nullptr;
                const bool checkLinks = (idx % std::max(1, p.linkStride) == 0);
                PR pr;
                if (!ev.eval(face, dirOf(nrm, uUp, axis, tilt, az), nullptr, seed, tilt, az,
                             checkLinks, pr))
                    return 1;                                // IK 失败
                if (!hardOk(pr, p)) return 2;                // 硬约束
                // ---- 腕部奇异防护（必须）
                //
                // |θ5| → 0 表示 5 轴伸直，这是 VS-087 的腕部奇异：J4 与 J6 轴共线，
                // 末端姿态对 θ4/θ6 的分配完全失去约束，实机在这一点附近会出现
                // 关节速度爆冲。hardOk **不查 cond**（cond 只报 P2 告警），
                // 所以一个 θ5=0 的奇异解能顺利通过 hardOk —— longwide 实测就踩到了：
                // i=788 从 θ5=−119.43° 被"救"成 θ5=0.00°、cond=0.0，margin 看着从
                // 0.57 涨到 1.29（假改善），代价是把一个贴限位的点位换成不可执行位姿。
                //
                // 判据取 |θ5| ≥ 6°：VS-087 腕部奇异的实际危险带在 θ5 小于 5~8°，
                // 取 6° 既能挡住奇异，又不会误伤真正需要小 θ5 的点（那些点根本
                // 不会是"J5 贴限位"的危机点，进不了这个分支）。
                constexpr double kJ5SingMin = 6.0;
                if (std::fabs(pr.theta[4]) < kJ5SingMin) return 4;   // 奇异
                // 容差 0.5°：目标 j5Tgt 是硬线，贴着它抖动（118.0→118.3）也算"没改善"，
                // 为 0.3° 去改姿态不划算，还会多一次姿态跳变。
                if (std::fabs(pr.theta[4]) >= std::fabs(results[idx].theta[4]) - 0.5) return 3;
                pr.idx = idx;
                pr.tiltDeg = tilt;
                pr.azDeg = az;
                pr.sprayAngle = angleBetweenDeg(-pr.gun, nrm);
                out = pr;
                return 0;
            };

            for (int i = 0; i < n; ++i) {
                if (!j5Hit[i]) continue;
                PR got;
                int rc = tryOne(i, j5T[i], j5A[i], got);
                if (rc == 0) {
                    got.reverted = std::fabs(j5T[i] - finalTilt[i]) > 1e-6 ||
                                   std::fabs(j5A[i] - azAt(i)) > 1e-6;
                    results[i] = got;
                    hasResult[i] = 1;
                    ++j5Rescued;
                    continue;
                }
                // ---- 解析方向过不去（硬约束）→ 在**附近**搜方位
                //
                // 根因：U* = −a_T 是**纯运动学**最优，它只知道"这样能让腕子伸直"，
                // 完全不知道几何。窄槽/深腔件上，"指向前臂反方向"往往正好是把枪管
                // 顶进槽壁的方向 —— longstripe 实测 27 个退回点里 clr=−20.9mm 的占一半，
                // 那是真撞，不是算错。
                //
                // 所以不立刻放弃：**倾角大致保持**（倾角决定 θ5），把方位扫一圈，
                // 找第一个既过 hardOk 又把 θ5 拉回来的。方位改了 θ5 照样能变小
                // （θ5 只关心枪轴与前臂的夹角，方位可以转），只是不如 U* 省倾角。
                // 这就是把运动学结论当**先验**、几何当**否决权**。
                if (rc == 2 || rc == 4) {
                    const double azTries[9] = {0, 40, -40, 75, -75, 110, -110, 145, -145};
                    bool hasAlt = false;
                    double bestJ5 = 1e9;
                    for (int t = 0; t < 9 && !hasAlt; ++t) {
                        double az2 = j5A[i] + azTries[t];
                        while (az2 <= -180.0) az2 += 360.0;
                        while (az2 > 180.0) az2 -= 360.0;
                        // 方位一转，同一方位上需要的倾角也会变，给三档
                        const double tilts[3] = {j5T[i], j5T[i] + 6.0,
                                                 std::max(2.0, j5T[i] - 6.0)};
                        for (int q = 0; q < 3; ++q) {
                            if (tilts[q] > p.tiltMaxDeg + 1e-9) continue;
                            if (tilts[q] > p.maxSprayAngleDeg + 1e-9) continue;
                            PR alt;
                            if (tryOne(i, tilts[q], az2, alt) != 0) continue;
                            const double j5v = std::fabs(alt.theta[4]);
                            if (j5v < bestJ5) { bestJ5 = j5v; got = alt; hasAlt = true; }
                        }
                    }
                    if (hasAlt) {
                        got.reverted = true;
                        results[i] = got;
                        hasResult[i] = 1;
                        ++j5Rescued;
                        ++azFb;
                        continue;
                    }
                }
                ++j5Rejected;
                if (rc == 1) { ++rejIk; }
                else if (rc == 2) { ++rejHard; }
                else if (rc == 4) { ++rejSing; }
                else { ++rejWorse; }
                if (dbgJ5) std::fprintf(stderr, "  [j5rej] i=%d rc=%d tilt=%.2f az=%.1f\n",
                                        i, rc, j5T[i], j5A[i]);
            }
            if (g_prof.on)
                std::fprintf(stderr,
                             "  [j5] 退回原因：IK失败 %d · 硬约束 %d · θ5未改善 %d · 腕部奇异 %d"
                             " · 方位搜索救回 %d\n",
                             rejIk, rejHard, rejWorse, rejSing, azFb);
            if (g_prof.on)
                std::fprintf(stderr, "  [j5] 退回原因：IK失败 %d · 硬约束 %d · θ5未改善 %d\n",
                             rejIk, rejHard, rejWorse);
        }
        if (g_prof.on) {
            int over0 = 0, over1 = 0;
            for (int i = 0; i < n; ++i) {
                if (!t.spray[i] || !hasResult[i]) continue;
                if (j5Base[i] > j5Tgt) ++over0;
                if (std::fabs(results[i].theta[4]) > j5Tgt) ++over1;
            }
            std::fprintf(stderr,
                         "  [j5] 目标|θ5|<=%.1f° · 危机点 %d · 采纳 %d · 退回基线 %d · 归 %d 段 · "
                         "仍超目标 %d→%d\n",
                         j5Tgt, over0, j5Rescued, j5Rejected, j5Segs, over0, over1);
        }
    }

    // ---- 过渡点抬升（spray=0 ⇒ 位置可改，只求不撞）
    g_prof.mark("阶段3 重解");
    int lifted = 0;
    if (p.fixTransit) {
        // 逐点独立抬升会把空移线抬成锯齿：开口盒子实测 A1 段 139 个点的 lift 在
        // 45~56mm 之间逐点跳，z 每步变 2~4mm，整段像波浪；相邻点 lift 差最大 57.5mm
        // （一步之内垂直跳 57mm）。空移是关枪过渡，同一段内应该用**同一个抬升量**：
        // 路径仍是原来那条折线，只是整体抬高，机器人不会上下抖。
        // 做法：先算出每点需要的 lift，再把"连续被抬的点串"统一到串内最大值。
        std::vector<double> need(n, 0.0);
        for (int i = 0; i < n; ++i) {
            PR& pr = results[i];
            if (!hasResult[i] || t.spray[i] != 0) continue;
            if (!std::isfinite(pr.clearance) || pr.clearance >= p.safety) continue;
            need[i] = std::min(p.transitLiftMax, p.safety - pr.clearance + 5.0);
        }
        for (int i = 0; i < n;) {
            if (need[i] <= 0.0) { ++i; continue; }
            int j = i;
            double mx = need[i];
            while (j + 1 < n && need[j + 1] > 0.0) { ++j; mx = std::max(mx, need[j]); }
            for (int k = i; k <= j; ++k) need[k] = mx;
            i = j + 1;
        }
        for (int i = 0; i < n; ++i) {
            if (need[i] <= 0.0) continue;
            PR& pr = results[i];
            // 抬升沿重力反向（远离地面/工件），不是沿点云 +Z —— 相机系里 +Z 是光轴，
            // 沿它抬等于把空移点横着推出工件，间隙数字变好了但实际没抬离。
            pr.pos = pr.pos + p.upW * need[i];
            pr.lifted = need[i];
            pr.clearance += need[i];
            // 分组距离必须同步抬升，否则总结论（已避开）与"连杆最小距离"（还是抬升前）
            // 会自相矛盾 —— curved 上就出现过 clearance=20.02 而 clrTool=11.39。
            if (pr.clrTool < 1e299) pr.clrTool += need[i];
            if (pr.clrLink < 1e299) pr.clrLink += need[i];
            pr.viol &= ~(RS_VIOL_COLLISION | RS_VIOL_CLEARANCE | RS_VIOL_LINK_HIT | RS_VIOL_TOOL_HIT);
            ++lifted;
        }
    }

    // ---- 姿态跳变 / 位置跳变：都是"平滑性"的代价，必须报出来（P2，速度/加速度风险）
    //
    // 两条都**不再限制在行内**。原来只在 `r.first+1..r.second` 里比，恰好漏掉了最容易
    // 出问题的那一对：某条喷涂线的最后一点 → 下一条空移线的第一点。垃圾桶实测
    // #309(喷涂末点, tilt 43.15°) → #310(空移首点, tilt 0)，raw 相距 4.00mm，
    // 修正后 74.90mm，全程 VERDICT=PASS —— 因为这一对跨行，根本没进循环。
    // 步长阈值要跟"本件自己的点距"比，不能写死：row80 稀疏件点距本身就 40mm+，
    // 按固定 12mm 报会把整条轨迹刷满假告警；反过来点距 3mm 的密件，
    // 12mm 已经是一次 4 倍的位移突变。
    //
    // 下限取 15 而不是 12：阶段 2.5 的限速本身就允许「每点 Δtilt ≤ dTiltCap」，
    // 折成位移正好是 12mm 量级，再加上本点自己的 Δpos 就是 13~18mm。
    // 这些是**设计内的斜坡步长**，不是跳点；把它们一起报出来，只会让
    // "修完之后告警更多了"这种误判发生（实测 12 阈值下 38 段，其中 30 段是斜坡）。
    // thrStep 放在 if 外面：下面"步长硬上限插点"用同一个阈值——报出来之前
    // 先 try 一次插点消掉，消不掉才报，两处口径必须一致。
    std::vector<double> gaps;
    gaps.reserve(n);
    for (int i = 1; i < n; ++i)
        if (hasResult[i] && hasResult[i - 1])
            gaps.push_back(norm(results[i].pos - results[i - 1].pos));
    double medGap = 0.0;
    if (!gaps.empty()) {
        std::vector<double> g2 = gaps;
        std::nth_element(g2.begin(), g2.begin() + g2.size() / 2, g2.end());
        medGap = g2[g2.size() / 2];
    }
    const double thrStep = std::max(15.0, 3.0 * medGap);

    if (p.tiltRateDeg > 0) {
        const double thrTilt = std::max(3.0 * p.tiltRateDeg, 10.0);
        for (int i = 1; i < n; ++i) {
            // 没有位姿的点（IK 失败）谈跳变没有意义，别让它再刷一条违规。
            // 注意判 hasResult 而不是 ikOk：不要求可达时 ikOk 恒为 0，但位姿是有的，
            // 跳变（换面/回退/换段造成）恰恰是那时候最该报的东西。
            if (!hasResult[i] || !hasResult[i - 1]) continue;
            // 姿态：用枪轴 3D 夹角。只比 tiltDeg 会漏掉整类问题 —— 逐点口径下
            // 方位翻转时倾角可以一个字都不变（垃圾桶桶底小坑那段全程 42.1~42.2°，
            // 方位绕满一圈，相邻枪轴夹角却有 18~58°）。
            if (angleBetweenDeg(results[i].gun, results[i - 1].gun) > thrTilt)
                results[i].viol |= RS_VIOL_TILT_JUMP;
            // 位置：直线步长。换段处不抬枪、不补偿，是这条判据的主战场。
            if (norm(results[i].pos - results[i - 1].pos) > thrStep)
                results[i].viol |= RS_VIOL_STEP_JUMP;
        }
    }

    // ---- L0 输入层自检（这一层不通，后面所有几何结论都是废的）
    int l0Issues = 0;
    if (p.requireReach) {
        // 只有要求可达时，"法兰→面点"的总长才必须落在臂展内
        double reachNeed = p.standoff + p.toolLen;
        double limit = kReach * 0.90;
        if (reachNeed > limit) ++l0Issues;
    }
    double standoffMedian = -1.0;
    {
        std::vector<int> idxSpray;
        for (int i = 0; i < n; ++i)
            if (t.spray[i]) idxSpray.push_back(i);
        int pick = std::min(60, (int)idxSpray.size());
        double tube = std::max(3.0 * spacing, 5.0);
        std::vector<double> ds;
        for (int k = 0; k < pick; ++k) {
            int i = idxSpray[(int)((long long)k * (idxSpray.size() - 1) / std::max(1, pick - 1))];
            double v;
            if (rayStandoff(cloud, t.pos[i], t.gun[i], tube, v)) ds.push_back(v);
        }
        if (!ds.empty()) {
            std::sort(ds.begin(), ds.end());
            standoffMedian = ds[ds.size() / 2];
            if (std::fabs(standoffMedian - p.standoff) > std::max(8.0, 0.15 * p.standoff)) ++l0Issues;
        }
    }

    // ---- 步长硬上限：相邻点位移超预算的段落插补过渡点（"千万不能跳点"）
    //
    // 倾斜修正会把「枪轴绕小特征扫一圈」放大成 TCP 的大回摆：垃圾桶桶底 Ø9mm 小坑，
    // raw 枪轴 7 点扫满 360°、TCP 每点只走 3mm；加 42° 倾角后，73.5mm 的倾斜矢量
    // 跟着枪轴一起绕圈，TCP 被甩到 Ø110mm 圆上，每点 22~55mm。
    // 倾角调度救不了它——喷这个坑就必然要扫这一圈（枪轴极角恒定、方位转满 360° 是
    // raw 的真实几何），倾斜矢量从 71.7mm 收到 12mm 要径向走 60mm，快速自旋区里
    // 这 60mm 无论怎么摊都至少 30mm/步。
    //
    // 所以唯一可靠的办法：**不改轨迹形状，把大步切成小步**——位置线性插值、
    // 枪轴球面插值（slerp），机器人沿真实回摆弧走，而不是一步甩过去。
    // 工艺（倾角/喷角/开关枪/速度）一格不改，插点 spray 只在两端都开枪时才开。
    //
    // 每个插补点过一遍枪体闸门（阶段 2.5 同款轻量评估器）：弦比弧贴轴，
    // 弦太长可能啃到特征本体。不过线就加倍细分——分得足够密时弦→弧、间隙→端点间隙，
    // 弧本身是端点之间的既有路径；6 次加倍仍不过就放弃插点，保留 STEP_JUMP 如实上报。
    // 插成功才撤掉这一段的 STEP_JUMP / TILT_JUMP（跳点已经不存在了）。
    std::vector<PR> emit;
    std::vector<int> emitSrc;
    std::vector<char> emitSpray;
    emit.reserve((size_t)n + 64);
    emitSrc.reserve((size_t)n + 64);
    emitSpray.reserve((size_t)n + 64);
    int insertedPts = 0, densifiedSegs = 0, densifyGiveUp = 0;
    {
        const double SwD = p.standoff > 1e-6 ? p.standoff : 100.0;
        auto toolClrOf = [&](const V3& pos, const V3& gunDir) {
            const V3 g = unit(gunDir);
            PR pr;
            evTool.eval(pos + g * SwD, g * -1.0, nullptr, nullptr, 0.0, 0.0, false, pr);
            return pr.clearance;   // requireReach=0 口径 = 枪体最小表面距离
        };
        for (int i = 0; i < n; ++i) {
            if (i > 0 && hasResult[i] && hasResult[i - 1]) {
                const PR& A = results[i - 1];
                const PR& B = results[i];
                const double step = norm(B.pos - A.pos);
                if (step > thrStep) {
                    const V3 a = unit(A.gun), b = unit(B.gun);
                    const double omega =
                        std::acos(std::max(-1.0, std::min(1.0, dot(a, b))));
                    auto slerp = [&](double s) -> V3 {
                        if (omega < 1e-9) return unit(a + (b - a) * s);
                        const double w1 = std::sin((1.0 - s) * omega) / std::sin(omega);
                        const double w2 = std::sin(s * omega) / std::sin(omega);
                        return unit(a * w1 + b * w2);
                    };
                    const double ca = toolClrOf(A.pos, A.gun);
                    const double cb = toolClrOf(B.pos, B.gun);
                    // 任一端没量出枪体间隙（不该发生：阶段 3 全程建了碰撞体）就不设闸
                    const bool verify = std::isfinite(ca) && std::isfinite(cb);
                    const double floorClr = std::min(ca, cb) - 2.0;
                    int k = (int)std::ceil(step / thrStep) - 1;
                    if (k < 1) k = 1;
                    std::vector<PR> ins;
                    bool ok = false;
                    while (true) {
                        ins.clear();
                        ok = true;
                        for (int j = 1; j <= k && ok; ++j) {
                            const double s = (double)j / (k + 1);
                            PR q;
                            q.pos = A.pos + (B.pos - A.pos) * s;
                            q.gun = slerp(s);
                            q.sprayAngle = A.sprayAngle + (B.sprayAngle - A.sprayAngle) * s;
                            q.tiltDeg = q.sprayAngle;
                            q.azDeg = B.azDeg;
                            if (verify) {
                                const double cc = toolClrOf(q.pos, q.gun);
                                if (!(cc >= floorClr)) { ok = false; break; }
                                q.clearance = cc;
                                q.clrTool = cc;
                            } else {
                                q.clearance = std::min(A.clearance, B.clearance);
                                q.clrTool = std::min(A.clrTool, B.clrTool);
                            }
                            // 其余报告字段取两端较保守值：插补点是几何中间态，
                            // 不是独立评估过的位姿，别伪装成比端点更安全。
                            q.clrLink = std::min(A.clrLink, B.clrLink);
                            q.selfGap = std::min(A.selfGap, B.selfGap);
                            q.margin = std::min(A.margin, B.margin);
                            q.ikOk = A.ikOk && B.ikOk;
                            q.standoffErr = 0.0;   // pos/gun 是端点凸组合，偏差有界不另测
                            q.viol = 0;
                            ins.push_back(q);
                        }
                        if (ok || k >= 96) break;
                        k *= 2;
                    }
                    if (ok) {
                        const char spr = (t.spray[i - 1] && t.spray[i]) ? 1 : 0;
                        for (const PR& q : ins) {
                            emit.push_back(q);
                            emitSrc.push_back(i);
                            emitSpray.push_back(spr);
                        }
                        // 跳点已经切小，相应的 P2 告警随之撤掉
                        results[i].viol &= ~(unsigned)(RS_VIOL_STEP_JUMP | RS_VIOL_TILT_JUMP);
                        ++densifiedSegs;
                        insertedPts += (int)ins.size();
                    } else {
                        ++densifyGiveUp;
                    }
                }
            }
            emit.push_back(results[i]);
            emitSrc.push_back(i);
            emitSpray.push_back(t.spray[i] ? 1 : 0);
        }
    }
    const int nOut = (int)emit.size();
    if (g_prof.on && (densifiedSegs || densifyGiveUp)) {
        std::fprintf(stderr,
                     "  [dens] 超步长(>%.1fmm)切段=%d 插点=%d 放弃=%d 段 · 输出 %d 点（原 %d）\n",
                     thrStep, densifiedSegs, insertedPts, densifyGiveUp, nOut, n);
    }

    // ---- 输出
    g_prof.mark("过渡+自检");
    RobotSafetyPoint* out = new RobotSafetyPoint[nOut];
    std::memset(out, 0, sizeof(RobotSafetyPoint) * (size_t)nOut);
    std::vector<double> angles;
    angles.reserve(n);
    double minClr = std::numeric_limits<double>::infinity();
    int tilted = 0;
    double minClrTool = std::numeric_limits<double>::infinity();
    double minClrLink = std::numeric_limits<double>::infinity();
    int reachableOut = 0;   // 按输出的点（含插补点）重算，见循环内
    // 未测到（IK 失败 / 本点没查连杆）统一写 nan：写 0 会被读成"贴脸 0mm"，
    // 那比"没测"更误导。
    auto measured = [](double v) {
        return (std::isfinite(v) && v < 1e299) ? v : std::numeric_limits<double>::quiet_NaN();
    };
    vc.total = 0;
    std::memset(vc.byRule, 0, sizeof(vc.byRule));

    for (int r = 0; r < nOut; ++r) {
        const PR& pr = emit[r];
        const int src = emitSrc[r];
        RobotSafetyPoint& q = out[r];
        q.x = pr.pos.x;
        q.y = pr.pos.y;
        q.z = pr.pos.z;
        q.i = pr.gun.x;
        q.j = pr.gun.y;
        q.k = pr.gun.z;
        q.speed = t.speed[src];
        q.spray = emitSpray[r];
        q.tiltDeg = pr.sprayAngle;
        q.clearanceMm = pr.clearance;
        q.selfGapMm = pr.selfGap;
        q.marginDeg = pr.margin;
        q.standoffErrMm = pr.standoffErr;
        q.liftedMm = pr.lifted;
        q.cond = pr.cond;
        // theta 只在真解过 IK 时才有意义；其余情况写 0（不是 0°，是"没解"）。
        for (int k = 0; k < 6; ++k) q.theta[k] = pr.ikOk ? pr.theta[k] : 0.0;
        q.ikOk = pr.ikOk ? 1 : 0;
        q.reverted = pr.reverted ? 1 : 0;
        q.violation = pr.viol;
        q.clearanceToolMm = measured(pr.clrTool);
        q.clearanceLinkMm = measured(pr.clrLink);
        // 最近部位按"有没有这个部位"决定，不再挂在 ikOk 上：
        // 不要求可达时没解 IK，但枪体的最近部位（握把/细管/喷嘴）是真量出来的，
        // 报告里"枪体撞在喷嘴"这句话必须留着。
        q.worstKind = pr.partKind;
        q.worstJ0 = pr.partJ0;
        q.worstJ1 = pr.partJ1;
        q.linkKind = pr.linkKind;
        q.linkJ0 = pr.linkJ0;
        q.linkJ1 = pr.linkJ1;
        if (pr.viol) {
            for (int b = 0; b < RS_VIOL_COUNT; ++b)
                if (pr.viol & (1u << b)) ++vc.byRule[b];
            ++vc.total;
        }
        angles.push_back(pr.sprayAngle);
        // 只统计**有限**的间隙：IK 失败的点 clearance 是 -inf，混进来会让
        // "最小表面距离" 变成一个没有意义的 -inf（真正要看的是不可达率）。
        if (std::isfinite(pr.clearance) && pr.clearance < minClr) minClr = pr.clearance;
        double ct = measured(pr.clrTool), cl = measured(pr.clrLink);
        if (ct < minClrTool) minClrTool = ct;
        if (cl < minClrLink) minClrLink = cl;
        if (pr.sprayAngle > 0.5) ++tilted;
        // 可达率按**输出的点**（含插补点）重算：插补点 ikOk = 两端与，
        // 分子分母才是同一个口径。
        if (!p.requireReach || pr.ikOk) ++reachableOut;
    }

    double sum = 0.0;
    for (double a : angles) sum += a;

    // 不要求可达时根本没有解过 IK（ikOk 恒 0）。按用户口径「默认全部通过」记成全部点
    // （循环里已逐点计入 reachableOut），但这只是记账口径：CSV 里 ikOk 仍是 0、
    // theta 仍是 0，界面/报告上必须显著标成「未检查」，不能把 100% 当成真测出来的可达率。
    // minClearanceLinkMm 同理：没建连杆碰撞体，它会是 nan，不是 0。
    R.points = out;
    R.count = nOut;
    R.reachable = reachableOut;
    R.tiltedPoints = tilted;
    R.poseReverted = reverted;
    R.violations = vc.total;
    for (int b = 0; b < RS_VIOL_COUNT; ++b) R.byRule[b] = vc.byRule[b];
    R.minClearanceMm = minClr;
    R.minClearanceToolMm = minClrTool;
    R.minClearanceLinkMm = minClrLink;
    R.tiltMeanDeg = nOut ? sum / nOut : 0.0;
    R.tiltP95Deg = percentile(angles, 0.95);
    R.tiltMaxDeg = angles.empty() ? 0.0 : *std::max_element(angles.begin(), angles.end());
    R.standoffMaxErrMm = standoffMaxErr;
    R.baseAuto = baseAuto ? 1 : 0;
    R.cloudPoints = cloud.size();
    R.cloudSpacingMm = spacing;
    R.azRows = (int)rows.size();
    R.azStanceRows = azStanceRows;
    R.l0Issues = l0Issues;
    R.standoffMedianMm = standoffMedian;

    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) R.basePose[r * 4 + c] = (r < 3 && c < 3) ? base.R.m[r][c] : (r == c ? 1.0 : 0.0);
    }
    for (int r = 0; r < 3; ++r) R.basePose[r * 4 + 3] = base.t[r];

    (void)hasPrevSol;
    (void)lifted;
    if (g_prof.on) {
        g_prof.mark("输出+组装");
        double total = std::chrono::duration<double, std::milli>(
                           std::chrono::steady_clock::now() - g_prof.t0).count();
        std::fprintf(stderr,
                     "  [prof] 总计 %.1f ms / %d 点 = %.3f ms/点 · IK %lld · eval %lld（%.1f/点）"
                     " · 段查询 %lld（%.1f/点）· 点云 %d 点\n",
                     total, n, n ? total / n : 0.0, g_prof.ik.load(), g_prof.evals.load(),
                     n ? (double)g_prof.evals.load() / n : 0.0, g_prof.seg.load(),
                     n ? (double)g_prof.seg / n : 0.0, cloud.size());
        // eval 内部分项（多线程累加，按 CPU 时间算，若 > 单墙钟是并行导致的正常现象）
        double evals = (double)std::max(1LL, (long long)g_prof.evals.load());
        std::fprintf(stderr,
                     "  [prof] eval 内部分项 CPU 占比（%/总）：IK %.0f%% · 建碰撞体 %.0f%% · "
                     "点云查询 %.0f%% · 自碰撞 %.0f%%（每次 eval：IK %.1f μs / 建 %.2f μs / "
                     "查询 %.1f μs / 自 %.2f μs）\n",
                     100.0 * (double)g_prof.nsIk / 1e6 / (total * 1.0) * 1.0,
                     100.0 * (double)g_prof.nsBuild / 1e6 / (total * 1.0),
                     100.0 * (double)g_prof.nsSeg / 1e6 / (total * 1.0),
                     100.0 * (double)g_prof.nsSelf / 1e6 / (total * 1.0),
                     (double)g_prof.nsIk / 1e3 / evals, (double)g_prof.nsBuild / 1e3 / evals,
                     (double)g_prof.nsSeg / 1e3 / evals, (double)g_prof.nsSelf / 1e3 / evals);
        double ikN = (double)std::max(1LL, (long long)g_prof.ik.load());
        std::fprintf(stderr,
                     "  [prof] ikOnce 总 %lld 次（每次 ik() 平均 %.1f 次）· 种子直接命中 %lld（%.0f%%）"
                     " · 种子几何收敛但超限 %lld（%.0f%%）· 种子失败 %lld（%.0f%%）\n",
                     (long long)g_prof.ikOnceN, (double)g_prof.ikOnceN / ikN,
                     (long long)g_prof.ikFast, 100.0 * (double)g_prof.ikFast / ikN,
                     (long long)g_prof.ikOver, 100.0 * (double)g_prof.ikOver / ikN,
                     (long long)g_prof.ikFail, 100.0 * (double)g_prof.ikFail / ikN);
        if (g_prof.ikOver > 0) {
            double ov = (double)g_prof.ikOver;            std::fprintf(stderr,
                         "  [prof] 超限分支：总 ikOnce %lld（平均 %.1f 次/点）· 换来合法解 %lld（%.0f%%）"
                         " · 白跑 %lld（%.0f%%）\n",
                         (long long)g_prof.ikOverOnce, (double)g_prof.ikOverOnce / ov,
                         (long long)g_prof.ikOverHit, 100.0 * (double)g_prof.ikOverHit / ov,
                         (long long)g_prof.ikOverMiss, 100.0 * (double)g_prof.ikOverMiss / ov);
            // 超限深度直方图：每一档里"换了种子之后有没有救回来"
            std::fprintf(stderr, "  [prof] 超限深度 → 救回率：");
            for (int i = 0; i < Prof::kOverBucketN; ++i) {
                long long c = g_prof.overCnt[i].load();
                long long h = g_prof.overHitB[i].load();
                std::fprintf(stderr, "%s%s %lld 次/救回 %lld",
                             i ? " · " : "", g_prof.kOverNames[i], c, h);
            }
            std::fprintf(stderr, "\n");
        }
        if (g_prof.sweepSeeds > 0)
            std::fprintf(stderr, "  [prof] 全种子搜 %lld 组 · 其中与前组逐位重复 %lld（%.1f%%）\n",
                         (long long)g_prof.sweepSeeds, (long long)g_prof.sweepDup,
                         100.0 * (double)g_prof.sweepDup / (double)g_prof.sweepSeeds);
        if (g_prof.ikOnceN > 0)
            std::fprintf(stderr, "  [prof] ikOnce 平均迭代 %.1f 次（上限 90）\n",
                         (double)g_prof.ikIters / (double)g_prof.ikOnceN);
    }
    return R;
}
