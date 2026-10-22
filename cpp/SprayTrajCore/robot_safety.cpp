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
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point last = t0;

    void reset() {
        ik = 0; seg = 0; evals = 0;
    }
    void mark(const char* what) {
        if (!on) return;
        auto now = std::chrono::steady_clock::now();
        double dt = std::chrono::duration<double, std::milli>(now - last).count();
        double all = std::chrono::duration<double, std::milli>(now - t0).count();
        std::fprintf(stderr, "  [prof] %-14s %9.1f ms（累计 %9.1f）\n", what, dt, all);
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

    for (int it = 0; it < 90; ++it) {
        FkOut f;
        fk(th, toolLen, f);
        V3 ep = pTarget - f.tcp.t;
        V3 er = matToRotvec(m3mul(RTarget, m3T(f.tcp.R)));
        double err[6] = {ep.x, ep.y, ep.z, wRot * er.x, wRot * er.y, wRot * er.z};
        double en = 0.0;
        for (int i = 0; i < 6; ++i) en += err[i] * err[i];
        en = std::sqrt(en);
        if (en < bestErr) {
            bestErr = en;
            lam = std::max(lam * 0.7, 0.5);
        } else {
            lam = std::min(lam * 0.5 + lam, 5e3);
        }
        if (norm(ep) < 0.02 && norm(er) * kRad < 0.02) break;

        double J[6][6];
        jacobian(f, J);
        for (int i = 3; i < 6; ++i)
            for (int j = 0; j < 6; ++j) J[i][j] *= wRot;
        double JJt[6][6];
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

    FkOut f;
    fk(th, toolLen, f);
    IkSol sol;
    for (int i = 0; i < 6; ++i) sol.theta[i] = th[i];
    sol.posErr = norm(pTarget - f.tcp.t);
    sol.rotErrDeg = norm(matToRotvec(m3mul(RTarget, m3T(f.tcp.R)))) * kRad;
    sol.ok = (sol.posErr <= std::max(0.02 * 20, 0.5) && sol.rotErrDeg <= std::max(0.02 * 20, 0.5));
    sol.over = overLimit(th);
    sol.margin = limitMarginMin(th);
    sol.cond = condNumber(f);
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
    ++g_prof.ik;
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

    // 先跑上一点种子：连续且合法就直接返回，保持主路径速度。
    if (seed) {
        IkSol sol = ikOnce(pTarget, RTarget, toolLen, seed);
        best = sol;
        hasBest = true;
        if (valid(sol)) {
            if (maxJointDelta(sol.theta, seed) <= kBranchContinuityTrigDeg) return sol;
            // 跳得大：先把它作为“最近的合法解”的兜底，再找有没有更近的解支。
            bestValid = sol;
            hasValid = true;
            bestValidDist = jointDist2(sol.theta, seed);
        }
    }

    // 解析种子 / 固定种子。有 seed 且上一点解跳得大时，选关节距离最近的合法解；
    // 没有 seed 时保持旧行为（第一个合法解直接返回）。
    for (int i = (seed ? 1 : 0); i < ns; ++i) {
        IkSol sol = ikOnce(pTarget, RTarget, toolLen, seeds[i]);
        if (valid(sol)) {
            if (!seed) return sol;
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
    if (hasValid) return bestValid;
    return best;
}

// ---------------------------------------------------------------- 碰撞体（胶囊）

struct Sizing {
    double jointR, joint1R, linkR, wristLinkR, baseColumnR, baseNeckR;
    double gripR, toolPipeR, nozzleR;
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
    s.nozzleR = s.gripR * 0.30;
    return s;
}

constexpr double kGripFraction = 0.72;    // 握把占「法兰→TCP」全长比例
constexpr double kNozzleFraction = 0.18;  // 喷嘴伸出 TCP 的比例

enum SegGroup { SG_TOOL = 0, SG_LINK = 1, SG_JOINT = 2, SG_BASE = 3 };

struct Seg {
    V3 a, b;
    double r;
    int group;
    int j0, j1;    // 关联关节号（自碰撞排除用；-1 表示无）
    const char* part;
    int kind;      // RS_SEG_*，用于把「最近部位」报成人能看懂的名字
};

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
    V3 flange = f.flange.t, tcp = f.tcp.t;
    V3 d = tcp - flange;
    double L = norm(d);
    if (L > 0.5) {
        V3 u = d * (1.0 / L);
        V3 neck = flange + u * (L * kGripFraction);
        V3 tip = tcp + u * (L * kNozzleFraction);
        out.push_back({flange, neck, sz.gripR, SG_TOOL, -1, -1, "工具握把", RS_SEG_GRIP});
        out.push_back({neck, tcp, sz.toolPipeR, SG_TOOL, -1, -1, "工具细管", RS_SEG_PIPE});
        out.push_back({tcp, tip, sz.nozzleR, SG_TOOL, -1, -1, "喷嘴", RS_SEG_NOZZLE});
    }
}

// 部位名（含关节号），给控制台/CSV/界面用
std::string segName(int kind, int j0, int j1) {
    char buf[48];
    switch (kind) {
        case RS_SEG_GRIP:      return "工具握把";
        case RS_SEG_PIPE:      return "工具细管";
        case RS_SEG_NOZZLE:    return "喷嘴";
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
        ++g_prof.seg;
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
// 方位角只能取一个（22 笔朝向不同的壁共用一个进刀侧）、倾角"只升不降"的平滑会跨折返
// 传递（上一笔的高倾角把下一笔整段抬起来）。
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
    double tiltNonUpDeg = 30.0;      // α=90°（竖直侧壁）时的倾角，实际按 α/90 线性缩放
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
inline int poseCmp(const PR& a, const PR& b) {
    int oa = (a.selfGap >= 0.0) ? 1 : 0, ob = (b.selfGap >= 0.0) ? 1 : 0;
    if (oa != ob) return oa < ob ? -1 : 1;
    if (a.clearance != b.clearance) return a.clearance < b.clearance ? -1 : 1;
    if (a.margin != b.margin) return a.margin < b.margin ? -1 : 1;
    return 0;
}

// 硬约束：IK 成功、不自碰撞、间隙不低于绝对安全距离、关节不超限。
// 注意 margin < jointMarginMin 只是 P1 软告警，不在这里判死；over（真正超限）才是硬约束。
inline bool hardOk(const PR& q, double hardClearance) {
    return q.ikOk
        && q.selfGap >= 0.0
        && q.clearance >= hardClearance
        && !q.over;
}

// 两个“都不满足硬约束”的位姿二选一：先比硬约束通过项，再比自间隙/间隙/关节余量。
// 与 poseCmp 的区别：selfGap<0 时还要比谁的负值更小；关节真正超限（over）
// 也单独判负，不能因为另一个位姿间隙大就选一个关节超限的位姿。
inline int hardCmp(const PR& a, const PR& b) {
    int ia = a.ikOk ? 1 : 0, ib = b.ikOk ? 1 : 0;
    if (ia != ib) return ia < ib ? -1 : 1;
    int sa = (a.selfGap >= 0.0) ? 1 : 0, sb = (b.selfGap >= 0.0) ? 1 : 0;
    if (sa != sb) return sa < sb ? -1 : 1;
    int oa = a.over ? 0 : 1, ob = b.over ? 0 : 1;
    if (oa != ob) return oa < ob ? -1 : 1;
    if (a.selfGap != b.selfGap) return a.selfGap < b.selfGap ? -1 : 1;
    if (a.clearance != b.clearance) return a.clearance < b.clearance ? -1 : 1;
    if (a.margin != b.margin) return a.margin < b.margin ? -1 : 1;
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
    M3 toBaseR(const M3& R) const { return m3mul(baseRT_, R); }
    V3 toCloudP(const V3& p) const { return m3mul(baseR_, p) + baseT_; }
    // 方向向量：基座系 → 点云系（只有旋转，不平移）。refX 这类「来自 FK 的轴」必须过这一步。
    V3 toCloudR(const V3& v) const { return m3mul(baseR_, v); }
    // 基座原点（点云系）。方位决策要用它判"机器人从哪一侧够过来"。
    V3 baseOrigin() const { return baseT_; }

    // 给定枪身朝外方向 d（倾斜后的法向 n'），求 IK + 碰撞 + 自碰撞。
    // **工具 +Z 是 −d**：d 是从工件指向枪的外法向，工具 +Z 是喷射方向（法兰→喷嘴）= −d。
    // 这条搞反过一整轮：位置全对、喷距守恒，但 IK 只能选"反腕 + 大臂下探"的解
    // （θ2=154°、J2 超限 19°、小臂从工件底下穿过去 24.7mm）。
    bool eval(const V3& face, const V3& d, const V3* refX, const double* seed,
              double tiltDeg, double azDeg, bool checkLinks, PR& out) const {
        ++g_prof.evals;
        V3 zTool = -unit(d);
        V3 pC = face + d * p_.standoff;
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
        IkSol sol = ik(toBaseP(pC), toBaseR(Rc), p_.toolLen, seed);
        if (!sol.ok) return false;

        FkOut f;
        fk(sol.theta, p_.toolLen, f);
        std::vector<Seg> segs;
        buildSegments(f, sz_, p_.includeBase, segs);

        double worstTool = 1e300, worstLink = 1e300;
        const Seg* segTool = nullptr;
        const Seg* segLink = nullptr;
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
        const Seg* worstSeg = (worstTool <= worstLink) ? segTool : segLink;
        double worst = std::min(worstTool, worstLink);
        double selfGap = 1e300;
        for (const auto& pr : pairs_) {
            const Seg& sa = segs[pr.first];
            const Seg& sb = segs[pr.second];
            double dax = segSegDist(sa.a, sa.b, sb.a, sb.b);
            selfGap = std::min(selfGap, dax - sa.r - sb.r);
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
        out.cond = sol.cond;
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
inline void tiltAxisAndDir(const V3& nIn, V3& axis, V3& uUp) {
    V3 n = unit(nIn);
    V3 up(0, 0, 1);
    V3 u = up - n * dot(up, n);
    if (norm(u) < std::sin(5.0 * kDeg)) {
        axis = V3(0, 0, 0);
        uUp = V3(0, 0, 0);
        return;
    }
    uUp = unit(u);
    axis = unit(cross(n, uUp));
}

// 该法向**需要多大倾角**（0 = 正朝上，α=90° 侧壁时给满 tiltNonUpDeg）。
// 抽成独立函数是因为行级方位打分也要用它 —— 两处若各写一份，
// 打分用的倾角和真正执行的倾角会悄悄分叉，A/B 结论就不可信了。
inline double tiltNeedDeg(const V3& nIn, const Params& p) {
    V3 n = unit(nIn);
    double alphaDeg = std::acos(std::max(-1.0, std::min(1.0, n.z))) * kRad;
    double t = p.tiltNonUpDeg * (alphaDeg / 90.0);
    if (t > p.tiltMaxDeg) t = p.tiltMaxDeg;
    // 朝下的悬垂面：往上扳会把枪顶进工件本体，这条路物理上不成立，不做线性外推
    if (n.z < -0.7) t = 0.0;
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
    //   α = 法向与 +Z 的夹角（0° = 正朝上，90° = 竖直侧壁，180° = 完全倒挂）
    //   tilt0 = tiltNonUpDeg · (α / 90°)，封顶 tiltMaxDeg
    // 取 90° 作满量程锚点是因为它物理上就是"面立起来了"的界限：此时给满基准倾角，
    // 比侧壁更倒挂的面（α>90°）继续线性增长，直到上限。
    // 阶梯做法在 45° 边界上产生 20° 的阶跃，反映到轨迹上就是相邻点倾角突跳，
    // 实机过不了；线性化之后法向连续 → α 连续 → 倾角需求连续。
    double tilt0 = tiltNeedDeg(n, p);

    V3 axis, uUp;
    tiltAxisAndDir(n, axis, uUp);
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

// 只升不降的限幅平滑：final[i] = max(need[i], final[邻点] - rate)
void smoothOnlyUp(std::vector<double>& v, const std::vector<std::pair<int, int>>& rows, double rate) {
    for (const auto& r : rows) {
        for (int i = r.first + 1; i <= r.second; ++i) v[i] = std::max(v[i], v[i - 1] - rate);
        for (int i = r.second - 1; i >= r.first; --i) v[i] = std::max(v[i], v[i + 1] - rate);
    }
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

RobotSafetyParams robot_safety_default_params(void) {
    RobotSafetyParams p;
    std::memset(&p, 0, sizeof(p));
    p.standoff = 100.0;
    p.safety = 20.0;
    p.hardClearance = 2.0;   // 绝对安全距离：只在这个硬下限被破坏时才回退
    p.search = 120.0;
    p.toolLen = kDefaultToolLen;
    p.tiltNonUpDeg = 30.0;   // α=90°（竖直侧壁）时的倾角 = 工艺上限。见 needTilt 的说明
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
    p.maxCloudPoints = 0;
    p.threads = 1;           // 默认单线程：结果可复现，便于对拍/回归
    p.azMode = 0;            // 默认「按站位逐行择优」；1 = 只用向上（旧口径）
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

RobotSafetyResult robot_safety_check(const SprayTrajPoint* pts, int n,
                                     const double* cloudXyz, int cloudCount,
                                     const double* basePose16,
                                     const RobotSafetyParams* paramsIn) {
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
        // 语义：0（含零初始化结构）= 默认＝按站位择优；只有显式填 1 才退回旧口径
        p.azStance = (s.azMode != 1);
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
        V3 nm(0, 0, 0);
        int cnt = 0;
        for (int i = 0; i < n; ++i)
            if (t.spray[i]) {
                nm = nm + t.normal(i);
                ++cnt;
            }
        V3 nMean = cnt ? unit(nm) : V3(0, 0, 1);
        V3 h(nMean.x, nMean.y, 0.0);
        if (norm(h) < 0.2) h = V3(1, 0, 0);
        h = unit(h);
        V3 origin = cloud.centroid - h * (0.50 * kReach);
        origin.z = cloud.centroid.z - 0.35 * kReach;
        V3 x = unit(cloud.centroid - origin);
        V3 z(0, 0, 1);
        V3 y = unit(cross(z, x));
        x = unit(cross(y, z));
        base.R.setCol(0, x);
        base.R.setCol(1, y);
        base.R.setCol(2, z);
        base.t = origin;
    }

    Evaluator ev(cloud, base, p);

    int nThreads = std::max(1, std::min(p.threads, 64));

    // 点 → 行 的归属在建行之后就定死，阶段 0/1/3 都按它取行级方位
    std::vector<int> rowOf(n, -1);
    for (size_t ri = 0; ri < rows.size(); ++ri)
        for (int i = rows[ri].first; i <= std::min(rows[ri].second, n - 1); ++i) rowOf[i] = (int)ri;

    std::vector<char> stanceRow(rows.size(), 0);   // 该行方位是否被改判为「站位侧」
    // ---- 进刀侧（方位角）逐行决策放在「阶段 2」之后 —— 必须先有真实倾角，
    // 否则判的是 tilt0 那个位姿，而实际执行倾角常常大 20~30°，P0 结论会反号（见阶段 2 末尾）。
    int azStanceRows = 0;

    // ---- 阶段 1：每点的倾角需求（带 IK 种子连续性）
    std::vector<double> need(n, 0.0), azOf(n, 0.0);
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

    // ---- 阶段 2：倾角平滑（只升不降）
    g_prof.mark("阶段1 求倾角");
    std::vector<double> finalTilt = need;
    smoothOnlyUp(finalTilt, rows, p.tiltRateDeg);

    // 方位：同一行内有任一点需要换方位 → 整行统一用该方位（否则倾斜方向逐点跳，姿态更乱）。
    // 基准值取阶段 0 的行级决策（azStance=false 时恒为 0 = 向上，与旧口径逐字节一致）。
    // 取该行出现次数最多的方位；并列时取较小的（与 numpy.unique+argmax 同口径）。
    std::vector<double> azRow(rows.size(), 0.0);
    for (size_t ri = 0; ri < rows.size(); ++ri) {
        int a = rows[ri].first, b = rows[ri].second;
        std::vector<double> vals;
        for (int i = a; i <= b; ++i) {
            // 旧口径只统计非零方位（0 是默认值，投它等于没投）；
            // 站位口岸位非零，0 就变成了一个**有意义的候选**，必须一起投票。
            if (std::fabs(azOf[i]) > 1e-9)
                vals.push_back(std::round(azOf[i] * 1000.0) / 1000.0);
        }
        if (!vals.empty()) {
            std::sort(vals.begin(), vals.end());
            double bestV = vals[0];
            int bestC = 0;
            for (size_t k = 0; k < vals.size();) {
                size_t k2 = k;
                while (k2 < vals.size() && vals[k2] == vals[k]) ++k2;
                int c = (int)(k2 - k);
                if (c > bestC) {
                    bestC = c;
                    bestV = vals[k];
                }
                k = k2;
            }
            azRow[ri] = bestV;
        }
    }

    // 说明：曾在这里加过"拿到真实倾角后复核方位"的一轮（阶段 2.5），实测是负收益，
    // 已删。原因写在阶段 0 的注释里：复核用的是抽样，判断与阶段 3 的逐点复核
    // 不一致时会整行被逐点回退 —— complex 上花了 15~20% 时间，
    // 输出与旧口径逐字节相同。要真正按"最差点"决策，就必须整行逐点评估两个方位，
    // 那是 2× 阶段 1 的开销（阶段 1 已占 84%），产线不接受。
    // ---- 阶段 2.5：进刀侧（方位角）逐行决策
    //
    // 为什么必须放在这里（倾角已定之后），而不是更早：
    // 早期版本在阶段 1 之前就定方位，只能用 tilt0 打分。实测 complex 上就是这么翻车的
    // —— 那行执行倾角是 40~55°（倾角搜索为了躲碰撞自己加上去的），而 tilt0 只有 20°，
    // 在 tilt0 上"向上撞、站位不撞"，到真倾角上完全不成立：换完方位冒出 110 点关节超限，
    // 间隙从 77mm 掉到 37mm。**判 P0 必须判真正要执行的那个位姿。**
    //
    // 方位角只能是行级常量（它描述"枪从哪一侧压过来"，逐点变等于沿轨迹反复拧腕，
    // 实机执行不了），所以决策也在行级：候选只有 0（向上，旧口径）与站位侧。
    // 两者都不是先验正确的 —— 向上只是"重力上方通常空旷"的代理，内腔里上方常常正对
    // 对面的壁；站位侧则可能隔着工件本体。所以按实测择优，并且：
    //   1. 只在 P0（穿透 / 自碰撞 / 不可达）上严格改善才换；
    //   2. 采样先稀后密，稀采样看不出优势就直接放弃，看出优势再密集确认。
    // 于是"旧口径本来安全的样本"逐字节不变，方位只在旧口径真的撞了的行上动手。
    {
        struct AzScore {
            int hit = 0;        // 枪体/连杆 穿透点数
            int selfHit = 0;    // 自碰撞点数
            int unreach = 0;    // IK 解不出
            int jointBad = 0;   // 超软限位 / 余量不足（只统计，不参与换方位决策）
            int sing = 0;       // 近奇异（只统计）
            double minClr = 1e300;
        };
        auto p0Wins = [](const AzScore& a, const AzScore& b) {
            if (a.hit != b.hit) return a.hit < b.hit;
            if (a.unreach != b.unreach) return a.unreach < b.unreach;
            return a.selfHit < b.selfHit;
        };
        // 按**该点实际执行倾角** finalTilt[i] 打分
        auto scoreAz = [&](int a, int b, double az, int S) {
            AzScore s;
            int len = b - a + 1;
            S = std::min(S, len);
            for (int k = 0; k < S; ++k) {
                int i = (S <= 1) ? a : a + (int)((long long)k * (len - 1) / (S - 1));
                if (t.spray[i] == 0) continue;          // 过渡点位置可改，不参与方位决策
                V3 n = t.normal(i);
                V3 face = t.pos[i] + t.gun[i] * p.standoff;
                V3 axis, uUp;
                tiltAxisAndDir(n, axis, uUp);
                if (norm(uUp) < 1e-6) continue;
                double tilt = finalTilt[i];
                V3 d = dirOf(n, uUp, axis, tilt, az);
                PR pr;
                if (!ev.eval(face, d, nullptr, nullptr, tilt, az, true, pr)) { ++s.unreach; continue; }
                if (pr.clrTool < 0.0 || pr.clrLink < 0.0) ++s.hit;
                if (pr.selfGap < 0.0) ++s.selfHit;
                if (pr.over || pr.margin < p.jointMarginMin) ++s.jointBad;
                if (pr.cond > p.condMax) ++s.sing;
                if (std::isfinite(pr.clearance) && pr.clearance < s.minClr) s.minClr = pr.clearance;
            }
            return s;
        };
        auto decideRow = [&](size_t ri) {
            int a = rows[ri].first, b = rows[ri].second;
            if (b < a) return;
            int im = a + (b - a) / 2;
            V3 n = t.normal(im);
            V3 face = t.pos[im] + t.gun[im] * p.standoff;
            V3 axis, uUp;
            tiltAxisAndDir(n, axis, uUp);
            double azs = stanceAzDeg(n, face, ev.baseOrigin(), uUp, axis);
            // 无倾向（基座在法向延长线上）或与当前方位几乎重合 → 不做无谓计算
            if (!std::isfinite(azs) || std::fabs(azs - azRow[ri]) < 5.0) return;
            // 候选必须包含**当前真正执行的方位**：它常常已经不是 0 —— 阶段 1 的兜底
            // （azPref±30/60）会把整行投成 ±30。早期版本拿"0"当基准、见到非零就跳过，
            // 结果最需要救的 longstripe 两行（azRow=±30）根本没进决策，阶段 2.5 只花 0.3ms。
            std::vector<double> cands{azRow[ri]};
            if (std::fabs(0.0 - azRow[ri]) >= 5.0) cands.push_back(0.0);
            cands.push_back(azs);
            int len = b - a + 1;
            // 稀采样先筛：只有出现了"不是当前方位"的赢家才值得上密集采样
            int best = 0;
            AzScore bestS = scoreAz(a, b, cands[0], 5);
            for (size_t k = 1; k < cands.size(); ++k) {
                AzScore s = scoreAz(a, b, cands[k], 5);
                if (p0Wins(s, bestS)) { bestS = s; best = (int)k; }
            }
            if (best == 0) return;
            int Sd = std::min(120, std::max(8, len / 20));
            bestS = scoreAz(a, b, cands[0], Sd);
            best = 0;
            for (size_t k = 1; k < cands.size(); ++k) {
                AzScore s = scoreAz(a, b, cands[k], Sd);
                if (p0Wins(s, bestS)) { bestS = s; best = (int)k; }
            }
            if (best != 0) {
                azRow[ri] = cands[best];
                // 只有把方位换成**站位侧**才记账；换成 0（向上）不算本层做了什么
                stanceRow[ri] = (std::fabs(cands[best] - azs) < 1e-6) ? 1 : 0;
            }
        };
        if (p.azStance) {
            if (nThreads <= 1) {
                for (size_t ri = 0; ri < rows.size(); ++ri) decideRow(ri);
            } else {
                std::atomic<size_t> nextRow{0};
                std::vector<std::thread> pool;
                int nt = std::max(1, std::min<int>(nThreads, (int)rows.size()));
                pool.reserve(nt);
                for (int k = 0; k < nt; ++k) {
                    pool.emplace_back([&] {
                        for (;;) {
                            size_t ri = nextRow.fetch_add(1);
                            if (ri >= rows.size()) break;
                            decideRow(ri);
                        }
                    });
                }
                for (auto& th : pool) th.join();
            }
        }
    }
    g_prof.mark("阶段2.5 进刀侧");
    int azChangedRows = 0;
    for (char c : stanceRow) if (c) ++azChangedRows;
    if (g_prof.on) {
        std::fprintf(stderr, "  [az] rows=%zu 非零方位行=%d :", rows.size(), azChangedRows);
        for (size_t t = 0; t < azRow.size(); ++t)
            if (std::fabs(azRow[t]) > 1e-9) std::fprintf(stderr, " r%zu=%.1f", t, azRow[t]);
        int nonZeroPts = 0;
        for (int i = 0; i < n; ++i) if (rowOf[i] >= 0 && std::fabs(azRow[rowOf[i]]) > 1e-9) ++nonZeroPts;
        std::fprintf(stderr, " · 受影响点=%d\n", nonZeroPts);
    }

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
            double az = rowOf[i] >= 0 ? azRow[rowOf[i]] : 0.0;

            V3 axis, uUp;
            tiltAxisAndDir(nrm, axis, uUp);
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
                    // 所有候选倾角都解不出逆解 → 不可达（P0）
                    PR bad;
                    bad.idx = i;
                    bad.tiltDeg = tilt;
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
                pr.idx = i;
                // ---- 平滑回退：只回退硬约束，不再“差一点就回退”
                // 用户现场口径：有凸点没关系，只要不碰撞就以平滑为主。
                // 所以这里把“平滑后间隙比本点最优小”从回退条件里拿掉：
                //   · 平滑位姿满足硬约束（IK / 自碰撞 / 间隙≥hardClearance / 关节不超限）
                //     → 保留平滑，哪怕间隙比本点最优小几毫米；
                //   · 只有硬约束被破坏，且本点最优更安全时，才退回本点最优。
                // 注意：margin < jointMarginMin 仍只是 P1 软告警，不触发回退；
                //       safety（默认 20mm）也只影响报告标橙，不再是回退阈值。
                if (hasCache[i]) {
                    bool smoothHard = hardOk(pr, p.hardClearance);
                    if (!smoothHard) {
                        bool cacheHard = hardOk(cache[i], p.hardClearance);
                        if (cacheHard || hardCmp(pr, cache[i]) < 0) {
                            pr = cache[i];
                            pr.idx = i;
                            pr.reverted = true;
                            ++acc.reverted;
                        }
                    }
                }
            }
            pr.sprayAngle = angleBetweenDeg(-pr.gun, nrm);
            acc.standoffMaxErr = std::max(acc.standoffMaxErr, pr.standoffErr);

            // 穿透按部位细分：连杆/关节撞 与 枪体/喷嘴撞 是两类完全不同的现场问题
            // （前者要改站位或姿态，后者只要抬高喷距/换角度），混在一个 COLLISION 里没法决策。
            // RS_VIOL_COLLISION 保留定义但不再置位，避免 byRule 重复计数。
            if (pr.clrTool < 0.0) pr.viol |= RS_VIOL_TOOL_HIT;
            if (pr.clrLink < 0.0) pr.viol |= RS_VIOL_LINK_HIT;
            if (pr.clearance >= 0.0 && pr.clearance < p.safety) pr.viol |= RS_VIOL_CLEARANCE;
            if (pr.selfGap < 0.0) pr.viol |= RS_VIOL_SELF_COLLISION;
            if (pr.over || pr.margin < p.jointMarginMin) pr.viol |= RS_VIOL_JOINT_LIMIT;
            if (pr.sprayAngle > p.maxSprayAngleDeg + 1e-9) pr.viol |= RS_VIOL_SPRAY_ANGLE;
            if (pr.cond > p.condMax) pr.viol |= RS_VIOL_NEAR_SINGULAR;

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
        for (const auto& a : accs) {
            recomputed += a.recomputed;
            reverted += a.reverted;
            reachable += a.reachable;
            standoffMaxErr = std::max(standoffMaxErr, a.standoffMaxErr);
        }
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
            if (!stanceRow[ri]) continue;
            if (!results[i].reverted && std::fabs(results[i].azDeg - azRow[ri]) < 1e-6) applied[ri] = 1;
        }
        azStanceRows = 0;
        for (char c : applied) if (c) ++azStanceRows;
    }

    // ---- 过渡点抬升（spray=0 ⇒ 位置可改，只求不撞）
    g_prof.mark("阶段3 重解");
    int lifted = 0;
    if (p.fixTransit) {
        for (int i = 0; i < n; ++i) {
            PR& pr = results[i];
            if (!hasResult[i] || t.spray[i] != 0) continue;
            if (!std::isfinite(pr.clearance) || pr.clearance >= p.safety) continue;
            double lift = std::min(p.transitLiftMax, p.safety - pr.clearance + 5.0);
            pr.pos = pr.pos + V3(0, 0, lift);
            pr.lifted = lift;
            pr.clearance += lift;
            // 分组距离必须同步抬升，否则总结论（已避开）与"连杆最小距离"（还是抬升前）
            // 会自相矛盾 —— curved 上就出现过 clearance=20.02 而 clrTool=11.39。
            if (pr.clrTool < 1e299) pr.clrTool += lift;
            if (pr.clrLink < 1e299) pr.clrLink += lift;
            pr.viol &= ~(RS_VIOL_COLLISION | RS_VIOL_CLEARANCE | RS_VIOL_LINK_HIT | RS_VIOL_TOOL_HIT);
            ++lifted;
        }
    }

    // ---- 姿态跳变：平滑回退的代价，必须报出来（P2，速度/加速度风险，不是碰撞）
    if (p.tiltRateDeg > 0) {
        double thr = std::max(3.0 * p.tiltRateDeg, 10.0);
        for (const auto& r : rows) {
            for (int i = r.first + 1; i <= r.second; ++i) {
                // 不可达点没有位姿，谈"倾角跳变"没有意义，别让它再刷一条违规
                if (!results[i].ikOk || !results[i - 1].ikOk) continue;
                double dlt = std::fabs(results[i].sprayAngle - results[i - 1].sprayAngle);
                if (dlt > thr) results[i].viol |= RS_VIOL_TILT_JUMP;
            }
        }
    }

    // ---- L0 输入层自检（这一层不通，后面所有几何结论都是废的）
    int l0Issues = 0;
    {
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

    // ---- 输出
    g_prof.mark("过渡+自检");
    RobotSafetyPoint* out = new RobotSafetyPoint[n];
    std::memset(out, 0, sizeof(RobotSafetyPoint) * (size_t)n);
    std::vector<double> angles;
    angles.reserve(n);
    double minClr = std::numeric_limits<double>::infinity();
    int tilted = 0;
    double minClrTool = std::numeric_limits<double>::infinity();
    double minClrLink = std::numeric_limits<double>::infinity();
    // 未测到（IK 失败 / 本点没查连杆）统一写 nan：写 0 会被读成"贴脸 0mm"，
    // 那比"没测"更误导。
    auto measured = [](double v) {
        return (std::isfinite(v) && v < 1e299) ? v : std::numeric_limits<double>::quiet_NaN();
    };
    vc.total = 0;
    std::memset(vc.byRule, 0, sizeof(vc.byRule));

    for (int i = 0; i < n; ++i) {
        const PR& pr = results[i];
        RobotSafetyPoint& q = out[i];
        q.x = pr.pos.x;
        q.y = pr.pos.y;
        q.z = pr.pos.z;
        q.i = pr.gun.x;
        q.j = pr.gun.y;
        q.k = pr.gun.z;
        q.speed = t.speed[i];
        q.spray = t.spray[i];
        q.tiltDeg = pr.sprayAngle;
        q.clearanceMm = pr.clearance;
        q.selfGapMm = pr.selfGap;
        q.marginDeg = pr.margin;
        q.standoffErrMm = pr.standoffErr;
        q.liftedMm = pr.lifted;
        q.cond = pr.cond;
        for (int k = 0; k < 6; ++k) q.theta[k] = pr.ikOk ? pr.theta[k] : 0.0;
        q.ikOk = pr.ikOk ? 1 : 0;
        q.reverted = pr.reverted ? 1 : 0;
        q.violation = pr.viol;
        q.clearanceToolMm = measured(pr.clrTool);
        q.clearanceLinkMm = measured(pr.clrLink);
        q.worstKind = pr.ikOk ? pr.partKind : -1;
        q.worstJ0 = pr.ikOk ? pr.partJ0 : -1;
        q.worstJ1 = pr.ikOk ? pr.partJ1 : -1;
        q.linkKind = pr.ikOk ? pr.linkKind : -1;
        q.linkJ0 = pr.ikOk ? pr.linkJ0 : -1;
        q.linkJ1 = pr.ikOk ? pr.linkJ1 : -1;
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
    }

    double sum = 0.0;
    for (double a : angles) sum += a;

    R.points = out;
    R.count = n;
    R.reachable = reachable;
    R.tiltedPoints = tilted;
    R.poseReverted = reverted;
    R.violations = vc.total;
    for (int b = 0; b < RS_VIOL_COUNT; ++b) R.byRule[b] = vc.byRule[b];
    R.minClearanceMm = minClr;
    R.minClearanceToolMm = minClrTool;
    R.minClearanceLinkMm = minClrLink;
    R.tiltMeanDeg = n ? sum / n : 0.0;
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
    }
    return R;
}
