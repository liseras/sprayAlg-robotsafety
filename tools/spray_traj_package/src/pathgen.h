// pathgen.h - 分区喷涂轨迹：平面之字 / 孔周整圆 / 斜面直线 + 连续串接 + 工具位姿
//
// 工艺铁律：进入工作距离后不允许关枪。所有路点都是 spray，仅用 zone 区分工艺段。
#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <queue>
#include <functional>
#include "vecmath.h"
#include "coverage.h"
#include "heightfield.h"
#include "regions.h"

static const double PI_P = 3.14159265358979323846;

struct SprayParams {
    double passSpacing = 100;   // 行距 mm（喷粉出粉量大, 取喷幅 50~70%）
    double standoff = 200;      // 喷距 mm（沿表面法向偏置）
    double margin = 12;         // 边界内缩 mm
    double step = 5;            // 轨迹点间距 mm
    double minRunLen = 25;      // 最短行程 mm（短于此不喷, 避免原地抖动）
    double tiltDeg = 0;         // 前倾角 deg
    double speed = 400;         // mm/s
    bool flat = false;          // true = 完全按平面（忽略起伏）
    double slopeAngle = 25;     // 陡峭判定角 deg（超过即视为独立斜面区）
    double slopeSpacing = 0;    // 斜面内线间距 mm（0 = 同 passSpacing）
    double slopeMargin = 6;     // 斜面两端内缩 mm
    double ringOffset = 30;     // 孔周环线半径 = 孔等效半径 + 该值
    int ringCount = 1;          // 孔周环线条数
    double ringSpacing = 0;     // 多环间距（0 = 同 passSpacing）
};

struct Waypoint {
    Vec3 pos;                    // 喷枪 TCP 位置（点云坐标系 = 机器人基座）
    Vec3 nref;                   // 期望表面法向（工具 Z = -nref）
    Vec3 toolZ;                  // 最终工具 Z 轴
    double qw = 1, qx = 0, qy = 0, qz = 0;
    std::string type = "spray";  // 恒为 spray（工作距离内不关枪）
    std::string zone = "plane";  // plane / ring / slope / travel
    int pass = 0;                // 行程/线段序号
    int group = 0;               // 区域组号：环线用孔序号, 斜面用线序号
    double speed = 400;
    double u = 0, v = 0, h = 0;  // 调试: 平面坐标与采样高度
};

// 一条连续的工艺线段（TCP 轨迹点，已含喷距偏置）
struct Run {
    std::string zone;
    int group = 0;
    Vec3 n;                      // 该段的表面法向（工具 Z = -n）
    std::vector<Vec3> p;
    std::vector<double> uu, vv, hh;   // 区域局部坐标（平面/环线=主平面uv, 斜面=斜面st）
    std::vector<double> ppu, ppv;     // 表面点在主平面的投影（串接/A* 统一用）
    bool closed = false;
};

// 工具位姿：Z=-法向，X=行进切向，Y=Z×X；可选前倾角
inline void makeOrientation(const Vec3& n, const Vec3& tangent, double tiltDeg,
                            Vec3& X, Vec3& Y, Vec3& Z) {
    Z = n * -1.0;
    Vec3 t = tangent - Z * dot(tangent, Z);
    X = (norm(t) > 1e-9) ? normalize(t) : Vec3(1, 0, 0);
    Y = cross(Z, X);
    if (tiltDeg != 0) {
        double a = tiltDeg * PI_P / 180.0;
        Vec3 Zn = normalize(Z * std::cos(a) + X * std::sin(a));
        Vec3 Xn = normalize(X * std::cos(a) - Z * std::sin(a));
        Z = Zn; X = Xn; Y = cross(Z, X);
    }
}

// A* 最短路（8 连通, 只走掩膜格）：过渡段绕开孔洞沿边界走
inline std::vector<std::pair<int,int>> astarPath(const Mask& mk, int si, int sj, int ei, int ej) {
    const int W = mk.W, H = mk.H;
    auto id = [&](int i, int j) { return (size_t)j * W + i; };
    auto nearestMask = [&](int i0, int j0) -> std::pair<int,int> {
        if (mk.inside(i0, j0) && mk.at(i0, j0)) return {i0, j0};
        for (int r = 1; r <= 8; ++r)
            for (int dj = -r; dj <= r; ++dj)
                for (int di = -r; di <= r; ++di) {
                    if (std::max(std::abs(di), std::abs(dj)) != r) continue;
                    int ni = i0 + di, nj = j0 + dj;
                    if (mk.inside(ni, nj) && mk.at(ni, nj)) return {ni, nj};
                }
        return {-1, -1};
    };
    auto s = nearestMask(si, sj), e = nearestMask(ei, ej);
    if (s.first < 0 || e.first < 0) return {};
    si = s.first; sj = s.second; ei = e.first; ej = e.second;
    if (si == ei && sj == ej) return {{si, sj}};

    const float INF = 1e18f;
    std::vector<float> g((size_t)W * H, INF);
    std::vector<int> came((size_t)W * H, -1);
    using QI = std::pair<float, int>;
    std::priority_queue<QI, std::vector<QI>, std::greater<QI>> open;
    auto heur = [&](int i, int j) {
        float dx = (float)std::abs(i - ei), dy = (float)std::abs(j - ej);
        return std::max(dx, dy) + 0.4142f * std::min(dx, dy);
    };
    g[id(si, sj)] = 0;
    open.push({heur(si, sj), (int)id(si, sj)});
    const int di[8] = {1, -1, 0, 0, 1, 1, -1, -1};
    const int dj[8] = {0, 0, 1, -1, 1, -1, 1, -1};
    const float dc[8] = {1, 1, 1, 1, 1.4142f, 1.4142f, 1.4142f, 1.4142f};
    size_t goal = id(ei, ej);
    while (!open.empty()) {
        auto [f, k] = open.top();
        open.pop();
        if ((size_t)k == goal) break;
        int ci = k % W, cj = k / W;
        if (f > g[k] + heur(ci, cj) + 1e-4) continue;
        for (int n = 0; n < 8; ++n) {
            int ni = ci + di[n], nj = cj + dj[n];
            if (!mk.inside(ni, nj) || !mk.at(ni, nj)) continue;
            size_t nk = id(ni, nj);
            float ng = g[k] + dc[n];
            if (ng < g[nk]) {
                g[nk] = ng;
                came[nk] = k;
                open.push({ng + heur(ni, nj), (int)nk});
            }
        }
    }
    if (came[goal] < 0) return {};
    std::vector<std::pair<int,int>> path;
    size_t k = goal;
    while (true) {
        path.push_back({(int)(k % W), (int)(k / W)});
        if (k == (size_t)id(si, sj)) break;
        if (came[k] < 0) return {};
        k = (size_t)came[k];
    }
    std::reverse(path.begin(), path.end());
    return path;
}

// ============================================================
// 1) 平面区：横向之字（沿 u 长边行走，v 方向进给）
// ============================================================
inline std::vector<Run> buildPlaneRuns(const Mask& mkMain, const Plane& pl,
                                       const Vec3& uAx, const Vec3& vAx,
                                       const HeightField& hf, const SprayParams& P) {
    std::vector<Run> runs;
    double vLo = 1e18, vHi = -1e18, uLo = 1e18, uHi = -1e18;
    for (int j = 0; j < mkMain.H; ++j)
        for (int i = 0; i < mkMain.W; ++i)
            if (mkMain.at(i, j)) {
                vLo = std::min(vLo, mkMain.v(j)); vHi = std::max(vHi, mkMain.v(j));
                uLo = std::min(uLo, mkMain.u(i)); uHi = std::max(uHi, mkMain.u(i));
            }
    if (vHi < vLo) return runs;

    std::vector<double> rowV;
    double span = vHi - vLo;
    if (span <= P.passSpacing) {
        rowV.push_back(0.5 * (vLo + vHi));
    } else {
        // 首末行贴边, 实际行距 <= 设定值（保证不漏喷）
        int nRows = (int)std::ceil(span / P.passSpacing) + 1;
        double d = span / (nRows - 1);
        for (int r = 0; r < nRows; ++r) rowV.push_back(vLo + r * d);
    }

    int passId = 0;
    for (double v : rowV) {
        int j = std::min(mkMain.H - 1, std::max(0, (int)((v - mkMain.vmin) / mkMain.cell)));
        int i = 0;
        bool any = false;
        while (i < mkMain.W) {
            if (!mkMain.at(i, j)) { i++; continue; }
            int i0 = i;
            while (i < mkMain.W && mkMain.at(i, j)) i++;
            int i1 = i - 1;
            double u0 = mkMain.umin + i0 * mkMain.cell;
            double u1 = mkMain.umin + (i1 + 1) * mkMain.cell;
            if (u1 - u0 < P.minRunLen) continue;
            // 两端各内缩半个栅格, 避免贴到掩膜锯齿边
            u0 += mkMain.cell * 0.5;
            u1 -= mkMain.cell * 0.5;
            if (u1 - u0 < P.minRunLen) continue;
            int n = std::max(2, (int)std::ceil((u1 - u0) / P.step) + 1);
            Run r;
            r.zone = "plane";
            r.group = passId;
            r.n = pl.n;
            for (int k = 0; k < n; ++k) {
                double u = u0 + (u1 - u0) * k / (n - 1);
                double h = P.flat ? 0.0 : hf.sample(u, v);
                Vec3 surf = pl.c + uAx * u + vAx * v + pl.n * h;
                r.p.push_back(surf + pl.n * P.standoff);
                r.uu.push_back(u); r.vv.push_back(v); r.hh.push_back(h);
                r.ppu.push_back(u); r.ppv.push_back(v);
            }
            runs.push_back(std::move(r));
            any = true;
        }
        if (any) passId++;
    }
    (void)uLo; (void)uHi;
    return runs;
}

// ============================================================
// 2) 孔周区：沿整圆走一圈（半径 = 孔等效半径 + ringOffset）
//    半径搜索: 优先在腐蚀掩膜(mkTrav)内放, 放不下逐毫米收缩,
//    仍放不下退到未腐蚀掩膜(允许贴近板边卷边), 再不行才跳过
// ============================================================
inline bool circleFits(const Mask& mk, double cu, double cv, double r, double step) {
    int N = std::max(72, (int)std::ceil(2 * PI_P * r / step));
    for (int k = 0; k < N; ++k) {
        double a = 2 * PI_P * k / N;
        double u = cu + r * std::cos(a), v = cv + r * std::sin(a);
        int ii = (int)((u - mk.umin) / mk.cell);
        int jj = (int)((v - mk.vmin) / mk.cell);
        if (!mk.inside(ii, jj) || !mk.at(ii, jj)) return false;
    }
    return true;
}

inline std::vector<Run> buildRingRuns(const std::vector<HoleInfo>& holes, const Mask& mkTrav,
                                      const Mask& mkFullMask,
                                      const Plane& pl, const Vec3& uAx, const Vec3& vAx,
                                      const HeightField& hf, const SprayParams& P,
                                      std::vector<std::string>* log = nullptr) {
    std::vector<Run> runs;
    double ringSp = P.ringSpacing > 0 ? P.ringSpacing : P.passSpacing;
    for (size_t hi = 0; hi < holes.size(); ++hi) {
        const HoleInfo& ho = holes[hi];
        for (int ri = 0; ri < P.ringCount; ++ri) {
            double r0 = ho.rEq + P.ringOffset - ri * ringSp;
            double rFloor = ho.rEq + 6.0;
            bool placed = false;
            // 两级掩膜: 先腐蚀后(安全), 再未腐蚀(贴边兜底)
            for (int lvl = 0; lvl < 2 && !placed; ++lvl) {
                const Mask& mk = (lvl == 0) ? mkTrav : mkFullMask;
                for (double r = r0; r >= rFloor && !placed; r -= 1.0) {
                    if (!circleFits(mk, ho.cu, ho.cv, r, P.step)) continue;
                    int N = std::max(32, (int)std::ceil(2 * PI_P * r / P.step));
                    Run run;
                    run.zone = "ring";
                    run.group = (int)hi;
                    run.n = pl.n;
                    run.closed = true;
                    for (int k = 0; k < N; ++k) {
                        double a = 2 * PI_P * k / N;
                        double u = ho.cu + r * std::cos(a), v = ho.cv + r * std::sin(a);
                        double h = P.flat ? 0.0 : hf.sample(u, v);
                        Vec3 surf = pl.c + uAx * u + vAx * v + pl.n * h;
                        run.p.push_back(surf + pl.n * P.standoff);
                        run.uu.push_back(u); run.vv.push_back(v); run.hh.push_back(h);
                        run.ppu.push_back(u); run.ppv.push_back(v);
                    }
                    runs.push_back(std::move(run));
                    placed = true;
                    if (log) {
                        char buf[256];
                        std::snprintf(buf, sizeof(buf),
                                      "孔#%zu 中心(%.0f,%.0f) rEq %.1f 环线 R=%.1f%s (%d 点)",
                                      hi, ho.cu, ho.cv, ho.rEq, r,
                                      lvl == 0 ? "" : " [贴边模式]", N);
                        log->push_back(buf);
                    }
                }
            }
            if (!placed && log) {
                // 诊断: 初始半径下圆周有多少点出界
                int bad = 0;
                int N = 360;
                for (int k = 0; k < N; ++k) {
                    double a = 2 * PI_P * k / N;
                    double u = ho.cu + r0 * std::cos(a), v = ho.cv + r0 * std::sin(a);
                    int ii = (int)((u - mkTrav.umin) / mkTrav.cell);
                    int jj = (int)((v - mkTrav.vmin) / mkTrav.cell);
                    if (!mkTrav.inside(ii, jj) || !mkTrav.at(ii, jj)) bad++;
                }
                char buf[256];
                std::snprintf(buf, sizeof(buf),
                              "孔#%zu 环线放不下: rEq %.1f, R=%.1f 时 %d/360 出界, 已收缩到 R>=%.1f",
                              hi, ho.rEq, r0, bad, rFloor);
                log->push_back(buf);
            }
        }
    }
    return runs;
}

// ============================================================
// 3) 斜面区：拟合斜面自身平面，沿长度方向走直线，工具垂直斜面
// ============================================================
inline bool buildSlopeZone(const std::vector<Vec3>& pts, const std::vector<uint8_t>& steepMask,
                           const Mask& mk, const Plane& mainPl, const Vec3& uAx, const Vec3& vAx,
                           SlopeZone& sz, int minPts = 400) {
    sz.ptIdx.clear();
    for (size_t i = 0; i < pts.size(); ++i) {
        Vec3 d = pts[i] - mainPl.c;
        int a = (int)((dot(d, uAx) - mk.umin) / mk.cell);
        int b = (int)((dot(d, vAx) - mk.vmin) / mk.cell);
        if (!mk.inside(a, b)) continue;
        if (steepMask[(size_t)b * mk.W + a]) sz.ptIdx.push_back((int)i);
    }
    if ((int)sz.ptIdx.size() < minPts) return false;
    sz.plane = fitPlaneTLS(pts, sz.ptIdx);
    if (dot(sz.plane.n, mainPl.n) < 0) sz.plane.n = sz.plane.n * -1.0;

    // 面内主轴：最大方差方向 = 斜面长度方向
    Vec3 c(0, 0, 0);
    for (int i : sz.ptIdx) c = c + pts[i];
    c = c / (double)sz.ptIdx.size();
    double C[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int i : sz.ptIdx) {
        Vec3 d = pts[i] - c;
        C[0][0] += d.x * d.x; C[0][1] += d.x * d.y; C[0][2] += d.x * d.z;
        C[1][1] += d.y * d.y; C[1][2] += d.y * d.z; C[2][2] += d.z * d.z;
    }
    C[1][0] = C[0][1]; C[2][0] = C[0][2]; C[2][1] = C[1][2];
    double eig[3]; Vec3 vec[3];
    jacobiEigen3(C, eig, vec);
    Vec3 lenAx = vec[2];
    lenAx = normalize(lenAx - sz.plane.n * dot(lenAx, sz.plane.n));
    sz.lenAxis = lenAx;
    sz.widAxis = normalize(cross(sz.plane.n, lenAx));
    sz.angleToMain = std::acos(std::min(1.0, std::fabs(dot(sz.plane.n, mainPl.n)))) * 180.0 / PI_P;

    sz.hf.build(pts, sz.ptIdx, sz.plane, sz.lenAxis, sz.widAxis, 4.0);
    sz.lenMin = 1e18; sz.lenMax = -1e18; sz.widMin = 1e18; sz.widMax = -1e18;
    sz.ss.clear(); sz.tt.clear();
    sz.ss.reserve(sz.ptIdx.size()); sz.tt.reserve(sz.ptIdx.size());
    for (int i : sz.ptIdx) {
        Vec3 d = pts[i] - sz.plane.c;
        double s = dot(d, sz.lenAxis), t = dot(d, sz.widAxis);
        sz.ss.push_back(s); sz.tt.push_back(t);
        sz.lenMin = std::min(sz.lenMin, s); sz.lenMax = std::max(sz.lenMax, s);
        sz.widMin = std::min(sz.widMin, t); sz.widMax = std::max(sz.widMax, t);
    }
    sz.valid = true;
    return true;
}

// 斜面直线: 沿长度方向走, 每条线按工件掩膜裁剪成若干段（不得横穿孔洞/越过板边）
inline std::vector<Run> buildSlopeRuns(const SlopeZone& sz, const Mask& clipMask,
                                       const Plane& mainPl, const Vec3& uAx, const Vec3& vAx,
                                       const SprayParams& P) {
    std::vector<Run> runs;
    if (!sz.valid) return runs;
    double sp = P.slopeSpacing > 0 ? P.slopeSpacing : P.passSpacing;
    double w0 = sz.widMin, w1 = sz.widMax;
    double width = w1 - w0;
    if (width < 3) return runs;

    std::vector<double> lineT;
    if (width <= sp) {
        lineT.push_back(0.5 * (w0 + w1));      // 一条线居中即可覆盖整条窄带
    } else {
        int n = (int)std::ceil(width / sp) + 1;
        double d = width / (n - 1);
        for (int k = 0; k < n; ++k) lineT.push_back(w0 + k * d);
    }

    int lineId = 0;
    for (double t : lineT) {
        // 该线附近的点决定长度范围（自然裁剪到斜面真实范围）
        double band = std::max(2.0, 0.4 * sp);
        double sLo = 1e18, sHi = -1e18;
        for (size_t k = 0; k < sz.ss.size(); ++k) {
            if (std::fabs(sz.tt[k] - t) > band) continue;
            sLo = std::min(sLo, sz.ss[k]);
            sHi = std::max(sHi, sz.ss[k]);
        }
        if (sHi < sLo) continue;
        sLo += P.slopeMargin; sHi -= P.slopeMargin;
        if (sHi - sLo < P.minRunLen) continue;
        int n = std::max(2, (int)std::ceil((sHi - sLo) / P.step) + 1);
        // 采样 + 按主平面掩膜分段（不得越界/穿孔）
        std::vector<uint8_t> on((size_t)n, 0);
        std::vector<double> su((size_t)n);
        for (int k = 0; k < n; ++k) {
            double s = sLo + (sHi - sLo) * k / (n - 1);
            su[(size_t)k] = s;
            Vec3 pos = sz.plane.c + sz.lenAxis * s + sz.widAxis * t;
            Vec3 rel = pos - mainPl.c;
            double uu = dot(rel, uAx), vv = dot(rel, vAx);
            int a = (int)((uu - clipMask.umin) / clipMask.cell);
            int b = (int)((vv - clipMask.vmin) / clipMask.cell);
            on[(size_t)k] = (clipMask.inside(a, b) && clipMask.at(a, b)) ? 1 : 0;
        }
        // 连续 on 段 -> 独立行程
        int k = 0;
        while (k < n) {
            if (!on[(size_t)k]) { k++; continue; }
            int k0 = k;
            while (k < n && on[(size_t)k]) k++;
            int k1 = k - 1;
            double a0 = su[(size_t)k0] + P.slopeMargin * 0.5;
            double a1 = su[(size_t)k1] - P.slopeMargin * 0.5;
            if (a1 - a0 < P.minRunLen) continue;
            int m = std::max(2, (int)std::ceil((a1 - a0) / P.step) + 1);
            Run r;
            r.zone = "slope";
            r.group = lineId;
            r.n = sz.plane.n;
            for (int q = 0; q < m; ++q) {
                double s = a0 + (a1 - a0) * q / (m - 1);
                double h = P.flat ? 0.0 : sz.hf.sample(s, t);
                Vec3 surf = sz.plane.c + sz.lenAxis * s + sz.widAxis * t + sz.plane.n * h;
                r.p.push_back(surf + sz.plane.n * P.standoff);
                r.uu.push_back(s); r.vv.push_back(t); r.hh.push_back(h);
                Vec3 srel = surf - mainPl.c;
                r.ppu.push_back(dot(srel, uAx)); r.ppv.push_back(dot(srel, vAx));
            }
            runs.push_back(std::move(r));
        }
        ++lineId;
    }
    return runs;
}

// ============================================================
// 4) 段串接：按给定顺序把各段接成一条连续轨迹（全程出粉）
//    closed 段（环线）自动旋转到离当前点最近的起点
// ============================================================
inline std::vector<Waypoint> sequenceRuns(std::vector<Run>& runs, const Mask& mkTrav,
                                          const Plane& pl, const Vec3& uAx, const Vec3& vAx,
                                          const HeightField& hf, const SprayParams& P) {
    std::vector<Waypoint> wp;
    if (runs.empty()) return wp;
    auto uv2cell = [&](double u, double v) {
        int i = (int)std::round((u - mkTrav.umin) / mkTrav.cell - 0.5);
        int j = (int)std::round((v - mkTrav.vmin) / mkTrav.cell - 0.5);
        return std::make_pair(std::min(mkTrav.W - 1, std::max(0, i)),
                              std::min(mkTrav.H - 1, std::max(0, j)));
    };
    Vec3 prevEnd;
    Vec3 prevN = pl.n;    // 上一段的表面法向（过渡段偏置方向从它渐变到下一段法向）
    double prevU = 0, prevV = 0;
    bool havePrev = false;

    for (auto& r : runs) {
        if (r.p.size() < 2) continue;
        size_t startIdx = 0;
        bool fwd = true;
        if (r.closed) {
            double bd = 1e18;
            for (size_t k = 0; k < r.p.size(); ++k) {
                double d = havePrev ? norm(r.p[k] - prevEnd) : norm(r.p[k]);
                if (d < bd) { bd = d; startIdx = k; }
            }
        } else if (havePrev) {
            fwd = norm(r.p.front() - prevEnd) <= norm(r.p.back() - prevEnd);
            startIdx = fwd ? 0 : r.p.size() - 1;
        }

        // --- 过渡（连续喷涂行走, 沿掩膜 A*, 遇孔绕行）---
        if (havePrev) {
            double su = r.ppu[startIdx], sv = r.ppv[startIdx];
            auto sc = uv2cell(prevU, prevV), ec = uv2cell(su, sv);
            auto cells = astarPath(mkTrav, sc.first, sc.second, ec.first, ec.second);
            std::vector<std::pair<double,double>> poly;
            poly.push_back({prevU, prevV});
            if (cells.size() >= 2) {
                for (size_t k = 1; k + 1 < cells.size(); ++k) {
                    int dx = cells[k].first - cells[k - 1].first;
                    int dy = cells[k].second - cells[k - 1].second;
                    int nx = cells[k + 1].first - cells[k].first;
                    int ny = cells[k + 1].second - cells[k].second;
                    if (dx != nx || dy != ny)
                        poly.push_back({mkTrav.u(cells[k].first), mkTrav.v(cells[k].second)});
                }
            }
            poly.push_back({su, sv});
            double total = 0;
            for (size_t k = 1; k < poly.size(); ++k)
                total += std::hypot(poly[k].first - poly[k - 1].first,
                                    poly[k].second - poly[k - 1].second);
            int n = std::max(1, (int)std::ceil(total / std::max(P.step, mkTrav.cell)));
            for (int k = 1; k <= n; ++k) {
                double target = total * k / n, acc = 0;
                double uu = poly.back().first, vv = poly.back().second;
                for (size_t m = 1; m < poly.size(); ++m) {
                    double L = std::hypot(poly[m].first - poly[m - 1].first,
                                          poly[m].second - poly[m - 1].second);
                    if (acc + L >= target - 1e-9 || m + 1 == poly.size()) {
                        double f = L > 1e-9 ? (target - acc) / L : 0;
                        f = std::max(0.0, std::min(1.0, f));
                        uu = poly[m - 1].first + (poly[m].first - poly[m - 1].first) * f;
                        vv = poly[m - 1].second + (poly[m].second - poly[m - 1].second) * f;
                        break;
                    }
                    acc += L;
                }
                double h = P.flat ? 0.0 : hf.sample(uu, vv);
                Vec3 surf = pl.c + uAx * uu + vAx * vv + pl.n * h;
                // 偏置方向从上段法向渐变到本段法向（换区时枪位置不跳变）
                double f = (double)k / n;
                Vec3 dir = normalize(prevN * (1 - f) + r.n * f);
                if (norm(dir) < 1e-9) dir = r.n;
                Waypoint w;
                w.pos = surf + dir * P.standoff;
                w.u = uu; w.v = vv; w.h = h;
                w.nref = dir;
                w.toolZ = dir * -1.0;
                w.zone = "travel";
                w.pass = r.group;
                w.group = r.group;
                w.type = "spray";
                wp.push_back(std::move(w));
            }
        }

        // --- 本段 ---
        size_t n = r.p.size();
        size_t emitLast = startIdx;
        for (size_t k = 0; k < n; ++k) {
            size_t idx = r.closed ? (startIdx + k) % n : (fwd ? k : n - 1 - k);
            emitLast = idx;
            Waypoint w;
            w.pos = r.p[idx];
            w.u = r.uu[idx]; w.v = r.vv[idx]; w.h = r.hh[idx];
            w.nref = r.n;
            w.toolZ = r.n * -1.0;
            w.zone = r.zone;
            w.pass = r.group;
            w.group = r.group;
            w.type = "spray";
            wp.push_back(std::move(w));
        }
        prevEnd = wp.back().pos;
        prevU = r.ppu[emitLast]; prevV = r.ppv[emitLast];
        prevN = r.n;
        havePrev = true;
    }
    return wp;
}

// ============================================================
// 5) 姿态后处理：法向插值/平滑 → 限幅求 Z → 限幅求 X → 四元数
//    解决分区之间（平面→斜面 58°）与 U 型弯处的姿态跳变
// ============================================================
inline void finalizeTrajectory(std::vector<Waypoint>& wp, const SprayParams& P) {
    const int n = (int)wp.size();
    if (n < 2) return;

    // 1) 过渡段法向：在前后工艺段法向之间线性插值
    for (int i = 0; i < n;) {
        if (wp[(size_t)i].zone != "travel") { ++i; continue; }
        int j = i;
        while (j < n && wp[(size_t)j].zone == "travel") ++j;
        Vec3 n0 = (i > 0) ? wp[(size_t)i - 1].nref : wp[(size_t)(j < n ? j : n - 1)].nref;
        Vec3 n1 = (j < n) ? wp[(size_t)j].nref : n0;
        for (int k = i; k < j; ++k) {
            double f = (double)(k - i + 1) / (j - i + 1);
            Vec3 v = n0 * (1 - f) + n1 * f;
            wp[(size_t)k].nref = norm(v) > 1e-9 ? normalize(v) : n0;
        }
        i = j;
    }
    // 2) 法向窗口平滑
    std::vector<Vec3> sm((size_t)n);
    for (int k = 0; k < n; ++k) {
        Vec3 acc(0, 0, 0); double ws = 0;
        for (int m = std::max(0, k - 4); m <= std::min(n - 1, k + 4); ++m) {
            double w = 1.0 - 0.18 * std::abs(m - k);
            acc = acc + wp[(size_t)m].nref * w;
            ws += w;
        }
        sm[(size_t)k] = normalize(acc * (1.0 / ws));
        if (norm(sm[(size_t)k]) < 1e-9) sm[(size_t)k] = wp[(size_t)k].nref;
    }
    for (int k = 0; k < n; ++k) wp[(size_t)k].nref = sm[(size_t)k];

    // 3) 平滑切向
    std::vector<Vec3> dirs((size_t)n - 1);
    for (int k = 0; k + 1 < n; ++k) {
        Vec3 d = wp[(size_t)k + 1].pos - wp[(size_t)k].pos;
        dirs[(size_t)k] = norm(d) > 1e-9 ? normalize(d)
                                        : (k > 0 ? dirs[(size_t)k - 1] : Vec3(1, 0, 0));
    }
    const int W = 8;
    std::vector<Vec3> tgt((size_t)n);
    for (int i = 0; i < n; ++i) {
        Vec3 acc(0, 0, 0); double ws = 0;
        for (int j = std::max(0, i - W); j <= std::min(n - 2, i + W); ++j) {
            double w = 1.0 - std::fabs((double)(j - i) + 0.5) / (W + 1);
            acc = acc + dirs[(size_t)j] * w;
            ws += w;
        }
        Vec3 t = acc * (1.0 / ws);
        Vec3 Z = wp[(size_t)i].nref * -1.0;
        Vec3 tp = t - Z * dot(t, Z);
        tgt[(size_t)i] = norm(tp) > 1e-6 ? normalize(tp)
                                        : (i > 0 ? tgt[(size_t)i - 1] : Vec3(1, 0, 0));
    }

    // 4) 逐点限幅：Z 每步 <= 8°, X 每步 <= 15°
    const double maxZ = 8.0 * PI_P / 180.0, maxX = 15.0 * PI_P / 180.0;
    auto rotToward = [](const Vec3& from, const Vec3& to, double maxAng) -> std::pair<Vec3, bool> {
        double c = std::max(-1.0, std::min(1.0, dot(from, to)));
        double ang = std::acos(c);
        if (ang <= maxAng) return {to, false};
        Vec3 ax = cross(from, to);
        if (norm(ax) < 1e-9) {
            // from 与 to 共线（180° 掉头）: 取一个与 from 垂直的轴作旋转轴, 逐步转
            Vec3 helper = (std::fabs(from.x) < 0.9) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
            ax = normalize(cross(from, helper));
        } else {
            ax = normalize(ax);
        }
        double cc = std::cos(maxAng), ss = std::sin(maxAng);
        Vec3 v = from * cc + cross(ax, from) * ss + ax * (dot(ax, from) * (1 - cc));
        return {normalize(v), true};
    };

    Vec3 Zprev = wp[0].nref * -1.0;
    Vec3 Xprev = normalize(tgt[0] - Zprev * dot(tgt[0], Zprev));
    if (norm(Xprev) < 1e-9) Xprev = Vec3(1, 0, 0);
    std::vector<uint8_t> fast((size_t)n, 0);
    for (int i = 0; i < n; ++i) {
        Vec3 Zt = wp[(size_t)i].nref * -1.0;
        auto rz = rotToward(Zprev, Zt, maxZ);
        Vec3 Z = rz.first;
        if (rz.second) fast[(size_t)i] = 1;
        Vec3 Xt = tgt[(size_t)i] - Z * dot(tgt[(size_t)i], Z);
        if (norm(Xt) < 1e-6) Xt = Xprev - Z * dot(Xprev, Z);
        Xt = normalize(Xt);
        auto rx = rotToward(Xprev, Xt, maxX);
        Vec3 X = rx.first;
        if (rx.second) fast[(size_t)i] = 1;
        X = normalize(X - Z * dot(X, Z));

        Vec3 Zf = Z, Xf = X;
        if (P.tiltDeg != 0) {   // 前倾角：绕工具 Y 轴向行进方向倾斜
            double a = P.tiltDeg * PI_P / 180.0;
            Zf = normalize(Z * std::cos(a) + X * std::sin(a));
            Xf = normalize(X * std::cos(a) - Z * std::sin(a));
        }
        Vec3 Y = cross(Zf, Xf);
        quatFromBasis(Xf, Y, Zf, wp[(size_t)i].qw, wp[(size_t)i].qx,
                      wp[(size_t)i].qy, wp[(size_t)i].qz);
        wp[(size_t)i].toolZ = Zf;
        Zprev = Z; Xprev = X;
    }
    // 5) 姿态变化快处降速（拐角/换区）
    for (int i = 0; i < n; ++i) {
        if (!fast[(size_t)i]) continue;
        for (int j = std::max(0, i - 3); j <= std::min(n - 1, i + 3); ++j)
            wp[(size_t)j].speed = std::min(wp[(size_t)j].speed, P.speed * 0.3);
    }
}
