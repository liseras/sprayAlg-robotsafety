// ordering.h - 工艺段顺序优化（开放路径 TSP：段可反向、起点自由、空移最短）
//
// 解决的问题
// ----------
// pathgen.h 的 sequenceRuns() 只对"即将接上的那一段"做就近选端（fwd / rev），
// **整套顺序沿用生成顺序**：平面之字行全在前、孔周环按单端贪心插在中间、斜面线
// 完全不排序。段一多就会出现"从工件这头跑到那头、再跑回来"的长空移
// （p029 实测过渡段占全长 41%）。
//
// 本模块把每一段看成一条"可反向的有向线段"（两端点任选其一作为进入点），
// 求一条经过所有段恰好一次的**开放路径**，使总空移最短：
//   · 段内工艺路径不变（之字行的行内走向、环线的圆周、斜面线的直线都照旧）；
//   · 每段的进入端自由 → 等价于"找最近的端点作为起点"；
//   · 下一段在所有未喷段里按代价最小挑选 → 等价于"空移选最近的一段"；
//   · 最后用 2-opt / Or-opt 消除贪心的短视（跨区往返、夹花）。
//
// 代价模型（三层，都可调）
// ------------------------
//   ① 行程长度：掩膜上的最短路（8 连通 Dijkstra 距离场），自动绕开大孔、沿板边走；
//      也可退化为直线距离（metric=0），用于超大栅格或对速度敏感的场景。
//   ② 换区惩罚：平面/孔周/斜面之间换区要换姿态，同距离下优先同区
//      （zonePenalty，默认 120mm，等于一个多行距）。
//   ③ 姿态转角：两段法向夹角折算成枪位弧长（angleWeight × standoff × 弧度），
//      让"平面 → 斜面 58°"这种换姿态的转移排在更后面。
//
// 求解策略
// --------
//   N <= 14 段：位压 DP（Held-Karp 变体，状态 = 已喷集合 + 当前段朝向），开放式
//               起点，精确最优；
//   N >  14 段：多起点最近邻 + 2-opt（子链反转）+ Or-opt（单段移位），
//               方向由 2 状态 Viterbi 精确给出。
//   两种结果都与原顺序比一次，**绝不会比原顺序更差**。
#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <queue>
#include <limits>
#include <cstdio>
#include "vecmath.h"
#include "coverage.h"
#include "pathgen.h"

struct OrderOpts {
    int    mode = 2;        // 0=关闭（保持输入顺序） 1=最近邻+2opt  2=自动（小规模精确 DP）
    int    metric = 1;      // 0=直线距离  1=掩膜绕行最短路（默认，绕孔/沿边）
    double zonePenalty = 120.0;   // 换区惩罚 mm
    double angleWeight = 0.25;    // 姿态转角折算权重（× standoff × 弧度）
    double standoff = 200.0;      // 喷距 mm（只用于②③的量纲换算）
    double startU = std::numeric_limits<double>::quiet_NaN();  // 固定起点（主平面坐标），NaN=自由
    double startV = std::numeric_limits<double>::quiet_NaN();
    bool   verbose = true;
};

// 掩膜上的单源距离场（8 连通 Dijkstra，单位 mm）
// 一次求解 = 从某个出口端点出发，到掩膜内任意格的最短行程；
// 之后查任意入口端点都是 O(1)，比逐对跑 A* 省一个数量级。
struct MaskField {
    const Mask* mk = nullptr;
    std::vector<float> d;

    void bind(const Mask& m) {
        mk = &m;
        d.assign((size_t)m.W * m.H, 1e18f);
    }
    // (u,v) 附近最近的掩膜格；找不到（离工件太远）返回 false
    bool nearestCell(double u, double v, int& oi, int& oj) const {
        int i0 = (int)std::floor((u - mk->umin) / mk->cell);
        int j0 = (int)std::floor((v - mk->vmin) / mk->cell);
        if (mk->inside(i0, j0) && mk->at(i0, j0)) { oi = i0; oj = j0; return true; }
        for (int r = 1; r <= 12; ++r)
            for (int dj = -r; dj <= r; ++dj)
                for (int di = -r; di <= r; ++di) {
                    if (std::max(std::abs(di), std::abs(dj)) != r) continue;
                    int ni = i0 + di, nj = j0 + dj;
                    if (mk->inside(ni, nj) && mk->at(ni, nj)) { oi = ni; oj = nj; return true; }
                }
        return false;
    }
    void runFrom(int si, int sj) {
        const int W = mk->W, H = mk->H;
        std::fill(d.begin(), d.end(), 1e18f);
        using QF = std::pair<float, int>;
        std::priority_queue<QF, std::vector<QF>, std::greater<QF>> q;
        d[(size_t)sj * W + si] = 0.0f;
        q.push({0.0f, sj * W + si});
        const float a = (float)mk->cell, b = a * 1.41421356f;
        const int di[8] = {1, -1, 0, 0, 1, 1, -1, -1};
        const int dj[8] = {0, 0, 1, -1, 1, -1, 1, -1};
        const float dc[8] = {a, a, a, a, b, b, b, b};
        while (!q.empty()) {
            auto [du, k] = q.top();
            q.pop();
            if (du > d[(size_t)k] + 1e-3f) continue;
            const int ci = k % W, cj = k / W;
            for (int n = 0; n < 8; ++n) {
                const int ni = ci + di[n], nj = cj + dj[n];
                if (!mk->inside(ni, nj) || !mk->at(ni, nj)) continue;
                const size_t nk = (size_t)nj * W + ni;
                const float nd = du + dc[n];
                if (nd < d[nk]) { d[nk] = nd; q.push({nd, (int)nk}); }
            }
        }
    }
    float atUV(double u, double v) const {
        int i = (int)std::floor((u - mk->umin) / mk->cell);
        int j = (int)std::floor((v - mk->vmin) / mk->cell);
        if (!mk->inside(i, j)) return 1e18f;
        if (!mk->at(i, j)) {                       // 端点落在孔/边界外：就近归位
            int ni, nj;
            if (!nearestCell(u, v, ni, nj)) return 1e18f;
            i = ni; j = nj;
        }
        return d[(size_t)j * mk->W + i];
    }
};

// 优化结果
struct RunOrder {
    std::vector<int>  order;      // 新顺序（runs 下标）
    std::vector<char> rev;        // 该段是否反向执行
    size_t            segCount = 0;     // 参与优化的段数
    double            costBefore = 0;   // 原顺序的空移代价
    double            costAfter = 0;    // 优化后
    bool              changed = false;
    bool              exact = false;    // 是否精确解（DP）
};

// 段顺序问题的代价表
struct OrderProblem {
    int N = 0;
    std::vector<double> trans;    // (2N)x(2N)：出口端口 → 入口端口 的行程 mm
    std::vector<int>    zone;     // 每段的区域编号（0 平面 / 1 孔周 / 2 斜面 / 3 其它）
    std::vector<double> angPen;   // N*N：姿态转角折算的等效行程 mm
    double zonePenalty = 120.0;

    // 段 k 的端口：A=front → 2k，B=back → 2k+1
    // dir=0 正序：进 A 出 B；dir=1 反序：进 B 出 A
    inline int outPort(int i, int d) const { return 2 * i + (d ? 0 : 1); }
    inline int inPort (int j, int e) const { return 2 * j + (e ? 1 : 0); }
    inline double cost(int i, int di, int j, int dj) const {
        const size_t P = (size_t)(2 * N);
        double c = trans[(size_t)outPort(i, di) * P + (size_t)inPort(j, dj)];
        c += angPen[(size_t)i * (size_t)N + (size_t)j];
        if (zone[(size_t)i] != zone[(size_t)j]) c += zonePenalty;
        return c;
    }
};

// 给定段顺序，用 2 状态 Viterbi 求最优方向与总代价（起点方向自由）
inline double evalSeq(const OrderProblem& Q, const std::vector<int>& order,
                      std::vector<char>* revOut = nullptr) {
    const int n = (int)order.size();
    if (revOut) revOut->assign((size_t)n, 0);
    if (n < 2) return 0.0;
    std::vector<char> src0((size_t)n, 0), src1((size_t)n, 0);
    double c0 = 0.0, c1 = 0.0;          // 上一段以 正/反 结束的最小代价
    for (int k = 1; k < n; ++k) {
        const int i = order[(size_t)k - 1], j = order[(size_t)k];
        const double a00 = c0 + Q.cost(i, 0, j, 0);
        const double a10 = c1 + Q.cost(i, 1, j, 0);
        const double a01 = c0 + Q.cost(i, 0, j, 1);
        const double a11 = c1 + Q.cost(i, 1, j, 1);
        const double n0 = std::min(a00, a10), n1 = std::min(a01, a11);
        src0[(size_t)k] = (a10 < a00) ? 1 : 0;
        src1[(size_t)k] = (a11 < a01) ? 1 : 0;
        c0 = n0; c1 = n1;
    }
    const char lastDir = (c1 < c0) ? 1 : 0;
    if (revOut) {
        char d = lastDir;
        for (int k = n - 1; k >= 1; --k) {
            (*revOut)[(size_t)k] = d;
            d = d ? src1[(size_t)k] : src0[(size_t)k];
        }
        (*revOut)[0] = d;
    }
    return std::min(c0, c1);
}

// 位压 DP：开放式最短路径（起点自由或指定起点、每段方向自由）
// startCost: 长度 2N，从"机器人起始接近点"到各端口的行程；nullptr = 起点自由
inline double solveExact(const OrderProblem& Q, std::vector<int>& order,
                         std::vector<char>& rev, const std::vector<double>* startCost = nullptr) {
    const int N = Q.N, P = 2 * N;
    const int FULL = 1 << N;
    const double INF = 1e30;
    std::vector<float> f((size_t)FULL * P, (float)INF);
    std::vector<uint8_t> pre((size_t)FULL * P, 255);

    // 首段：起点自由时为 0；指定起点时为首段入口到起点的行程
    for (int i = 0; i < N; ++i) {
        const double s0 = startCost ? (*startCost)[(size_t)2 * i] : 0.0;
        const double s1 = startCost ? (*startCost)[(size_t)2 * i + 1] : 0.0;
        f[(size_t)(1 << i) * P + 2 * i] = (float)s0;
        f[(size_t)(1 << i) * P + 2 * i + 1] = (float)s1;
    }
    for (int mask = 1; mask < FULL; ++mask) {
        const size_t base = (size_t)mask * P;
        for (int s = 0; s < P; ++s) {
            const float cur = f[base + s];
            if (cur >= INF * 0.5f) continue;
            const int i = s >> 1, di = s & 1;
            for (int j = 0; j < N; ++j) {
                if (mask & (1 << j)) continue;
                const int nm = mask | (1 << j);
                const size_t nbase = (size_t)nm * P;
                for (int dj = 0; dj < 2; ++dj) {
                    const float nv = cur + (float)Q.cost(i, di, j, dj);
                    if (nv < f[nbase + 2 * j + dj]) {
                        f[nbase + 2 * j + dj] = nv;
                        pre[nbase + 2 * j + dj] = (uint8_t)s;
                    }
                }
            }
        }
    }
    // 收尾：全访问集合里最优点
    const size_t base = (size_t)(FULL - 1) * P;
    int bestS = 0;
    double bestC = INF;
    for (int s = 0; s < P; ++s) {
        const double c = f[base + s];
        if (c < bestC) { bestC = c; bestS = s; }
    }
    if (bestC >= INF * 0.5) { order.clear(); rev.clear(); return INF; }

    std::vector<int> rst;
    std::vector<char> rdir;
    int mask = FULL - 1, s = bestS;
    while (mask) {
        const int i = s >> 1, di = s & 1;
        rst.push_back(i);
        rdir.push_back((char)di);
        const uint8_t ps = pre[(size_t)mask * P + s];
        mask &= ~(1 << i);
        if (ps == 255) break;
        s = ps;
    }
    std::reverse(rst.begin(), rst.end());
    std::reverse(rdir.begin(), rdir.end());
    order = std::move(rst);
    rev = std::move(rdir);
    return bestC;
}

// 最近邻 + 2-opt + Or-opt（大 N 用）
inline double solveHeuristic(const OrderProblem& Q, std::vector<int>& order,
                             std::vector<char>& rev, const std::vector<double>* startCost = nullptr) {
    const int N = Q.N;
    const double INF = 1e30;
    std::vector<int> bestOrd;
    double bestC = INF;

    // 顺序代价 = 段间空移 + （可选）从起点进首段的空移
    auto evalWithStart = [&](const std::vector<int>& o, std::vector<char>& r2) {
        double c = evalSeq(Q, o, &r2);
        if (startCost && !o.empty())
            c += (*startCost)[(size_t)2 * o.front() + (r2.empty() ? 0 : (size_t)r2[0])];
        return c;
    };

    // 多起点 NN：每个段 × 每个朝向都当一次首段
    for (int s0 = 0; s0 < 2 * N; ++s0) {
        std::vector<char> used((size_t)N, 0);
        std::vector<int> cur;
        std::vector<char> dir;
        int i = s0 >> 1, di = s0 & 1;
        used[(size_t)i] = 1;
        cur.push_back(i); dir.push_back((char)di);
        for (int step = 1; step < N; ++step) {
            double bd = INF;
            int bj = -1; char be = 0;
            for (int j = 0; j < N; ++j) {
                if (used[(size_t)j]) continue;
                for (int e = 0; e < 2; ++e) {
                    const double c = Q.cost(i, di, j, e);
                    if (c < bd) { bd = c; bj = j; be = (char)e; }
                }
            }
            if (bj < 0) break;
            used[(size_t)bj] = 1;
            cur.push_back(bj); dir.push_back(be);
            i = bj; di = be;
        }
        // 顺序定了，方向交给 Viterbi 重新求（NN 里的方向是逐步贪心的，未必最优）
        std::vector<char> r2;
        const double c = evalWithStart(cur, r2);
        if (c < bestC) { bestC = c; bestOrd = cur; rev = r2; }
    }
    if (bestOrd.empty()) return INF;

    // 2-opt：反转子链 / Or-opt：单段移位
    std::vector<int> cand = bestOrd;
    bool improved = true;
    int rounds = 0;
    while (improved && rounds++ < 60) {
        improved = false;
        // 2-opt
        for (int a = 0; a + 1 < N; ++a) {
            for (int b = a + 1; b < N; ++b) {
                std::vector<int> t = cand;
                std::reverse(t.begin() + a, t.begin() + b + 1);
                std::vector<char> r2;
                const double c = evalWithStart(t, r2);
                if (c + 1e-6 < bestC) { bestC = c; cand.swap(t); rev = r2; improved = true; }
            }
        }
        // Or-opt（单段前移/后移一位的邻域）
        for (int a = 0; a < N; ++a) {
            for (int b = 0; b < N; ++b) {
                if (a == b) continue;
                std::vector<int> t = cand;
                const int v = t[(size_t)a];
                t.erase(t.begin() + a);
                t.insert(t.begin() + b, v);
                std::vector<char> r2;
                const double c = evalWithStart(t, r2);
                if (c + 1e-6 < bestC) { bestC = c; cand.swap(t); rev = r2; improved = true; }
            }
        }
    }
    order = std::move(cand);
    return bestC;
}

// 主入口：优化 runs 的执行顺序与方向（不修改 runs，只给出方案）
inline RunOrder optimizeRunOrder(const std::vector<Run>& runs, const Mask& mkTrav,
                                 const Plane& pl, const Vec3& uAx, const Vec3& vAx,
                                 const OrderOpts& O) {
    RunOrder res;
    const int N = (int)runs.size();
    res.segCount = (size_t)N;
    if (N < 2 || O.mode == 0) {
        res.order.resize((size_t)N);
        res.rev.assign((size_t)N, 0);
        for (int i = 0; i < N; ++i) res.order[(size_t)i] = i;
        return res;
    }

    OrderProblem Q;
    Q.N = N;
    Q.zonePenalty = O.zonePenalty;
    Q.zone.resize((size_t)N);
    Q.angPen.assign((size_t)N * (size_t)N, 0.0);
    for (int i = 0; i < N; ++i) {
        const std::string& z = runs[(size_t)i].zone;
        Q.zone[(size_t)i] = (z == "plane") ? 0 : (z == "ring") ? 1 : (z == "slope") ? 2 : 3;
    }
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            const double c = std::max(-1.0, std::min(1.0, dot(runs[(size_t)i].n, runs[(size_t)j].n)));
            Q.angPen[(size_t)i * (size_t)N + (size_t)j] =
                O.angleWeight * O.standoff * std::acos(c);
        }

    // ---- 端点与面内投影 ----
    // 优先用段自带的表面点投影 ppu/ppv；缺失时退回枪位在主平面的投影。
    auto uvOf = [&](const Run& r, size_t k) -> std::pair<double, double> {
        if (k < r.ppu.size() && k < r.ppv.size()) return {r.ppu[k], r.ppv[k]};
        const Vec3 d = r.p[k] - pl.c;
        return {dot(d, uAx), dot(d, vAx)};
    };
    std::vector<Vec3> portP((size_t)2 * N);
    std::vector<double> portU((size_t)2 * N), portV((size_t)2 * N);
    for (int i = 0; i < N; ++i) {
        const Run& r = runs[(size_t)i];
        const size_t last = r.p.size() - 1;
        portP[(size_t)2 * i] = r.p.front();
        auto uv0 = uvOf(r, 0);
        portU[(size_t)2 * i] = uv0.first; portV[(size_t)2 * i] = uv0.second;
        portP[(size_t)2 * i + 1] = r.p.back();
        auto uv1 = uvOf(r, last);
        portU[(size_t)2 * i + 1] = uv1.first; portV[(size_t)2 * i + 1] = uv1.second;
    }

    // ---- 行程代价表 ----
    // 每个出口端口跑一次距离场，再查所有入口：O(2N) 次 Dijkstra，而不是 O(4N²) 次 A*。
    const size_t P = (size_t)2 * N;
    Q.trans.assign(P * P, 0.0);
    const bool useMask = (O.metric == 1) && mkTrav.W > 1 && mkTrav.H > 1
                         && (double)mkTrav.W * (double)mkTrav.H * (double)P < 4.0e8;
    if (useMask) {
        MaskField f;
        f.bind(mkTrav);
        for (size_t a = 0; a < P; ++a) {
            int si, sj;
            if (!f.nearestCell(portU[a], portV[a], si, sj)) {
                // 端点离行程区太远：这一行退回直线距离
                for (size_t b = 0; b < P; ++b)
                    Q.trans[a * P + b] = norm(portP[b] - portP[a]);
                continue;
            }
            f.runFrom(si, sj);
            for (size_t b = 0; b < P; ++b) {
                float d = f.atUV(portU[b], portV[b]);
                Q.trans[a * P + b] = (d < 1e17f) ? (double)d : norm(portP[b] - portP[a]);
            }
            // 闭环段（孔周整圆）的进入点可以自由旋转：取环上离出口最近的那一点
            for (int j = 0; j < N; ++j) {
                const Run& rj = runs[(size_t)j];
                if (!rj.closed) continue;
                double bd = Q.trans[a * P + (size_t)2 * j];
                for (size_t k = 0; k < rj.p.size(); ++k) {
                    auto uv = uvOf(rj, k);
                    float d = f.atUV(uv.first, uv.second);
                    if (d < 1e17f && (double)d < bd) bd = (double)d;
                }
                Q.trans[a * P + (size_t)2 * j] = bd;
                Q.trans[a * P + (size_t)2 * j + 1] = bd;
            }
        }
    } else {
        for (size_t a = 0; a < P; ++a)
            for (size_t b = 0; b < P; ++b)
                Q.trans[a * P + b] = norm(portP[b] - portP[a]);
        // 闭环段：入口取离出口最近的环上点（直线口径）
        for (size_t a = 0; a < P; ++a)
            for (int j = 0; j < N; ++j) {
                const Run& rj = runs[(size_t)j];
                if (!rj.closed) continue;
                double bd = Q.trans[a * P + (size_t)2 * j];
                for (const Vec3& p : rj.p) bd = std::min(bd, norm(p - portP[a]));
                Q.trans[a * P + (size_t)2 * j] = bd;
                Q.trans[a * P + (size_t)2 * j + 1] = bd;
            }
    }
    // ---- 指定起点（可选）----
    // 机器人从工件外的某个接近点出发时，首段的入口端应按"离该点最近"来选。
    // 用同一套掩膜距离场算出 startCost[2N]，作为求解器首段的初值代价。
    std::vector<double> startCost;
    if (std::isfinite(O.startU) && std::isfinite(O.startV)) {
        startCost.assign(P, 0.0);
        bool ok = false;
        if (useMask) {
            MaskField fs;
            fs.bind(mkTrav);
            int si, sj;
            if (fs.nearestCell(O.startU, O.startV, si, sj)) {
                fs.runFrom(si, sj);
                for (size_t b = 0; b < P; ++b) {
                    const float d = fs.atUV(portU[b], portV[b]);
                    startCost[b] = (d < 1e17f) ? (double)d : 1e6;
                }
                ok = true;
            }
        }
        if (!ok)
            for (size_t b = 0; b < P; ++b)
                startCost[b] = std::hypot(portU[b] - O.startU, portV[b] - O.startV);
    }
    const std::vector<double>* startP = startCost.empty() ? nullptr : &startCost;

    // ---- 原顺序代价 ----
    std::vector<int> orig((size_t)N);
    for (int i = 0; i < N; ++i) orig[(size_t)i] = i;
    double cBefore = evalSeq(Q, orig, nullptr);
    if (startP) {
        // 原顺序从它自己的第一段进入，入口端取更近的一端
        const double d0 = (*startP)[(size_t)2 * orig.front()];
        const double d1 = (*startP)[(size_t)2 * orig.front() + 1];
        cBefore += std::min(d0, d1);
    }

    // ---- 求解 ----
    std::vector<int> best;
    std::vector<char> bestRev;
    double bestC;
    bool exact = false;
    if (O.mode == 2 && N <= 14) {
        bestC = solveExact(Q, best, bestRev, startP);
        exact = true;
        if (best.empty()) { bestC = solveHeuristic(Q, best, bestRev, startP); exact = false; }
    } else {
        bestC = solveHeuristic(Q, best, bestRev, startP);
    }
    if (best.size() != (size_t)N) {          // 求解异常 → 保持原顺序
        best = orig;
        bestRev.assign((size_t)N, 0);
        bestC = cBefore;
    }

    res.order = best;
    res.rev = bestRev;
    res.costBefore = cBefore;
    res.costAfter = bestC;
    res.exact = exact;
    res.changed = (bestC + 1e-6 < cBefore);
    if (!res.changed) {                       // 优化无收益 → 原样返回，避免引入无谓扰动
        res.order = orig;
        res.rev.assign((size_t)N, 0);
        res.costAfter = cBefore;
    }
    if (O.verbose) {
        fprintf(stderr, "[order] 段顺序优化: %d 段，空移 %.0f → %.0f mm（%+.1f%%，%s%s）\n",
                N, cBefore, res.costAfter,
                100.0 * (res.costAfter - cBefore) / std::max(1.0, cBefore),
                exact ? "精确DP" : "NN+2opt",
                res.changed ? "" : "，无改善→保持原序");
    }
    return res;
}

// 按优化结果重排 runs（就地反转需要反向的段）
inline void applyRunOrder(std::vector<Run>& runs, const RunOrder& o) {
    if (o.order.size() != runs.size()) return;
    std::vector<Run> out;
    out.reserve(runs.size());
    for (size_t k = 0; k < o.order.size(); ++k) {
        Run r = std::move(runs[(size_t)o.order[k]]);
        if (k < o.rev.size() && o.rev[k] && !r.closed) {
            std::reverse(r.p.begin(), r.p.end());
            std::reverse(r.uu.begin(), r.uu.end());
            std::reverse(r.vv.begin(), r.vv.end());
            std::reverse(r.hh.begin(), r.hh.end());
            std::reverse(r.ppu.begin(), r.ppu.end());
            std::reverse(r.ppv.begin(), r.ppv.end());
        }
        out.push_back(std::move(r));
    }
    runs.swap(out);
}
