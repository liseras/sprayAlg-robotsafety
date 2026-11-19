// coverage.h - 平面投影占据栅格、连通域、距离变换腐蚀（确定喷涂覆盖区域）
#pragma once
#include <vector>
#include <cstdint>
#include <cmath>
#include <queue>
#include <algorithm>
#include "vecmath.h"
#include "plane_fit.h"

struct Mask {
    int W = 0, H = 0;
    double cell = 3.0;   // mm
    double umin = 0, vmin = 0;
    std::vector<uint8_t> m;
    uint8_t& at(int i, int j) { return m[(size_t)j * W + i]; }
    uint8_t at(int i, int j) const { return m[(size_t)j * W + i]; }
    bool inside(int i, int j) const { return i >= 0 && i < W && j >= 0 && j < H; }
    double u(int i) const { return umin + (i + 0.5) * cell; }
    double v(int j) const { return vmin + (j + 0.5) * cell; }
};

// 把点云投影到平面坐标系 (u,v)，占据栅格化（每格至少 minPts 个点才算覆盖）
inline Mask buildMask(const std::vector<Vec3>& pts, const std::vector<int>& idx,
                      const Plane& pl, const Vec3& uAxis, const Vec3& vAxis,
                      double cell, int minPtsPerCell) {
    std::vector<double> us, vs;
    us.reserve(idx.size()); vs.reserve(idx.size());
    double umin = 1e18, vmin = 1e18, umax = -1e18, vmax = -1e18;
    for (int i : idx) {
        Vec3 d = pts[i] - pl.c;
        double u = dot(d, uAxis), v = dot(d, vAxis);
        us.push_back(u); vs.push_back(v);
        umin = std::min(umin, u); umax = std::max(umax, u);
        vmin = std::min(vmin, v); vmax = std::max(vmax, v);
    }
    Mask mk;
    mk.cell = cell; mk.umin = umin; mk.vmin = vmin;
    mk.W = std::max(1, (int)std::ceil((umax - umin) / cell));
    mk.H = std::max(1, (int)std::ceil((vmax - vmin) / cell));
    mk.m.assign((size_t)mk.W * mk.H, 0);
    std::vector<int> cnt((size_t)mk.W * mk.H, 0);
    for (size_t k = 0; k < us.size(); ++k) {
        int i = (int)((us[k] - umin) / cell);
        int j = (int)((vs[k] - vmin) / cell);
        if (i >= 0 && i < mk.W && j >= 0 && j < mk.H) cnt[(size_t)j * mk.W + i]++;
    }
    for (size_t k = 0; k < cnt.size(); ++k) mk.m[k] = cnt[k] >= minPtsPerCell ? 1 : 0;
    return mk;
}

// 最大 4 连通域（孤立杂散点/悬挂杆会被丢弃）
inline void keepLargestComponent(Mask& mk) {
    std::vector<int> label((size_t)mk.W * mk.H, -1);
    int bestId = -1, bestSize = 0, cur = 0;
    for (int j = 0; j < mk.H; ++j) {
        for (int i = 0; i < mk.W; ++i) {
            if (!mk.at(i, j) || label[(size_t)j * mk.W + i] >= 0) continue;
            int size = 0;
            std::queue<std::pair<int,int>> q;
            q.push({i, j});
            label[(size_t)j * mk.W + i] = cur;
            while (!q.empty()) {
                auto [ci, cj] = q.front(); q.pop();
                size++;
                const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
                for (int k = 0; k < 4; ++k) {
                    int ni = ci + di[k], nj = cj + dj[k];
                    if (mk.inside(ni, nj) && mk.at(ni, nj) &&
                        label[(size_t)nj * mk.W + ni] < 0) {
                        label[(size_t)nj * mk.W + ni] = cur;
                        q.push({ni, nj});
                    }
                }
            }
            if (size > bestSize) { bestSize = size; bestId = cur; }
            cur++;
        }
    }
    for (size_t k = 0; k < mk.m.size(); ++k)
        mk.m[k] = (label[k] == bestId) ? 1 : 0;
}

// chamfer 距离变换：每个栅格到最近背景格的近似距离（单位 mm）
// 注：距离针对"背景"计算，即掩膜格的值为到边界的距离
inline std::vector<float> distanceTransformMM(const Mask& mk) {
    const float INF = 1e18f;
    const float a = mk.cell, b = mk.cell * 1.4142f;  // 正交/对角权重
    std::vector<float> d((size_t)mk.W * mk.H, INF);
    for (int j = 0; j < mk.H; ++j)
        for (int i = 0; i < mk.W; ++i)
            if (!mk.at(i, j)) d[(size_t)j * mk.W + i] = 0;
    auto get = [&](int i, int j) -> float {
        if (!mk.inside(i, j)) return INF;
        return d[(size_t)j * mk.W + i];
    };
    for (int j = 0; j < mk.H; ++j) {
        for (int i = 0; i < mk.W; ++i) {
            float& cur = d[(size_t)j * mk.W + i];
            cur = std::min(cur, get(i - 1, j) + a);
            cur = std::min(cur, get(i, j - 1) + a);
            cur = std::min(cur, get(i - 1, j - 1) + b);
            cur = std::min(cur, get(i + 1, j - 1) + b);
        }
    }
    for (int j = mk.H - 1; j >= 0; --j) {
        for (int i = mk.W - 1; i >= 0; --i) {
            float& cur = d[(size_t)j * mk.W + i];
            cur = std::min(cur, get(i + 1, j) + a);
            cur = std::min(cur, get(i, j + 1) + a);
            cur = std::min(cur, get(i + 1, j + 1) + b);
            cur = std::min(cur, get(i - 1, j + 1) + b);
        }
    }
    return d;
}

// 腐蚀：只保留到边界距离 >= margin mm 的区域
inline void erodeMask(Mask& mk, double marginMM) {
    if (marginMM <= 0) return;
    auto d = distanceTransformMM(mk);
    for (size_t k = 0; k < mk.m.size(); ++k)
        mk.m[k] = (d[k] >= marginMM) ? mk.m[k] : 0;
}

// 填充封闭的背景区域（小孔/槽）：喷涂时忽略, 面积 <= maxAreaMM2 才填。
// 返回被填充格子的标记（用于预览高亮），填充个数见 filledCount。
inline std::vector<uint8_t> fillEnclosedHoles(Mask& mk, double maxAreaMM2, int* filledCount = nullptr) {
    std::vector<uint8_t> filled(mk.m.size(), 0);
    const int W = mk.W, H = mk.H;
    std::vector<uint8_t> visited(mk.m.size(), 0);
    std::vector<int> comp;
    int count = 0;
    for (int j = 0; j < H; ++j) {
        for (int i = 0; i < W; ++i) {
            size_t k0 = (size_t)j * W + i;
            if (mk.at(i, j) || visited[k0]) continue;
            // BFS 收集该背景连通域
            comp.clear();
            std::queue<int> q;
            q.push((int)k0);
            visited[k0] = 1;
            bool touchBorder = false;
            while (!q.empty()) {
                int k = q.front(); q.pop();
                comp.push_back(k);
                int ci = k % W, cj = k / W;
                if (ci == 0 || ci == W - 1 || cj == 0 || cj == H - 1) touchBorder = true;
                const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
                for (int n = 0; n < 4; ++n) {
                    int ni = ci + di[n], nj = cj + dj[n];
                    if (ni < 0 || ni >= W || nj < 0 || nj >= H) continue;
                    size_t nk = (size_t)nj * W + ni;
                    if (!mk.at(ni, nj) && !visited[nk]) {
                        visited[nk] = 1;
                        q.push((int)nk);
                    }
                }
            }
            double areaMM2 = (double)comp.size() * mk.cell * mk.cell;
            if (!touchBorder && areaMM2 <= maxAreaMM2) {
                for (int k : comp) {
                    mk.m[(size_t)k] = 1;
                    filled[(size_t)k] = 1;
                }
                count++;
            }
        }
    }
    if (filledCount) *filledCount = count;
    return filled;
}
