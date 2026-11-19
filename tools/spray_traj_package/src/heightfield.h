// heightfield.h - 表面高度场：平面残差的网格化、平滑与采样（轨迹贴合表面起伏）
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "vecmath.h"
#include "plane_fit.h"

struct HeightField {
    int W = 0, H = 0;
    double cell = 5.0, umin = 0, vmin = 0;
    std::vector<float> h;   // 平均残差
    std::vector<uint8_t> known;

    void build(const std::vector<Vec3>& pts, const std::vector<int>& idx,
               const Plane& pl, const Vec3& uAxis, const Vec3& vAxis, double cellMM) {
        cell = cellMM;
        double umin_ = 1e18, vmin_ = 1e18, umax_ = -1e18, vmax_ = -1e18;
        for (int i : idx) {
            Vec3 d = pts[i] - pl.c;
            double u = dot(d, uAxis), v = dot(d, vAxis);
            umin_ = std::min(umin_, u); umax_ = std::max(umax_, u);
            vmin_ = std::min(vmin_, v); vmax_ = std::max(vmax_, v);
        }
        if (umax_ < umin_) { W = H = 1; umin = vmin = 0; h.assign(1, 0.f); known.assign(1, 0); return; }
        W = std::max(1, (int)std::ceil((umax_ - umin_) / cell));
        H = std::max(1, (int)std::ceil((vmax_ - vmin_) / cell));
        umin = umin_; vmin = vmin_;
        std::vector<double> sum((size_t)W * H, 0);
        std::vector<int> cnt((size_t)W * H, 0);
        for (int i : idx) {
            Vec3 d = pts[i] - pl.c;
            int a = (int)((dot(d, uAxis) - umin) / cell);
            int b = (int)((dot(d, vAxis) - vmin) / cell);
            if (a < 0 || a >= W || b < 0 || b >= H) continue;
            sum[(size_t)b * W + a] += pl.d(pts[i]);
            cnt[(size_t)b * W + a]++;
        }
        h.assign((size_t)W * H, 0.f);
        known.assign((size_t)W * H, 0);
        for (size_t k = 0; k < h.size(); ++k)
            if (cnt[k]) { h[k] = (float)(sum[k] / cnt[k]); known[k] = 1; }
        // 空格子（孔洞/缺口）从邻域扩散填充, 避免轨迹过孔时高度采样跳变
        for (int it = 0; it < 400; ++it) {
            bool changed = false;
            for (int j = 0; j < H; ++j)
                for (int i = 0; i < W; ++i) {
                    size_t k = (size_t)j * W + i;
                    if (known[k]) continue;
                    double s = 0; int n = 0;
                    const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
                    for (int t = 0; t < 4; ++t) {
                        int ni = i + di[t], nj = j + dj[t];
                        if (ni < 0 || ni >= W || nj < 0 || nj >= H) continue;
                        size_t nk = (size_t)nj * W + ni;
                        if (known[nk]) { s += h[nk]; n++; }
                    }
                    if (n) { h[k] = (float)(s / n); changed = true; }
                }
            for (size_t k = 0; k < h.size(); ++k)
                if (!known[k] && h[k] != 0.f) known[k] = 1;
            if (!changed) break;
        }
        boxBlur(); boxBlur();
    }

    void boxBlur() {
        std::vector<float> out((size_t)W * H, 0.f);
        for (int j = 0; j < H; ++j)
            for (int i = 0; i < W; ++i) {
                double s = 0; int n = 0;
                for (int dj = -1; dj <= 1; ++dj)
                    for (int di = -1; di <= 1; ++di) {
                        int ni = i + di, nj = j + dj;
                        if (ni < 0 || ni >= W || nj < 0 || nj >= H) continue;
                        s += h[(size_t)nj * W + ni]; n++;
                    }
                out[(size_t)j * W + i] = (float)(s / n);
            }
        h.swap(out);
    }

    // 双线性插值采样
    double sample(double u, double v) const {
        if (W < 2 || H < 2) return 0;
        double x = (u - umin) / cell, y = (v - vmin) / cell;
        int i = (int)std::floor(x), j = (int)std::floor(y);
        if (i < 0 || i >= W - 1 || j < 0 || j >= H - 1) return 0;
        double fx = x - i, fy = y - j;
        double h00 = h[(size_t)j * W + i],       h10 = h[(size_t)j * W + i + 1];
        double h01 = h[(size_t)(j + 1) * W + i], h11 = h[(size_t)(j + 1) * W + i + 1];
        return h00 * (1 - fx) * (1 - fy) + h10 * fx * (1 - fy) +
               h01 * (1 - fx) * fy + h11 * fx * fy;
    }

    // 该处是否有真实点云数据（非扩散填充）
    bool knownAt(double u, double v) const {
        if (W < 1 || H < 1) return false;
        int i = (int)((u - umin) / cell), j = (int)((v - vmin) / cell);
        if (i < 0 || i >= W || j < 0 || j >= H) return false;
        return known[(size_t)j * W + i] != 0;
    }
};
