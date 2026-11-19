// regions.h - 区域分割：大孔识别 / 陡峭(斜面)区域识别 / 掩膜形态学
#pragma once
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include "vecmath.h"
#include "coverage.h"
#include "heightfield.h"

static const double PI_R = 3.14159265358979323846;

// ---- 大孔：封闭背景连通域中面积超过阈值的（真窗口，必须绕行） ----
struct HoleInfo {
    double cu = 0, cv = 0;   // 孔中心（平面坐标 mm）
    double rEq = 0;          // 等效半径 sqrt(area/pi)
    double rMax = 0;         // 到孔心最远距离（外接半径）
    double area = 0;
};

inline std::vector<HoleInfo> findOpenHoles(const Mask& mk, double minAreaMM2) {
    std::vector<HoleInfo> out;
    const int W = mk.W, H = mk.H;
    std::vector<uint8_t> vis((size_t)W * H, 0);
    std::vector<int> comp;
    for (int j = 0; j < H; ++j)
        for (int i = 0; i < W; ++i) {
            if (mk.at(i, j) || vis[(size_t)j * W + i]) continue;
            comp.clear();
            std::queue<int> q;
            q.push(j * W + i);
            vis[(size_t)j * W + i] = 1;
            bool border = false;
            while (!q.empty()) {
                int k = q.front(); q.pop();
                comp.push_back(k);
                int ci = k % W, cj = k / W;
                if (ci == 0 || ci == W - 1 || cj == 0 || cj == H - 1) border = true;
                const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
                for (int n = 0; n < 4; ++n) {
                    int ni = ci + di[n], nj = cj + dj[n];
                    if (ni < 0 || ni >= W || nj < 0 || nj >= H) continue;
                    size_t nk = (size_t)nj * W + ni;
                    if (!mk.at(ni, nj) && !vis[nk]) { vis[nk] = 1; q.push((int)nk); }
                }
            }
            double area = (double)comp.size() * mk.cell * mk.cell;
            if (border || area < minAreaMM2) continue;
            HoleInfo h;
            h.area = area;
            double su = 0, sv = 0;
            for (int k : comp) { su += mk.u(k % W); sv += mk.v(k / W); }
            h.cu = su / comp.size();
            h.cv = sv / comp.size();
            h.rEq = std::sqrt(area / PI_R);
            double rmax = 0;
            for (int k : comp)
                rmax = std::max(rmax, std::hypot(mk.u(k % W) - h.cu, mk.v(k / W) - h.cv));
            h.rMax = rmax;
            out.push_back(h);
        }
    return out;
}

// ---- 3x3 多数表决：去掉孤立噪点 / 补单格缺口 ----
inline void maskOpen(std::vector<uint8_t>& m, int W, int H) {
    std::vector<uint8_t> out = m;
    for (int j = 0; j < H; ++j)
        for (int i = 0; i < W; ++i) {
            int cnt = 0;
            for (int dj = -1; dj <= 1; ++dj)
                for (int di = -1; di <= 1; ++di) {
                    int ni = i + di, nj = j + dj;
                    if (ni < 0 || ni >= W || nj < 0 || nj >= H) continue;
                    cnt += m[(size_t)nj * W + ni];
                }
            size_t k = (size_t)j * W + i;
            out[k] = m[k] ? (cnt >= 3) : (cnt >= 6);
        }
    m.swap(out);
}

inline void dilateMask(std::vector<uint8_t>& m, int W, int H, int iters) {
    for (int it = 0; it < iters; ++it) {
        std::vector<uint8_t> out = m;
        for (int j = 0; j < H; ++j)
            for (int i = 0; i < W; ++i) {
                if (m[(size_t)j * W + i]) continue;
                for (int dj = -1; dj <= 1 && !out[(size_t)j * W + i]; ++dj)
                    for (int di = -1; di <= 1; ++di) {
                        int ni = i + di, nj = j + dj;
                        if (ni < 0 || ni >= W || nj < 0 || nj >= H) continue;
                        if (m[(size_t)nj * W + ni]) { out[(size_t)j * W + i] = 1; break; }
                    }
            }
        m.swap(out);
    }
}

inline int countMask(const std::vector<uint8_t>& m) {
    int c = 0;
    for (uint8_t v : m) c += v ? 1 : 0;
    return c;
}

// ---- 陡峭区：粗尺度高度场梯度判定 ----
// 小窗口 PCA 会被冲压波纹(±几mm, 20~40mm 周期)淹没, 必须在大尺度上求梯度:
// 衰减波纹、保留真正的折边/斜面（几十 mm 内抬起几十 mm）
inline std::vector<uint8_t> steepFromHeight(const HeightField& hc, const Mask& mk,
                                            double angleThrDeg) {
    std::vector<uint8_t> out(mk.m.size(), 0);
    double d = std::max(hc.cell, 6.0);
    for (int j = 0; j < mk.H; ++j)
        for (int i = 0; i < mk.W; ++i) {
            if (!mk.at(i, j)) continue;
            double u = mk.u(i), v = mk.v(j);
            // 四个采样点必须都有真实数据——孔内是扩散填充的假高度, 会产生假梯度
            if (!hc.knownAt(u + d, v) || !hc.knownAt(u - d, v) ||
                !hc.knownAt(u, v + d) || !hc.knownAt(u, v - d)) continue;
            double gu = (hc.sample(u + d, v) - hc.sample(u - d, v)) / (2 * d);
            double gv = (hc.sample(u, v + d) - hc.sample(u, v - d)) / (2 * d);
            double ang = std::atan(std::hypot(gu, gv)) * 180.0 / PI_R;
            if (ang > angleThrDeg) out[(size_t)j * mk.W + i] = 1;
        }
    return out;
}

// ---- 斜面区域：拟合自身平面 + 长度/宽度方向 + 自身高度场 ----
struct SlopeZone {
    bool valid = false;
    Plane plane;                 // 斜面自身平面
    Vec3 lenAxis, widAxis;       // 长度方向(沿工件长边) / 宽度方向
    std::vector<int> ptIdx;
    std::vector<double> ss, tt;  // 斜面坐标（长度/宽度方向投影）
    double lenMin = 0, lenMax = 0, widMin = 0, widMax = 0;
    double angleToMain = 0;      // 与主面夹角(deg)
    HeightField hf;
};
