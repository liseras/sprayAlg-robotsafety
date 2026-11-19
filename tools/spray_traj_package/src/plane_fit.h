// plane_fit.h - 稳健平面拟合（截断最小二乘，自动剔除悬挂杆等杂散点）
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include "vecmath.h"

struct Plane {
    Vec3 n;         // 单位法向
    Vec3 c;         // 平面上参考点（质心）
    double d(const Vec3& p) const { return dot(p - c, n); }  // 点到平面有符号距离
};

// 总体最小二乘平面拟合（PCA 最小特征向量方向为法向）
inline Plane fitPlaneTLS(const std::vector<Vec3>& pts, const std::vector<int>& idx) {
    Vec3 c(0, 0, 0);
    for (int i : idx) c = c + pts[i];
    c = c / (double)idx.size();
    double C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i : idx) {
        Vec3 d = pts[i] - c;
        C[0][0] += d.x*d.x; C[0][1] += d.x*d.y; C[0][2] += d.x*d.z;
        C[1][1] += d.y*d.y; C[1][2] += d.y*d.z; C[2][2] += d.z*d.z;
    }
    C[1][0] = C[0][1]; C[2][0] = C[0][2]; C[2][1] = C[1][2];
    double eig[3]; Vec3 vec[3];
    jacobiEigen3(C, eig, vec);
    Plane pl; pl.n = vec[0]; pl.c = c;
    return pl;
}

// RANSAC 平面拟合：在主面上取最大一致集（阈值 thresh），
// 避免卷边/折边把参考平面拉歪（平面一旦倾斜, 卷边残差被吸收, 枪就不再垂直主面）
inline Plane fitPlaneRANSAC(const std::vector<Vec3>& pts, const std::vector<int>& idx,
                            double thresh = 6.0, int iters = 300, unsigned seed = 12345) {
    auto rnd = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return (seed >> 16) & 0x7fff;
    };
    const size_t n = idx.size();
    if (n < 10) return fitPlaneTLS(pts, idx);
    // 打分时抽样（每 3 个取 1 个），加速
    std::vector<int> sample;
    sample.reserve(n / 3 + 1);
    for (size_t k = 0; k < n; k += 3) sample.push_back(idx[k]);

    Plane best;
    size_t bestCount = 0;
    for (int it = 0; it < iters; ++it) {
        int a = (int)(rnd() % n), b = (int)(rnd() % n), c = (int)(rnd() % n);
        if (a == b || b == c || a == c) continue;
        Vec3 p0 = pts[idx[a]], p1 = pts[idx[b]], p2 = pts[idx[c]];
        Vec3 nn = cross(p1 - p0, p2 - p0);
        double L = norm(nn);
        if (L < 1e-9) continue;
        nn = nn / L;
        size_t cnt = 0;
        for (int k : sample)
            if (std::fabs(dot(pts[k] - p0, nn)) <= thresh) cnt++;
        if (cnt > bestCount) {
            bestCount = cnt;
            best.n = nn; best.c = p0;
        }
    }
    if (bestCount == 0) return fitPlaneTLS(pts, idx);
    // 用最佳一致集精修两轮
    std::vector<int> in;
    for (int rep = 0; rep < 2; ++rep) {
        in.clear();
        for (int i : idx)
            if (std::fabs(best.d(pts[i])) <= thresh) in.push_back(i);
        if (in.size() < 10) break;
        best = fitPlaneTLS(pts, in);
    }
    // 方向统一：与整体法向同侧（朝向扫描仪）
    Vec3 ref = fitPlaneTLS(pts, idx).n;
    if (dot(best.n, ref) < 0) best.n = best.n * -1.0;
    return best;
}

// 截断迭代拟合：阈值在首轮确定（基于 MAD），避免迭代过紧削掉缓曲面起伏
inline Plane fitPlaneRobust(const std::vector<Vec3>& pts, int iters = 4,
                            double absFloor = 2.0, std::vector<int>* inlierIdx = nullptr) {
    std::vector<int> idx(pts.size());
    for (size_t i = 0; i < pts.size(); ++i) idx[i] = (int)i;
    Plane pl = fitPlaneTLS(pts, idx);
    // 首轮确定剔除阈值
    std::vector<double> res0(idx.size());
    for (size_t k = 0; k < idx.size(); ++k) res0[k] = std::fabs(pl.d(pts[idx[k]]));
    std::vector<double> s0 = res0;
    std::sort(s0.begin(), s0.end());
    double med0 = s0[s0.size() / 2];
    std::vector<double> dev0;
    dev0.reserve(s0.size());
    for (double r : s0) dev0.push_back(std::fabs(r - med0));
    std::sort(dev0.begin(), dev0.end());
    double th = std::max(3.0 * dev0[dev0.size() / 2] * 1.4826, absFloor);
    for (int it = 0; it < iters; ++it) {
        std::vector<double> res(idx.size());
        for (size_t k = 0; k < idx.size(); ++k) res[k] = std::fabs(pl.d(pts[idx[k]]));
        std::vector<int> keep;
        keep.reserve(idx.size());
        for (size_t k = 0; k < idx.size(); ++k)
            if (res[k] <= th) keep.push_back(idx[k]);
        if ((int)keep.size() < 100 || keep.size() == idx.size()) { idx = keep; break; }
        idx = keep;
        pl = fitPlaneTLS(pts, idx);
    }
    if (inlierIdx) *inlierIdx = idx;
    return pl;
}
