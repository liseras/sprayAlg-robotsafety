// vecmath.h - 基础向量/矩阵/四元数运算（无第三方依赖）
#pragma once
#include <cmath>
#include <vector>
#include <algorithm>

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() {}
    Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
};

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalize(const Vec3& a) {
    double n = norm(a);
    return n > 1e-12 ? a / n : Vec3(0, 0, 0);
}

// 3x3 对称矩阵 Jacobi 特征分解，返回特征值(升序)与对应特征向量（列）
inline void jacobiEigen3(const double A[3][3], double eig[3], Vec3 vec[3]) {
    double a[3][3] = {{A[0][0], A[0][1], A[0][2]},
                      {A[1][0], A[1][1], A[1][2]},
                      {A[2][0], A[2][1], A[2][2]}};
    double V[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int it = 0; it < 200; ++it) {
        double off = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
        if (off < 1e-15) break;
        for (int p = 0; p < 3; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (std::fabs(a[p][q]) < 1e-18) continue;
                double theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
                double t = (theta >= 0 ? 1.0 : -1.0) /
                           (std::fabs(theta) + std::sqrt(theta * theta + 1));
                double c = 1 / std::sqrt(t * t + 1), s = t * c;
                for (int k = 0; k < 3; ++k) {
                    double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k) {
                    double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = c * vkp - s * vkq;
                    V[k][q] = s * vkp + c * vkq;
                }
                a[p][p] -= 0; // 已在对角化中更新
            }
        }
    }
    double vals[3] = {a[0][0], a[1][1], a[2][2]};
    int idx[3] = {0, 1, 2};
    std::sort(idx, idx + 3, [&](int i, int j) { return vals[i] < vals[j]; });
    for (int k = 0; k < 3; ++k) {
        eig[k] = vals[idx[k]];
        vec[k] = normalize(Vec3(V[0][idx[k]], V[1][idx[k]], V[2][idx[k]]));
    }
}

// 由工具坐标系基向量（基坐标系下的 X/Y/Z 轴）计算四元数 (qw,qx,qy,qz)
inline void quatFromBasis(const Vec3& X, const Vec3& Y, const Vec3& Z,
                          double& qw, double& qx, double& qy, double& qz) {
    double m00 = X.x, m01 = Y.x, m02 = Z.x;
    double m10 = X.y, m11 = Y.y, m12 = Z.y;
    double m20 = X.z, m21 = Y.z, m22 = Z.z;
    double tr = m00 + m11 + m22;
    if (tr > 0) {
        double s = std::sqrt(tr + 1.0) * 2;
        qw = 0.25 * s;
        qx = (m21 - m12) / s;
        qy = (m02 - m20) / s;
        qz = (m10 - m01) / s;
    } else if (m00 > m11 && m00 > m22) {
        double s = std::sqrt(1.0 + m00 - m11 - m22) * 2;
        qw = (m21 - m12) / s; qx = 0.25 * s;
        qy = (m01 + m10) / s; qz = (m02 + m20) / s;
    } else if (m11 > m22) {
        double s = std::sqrt(1.0 + m11 - m00 - m22) * 2;
        qw = (m02 - m20) / s; qx = (m01 + m10) / s;
        qy = 0.25 * s;        qz = (m12 + m21) / s;
    } else {
        double s = std::sqrt(1.0 + m22 - m00 - m11) * 2;
        qw = (m10 - m01) / s; qx = (m02 + m20) / s;
        qy = (m12 + m21) / s; qz = 0.25 * s;
    }
    double n = std::sqrt(qw * qw + qx * qx + qy * qy + qz * qz);
    qw /= n; qx /= n; qy /= n; qz /= n;
}
