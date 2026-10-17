// spray_traj.h —— 喷涂轨迹生成核心库 C ABI
// 平面/曲面点云 -> 连续S形喷涂轨迹（等距、全程喷涂）
#pragma once

#ifdef _WIN32
  #ifdef SPRAYTRAJ_EXPORTS
    #define SPRAYTRAJ_API __declspec(dllexport)
  #else
    #define SPRAYTRAJ_API __declspec(dllimport)
  #endif
#else
  #define SPRAYTRAJ_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

// 工艺参数（全部可调，WPF 侧传入）
typedef struct SprayTrajParams {
    double scanRowSpacing;   // 行距 mm
    double pointSpacing;     // 点距 mm
    double standoff;         // 离面距离 mm（轨迹点沿枪轴抬升）
    double speed;            // 行走速度 mm/s
    double turnRadiusFactor; // 折返圆弧半径 = 该系数 × 行距
    double flatnessTol;      // 平面重建内点阈值 mm（点到平面距离 < 此值即归入该平面）
    double curvatureTol;     // 曲面升级阈值 mm（拟合残差 < 此值则浅曲面当作平面处理）
    int    sorK;             // 统计离群滤波近邻数
    double sorZ;             // 统计离群滤波 Z 倍数
    int    normalFlip;       // 法向翻转：0=自动判向，1=强制翻转（法向取反）
} SprayTrajParams;

// 默认参数
SPRAYTRAJ_API SprayTrajParams spray_traj_default_params(void);

// 单个轨迹点
typedef struct SprayTrajPoint {
    double x, y, z;    // 位置 mm（基座系，已按离面距离抬升）
    double i, j, k;    // 喷枪轴线单位向量（指向表面）
    double speed;      // mm/s
    int    spray;      // 1=喷涂（全程恒为1）
} SprayTrajPoint;

// 单个重建平面（供显示/调试）
typedef struct SprayTrajPlane {
    double nx, ny, nz;   // 法向（定向后）
    double cx, cy, cz;   // 质心
    int    count;        // 该面内点数
    // 该面在自身局部基(u,v)下的轴对齐包围盒（供 WPF 直接画矩形，不做 PCA 近似）
    double ux, uy, uz;   // 局部基 u（主轴方向，单位向量）
    double vx, vy, vz;   // 局部基 v
    double minU, maxU;   // 内点在 u 方向的范围
    double minV, maxV;   // 内点在 v 方向的范围
} SprayTrajPlane;

// 调用结果
typedef struct SprayTrajResult {
    SprayTrajPoint* points;   // 轨迹点数组（调用方用 spray_traj_free 释放）
    int             count;    // 点数
    // 自检
    double          plane_rms;      // 平面拟合残差 RMS mm
    double          normal_x, normal_y, normal_z; // 重建面法向（定向后，指向外表面）
    double          plane_cx, plane_cy, plane_cz; // 重建面质心（内点质心）
    double          seg_len_mean;   // 相邻点距均值
    double          seg_len_max;    // 相邻点距最大值
    int             rows;           // 行数
    int             input_points;   // 输入点数
    int             removed_points; // 去噪剔除点数
    int             inlier_points;  // 平面内点数量（重建平面覆盖的点）
    int             plane_count;    // 分割出的平面数
    SprayTrajPlane  planes[8];      // 各分割平面（最多 8 个）
    const char*     error;          // 失败原因（成功为 NULL）
} SprayTrajResult;

// 主入口：生成连续S形轨迹
// 输入：点云数组（xyz），点数，参数
// 返回：轨迹结果；失败时 error 非空、count=0
SPRAYTRAJ_API SprayTrajResult spray_traj_generate(
    const double* xyz, int n, const SprayTrajParams* params);

// 释放结果（内部分配的点数组与错误串）
SPRAYTRAJ_API void spray_traj_free(SprayTrajResult* res);

#ifdef __cplusplus
}
#endif
