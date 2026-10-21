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
    double fillHoleSize;     // 孔洞填补阈值 mm：点云间隙/孔洞跨距 < 此值 → 覆盖重建（行连续跨越、重建面补洞）；>= 此值 → 真实孔/槽保留
    int    scanDirMode;      // 轨迹方向：0=自动(最平滑)，1=沿长边喷(行平行长边,换行少)，2=沿短边喷(行平行短边,换行多)
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

// 可喷面角度窗口：面法向与 +Z 的夹角（0~180°）≤ 此值才列入可喷清单。
// 130° 的物理含义：朝上（0°）、侧向（90°）、以及最多向下倾 40° 的面都算可喷；
// 开口盒的底面（180°，整体朝下）被排除。
// 注意：这是"可喷侧朝向"的窗口，不是内外壁判据——内外壁由 isInner 单独给出。
#define SPRAYTRAJ_NORMAL_WIN_DEG 130.0

// 单个重建平面（供显示/调试）
#define SPRAYTRAJ_PLANE_MAX_POLY_VERTS 384
#define SPRAYTRAJ_PLANE_MAX_POLY_TRIS  768
typedef struct SprayTrajPlane {
    double nx, ny, nz;   // 法向（定向后）
    double cx, cy, cz;   // 质心
    int    count;        // 该面内点数
    // ---- 面级法向诊断（内壁 / 外壁）----
    // nx,ny,nz      = 喷距法向：沿用原有"|nz|>0.05 归 +Z 半球 / 近竖直背离质心"的定向，
    //                 轨迹沿它抬升；
    // n0x,n0y,n0z   = 真实法向：射线自由空间判出的可喷侧（材料侧的反面）。
    // 两者点积为负说明该面是内壁（或单层壳二义），界面需要提示。
    // 本组字段只上报，不参与轨迹生成。
    double n0x, n0y, n0z;   // 射线判定后的真实法向（指向可喷侧）
    double angZ;            // 与 +Z 夹角，0~180°
    double tiltDeg;         // 相对水平面仰角 = 90-angZ，-90~90（0=竖直，+90=朝上，-90=朝下）
    double freePlus;        // 沿 +n 侧自由距离 mm（-1 = 量程内无命中）
    double freeMinus;       // 沿 -n 侧自由距离 mm（-1 = 量程内无命中）
    int    side;            // 角度分类：+1 正向(angZ<87°) / -1 反向(angZ>93°) / 0 竖直
    int    isInner;         // 1 = 内壁（可喷侧朝工件内部/腔体）/ 0 = 外壁 / -1 = 薄壁双层未拆
    int    sideSrc;         // 1=射线 2=两侧都自由(兜底) 3=强制翻转 4=薄壁双层(已按层拆行)
    // 该面在自身局部基(u,v)下的轴对齐包围盒
    double ux, uy, uz;   // 局部基 u（主轴方向，单位向量）
    double vx, vy, vz;   // 局部基 v
    double minU, maxU;   // 内点在 u 方向的范围
    double minV, maxV;   // 内点在 v 方向的范围
    // 规则化重建面片：轮廓提取 + 直线简化 + 矩形吸附 + 三角化后的网格。
    // polyVerts 为世界坐标顶点（每 3 个一组 xyz），polyTris 为顶点索引三元组。
    // 边界是直边（矩形面直接给 4 个角点），不再是逐格小方片拼出的锯齿。
    int    polyVertCount;                             // 顶点数
    int    polyTriCount;                              // 三角形数
    double polyVerts[SPRAYTRAJ_PLANE_MAX_POLY_VERTS*3];  // 顶点坐标（x,y,z 连续存放）
    int    polyTris[SPRAYTRAJ_PLANE_MAX_POLY_TRIS*3];    // 三角形索引
} SprayTrajPlane;

// 调用结果
// 槽位数从 8 提到 16：薄壁件一个面块要拆成"内壁 + 外壁"两行，
// 6 个面块拆完是 12 行，8 个槽会截断掉一半内壁。
#define SPRAYTRAJ_PLANE_SLOTS 16
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
    SprayTrajPlane  planes[SPRAYTRAJ_PLANE_SLOTS];  // 各分割平面 / 拆层后的面
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
