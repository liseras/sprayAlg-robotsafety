// robot_safety.h —— 机器人可执行性校验 C ABI（与 spray_traj.h 同一套约定）
//
// 职责边界：
//   spray_traj.cpp  负责「喷得对不对」（面分割、行铺满、姿态连续）—— 本文件不动它一行；
//   本文件          负责「机器人能不能安全执行」：可达、关节限位、枪体/连杆碰撞、自碰撞、喷角。
//
// 三条硬要求（用户原话）在这里落地：
//   1. 内壁法向不朝上 → 加倾斜角；
//   2. 倾角改变时，机器人末端到该面点的距离**严格不变**（球面旋转，不是平移）；
//   3. 要考虑碰撞：轨迹（枪体+连杆）会不会撞原始点云工件、关节/连杆会不会自碰撞。
//
// 为什么是「只新增文件」：把安全约束塞进 spray_traj.cpp 会立刻改变所有样本的 CSV，
// 让 master 的「9 个非曲面样本逐字节一致」这条回归红线失效。安全层做成独立入口，
// 下游（WPF / 上位机）自己决定要不要调用、以及调用后是否采纳修正后的轨迹。
#pragma once

#include "spray_traj.h"

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------- 参数
typedef struct RobotSafetyParams {
    double standoff;         // 喷距 mm（必须与生成轨迹时用的值一致；轨迹点 = 面点 + standoff·n）
    double safety;           // 间隙告警阈值 mm：表面距离小于它 → 告警（P1，软约束）
    double hardClearance;    // 绝对安全距离 mm：平滑姿态低于它才回退（硬约束）。
                             // 0 = 只要不穿透；默认 2。与上面的 safety 告警阈值分开：
                             // safety 只决定报告是否标橙，hardClearance 才决定是否回退。
    double search;           // 判定带 mm：超过它只报下界，不当实测值。默认 120（≈ DensoArmKine 口径）
    double toolLen;          // 法兰面 → 喷嘴出口 mm（VS-087 默认 150）
    double tiltNonUpDeg;     // α=90°（竖直侧壁）时的倾角；实际倾角 = 它 × α/90°，α 为法向与 +Z 夹角（默认 20）
    double upThresholdDeg;   // 已不参与倾角决策（倾角改为按法向线性给），保留字段仅为 ABI 稳定
    double tiltMaxDeg;       // 倾角上限（默认 60）
    double tiltStepDeg;      // 倾角搜索步长（默认 5）
    double tiltRateDeg;      // 已不参与倾角决策（姿态改为面级常量）：仅用作姿态跳变告警阈值
                             // （= max(3×它, 10°)）。保留字段仅为 ABI 稳定。
    double maxSprayAngleDeg; // 工艺允许的最大喷角（默认 30；超出报工艺违规，不是碰撞违规）
    double jointMarginMin;   // 关节限位余量下限（默认 1.0°）
    double condMax;          // 归一化雅可比条件数上限（默认 30，只告警不阻断）
    double transitLiftMax;   // 过渡点（spray=0）碰撞时沿 +Z 最大抬升量 mm
    double voxelMm;          // 点云体素降采样边长 mm（0 = 不降采样）。碰撞精度只需远小于 safety
    int    linkStride;       // 连杆碰撞每隔几个点查一次（工具段每点都查）
    int    includeBase;      // 点云里含地面/底座时打开（1=查底座段）
    int    fixTransit;       // 1=过渡点碰撞时允许抬升
    int    azFallback;       // 1=允许把「上翻方位」左右偏 30/60° 兜底
    int    maxCloudPoints;   // 保留字段（<=0 = 不抽稀，用全量点云；KD-tree 支撑全量查询）
    int    threads;          // 并行线程数（<=1 = 单线程，结果可复现；生产建议 = 物理核数）
    int    azMode;           // 进刀侧策略：0=按机器人站位逐行择优（默认）；1=只用「向上」（旧口径，供对拍）
} RobotSafetyParams;

// ---------------------------------------------------------------- 违规位
#define RS_VIOL_IK_FAIL        (1u << 0)   // P0 不可达
#define RS_VIOL_COLLISION      (1u << 1)   // P0 穿透点云
#define RS_VIOL_CLEARANCE      (1u << 2)   // P1 表面距离 < safety
#define RS_VIOL_SELF_COLLISION (1u << 3)   // P0 连杆自碰撞
#define RS_VIOL_JOINT_LIMIT    (1u << 4)   // P1 超出软限位 / 余量不足
#define RS_VIOL_SPRAY_ANGLE    (1u << 5)   // P1 喷角超工艺上限
#define RS_VIOL_NEAR_SINGULAR  (1u << 6)   // P2 接近奇异
#define RS_VIOL_TILT_JUMP      (1u << 7)   // P2 倾角相对上一点跳变过大
#define RS_VIOL_LINK_HIT       (1u << 8)   // P0 连杆/关节/底座撞点云（不含枪体）
#define RS_VIOL_TOOL_HIT       (1u << 9)   // P0 枪体/细管/喷嘴撞点云
#define RS_VIOL_COUNT          10

// 最近碰撞部位编码（Point.worstKind）。工具段报 kind 0/1/2，
// 连杆段报 kind 3 并用 j0/j1 标出两端关节号（如 3 + (2,3) = 「连杆 J2-J3」）。
#define RS_SEG_GRIP     0   // 工具握把（法兰→0.72 工具长）
#define RS_SEG_PIPE     1   // 工具细管（0.72→1.0 工具长）
#define RS_SEG_NOZZLE   2   // 喷嘴（超出 TCP 的 0.18 工具长）
#define RS_SEG_LINK     3   // 连杆 J{j0}-J{j1}
#define RS_SEG_JOINT    4   // 关节球 J{j0}
#define RS_SEG_BASE_COL 5
#define RS_SEG_BASE_NECK 6

// 单个点：既有输入（修正后），也有结论
// 布局刻意「doubles 在前、ints 在后」：P/Invoke 的 LayoutKind.Sequential 按 8 字节对齐，
// 字段混排时容易在 C# 侧算出不同的 padding，这里从根上避免。
typedef struct RobotSafetyPoint {
    double x, y, z;          // 修正后 TCP（倾斜 = 球面旋转，喷距不变；未修改的点与原值一致）
    double i, j, k;          // 修正后枪轴朝向（指向工件，= -n'）
    double speed;

    double tiltDeg;          // 实际喷角（枪轴与**原始**面法向的夹角）
    double clearanceMm;      // 最小表面距离 mm（负 = 穿透）
    double selfGapMm;        // 自碰撞最小间隙 mm（负 = 干涉）
    double marginDeg;        // 关节限位最小余量（度，负 = 超限）
    double standoffErrMm;    // |TCP-面点| 与 standoff 的偏差（球面旋转的守恒大验证）
    double liftedMm;         // 过渡点抬升量
    double cond;             // 归一化雅可比条件数
    double clearanceToolMm;  // 仅枪体（握把+细管+喷嘴）的最小表面距离，负=穿透
    double clearanceLinkMm;  // 仅连杆+关节+底座的最小表面距离，负=穿透
    double theta[6];         // 关节角（度）；ikOk=0 时为 0

    int    spray;            // 1=喷涂
    int    ikOk;
    int    reverted;         // 1 = 平滑位姿违反硬约束（低于绝对安全距离/碰撞/自碰撞/关节超限），
                             //     已退回本点最优；硬约束通过时即使间隙略小也保留平滑
    unsigned int violation;  // RS_VIOL_* 位或
    int    worstKind;        // 全局最近部位（RS_SEG_*）；-1 = 无位姿
    int    worstJ0, worstJ1; // 该部位两端关节号（连杆用；其余为 -1）
    int    linkKind;         // 连杆/关节/底座里最近的那一段（-1 = 没查连杆）
    int    linkJ0, linkJ1;
} RobotSafetyPoint;

// 整体结果
typedef struct RobotSafetyResult {
    RobotSafetyPoint* points;
    int    count;

    int    reachable;        // IK 成功的点数
    int    tiltedPoints;     // 倾角 > 0.5° 的点数
    int    poseReverted;     // 平滑回退点数
    int    violations;       // 违规条目总数
    int    byRule[RS_VIOL_COUNT];

    double minClearanceMm;      // 全部碰撞体（枪体 + 连杆）的最小表面距离
    double minClearanceToolMm;  // 仅枪体的最小表面距离（喷嘴/细管/握把）
    double minClearanceLinkMm;  // 仅连杆/关节/底座的最小表面距离
    double tiltMeanDeg;
    double tiltP95Deg;
    double tiltMaxDeg;
    double standoffMaxErrMm;

    double basePose[16];     // T_base_in_cloud（4x4，行主序）。NULL 输入时是自动假设值
    int    baseAuto;         // 1 = 基座位姿是猜的（不是标定值）
    int    cloudPoints;      // 参与判定的点数（抽稀后）
    double cloudSpacingMm;   // 点云平均间距

    int    azRows;           // 方位角决策单元数（= 笔画数：按折返切，不是"连续 spray 段"）
    int    azStanceRows;     // 其中**实际生效**地改用「站位侧」的笔画数（被逐点回退的不算）

    int    l0Issues;         // L0 输入层自检问题数（吊臂余量 / 反推喷距不一致），非 0 时结论不可信
    double standoffMedianMm; // 沿枪轴反推的喷距中位数（<0 = 未测到）

    const char* error;       // 失败原因（成功为 NULL）
} RobotSafetyResult;

// 默认参数（standoff=100 safety=20 search=120 toolLen=150 ...）
SPRAYTRAJ_API RobotSafetyParams robot_safety_default_params(void);

// 主入口：对一条轨迹做可执行性校验 + 倾角修正
//   pts       : 轨迹点（spray_traj_generate 的输出，8 个字段）
//   cloudXyz  : 原始点云 xyz（3*n 个 double），必须与轨迹同坐标系
//   basePose16: T_base_in_cloud（4x4 行主序），**强烈建议给真值**；NULL = 按工件质心自动假设
// 返回的 points 里位置/朝向已按倾角修正，调用方负责 free。
SPRAYTRAJ_API RobotSafetyResult robot_safety_check(
    const SprayTrajPoint* pts, int n,
    const double* cloudXyz, int cloudCount,
    const double* basePose16,
    const RobotSafetyParams* params);

SPRAYTRAJ_API void robot_safety_free(RobotSafetyResult* res);

// 违规位 -> 规则名（"COLLISION" 等），未知返回 NULL
SPRAYTRAJ_API const char* robot_safety_rule_name(unsigned int bit);

#ifdef __cplusplus
}
#endif
