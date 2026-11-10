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
// 两层口径（见 RobotSafetyParams.requireReach）：
//   · 「几何层」= 喷距 + 倾角/喷角 + 枪体不撞点云 —— 只跟点云和轨迹有关，任何机器人都成立；
//   · 「运动学层」= 可达 + 关节限位 + 连杆碰撞 + 自碰撞 —— 只对**指定型号 + 指定站位**成立。
// 默认两层都查（requireReach=1）；点云不是给这台机器人用时可以关掉运动学层，
// 只留几何层，这样"轨迹铺得对不对"不会被"905mm 臂展够不着"掩盖。
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
    double tiltNonUpDeg;     // δ=90°（竖直侧壁）时的倾角；实际倾角 = 它 × δ/90°，
                             // δ 为法向与**重力向上**的夹角。重力向上由 robot_safety_check_up
                             // 的 upCloud 给；老入口 robot_safety_check 里它就是点云 +Z。
                             // 默认 60（2026-10-26 由 30 改为 60，用户口径「60°×δ/90」：侧壁直接给满）。
                             // 注意它同时是喷角软告警的比较对象：倾角 > maxSprayAngleDeg 会报 SPRAY_ANGLE
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
    // 是否把最终轨迹当成「这台机器人要执行」。1=要求可达（默认，完整口径）；0=不要求。
    //
    // 为什么需要这个开关：可达率（IK / 关节限位 / 连杆自碰撞）是**机器人型号 + 站位**的
    // 属性，不是「点云 + 轨迹」的属性。同一份点云换个臂展/换条站位，结论可以从 REJECT
    // 变成 PASS；点云本身没有对错。拿它当点云的验收口径，会把「这条轨迹铺得对不对」
    // 这类真问题淹掉（开口盒子实测：1113 个 P0 里 887 个是 J5 超限 ±120°，纯粹是 905mm
    // 臂展够不到 292mm 深的盒底——把基座从"货架下"挪到地面，P0 立刻 1115 → 499）。
    //
    // 0 时实际做了什么（**只关运动学，不关几何**）：
    //   · 不解 IK、不建连杆/关节/底座碰撞体、不查自碰撞、不判关节限位；
    //   · 仍然逐点检验：喷距守恒（球面旋转）、倾角/工艺喷角、**枪体（握把+细管+喷嘴）
    //     与点云的间隙**。枪体位置由「目标 TCP + 枪轴」解析确定，与机器人无关，
    //     所以"枪会不会戳进工件"这一条在两种模式下是同一个判据、同一份数据。
    //   · R.reachable 直接记成 n（"按全部通过"），CSV 里 theta 全 0、ikOk 全 0，
    //     界面/报告上必须显著标成「未检查」，不能把 100% 当成真测出来的。
    int    requireReach;     // 0=不要求（只查几何）；1=要求（默认）。零初始化结构=0，
                             // 但所有调用方都走 robot_safety_default_params()，那里给 1。
    int    threads;          // 并行线程数（<=1 = 单线程，结果可复现；生产建议 = 物理核数）
    int    azMode;           // 进刀侧策略：0=按机器人站位逐行择优（默认）；1=只用「向上」（旧口径，供对拍）
    // 姿态求解策略。**按位解释**（旧调用方传 0/1，行为逐字节不变）：
    //   bit0 = 单次逆解（fast）：姿态纯解析（倾角=法向线性需求，方位=向上），每点只解 1 次 IK，
    //          解不出/过不了硬约束就照实报，不做大步回退（只保留 ±1 步救援）。
    //   bit1 = 倾角逐点连续：倾角按该点法向**逐点**线性给（tilt = tiltNonUpDeg × δ/90，
    //          δ = 法向与「正位置」+Z 的夹角），**不做面级常量**。
    //
    // 为什么 bit0 需要：搜索模式下平均每点 9.5 次 IK（朝上的面 27 次且全失败），开口盒子 19082 点
    // 实测 30.8s / 8 线程，其中 IK 占 95% CPU。fast 模式把 IK 降到 1 次/点，实测同件 1.6s，
    // 而**碰撞、可达、关节、喷角这四类判据一条不少**——省掉的只是"为了挑一个更容易解的倾角
    // 而反复试解"，不是安全判定。
    //
    // 为什么 bit1 需要：面级常量（一个面组锁一个倾角）会把**连续曲面切成台阶**。垃圾桶实测
    // 同一圈桶壁被切成 0°（306 点）和 25°（95 点）两块，中间没有过渡，相邻点枪轴夹角最大 48.7°；
    // 槽型盘更极端，315°~360° 那 348 点全是 0°、其余侧 23~35°。工艺上要的是
    // 「法向偏离正位置越大、倾得越多」的连续关系，而不是面块常量。
    int    fastPose;
    // 平放面（法向 ∥ 竖直）时，倾角往哪一侧扳。安全层工作在**基座系**，基座 +Z 竖直向上，
    // 所以这是一个基座系里的选择：
    //    0  = 不指定（旧口径：平放面倾角恒 0，"往上扳"没意义）
    //   +1  = 往基座 +X 侧倾（朝基座前方）
    //   −1  = 往基座 −X 侧倾（朝基座身后）
    // 只有基座位姿已知（手眼 / 恒等 / 手动）时才生效 —— 自动假设的基座是猜的，
    // 拿它的 X 轴当倾角方向等于把猜测当标定。竖直侧壁等非法向∥up 的面仍按"向上"扳。
    double tiltSideX;
} RobotSafetyParams;
#define RS_POSE_FAST        (1 << 0)   // 单次逆解
#define RS_POSE_TILT_PERP   (1 << 1)   // 倾角逐点连续（取消面级常量）
// 注意：requireReach 占用的就是原来那个没用上的 `maxCloudPoints` 保留槽位，
// fastPose 占的是结构体尾部原本的 4 字节对齐填充，所以 sizeof(RobotSafetyParams) 仍是 152、
// ABI 不变，旧调用方（不碰这两个字段）拿到的仍是 default_params 给的 1/0 —— 行为与以前逐字节一致。

// ---------------------------------------------------------------- 违规位
#define RS_VIOL_IK_FAIL        (1u << 0)   // P0 不可达
#define RS_VIOL_COLLISION      (1u << 1)   // P0 穿透点云
#define RS_VIOL_CLEARANCE      (1u << 2)   // P1 表面距离 < safety
#define RS_VIOL_SELF_COLLISION (1u << 3)   // P0 连杆自碰撞
#define RS_VIOL_JOINT_LIMIT    (1u << 4)   // P1 超出软限位 / 余量不足
#define RS_VIOL_SPRAY_ANGLE    (1u << 5)   // P1 喷角超工艺上限
#define RS_VIOL_NEAR_SINGULAR  (1u << 6)   // P2 接近奇异
#define RS_VIOL_TILT_JUMP      (1u << 7)   // P2 相邻点枪轴夹角跳变过大（不只比倾角，见下）
#define RS_VIOL_LINK_HIT       (1u << 8)   // P0 连杆/关节/底座撞点云（不含枪体）
#define RS_VIOL_TOOL_HIT       (1u << 9)   // P0 枪体/细管/喷嘴撞点云
#define RS_VIOL_STEP_JUMP      (1u << 10)  // P2 相邻点 TCP 直线步长过大（跨行/跨换段也查）
#define RS_VIOL_COUNT          11

// TILT_JUMP 与 STEP_JUMP 的分工（两个都报过，别合并）：
//   TILT_JUMP 看**姿态**——相邻枪轴夹角超 tiltRate 的 3 倍（最低 10°）。用 3D 夹角而不是
//     只比 tiltDeg：逐点口径下方位翻转时倾角可以恒定不变（垃圾桶桶底小坑那段全程
//     42.1~42.2°，方位绕了一整圈），只比倾角等于看不见。
//   STEP_JUMP 看**位置**——相邻 TCP 直线距离超限。它查的是"机器人要直接跳过去多远"，
//     和姿态是不是平滑无关：喷涂末点被倾斜推出 73.5mm、而紧邻的空移点按设计不动，
//     倾角跳变和步长跳变在这里同时成立，但别的地方（纯平移）只有步长跳变。

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

    int    reachable;        // IK 成功的点数。requireReach=0 时恒等于 count（未检查，按全部通过）
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

// 安全层的构建时间戳（如 "Oct 23 2026 13:20:05"），静态字符串、无需释放。
// 存在理由：WPF 程序用的是 bin 目录里的那份 SprayTraj.dll，C++ 改完后若忘了
// 同步（或程序开着导致 PostBuild 拷贝失败），界面会安静地跑旧算法 —— 现象是
// "明明改了，怎么还是老样子"。界面上把这一行显示出来，一眼就能判断用的是哪一版。
SPRAYTRAJ_API const char* robot_safety_build_stamp(void);

// 主入口：对一条轨迹做可执行性校验 + 倾角修正
//   pts       : 轨迹点（spray_traj_generate 的输出，8 个字段）
//   cloudXyz  : 原始点云 xyz（3*n 个 double），必须与轨迹同坐标系
//   basePose16: T_base_in_cloud（4x4 行主序），**强烈建议给真值**；NULL = 按工件质心自动假设
// 返回的 points 里位置/朝向已按倾角修正，调用方负责 free。
//
// ⚠ 本入口隐含假设「点云系 +Z = 重力向上」。点云若是相机系（+Z 是光轴），必须改用
//   robot_safety_check_up 并把真实重力向上传进去，否则倾角整条口径都错。
SPRAYTRAJ_API RobotSafetyResult robot_safety_check(
    const SprayTrajPoint* pts, int n,
    const double* cloudXyz, int cloudCount,
    const double* basePose16,
    const RobotSafetyParams* params);

// 同上，但显式给「重力向上」在**点云系**里的方向。这是倾角口径的正解。
//
// 为什么另开一个入口而不是给 RobotSafetyParams 加 upCloud[3]：
//   结构体现在 sizeof==152、fastPose 就塞在尾部填充里，再塞 3 个 double 必然变尺寸。
//   旧调用方（含已编译的 exe）按 152 字节分配、DLL 按 176 读，会读到结构体外的栈内存，
//   表现是「参数随机错乱、不崩」—— 比直接崩更难查。改签名只会让旧的调用方链接失败，
//   那是显式报错，可控。
//
// upCloud : 点云系里的重力向上方向，**不必归一化**（函数内归一）。
//           NULL 或长度 < 1e-6 ⇒ 退回 (0,0,1)，即 robot_safety_check 的行为。
//
// 怎么取这个方向 —— 别对矩阵的性质想当然：
//   已知 p_base = R·p_cam + t，则基座 +Z（重力向上）在点云系里 = **Rᵀ·ẑ**，
//   写成行向量正好是 **R 的第 3 行** (R[2][0], R[2][1], R[2][2])。
//   常见误用是取 R 的第 3 **列**，那是「相机光轴在基座系里的朝向」，物理上完全是另一回事
//   （本机实测两者夹角 149.4°，方向近乎相反）。
//
// 但更要紧的一条：手眼矩阵的 R 往往是手工对出来的，不一定可信；而「工件怎么摆」
// 是工艺事实。实测煤气灶背板（立着挂）：主面法向与点云 +Z 只差 5°，可手眼 R 给出的
// up 会让主面偏离竖直 16.5°，只取 "up=±点云 y 轴" 才得到 3.7° —— 立式背板 860 高 × 400 宽。
// 所以**优先用工艺事实定 up，矩阵只作交叉校验**。
SPRAYTRAJ_API RobotSafetyResult robot_safety_check_up(
    const SprayTrajPoint* pts, int n,
    const double* cloudXyz, int cloudCount,
    const double* basePose16,
    const RobotSafetyParams* params,
    const double upCloud[3]);

SPRAYTRAJ_API void robot_safety_free(RobotSafetyResult* res);

// 违规位 -> 规则名（"COLLISION" 等），未知返回 NULL
SPRAYTRAJ_API const char* robot_safety_rule_name(unsigned int bit);

#ifdef __cplusplus
}
#endif
