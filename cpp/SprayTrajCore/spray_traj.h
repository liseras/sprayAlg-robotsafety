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
    // 面吸收/覆盖阈值 mm：离这张面多远以内的点，仍算"这张面覆盖到"（进轮廓、覆盖率、
    // 残点归并），但**不参与解平面方程**。拟合精度另由点云噪声 3σ 自动决定，不受此值影响。
    // 生效值 = min(此值, 0.6×standoff)：吸收高度 h 处的真实净空是 standoff−h，
    // 取 h ≤ 0.6×standoff 才能保证最坏仍有 0.4×standoff 的余量（否则枪会顶进材料）。
    double flatnessTol;
    double curvatureTol;     // 曲面升级阈值 mm（拟合残差 < 此值则浅曲面当作平面处理）
    int    sorK;             // 统计离群滤波近邻数
    double sorZ;             // 统计离群滤波 Z 倍数
    int    normalFlip;       // 法向翻转：0=自动判向，1=强制翻转（法向取反）
    double fillHoleSize;     // 孔洞填补阈值 mm：点云间隙/孔洞跨距 < 此值 → 覆盖重建（行连续跨越、重建面补洞）；>= 此值 → 真实孔/槽保留
    int    scanDirMode;      // 轨迹方向：0=自动(最平滑)，1=沿长边喷(行平行长边,换行少)，2=沿短边喷(行平行短边,换行多)
    // 枪体安全半径 mm：喷头（枪帽/旋杯）端面到障碍物的最小允许距离，= 喷头外径/2 + 安全余量。
    // 内壁/内底布行时，枪位（面点沿法向抬 standoff 处）到工件点云距离 < 此值即不布行，
    // 避免枪头在角落顶进相邻壁。<=0 时回落到 max(12, 0.2*standoff)。
    double gunClear;
    // 行端点外扩 mm（布行边界外扩策略）。
    // 点云的边缘常常"缺一截"——扫描遮挡，或者边缘点被判给了相邻面——行端点于是
    // 落在比真实边沿更靠里的位置，边沿那一条带就漏喷（实测薄壁开口盒内壁顶部漏 9mm，
    // 与 gunClear 无关，纯粹是布行边界没有外扩）。
    // 每条行的两端沿行方向各外扩此值：<0 = 自动（= 1×行距）；0 = 不外扩。
    // 外扩逐毫米试探，每步都用 gunClear 校验枪位——会把枪头顶进相邻壁的方向自动停住。
    double edgeExtend;
    // 跨行方向外扩 mm：首末行相对点云边界各向外跨多少。
    // 和 edgeExtend 是两个方向：edgeExtend 管"每条行该走多长"（行方向两端），
    // 这个管"第一行、最后一行摆在哪"（跨行方向）。<0 = 自动（= 0.25×行距）；0 = 首末行压在点云边界上。
    // 原实现按固定网格从边界内 0.25×行距 起步，末行会被网格截断、最多少掉近一整个行距
    // （实测薄壁盒内壁顶部末行停在 281，真实边沿 300，漏 19mm）。
    double edgeExtendAcross;
    // 薄壳双面壁厚 mm：一层点代表一张板时，板的两面都要喷；此值 >0 才启用。
    // 上面所有内外壁判据都要求"点云里有两片点"，但薄壁件（1~4mm）激光通常只打到
    // 可及的那一面，只有一层点。这时"内外"不是两片点，而是**同一张板的两个朝向**：
    //   · 射线两侧自由距离都够伸枪 → 板的两面都是空的、都能喷，按此壁厚沿法向
    //     偏移复制一个反向面（sideSrc=5，界面标"薄壳双面"）；
    //   · 只有一侧能伸枪 → 实心件的外表面，另一面在材料里，不复制。
    // 0（默认）= 不启用，行为与之前完全一致。
    double shellThick;
    // ---- 空移让开方向（**点云坐标系**里的方向，不必归一化）----
    //
    // 语义：空移段抬枪绕行时，"离开工件"朝哪一侧走 —— 抬到该方向上的安全面，
    // 再在安全面上平移、最后降回下一面的枪位。整条避障路径的"高度"都按这个轴算。
    //
    // 为什么不复用 SPRAYTRAJ_UP：那是"可喷侧方向"，语义完全不同。
    //   可喷侧 = 相机看得见的那一面，随工件怎么摆而变（侧装时是水平的）；
    //   让开方向 = 工件正面朝哪边，只由工件相对机器人怎么放决定。
    // 侧装/斜装工件上两者近乎垂直，混用一个方向会出大问题：
    // 背板_X+800 实测可喷侧≈基座 +X（水平）、重力=基座 +Y，算法原来按重力 +Y 让开，
    // 结果空移从工件顶上翻过去，一次面间过渡走出 1.8 米（净位移 ΔY=336 / ΔZ=523），
    // 而工件本身只有 866×539×862 —— 全程挂在工件上沿外面绕。
    //
    // 现场口径（2026-11-11 用户拍板）：空移要**朝机器人正面让开**，即基座系 +X。
    //   这既是"机器人伸臂的方向"，也是这一侧没有夹具/背架、机器人一定探得到的方向；
    //   沿 +X 退让 = 枪先撤到工件前面，再从前面平移过去，物理上最自然。
    //
    // 三轴全 0（默认）= 退回旧口径（沿内部 +Z 让开，等价于"工件朝上摆、相机俯拍"）。
    // 调用方应把基座 +X 在**点云系**里的方向填进来（= 手眼矩阵 R 的第 1 行）。
    double liftDirX, liftDirY, liftDirZ;
    // 最小可喷面宽度 mm：重建面片的最小外接矩形短边 < 此值时，判为线状/退化面，
    // 不输出重建面、也不生成轨迹。<=0 关闭（保留旧行为）。默认 18。
    // 注意：平面拟合阈值 flatnessTol 是"点到平面的距离容差"，抓不到共线点集——
    // 一条线对任意包含它的平面残差都≈0，必须用面内宽度/二维秩单独判。
    // ⚠ 实测口径：报出的是**重建网格**宽度，含形态学沉降（+3 格 ≈ 6~18mm，随格距
    //    cs=clamp(sqrt(面宽×面高/8000),2,6) 变），不是材料真实宽度。挂钩细丝真实
    //    10.0mm 会报成 18.0mm，正好压在默认值上被放行 —— 填参时要按"真实宽度 +
    //    膨胀量"折算，或直接把此值当"网格宽度阈值"用（当前现场取 19 即可拦住）。
    double minFaceWidth;
    // 换面不关枪的长度上限 mm（**枪位直线距离**，不是面点距离）。
    //
    // 现场口径（2026-11-19 用户拍板）：换面时如果枪不需要"位移出来"（不用退枪、不用
    // 抬到安全面），直接沿枪位连过去就行，**没必要停枪** —— 停一次枪要减速、断喷、
    // 再起喷，既费节拍，接缝处还容易补枪过量。只有换面空移**太长**，或者必须退枪
    // 抬到安全面绕行时，才真的关枪。
    //
    // 判据（四条同时满足才算"短换面"）：
    //   ① 枪位直线距离 <= 本值；
    //   ② 两端枪轴夹角 <= crossFaceMaxTurnDeg（默认 60°）；
    //   ③ 直线全程离件 >= gunClearEff（枪体安全半径，不是 moveClear 的空移余量）；
    //   ④ 直线不穿过任何重建面（segCrossFace）。
    // 满足时这一段按 pointSpacing 铺满、spray=1、kind=1（查看器画橙色）；
    // 不满足时走原来的"关枪 + 抬枪绕行"（spray=0，蓝色）。
    // 0 = 关闭本行为（一律关枪，与本参数引入前完全一致）。
    // 拒掉时若 SPRAYTRAJ_TRANS_DBG=1，会打出是哪一条拒的、实测值多少。
    double crossFaceGunOnLen;
    // 换面不关枪的**枪轴转角上限**（度）。默认 60。
    //
    // 为什么必须单独给一个参数（2026-11-20 用户报障「我填了 800 阈值，但是还是关枪了」）：
    // 上面 crossFaceGunOnLen 只判**距离**，但短换面还有一条夹角硬闸——
    // 两端枪轴夹角超过 60°（dot<=0.5）时，插值出来的枪轴会快速旋转，真机上等于
    // "一边喷一边猛拧腕"。实测 p037 两块面夹角 82.1°、枪位只隔 190mm，
    // 距离门早就放行了，却在这里被静默拒掉。用户只看到"800 没生效"。
    //
    // 现在把它显式化：想让它放行就把本值调大（现场常见折角件 80~100°）。
    // 转角本身不是位置问题（法向已归一化，枪位不会陷进工件），是**腕部转速**问题：
    // 82° 在 190mm/200mm·s⁻¹ ≈ 0.95s 内转完 ≈ 86°/s，工业腕通常 150~300°/s 跟得上。
    // 90° 以上才真的需要停下重新定位，所以默认仍是 60（保守），按件调。
    // <=0 = 不判夹角（只看距离/净空/穿面），0 = 关闭需用 crossFaceGunOnLen<=0。
    double crossFaceMaxTurnDeg;
    // ---- 喷幅宽度 mm（**独立工艺量，不由行距反推**）----
    //
    // 用户 2026-11-20 16:19 拍板的口径：
    //   「喷幅应该跟行距不是一回事吧，行距设多少，喷幅只跟出粉量有关，
    //     还有关的是里面距离，毕竟远了喷不到」
    //
    // 三者是三个独立的物理量，关系是**行距由喷幅选出来**，不是喷幅由行距反推：
    //   喷幅 W    ← 出粉量、雾化气压、喷距（涂装参数，本参数）
    //   行距 L    ← 由 W 和搭接率 o 选：L = W(1−o)。要改覆盖就调 L，不调 W
    //   有效喷距  ← 喷距（= standoff）超过它喷不到，见下面 maxSprayDist
    //
    // 历史实现把 W 从 L 反推（W = L/(1−o)），于是"行距调大 → 喷幅跟着变宽"，
    // 物理上不成立：行距 120 时反推出 W=171mm、R=85.7mm，等于假设一把枪能喷 171mm 宽。
    // 现在 W 由调用方按实际喷枪给。<=0 = 未提供 → 退回旧的"由行距反推"口径，
    // 保持既有回归逐位不变（填了才走新口径）。
    double fanWidth;
    // ---- 有效喷距 mm：喷距超过此值就喷不到（锥面扩散、雾滴干涸）----
    //
    // 用户原话「还有关的是里面距离，毕竟远了喷不到」。约束的是**沿枪轴的喷距**，
    // 与上面"横向能不能扫到"（喷幅半径）是两件事，必须分开判：
    //   · 横向：点在喷锥的横截面内 → |偏移_⊥| <= R(喷距)
    //   · 轴向：喷距 d 必须在 [0.5·standoff, maxSprayDist] 内 → 太近会冲稀、太远喷不到
    // 现在代码只用 3D 欧氏距离判"已覆盖"，把一个高出主面 150mm 的凸台也算成
    // "距最近行点 85mm 所以喷得到"——而它的真实喷距是 100+150=250mm，早喷不到了。
    // <=0 = 不判（退回旧口径）。
    double maxSprayDist;
} SprayTrajParams;

// 默认参数
SPRAYTRAJ_API SprayTrajParams spray_traj_default_params(void);

// 单个轨迹点
typedef struct SprayTrajPoint {
    double x, y, z;    // 位置 mm（基座系，已按离面距离抬升）
    double i, j, k;    // 喷枪轴线单位向量（指向表面）
    double speed;      // mm/s
    int    spray;      // 1=喷涂段；0=空移段（换面/抬枪/行间跳转，枪关）
                       //   空移段**同样按点距铺满点**（每 pointSpacing 一个），
                       //   不是"只有拐点、两点间靠控制器插补"—— 上位机/DensoArmKine
                       //   必须逐点走到，不能因为 spray=0 就跳过。导出到 plan.json 时
                       //   这一位就是 GunOpen（0 → false）。
    int    kind;       // 换面标记（**只影响查看器配色，不参与开关枪判断**）：
                       //   0 = 普通点。**换行也算 0** —— 同一张面内行末掉头是正常
                       //       连续喷涂，枪没关、点也在覆盖范围内，不该标色
                       //       （用户 2026-11-20 明确：「换面才需要换颜色，换行不需要」）。
                       //   1 = **换面**短接：面与面之间的连接段，且枪未关（spray=1）
                       //       ⇒ 查看器画橙色："这是换面过渡，不是布行，但没停枪"。
                       //       换面才有关枪/不关枪这个需要区分的语义；关枪的那种
                       //       走 spray=0，查看器画蓝色。
                       //   占的是原结构体末尾的 4 字节对齐填充，sizeof 仍是 64，
                       //   老调用方按字段偏移读 x..spray 不受影响。
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
    // 该面在点云里是"单层壳"还是"薄壁双层"（沿法向投影后找最大间隙得出）：
    //   1 = 单层壳：只有一层点，没有壁厚 —— 这类面谈不上"有外壁"，界面要说明，
    //       否则"内壁 4 / 外壁 0"会被误读成算法漏了 4 个外壁；
    //   2 = 双层：拆层前是内外壁合在一起的一整块（拆完每段仍是 1）；
    //   0 = 未判定（曲面块、点数太少、旧 DLL）。
    int    layerCnt;
    double layerThick;   // 沿法向 2%~98% 分位跨度 mm（单层壳 ≈ 噪声量级）
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

// C# 友好入口：结果写入调用方提供的缓冲区（SprayTrajResult 约 145KB，
// 直接用"返回结构体"的 ABI 会在托管侧栈上分配过大，故提供 out 参数版本）。
// 返回 0 = 成功；1 = 失败（out->error 非空）；-1 = 参数空。
SPRAYTRAJ_API int spray_traj_generate_ex(
    const double* xyz, int n, const SprayTrajParams* params, SprayTrajResult* out);

// ABI 尺寸对账：调用方 (C# Marshal.SizeOf) 必须与本函数返回值一致，否则禁止调用，
// 避免结构体布局不一致时原生按自身 sizeof 写越界（曾导致 0xC0000005 闪退）。
SPRAYTRAJ_API void spray_traj_abi_sizes(int* paramsSize, int* planeSize, int* resultSize);

// 结果标量（C# 用固定偏移读取会依赖结果结构体整体布局，故单独给一个扁平信息结构）。
typedef struct SprayTrajInfo {
    int    count;          // 轨迹点数
    int    rows;           // 行数
    int    input_points;
    int    removed_points;
    int    inlier_points;
    int    plane_count;
    double plane_rms;
    double normal_x, normal_y, normal_z;
    double plane_cx, plane_cy, plane_cz;
    double seg_len_mean, seg_len_max;
} SprayTrajInfo;

SPRAYTRAJ_API void spray_traj_result_info(const SprayTrajResult* res, SprayTrajInfo* out);
// 取第 idx 个重建面（0 <= idx < plane_count）；成功返回 0，越界返回 -1。
SPRAYTRAJ_API int  spray_traj_get_plane(const SprayTrajResult* res, int idx, SprayTrajPlane* out);
// 失败原因（成功时为 NULL）。
SPRAYTRAJ_API const char* spray_traj_result_error(const SprayTrajResult* res);

#ifdef __cplusplus
}
#endif
