# 机器人安全层（robot-safety 分支）· C++ 实现 + 上位机集成

> 分支 `robot-safety`。**既有算法源码一行未改**：`spray_traj.cpp` / `spray_traj.h` 完全没动，
> 只在 `CMakeLists.txt` 里把一个新源文件加进同一个 DLL，并在 `Native.cs` 里追加 P/Invoke 声明。
> 结果：10 个回归样本中 **9 个非曲面样本逐字节一致（9/9 IDENTICAL）**，红线守住。
> 上位机（`wpf/SprayTrajViewer`）已经接上：两个页签 + 结论卡 + 风险点标注 + 修正轨迹导出。

---

## 0. 一句话

`spray_traj.cpp` 负责**「喷得对不对」**（面分割、行铺满、姿态连续）；
本层负责**「机器人能不能安全执行」**——可达、关节限位、枪体/连杆撞不撞原始点云、连杆会不会自碰撞、喷角在不在工艺窗口内。
两者编在同一个 `SprayTraj.dll` 里，**互不依赖**：不调用 `robot_safety_check` 时，行为与以前完全一致。

---

## 1. 新增了什么

| 文件 | 行数 | 内容 |
|---|---|---|
| `cpp/SprayTrajCore/robot_safety.h` | 178 | C ABI：参数/逐点结果/整体结果 + 违规位定义 |
| `cpp/SprayTrajCore/robot_safety.cpp` | 1190 | 几何、VS-087 FK/IK/雅可比、KD-tree 点云、胶囊碰撞、倾角规划 |
| `cpp/SprayTrajCore/safety_driver.cpp` | 470 | 命令行驱动 `spray_safety.exe`（轨迹 CSV + 点云 PLY → 校验后 CSV） |
| `wpf/SprayTrajViewer/Native.cs` | +105 | P/Invoke 声明（结构体 / 三个出口函数 / 违规位常量） |
| `wpf/SprayTrajViewer/MainWindow.xaml(.cs)` | ~450 | 双页签界面、安全评估调用、风险点分层渲染、修正轨迹导出 |

**为什么单独一个 exe 而不是改 `driver.cpp`**：`driver.cpp` 的输出决定 9 个样本的 CSV 基线，
安全层不能有一丝机会碰到那条链。`spray_safety` 是并行的独立入口。

### C ABI

```c
RobotSafetyParams robot_safety_default_params(void);

RobotSafetyResult robot_safety_check(
    const SprayTrajPoint* pts, int n,          // spray_traj_generate 的输出
    const double* cloudXyz, int cloudCount,    // 原始点云，与轨迹同坐标系
    const double* basePose16,                  // T_base_in_cloud（4x4 行主序）；NULL = 自动假设
    const RobotSafetyParams* params);

void robot_safety_free(RobotSafetyResult* res);
```

`RobotSafetyPoint` 里除了修正后的 `x,y,z,i,j,k`，还带 `tiltDeg / clearanceMm / selfGapMm /
marginDeg / standoffErrMm / theta[6] / violation`（位掩码）。结构体布局刻意**doubles 在前、ints 在后**，
C 侧 `sizeof` 与 C# 侧 `Marshal.SizeOf` 实测一致：`Params=128 / Point=176 / Result=272`，字段偏移逐项相同。

---

## 2. 三条硬要求怎么落地

### 2.1 内壁法向不朝上 → 加倾斜角（按法向线性给）

倾角不是一个开关，而是**法向偏离朝上程度的线性函数**：

```
α     = 法向与 +Z 的夹角（0° = 正朝上，90° = 竖直侧壁，180° = 完全倒挂）
tilt0 = tiltNonUpDeg × (α / 90°)，封顶 tiltMaxDeg
```

| α | 面的形态 | 基准倾角 |
|---|---|---|
| 0° | 正朝上 | 0°（垂直喷，工艺最优） |
| 30° | 缓坡 | 10° |
| 45° | 斜坡 | 15° |
| 60° | 陡坡 | 20° |
| 90° | 竖直侧壁 | **30°（= `tiltNonUpDeg`，即工艺喷角上限）** |
| 120° | 倒挂内壁 | 40° |
| 150° | 接近倒扣 | 50° |
| `n_z < −0.7` | 朝下悬垂 | 0°（往上扳会顶穿工件，不做线性外推） |

取 90° 作满量程锚点，是因为它物理上就是"面立起来了"的界限：到这里给满工艺允许的倾角，
比侧壁更倒挂的面继续线性增长（必然超工艺，由 `SPRAY_ANGLE` 如实报出来）。

> **为什么 90° 给到最大是 30° 而不是随便一个数**：侧壁是机器人最难应付的形态——末端要
> 绕到面法向附近，腕部只能从侧面斜插进来，翻不过去就撞工件或接近奇异。倾角在这里必须
> 给到工艺允许的上限；再往上（>30°）膜厚就开始明显变薄，那是工艺的账，不是安全的账，
> 所以交给 `SPRAY_ANGLE` 报 P1，而不是偷偷加大。

> 旧实现是阶梯：`n_z < cos45°` 一律给 20°。45° 边界上直接跳 20°，落到轨迹上就是相邻点
> 倾角突跳——这正是实机跟不上的那种"看起来连续、实际一帧一个台阶"。

`tilt0` 只是起点。间隙不够或够不着时按 **固定 5° 网格**往上加到上限 60°；还不行就把
"往上扳"的方位左右偏 30°/60° 再试（内壁往上扳常常顶到对面的壁）。

> 网格必须钉在 5° 整数倍上，不能写成 `tilt0 + k·step`。`tilt0` 随法向连续变化，候选相位
> 会跟着漂移：`curved` idx 2990 在网格的 40° 刚好有 +1.18mm 自碰撞间隙，漂移到 36.48°
> 就变成 −1.16mm 的真干涉——同一个物理点，纯数值抖动决定成败。

**倾斜轴取"指向上方"**——这是唯一既能把枪身从内壁里翻出来、又不改变末端到面点距离的方向。

> 配套修掉一个早停判据的错误：原来"连续 4 档 clearance 变化 < 0.5mm 就停"。可间隙常常早早
> 饱和到 `search` 上限（120mm），饱和区里自碰撞和关节余量仍在随倾角剧烈变化。现在要三个
> 指标都不动才计数，否则会漏掉后面能过关的档位（同样是 2990 那个点）。

### 2.2 倾角改变，距离严格不变（球面旋转）

```cpp
V3 zTool = -unit(d);                     // 工具 +Z = 喷射方向 = −n'
V3 pC    = face + d * standoff;          // d(θ) = Rot(axis,θ)·n
M3 Rc    = frameFromZ(zTool);
```

`|TCP − face| ≡ standoff`，θ 只改喷角。逐点实测 `standoff_err`，**10 个样本全部 0.000000mm**。

> ⚠ 这里最容易搞反的一条：工具 **+Z 是 −d，不是 +d**（`d` 是外法向，`(I,J,K)` 是 −n）。
> 搞反的症状极具欺骗性——位置全对、喷距守恒、数字都"正常"，但 IK 只能选"反腕 + 大臂下探"的解
> （θ2=154°、J2 超限 19°、小臂从工件底下穿过去 24.7mm）。当时差点误判成"工位摆位不合理"去改基座。

### 2.3 碰撞：判的是**整条臂**，不是只判末端

* **工件碰撞**：关节球 + 6 根连杆 + 握把 + 细管 + 喷嘴，共 20 个胶囊（A-B 轴 + 半径），
  每段都对 KD-tree 点云求线段最短距离，`surface = 距离 − 半径`，`<0` 即穿透。
* **按部位分开报**：枪体（握把/细管/喷嘴）与连杆/关节/底座分别给最小距离，并报出**最近的是哪一段**
  （`连杆J3-J4`、`关节J4`…）。
* **自碰撞**：非相邻连杆胶囊两两之间，Ericson 5.1.9 解析解（排除"同一关节球 + 它自己的两根连杆"）。
* **可达/限位/奇异**：阻尼最小二乘 IK，种子顺序 = 上一点解 → 解析种子 → 固定种子表；
  关节软限位余量；归一化雅可比（位置行 ÷ 500mm）条件数。

> **为什么要分开报**：这两类撞车的处置完全不同——枪体撞改喷距/喷角就行，连杆撞必须改站位或
> 整条轨迹的姿态。混成一个"最小距离"，现场拿到 −8.57mm 也不知道该动哪儿。实测 `longstripe`：
> **枪体 21.50mm（安全），连杆 `关节J4` −8.57mm（穿透）**——只看喷嘴会把它判成通过。
> `big` 同理：枪体 47.14mm 很宽松，连杆 `J3-J4` 29.68mm 才是真正压线的那个。

> **修掉一处真实漏判**：原来为了省时间，`linkStride=3`——每 3 个点才判一次连杆，其余 2/3 的点
> 只判枪体。也就是说"连杆会不会撞"这件事，**有三分之二的点根本没查**。现在逐点全判；
> 性能靠剪枝拿回来：某段连"当前最差表面距离 + 半径"都够不着时直接跳过 KD 查询，
> 全量 10 样本仍只要 22 秒。

半径唯一来源 `ArmSizing.From`（reach 905）：关节 34.39、J1 42.99、连杆 21.32、腕 24.76、
握把 21.32、细管 9.59、喷嘴 6.39。

---

## 3. 一个反直觉但必须守住的优先级

**倾角与间隙不是单调的。** `complex` 的 idx 3163 实测：

| 倾角 | 0° | 20° | 45° | 55° |
|---|---|---|---|---|
| 间隙 | 14.2mm | **18.3mm** | 5.6mm | **−2.4mm** |

加大倾角在这里**不是变好，而是变坏**。而"只升不降"的限幅平滑隐含假设了"倾角越大越安全"，
"整行统一倾斜方位"又让个别点被统一方位顶进工件——两处都在**拿 P0 碰撞换 P2 姿态平滑**。

改法：平滑后的位姿只要比本点最优更差（自碰撞通过 / 间隙 / 关节余量逐项比较），就退回本点最优，
记 `reverted`。代价是姿态可能跳一下，由 `TILT_JUMP`（P2）报出来。

效果：`complex` 由穿透 **−5.81mm** → 最小间隙 **+19.87mm**，回退 853 点。

---

## 4. 验证

### 4.1 回归红线：9/9 逐字节一致

用新构建的 `spray_driver.exe`（默认参数）重跑 10 个样本，与基线 `E:/tmp_build/final_*.csv` 比 MD5：

```
synth IDENTICAL   panel IDENTICAL   tray IDENTICAL   zhe IDENTICAL   cao IDENTICAL
complex IDENTICAL longstripe IDENTICAL longwide IDENTICAL big IDENTICAL
curved DIFFER
```

`curved` 的差异**与本改动无关**，做了对照实验：用 `git show HEAD:CMakeLists.txt`（不含 `robot_safety.cpp`）
单独编一份 DLL，输出与含 `robot_safety.cpp` 的 DLL **MD5 完全相同**（`f126e595…`），
而基线 `final_curved.csv` 是 `b68ac2fa…` —— 属第十一轮基线与当前 HEAD 之间的既有差异。

### 4.2 与 DensoArmKine（C# 原版）独立实现对拍

用本层输出的关节角喂 `DensoArmKine.exe --cloud-check`，比最小表面距离：

| 样本 | 抽样 | 最大偏差 | 状态码 |
|---|---|---|---|
| complex | 8 点 | **0.01mm** | 8/8 一致（含 19.87 NEAR ↔ 20.01 CLEAR 的 `safety=20` 切换） |
| tray | 8 点 | **0.01mm** | 8/8 一致 |
| big | 8 点 | **0.01mm** | 8/8 一致（含连杆瓶颈点 idx 7 = 29.68mm） |

FK 公式、连杆表、胶囊半径、IK 收敛解支全部对上。`big` 那一组尤其说明问题：那 8 个点的
最小距离来自**连杆**而不是喷嘴，两边仍然逐点相同——连杆几何的移植是对的。

> 注意 `DensoArmKine` 只报"全局最小表面距离"，不区分部位；分部位是本层新增的能力。

### 4.3 10 个样本（`spray_safety.exe`，默认参数，自动基座）

| 样本 | 点 | 可达 | 倾角 均/p95/最大 | **枪体最小** | **连杆最小** | 最近连杆段 | 回退 | 结论 |
|---|---|---|---|---|---|---|---|---|
| synth | 723 | 100% | 0.0/0.0/0.0° | 64.70 | 77.01 | 关节J1 | 0 | PASS |
| panel | 2458 | 100% | 0.0/0.0/0.0° | 59.04 | 77.01 | 关节J1 | 0 | PASS |
| tray | 1323 | 100% | 4.3/29.2/29.2° | 54.87 | 77.01 | 关节J1 | 0 | PASS |
| zhe | 1863 | 100% | 8.9/21.1/30.0° | 26.32 | 77.01 | 关节J1 | 0 | PASS |
| cao | 1647 | 100% | 5.4/30.0/55.0° | 23.65 | 77.01 | 关节J1 | 1 | PASS_WITH_WARNINGS |
| complex | 4989 | 100% | 14.2/45.0/60.0° | 19.87 | 77.01 | 关节J1 | 852 | PASS_WITH_WARNINGS |
| curved | 3457 | 100% | 22.0/40.0/60.0° | 20.02 | 20.37 | 关节J6 | 88 | PASS_WITH_WARNINGS |
| longstripe | 11558 | **99.4%** | 34.9/55.0/60.0° | 21.50 | **−8.56** | **连杆J3-J4** | 6989 | **REJECT** |
| longwide | 12246 | 100% | 2.8/2.8/20.0° | 35.14 | 38.95 | 连杆J2-J3 | 6 | PASS_WITH_WARNINGS |
| big | 9761 | 100% | 0.0/0.0/0.0° | 47.14 | 29.68 | 连杆J3-J4 | 0 | PASS_WITH_WARNINGS |

全量 10 个样本、约 5.8 万点，**总耗时 24 秒**（含 70 万点点云的 KD-tree 构建；连杆逐点全判之后
单样本最多 15.9s，是 `longstripe` 22 万点）。

**看这张表最该注意的两行**：

* `longstripe` 枪体 21.50mm（安全）但连杆 **−8.56mm**（穿透）——只看喷嘴会判"通过"。
* `big` 枪体 47.14mm 很宽松，连杆 **29.68mm** 才是压线的那个——枪越短、件越大，这种情况越常见。

> 「倾斜点」不再单列：线性化之后只要面有 8° 的坡就有 1.8° 倾角，这个计数已经等价于
> "非正朝上的点数"，看倾角**分布**（均/p95/最大）才有信息量。

**`longstripe` 的 REJECT 是真实结论，不是 bug**：11.5k 点、跨距远超 VS-087 单站臂展（905mm），
自动假设的基座只能覆盖 99.4%，71 点不可达；另有 40 点连杆穿透（最深 −8.56mm，`连杆J3-J4`）。
单站作业覆盖不了这个件，需要导轨/变位机——这类结论正是安全层该给出的。

---

## 5. 已知边界

1. **基座位姿是最大不确定项。** 不给 `-base` 时按"工件质心前方 0.50·reach、下沉 0.35·reach"猜一个，
   报告里会显式标注 `baseAuto=1`。可达率低于 90% 基本就是猜错了，必须用真标定值。
   `longstripe` 就是这种情况。
2. **喷角上限 30° 是占位假设。** 交接文档第 184 行把它列为待用户拍板的开放问题。
   它只影响 `SPRAY_ANGLE`（P1，工艺告警），不影响碰撞/可达结论。
3. **点云遮挡**：KD-tree 查的是"点云上的最近点"，看不到被遮挡的背面。这条对**连杆**判定影响很大——
   单视角扫描通常只有正面点，连杆从背面穿过时不会被发现。凹腔深处若背面缺失，结论会偏乐观；
   这一点与 `DensoArmKine` 同源，不额外引入差异。
4. **碰撞是逐点判，不是连续扫掠**：只在每个轨迹点上摆一次臂，不检查两点之间连杆扫过的空间。
   点距 4mm、倾角跳变受 `tiltRateDeg` 限制时足够，但若有人把点距放大到几十毫米，需要另做扫掠检查。
5. **速度/加速度未校验**：本层只判几何与限位，不判关节角速度/加速度是否超限。

---

## 6. 怎么用

```bash
# 编译
cd cpp/SprayTrajCore
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release

# 校验（不给 -base 时用自动假设）
build/bin/Release/spray_safety.exe traj.csv cloud.ply checked.csv

# 给真基座位姿（mm / 度，ZYX 欧拉）
build/bin/Release/spray_safety.exe traj.csv cloud.ply checked.csv \
    -base 1200 0 600 180 0 0 -standoff 100 -safety 20 -tool 150

# 看结构体 ABI 尺寸（与 C# 侧 Marshal.SizeOf 对照）
build/bin/Release/spray_safety.exe --abi
```

输出 CSV：前 8 列按 master 格式原样，中间额外列原样带出，追加
`tilt_deg, clearance_mm, clearance_tool_mm, clearance_link_mm, worst_part, link_part,
self_gap_mm, margin_deg, standoff_err_mm, lifted_mm, cond, reverted, J1_deg..J6_deg, violation`。

* `clearance_tool_mm` 只看枪体（握把/细管/喷嘴），`clearance_link_mm` 只看连杆/关节/底座；
* `worst_part` 是全局最近的那一段，`link_part` 是连杆侧最近的一段（`连杆J3-J4` / `关节J4`…）；
* `violation` 里 `LINK_HIT` 与 `TOOL_HIT` 分开，不再混成 `COLLISION`。

未测到的值写 `nan`（不是 0）：IK 失败的点没有位姿，`0` 会被读成"贴脸 0mm"。
同目录另写 `trajectory_checked.csv.base.json`（`T_base_in_cloud` + 参数），供对拍时复用同一个基座位姿。

退出码：`0` 无 P0（可执行）、`1` 有 P0（撞了/不可达）、`2` 参数或读文件错误。

### C# 侧（上位机已接好）

`Native.cs` 已加好声明，调用方式：

```csharp
var sp = Native.robot_safety_default_params();
sp.standoff = 100; sp.safety = 20; sp.toolLen = 150;
var rs = Native.robot_safety_check(pts, pts.Length, cloudXyz, cloudXyz.Length / 3, null, ref sp);
try {
    // rs.points / rs.byRule / rs.minClearanceMm …
} finally {
    Native.robot_safety_free(ref rs);   // ⚠ 注意：free 会把 rs.count 清 0，标量要先留副本
}
```

### 6.1 上位机集成（`wpf/SprayTrajViewer`）

界面上线了，不再只是代码可调用。左侧拆成两个页签，两个主按钮都在首屏：

| 页签 | 内容 |
|---|---|
| **轨迹生成** | 原有点云/工艺参数/生成按钮 + 显示控制（新增「标记风险点」「显示修正后轨迹」） |
| **机器人安全** | 基座位姿（自动假设 / 手填标定值）、工具长、安全间隙、工艺喷角上限、倾角搜索上限、关节余量；「机器人安全评估」与「导出修正后轨迹 (CSV)」两个按钮；本页内嵌一份结论 |

右侧 3D 视图右上角固定一张结论卡：结论（PASS / PASS WITH WARNINGS / REJECT）、耗时、
可达率、最小表面距离、**枪体最小（含部位）/ 连杆最小（含 `连杆J3-J4` 这种段名）**、喷角均值/p95/最大、
喷距守恒偏差、平滑回退点数、逐规则违规计数（连杆撞 / 枪体撞分开）、基座位姿（并标明是否为假设值）。

**风险点分三层画**，否则几千个喷角超限点会把工件糊成一片、反而看不出真正的坑：

| 颜色 | 含义 | 画法 |
|---|---|---|
| 红（大点） | P0 碰撞 / 不可达 / 自碰撞 | 逐点全画 |
| 橙（中点） | 间隙不足 / 关节余量不足 | 逐点全画 |
| 黄（小点） | 喷角 / 接近奇异 / 倾角跳变 | 抽稀到 ≤400 个，只指示区段 |

青色叠加层是**倾斜修正后**的轨迹，默认关闭（多数样本与红色原轨迹重合，开着只会遮挡）。
连线判据用**原轨迹**相邻点距，不用修正后距离——倾斜会让点整体挪开 2·standoff·sin(θ/2)，
按修正后距离判会把每一行都判成断开。

命令行（无人值守验收 / 截图取证）：

```bash
SprayTrajViewer.exe --ply xxx.ply --auto --safety 1 --tab 1 --view iso \
    --shot 3d.png --shot-win win.png --export-fixed safe.csv --exit
```

`--safety` 不给就完全不调安全层，行为与以前一模一样（实测日志里 0 条 `SAFETY`）。
`--shot` 只截 3D 视图，`--shot-win` 截整窗（含结论卡）。

**已验证**（截图在 `wpf_shots/`，日志见 `%TEMP%\SprayTrajViewer.log`）：

| 样本 | 结论 | 关键数字 | 耗时 | 截图 |
|---|---|---|---|---|
| `tray` | PASS | 1323/1323 · 枪体 54.87mm · 连杆 77.01mm · 喷角 4.3/29.2/29.2° | 174 ms | `tray_pass.png` / `tray_pass_safety_tab.png` |
| `big` | PASS WITH WARNINGS | 9761/9761 · 枪体 47.14mm · **连杆 J3-J4 29.68mm** · 关节 11 | 1.2 s | `big_link_bottleneck.png` |
| `complex` | PASS WITH WARNINGS | 4989/4989 · 枪体 19.87mm · 连杆 77.01mm · 喷角 14.2/45/60° · 回退 852 | 880 ms | `complex_warnings.png` |
| `longstripe` | REJECT | 11487/11558 · **枪体 21.50mm（安全）/ 连杆 J4 −8.57mm（穿透）** · P0=111 | 15.6 s（22 万点） | `longstripe_reject.png` |

`big` 与 `longstripe` 这两张是这次改动的直接证据：**枪体和连杆分开之后，问题从"最小距离不够"
变成了"到底哪一段不够"**——`big` 是连杆压线而枪体很宽松，`longstripe` 是喷嘴安全但连杆穿进去了。

三组读数与 `spray_safety.exe` 命令行版逐项一致；界面导出的 CSV 与命令行输出逐点对比，
间隙最大差 0.01mm、倾角 0°、TCP 0.0024mm（都是 CSV 小数位舍入）。

**踩到的坑（会再犯）**：`robot_safety_free` 会把 `res->count` 清 0，若在 free 之后再用
`count`/遍历点数，就会得到"可达 1323/0"这种自相矛盾的结论。要么先复制标量，要么先读完再 free。

---

## 7. 附：Python 原型

`robot_safety/`（2526 行）是最初的原型实现，用来快速试算法与做交叉验证，
**不是产品路径**。倾角策略改成线性后两边同源（`_need_tilt` 与 `needTilt` 用同一套候选：
线性基准 + 固定 5° 网格 + 三项早停），逐点对数的差异只来自点云抽稀（C++ 全量、Python 抽稀到
25000 点）与 IK 收敛路径，量级都在千分位。
两者都保留，改算法时可以先在 Python 里试，再落到 C++。
