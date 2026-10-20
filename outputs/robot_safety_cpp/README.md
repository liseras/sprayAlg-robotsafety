# 机器人安全层（robot-safety 分支）· C++ 实现与验证

> 分支 `robot-safety`。**既有源码一行未改**：`spray_traj.cpp` / `spray_traj.h` 完全没动，
> 只在 `CMakeLists.txt` 里把一个新源文件加进同一个 DLL，并在 `Native.cs` 里追加 P/Invoke 声明。
> 结果：10 个回归样本中 **9 个非曲面样本逐字节一致（9/9 IDENTICAL）**，红线守住。

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

### 2.1 内壁法向不朝上 → 加倾斜角

| 面 | 判据 | 基准倾角 |
|---|---|---|
| 朝上 | `n_z ≥ cos45°` | **0°**（垂直喷，工艺最优） |
| 侧壁 / 内壁 | `n_z < cos45°` | **20°** |
| 朝下悬垂 | `n_z < −0.7` | **0°**（往上扳会顶穿工件） |

20° 只是起点。间隙不够或够不着时按 5° 步进加到上限 60°；还不行就把"往上扳"的方位左右偏 30°/60°
再试（内壁往上扳常常顶到对面的壁）。**倾斜轴取"指向上方"**——这是唯一既能把枪身从内壁里翻出来、
又不改变末端到面点距离的方向。

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

### 2.3 碰撞

* **工件碰撞**：关节球 + 6 根连杆 + 握把 + 细管 + 喷嘴，每段一个胶囊（A-B 轴 + 半径），
  对 KD-tree 点云求线段最短距离，`surface = 距离 − 半径`，`<0` 即穿透。
* **自碰撞**：非相邻连杆胶囊两两之间，Ericson 5.1.9 解析解（排除"同一关节球 + 它自己的两根连杆"）。
* **可达/限位/奇异**：阻尼最小二乘 IK，种子顺序 = 上一点解 → 解析种子 → 固定种子表；
  关节软限位余量；归一化雅可比（位置行 ÷ 500mm）条件数。

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

效果：`complex` 由穿透 **−5.81mm** → 最小间隙 **+19.87mm**，回退 116 点。

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

FK 公式、连杆表、胶囊半径、IK 收敛解支全部对上。

### 4.3 10 个样本（`spray_safety.exe`，默认参数，自动基座）

| 样本 | 点 | 可达 | 倾斜点 | 倾角 p95/最大 | 最小间隙 | 回退 | 结论 |
|---|---|---|---|---|---|---|---|
| synth | 723 | 100% | 0 | 0.0°/0.0° | 64.70 | 0 | PASS |
| panel | 2458 | 100% | 0 | 0.0°/0.0° | 59.04 | 0 | PASS |
| tray | 1323 | 100% | 194 | 20.0°/20.0° | 54.87 | 0 | PASS |
| zhe | 1863 | 100% | 562 | 20.0°/30.0° | 25.35 | 0 | PASS |
| cao | 1647 | 100% | 341 | 30.0°/40.0° | 23.65 | 1 | PASS_WITH_WARNINGS |
| complex | 4989 | 100% | 2129 | 45.0°/60.0° | 19.87 | 116 | PASS_WITH_WARNINGS |
| curved | 3457 | 100% | 3217 | 40.0°/60.0° | 11.53 | 20 | PASS_WITH_WARNINGS |
| longstripe | 11558 | **99.5%** | 11496 | 55.0°/60.0° | **−20.87** | 7536 | **REJECT** |
| longwide | 12246 | 100% | 120 | 0.0°/20.0° | 35.23 | 7 | PASS_WITH_WARNINGS |
| big | 9761 | 100% | 0 | 0.0°/0.0° | 32.44 | 0 | PASS_WITH_WARNINGS |

全量 10 个样本、约 5.8 万点，**总耗时 22 秒**（含 70 万点点云的 KD-tree 构建）。

**`longstripe` 的 REJECT 是真实结论，不是 bug**：11.5k 点、跨距远超 VS-087 单站臂展（905mm），
自动假设的基座只能覆盖 99.5%，两端 61 点不可达；另有 21 点真实穿透（最深 −20.87mm）。
单站作业覆盖不了这个件，需要导轨/变位机——这类结论正是安全层该给出的。

---

## 5. 已知边界

1. **基座位姿是最大不确定项。** 不给 `-base` 时按"工件质心前方 0.50·reach、下沉 0.35·reach"猜一个，
   报告里会显式标注 `baseAuto=1`。可达率低于 90% 基本就是猜错了，必须用真标定值。
   `longstripe` 就是这种情况。
2. **喷角上限 30° 是占位假设。** 交接文档第 184 行把它列为待用户拍板的开放问题。
   它只影响 `SPRAY_ANGLE`（P1，工艺告警），不影响碰撞/可达结论。
3. **点云遮挡**：KD-tree 查的是"点云上的最近点"，看不到被遮挡的背面。凹腔深处若背面有点云缺失，
   结论会偏乐观——这一点与 `DensoArmKine` 同源，不额外引入差异。
4. **速度/加速度未校验**：本层只判几何与限位，不判关节角速度/加速度是否超限。

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
`tilt_deg, clearance_mm, self_gap_mm, margin_deg, standoff_err_mm, lifted_mm, cond, reverted,
J1_deg..J6_deg, violation`。
同目录另写 `trajectory_checked.csv.base.json`（`T_base_in_cloud` + 参数），供对拍时复用同一个基座位姿。

退出码：`0` 无 P0（可执行）、`1` 有 P0（撞了/不可达）、`2` 参数或读文件错误。

### C# 侧

`Native.cs` 已加好声明，调用方式：

```csharp
var sp = Native.robot_safety_default_params();
sp.standoff = 100; sp.safety = 20; sp.toolLen = 150;
var rs = Native.robot_safety_check(pts, pts.Length, cloudXyz, cloudXyz.Length / 3, null, ref sp);
try {
    // rs.points / rs.byRule / rs.minClearanceMm …
} finally {
    Native.robot_safety_free(ref rs);
}
```

UI 上的按钮与结果面板尚未接（属于下一步，不在本阶段范围）。

---

## 7. 附：Python 原型

`robot_safety/`（2526 行）是最初的原型实现，用来快速试算法与做交叉验证，
**不是产品路径**。C++ 版与它逐项对数一致（`synth` 全等；`complex` 倾斜点 2129 vs 2128、
回退 116 vs 114，差异来自 C++ 用全量点云、Python 抽稀到 25000 点）。
两者都保留，改算法时可以先在 Python 里试，再落到 C++。
