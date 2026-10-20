# robot_safety —— 机器人安全层

> 分支 `robot-safety`。**只新增文件**：`robot_safety/` 与 `outputs/robot_safety/`，
> 不碰 `spray_traj.cpp` / `spray_traj.h` / `Native.cs`。
> `git status` 应只有这两处 `??`；master 的 9 个非曲面样本 CSV 基线因此逐字节不变。

---

## 0. 一句话

把 master 生成的那条轨迹（"喷得对不对"已经由 master 保证）**再过一遍机器人的账**：
这条轨迹 VS-087 能不能走到、枪体和连杆会不会撞工件、会不会自碰撞、非朝上的内壁有没有
按工艺加倾角而且**加倾角的同时喷距严格不变**。

## 1. 三条硬要求怎么落地的

### 1.1 内壁法向非朝上 → 加倾斜角

判据写在 `planner._need_tilt`：

| 面 | 判据 | 基准倾角 `tilt0` |
|---|---|---|
| 朝上的面 | `n_z ≥ cos(45°)` | **0°**（垂直喷，工艺最优） |
| 侧壁 / 内壁 | `n_z < cos(45°)` | **20°**（`--tilt-nonup`） |
| 朝下的悬垂面 | `n_z < −0.7` | **0°**（往上扳会顶穿工件，所以不扳） |

`tilt0` 只是**起点**，不是结论。起点下若间隙不够或够不着，按 5° 步进加大到 `--tilt-max`（默认 60°）；
还不行就把"往上扳"的方位左右偏 30°/60° 再试（内壁往上扳常常会顶到对面的壁）。

**为什么倾斜轴要"指向上方"**：这是唯一能把枪身从内壁里"翻出来"、同时又不改变末端到面点距离的方向
（见 1.2）。轴 = `up` 在法向切平面上的投影，法向与竖直夹角 < 5° 时无切向可扳，直接不倾斜。

### 1.2 倾角改变，但 TCP 到面点的距离不变

这一条是整个模块的几何地基，写死在 `planner.Evaluator.eval_pose` 里：

```
d(θ) = Rot(axis, θ) · n              # 枪身朝外方向，轴在 n 的切平面内
TCP  = face + standoff · d(θ)        # face 是**原始面点**，不随倾斜移动
⇒ |TCP − face| ≡ standoff            恒等式，不是近似
```

即 TCP 是**绕面点做球面旋转**，不是平移。所以：

* 变的只有喷角 θ（= 枪轴与法向夹角）；
* 喷距一分不变，**逐点实测** `standoff_err = ||TCP−face| − standoff|` 并写进报告的"喷距保持最大偏差"。

实测：所有样本 `standoff_max_err = 0.000000 mm`（浮点级）；FK 反算 TCP 与写回 CSV 的坐标差 0.05mm。

> ⚠ **踩过的坑（最严重的一个）**：工具坐标系 `+Z` 必须取 **−n**（= 喷射方向，指向工件），
> 一度取成 `+n`。症状极具迷惑性：位置全对、喷距也守恒，但 IK 只能选出"反腕 + 大臂下探"的解
> （θ2=154°、J2 超限 19°、J5-J6 连杆穿工件 24.7mm）。修正后同一批点从"全超限"变成余量 +46°。
> 这个错误还连累过一个**错误的 L0 结论**："tool_len(150) > standoff(100) ⇒ 法兰陷进工件"。
> 真相反：工具越长法兰离工件越远、越安全。`precheck` 里留了记录。

### 1.3 碰撞：轨迹 vs 工件、连杆 vs 连杆

| 项 | 口径 |
|---|---|
| 碰撞体 | 底座(可选)/J1..J6 关节球/6 根连杆/工具握把/细管/喷嘴，**全部按胶囊**（A-B 轴 + 半径 R） |
| 表面距离 | 点到轴最短距离 − R，**< 0 即穿透** |
| 半径来源 | `ArmSizing.From`：`jointR=clamp(reach·0.038,12,46)=34.39`、`Joint1R=42.99`、`LinkR=21.32`、`WristLinkR=24.76`、`GripR=21.32`、`ToolPipeR=9.59`、`NozzleR=6.39`、`BaseColumnR=89.4`、`BaseNeckR=39.5` |
| 判定带 | `--search`（默认 120mm，与 `CloudCheck.EffectiveSearchMm = max(120, safety·4)` 同口径）。带外**只报下界**，绝不当实测值 |
| 自碰撞 | 非相邻连杆胶囊两两最小距离 − 半径和（Ericson 5.1.9 解析解，0.015ms/次，无采样误差） |
| 连杆表 | 官方 `vs087.urdf`：J1(0,0,197.5)+Z±170 / J2(30,0,197.5)+Y−100~135 / J3(0,0,445)+Y−133~150 / J4(−20,0,287)+Z±270 / J5(0,0,143)+Y±120 / J6(0,0,80)+Z±360；reach 905；工具 150mm |

**FK 口径**与 `DensoArmKine/Kinematics/FkSolver.cs` 一致：`T_i = T_{i−1}·Trans(offset)·Rot(axis, θ)`，
转轴取 `R_{i−1}·axis`。

## 2. 最重要的一个发现：平滑不许把好点变坏

倾角和间隙**不是单调关系**。在凹槽/内壁里，倾角越大、喷嘴反而越往死角里扎：

```
complex 样本 idx 3163 的实测曲线（同一面点，只改倾角）
  tilt   0°    5°   10°   15°   20°   25°   30°   35°   40°   45°   50°   55°   60°
  clear 14.2  13.8  14.8  17.1  18.3  17.9  14.7  10.5   7.5   5.6   1.8  -2.4  -2.0   (mm)
```

而"姿态连续"这件事逼出了两处平滑，两处都在**拿 P0 碰撞换 P2 姿态平滑**：

1. `_smooth_only_up`：只升不降的限幅平滑，隐含假设"倾角越大越安全"——上表直接否掉它。
2. `az_row`：为了让倾斜方位不逐点跳，整行统一成一个方位——个别点会被统一方位顶进工件。

实测后果：`complex` 原本 `REJECT`（喷嘴穿透 **−5.81mm**，64 个点 COLLISION/CLEARANCE）。
修复方式：**平滑只在不使任何点变差时才生效**。阶段 3 重新求解后，若平滑位姿的
`(自碰撞是否通过, 间隙, 关节余量)` 劣于该点自身最优位姿，就退回后者，并记 `pose_reverted`。

| | 修复前 | 修复后 |
|---|---|---|
| 结论 | REJECT | **P0 = 0** |
| 最小表面距离 | −5.81 mm（穿透） | **+19.87 mm** |
| COLLISION 点数 | 6 | **0** |
| CLEARANCE 点数 | 58 | **2**（19.87 / 19.99，就卡在 20mm 线上） |
| 平滑回退点 | — | 114 |

代价是姿态可能在个别点跳一下 —— 这正是 CSV 里 `tilt_jump_deg` 列存在的意义。
取舍原则：**P0（碰撞）永远压过 P2（姿态跳变）**。

## 3. 用法

```bash
PY="C:/Users/DYF/.workbuddy/binaries/python/envs/default/Scripts/python.exe"

# 基本：给轨迹 + 原始点云
$PY -m robot_safety.run_safety \
    --traj E:/tmp_build/final_tray.csv \
    --cloud cpp/_regress/tray.ply \
    --out-dir outputs/robot_safety/tray

# 有真实基座位姿时（点云系 -> 基座系），两种给法
$PY -m robot_safety.run_safety ... --base-pose 0 0 0 0 0 0
$PY -m robot_safety.run_safety ... --base-json T_base_in_cloud.json

# 只量现状、不倾斜（给"当前工艺"出一份基线）
$PY -m robot_safety.run_safety ... --no-tilt

# 独立实现交叉验证（C# 侧 DensoArmKine 对拍）
$PY -m robot_safety.verify_against_denso --out-dir outputs/robot_safety/tray \
    --cloud cpp/_regress/tray.ply --samples 8
```

退出码：`0` 无 P0 / `1` 有 P0 / `2` 输入或运行错误。

主要参数：`--standoff`(100) `--safety`(20) `--search`(120) `--tilt-nonup`(20) `--up-threshold`(45)
`--tilt-max`(60) `--tilt-rate`(3) `--max-spray-angle`(30) `--tool-len`(150) `--link-stride`(3)
`--max-points`(25000) `--include-base` `--no-tilt` `--no-fix-transit`。

## 4. 输出物

| 文件 | 内容 |
|---|---|
| `trajectory_checked.csv` | **只加列不改列**：原 8 列原样保留，追加 `tilt_deg, tilt_jump_deg, clearance_mm, self_gap_mm, margin_deg, lifted_mm, violation, J1_deg..J6_deg` |
| `verdict.json` | 机器可读：`pass / verdict / counts{P0,P1,P2} / stats / cloud / params / violations[]` |
| `report.md` | 人读报告：关键数字、违规按**连续区段**聚合、方法与口径、坐标声明 |
| `report.png` | 三视图（X-Z / Y-Z / X-Y）：点云灰、安全绿、告警黄、穿透红 |
| `base_pose.json` | `T_base_in_cloud` + 口径参数；对拍脚本与"换工位重算"都用它 |

`verdict` 三档：`PASS`（无违规）/ `PASS_WITH_WARNINGS`（有 P1：间隙偏小、喷角超工艺窗口、贴限位）/
`REJECT`（有 P0：碰撞、自碰撞、不可逆解）。

## 5. 已经验证到什么程度

### 5.1 运动学自检

* 60 组随机姿态 IK→FK 闭环：`worst_pos = 0.0053 mm`、`worst_rot = 0.0006°`；
* 连续 200 点圆弧：`mean 1.70ms / p95 2.89ms`，接近 100% 收敛；
* 加速手段（上一点解做种子 → 解析种子 → 固定种子表）**不得牺牲正确性**：
  只采纳"收敛 + 残差好 + 不超限"的解，否则继续试下一个种子。

### 5.2 与 `DensoArmKine`（C# 原版）对拍

同一组 `J1..J6` 喂给两边，比较各碰撞段的最小表面距离（刚体变换不变量，可直接比）：

| 样本 | 点云 | 抽样 | 最大偏差 | 碰撞判定 |
|---|---|---|---|---|
| synth | 30,351 | 8 | **0.05 mm** | 全部一致（CLEAR） |
| tray | 13,615 | 8 | **0.04 mm** | 全部一致（CLEAR） |
| complex | 83,680 | 8 | **0.07 mm** | 全部一致（含 NEAR↔CLEAR 临界切换） |

连 19.83 / 19.87 / 20.02 这种卡在 safety=20mm 线上的点，两边的 `NEAR`/`CLEAR` 状态都同步。
**结论：FK/IK/连杆表/胶囊半径的移植没有偏差。**

### 5.3 成绩与性能

见 `outputs/robot_safety/SUMMARY.md`（每次跑完由 `run_all` 重写）。

## 6. 已知边界（必须知道的 5 件事）

1. **坐标闭环（P0）尚未完成**。点云系 → 基座系必须来自手眼标定或工位图纸。
   没给 `--base-pose/--base-json` 时，本层按"工件质心落在基座系 `(0.50·reach, 0, 0.35·reach)` ≈
   `(452, 0, 317)`、+X 指向工件"猜一个基座，并在报告里显式标注。
   **「可达率」与「连杆碰撞」结论依赖这个假设；枪体碰撞与喷距/喷角结论不依赖它。**
   实测手眼矩阵 `handeye_real_20260911.json`（绕 Z 转 32°、平移 (520,−180,940)）与点云 AABB
   （synth z 1127~1142、panel z 1343~1374）不自洽，说明轨迹坐标其实已在基座系 ——
   这个矛盾需要你确认后才能真正闭环。
2. **喷角上限 30° 是占位假设**。`ROBOT_SAFETY_交接说明` §提问清单明确把"喷角上限"列为待你拍板项。
   报告里会把"为保安全间隙需要的喷角分布"和这个假设上限并列，冲突点数一语道破，改 `--max-spray-angle` 即可。
3. **不含外部轴 / 变位机 / 导轨**。基座按固定安装处理。
4. **点云抽稀影响精度**。默认 `--max-points 25000`（按采样密度自适应体素）。抽稀会略降精度，
   要更高置信度就把 `--max-points` 调大（代价是线性变慢）。
5. **碰撞体是胶囊近似**，不含气管走向、夹具、料架、围栏。`--include-base` 只在点云含地面/底座时开。

## 7. 模块结构

```
robot_safety/
  geom.py        向量/4x4/罗德里格斯/点线距离/胶囊-胶囊距离(Ericson 解析解)/倾斜轴
  robot.py       VS-087 连杆表、fk、jacobian、ik(阻尼最小二乘+解析种子)、arm_segments、ArmSizing
  cloud.py       PLY 读写(ascii/binary LE/BE)、KD-tree 胶囊段查询(全向量化)、按密度抽稀
  traj.py        Traj 容器、CSV 读写(只加列不改列，支持字符串列)、行(行=连续 spray 段)
  planner.py     Params / Evaluator / _need_tilt / _smooth_only_up / precheck(L0) / plan / apply_result
  report.py      verdict.json / report.md / report.png
  run_safety.py  CLI
  verify_against_denso.py   与 DensoArmKine.exe 对拍
```

## 8. 留给下一个人的坑位提示

* `geom.tilt_axis_and_dir` 的退化判据**必须在归一化之前**取 `|u| = sin∠(n, up)`。
  先 `unit()` 再比阈值恒成立 → 平整朝上的面会被判成"非朝上"，喷角被加到 40°，
  而喷距仍然守恒、数字"看着正常"，极难发现。默认阈值 5°。
* 归一化雅可比条件数才有意义（位置行 / 500mm）→ 正常姿态 ~30；原始雅可比天然上千，
  正常姿态就 1118，拿它当奇异判据会误报一片。
* `search` 必须 ≥ `standoff`，否则所有段都饱和、报告里只剩"≥N"。
* 点云抽稀要按**采样密度**定体素（`s' = s·√(N/target)`），按 AABB 对角线会把 111k 点抽成 190 点。
* 报告里的表如果按逐点展开，几百行会淹掉结论 —— 违规一律按连续区段聚合，逐点清单留在 `verdict.json`。
