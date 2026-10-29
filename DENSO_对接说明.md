# 本仓轨迹 → DensoArmKine 对接说明

> 目标：本仓（喷涂算法自研）生成的轨迹，能在 `E:\Code\新喷涂框架\SprayVision\docs\SprayVisionSimple\tools\DensoArmKine` 里直接跑起来（播放 / 逆解 / 碰撞）。

## 0. 最省事的一条路：界面自动导出（2026-10-26 起）

界面（`wpf/SprayTrajViewer`）勾着「自动导出（点云同名 + 模拟平台 JSON）」时，
**生成轨迹就到盘、安全评估完再覆盖一次**，全程不用点按钮、不用选目录：

```
<点云所在目录>\traj_out\
    面板_traj.csv           原始轨迹（X,Y,Z,I,J,K,speed,spray）
    面板_safe.csv           修正后轨迹（含 tilt/clearance/关节角/违规位）
    面板.plan.json          DensoArmKine 直接导入的这一份（有修正结果时用修正后的位姿）
    导入说明.txt            加载口径速查（只写一次）
```

* 文件名 = **点云文件名**（去扩展名），不同工件不会互相覆盖。
* 重新生成轨迹时会**先删掉上一轮的 `_safe.csv`**：那份对应的是另一条轨迹，
  留着会造成"新轨迹 + 旧姿态"的组合，而这件事从文件名上看不出来。
* 无人值守：`SprayTrajViewer.exe --ply x.ply --auto --safety --fast 1 --shot-win s.png --exit`
  （`--exit` 只在带 `--shot`/`--shot-win` 时才生效，这是既有行为）。
* `--fast 0|1` 对应界面「单次逆解」：19082 点的开口盒子 30.8s → 1.59s（19×），
  代价与适用边界见 `robot_safety/README.md` 3.2。

JSON 是**在 C# 里直接写的**（`MainWindow.WriteDensoPlanJson`），不再依赖 Python 环境。
与 `traj_to_denso_plan.py` 逐点对拍：位置最大偏差 0.0005mm（CSV 三位小数）、
姿态旋转矩阵最大偏差 0.0006（≈0.03°，来自 CSV 四位小数的轴向），GunOpen/Label 全一致。

## 1. 三段链路

```
点云 ply ──(spray_driver.exe)──► 轨迹 csv ──(traj_to_denso_plan.py)──► spray_plan json ──► DensoArmKine --traj-run
          X,Y,Z,I,J,K,speed,spray      Points[].Cartesian       逐点逆解统计
```

一条命令跑完（推荐）：

```bash
# 已经平移到基座可达区的点云（C:\Users\DYF\Desktop\点云\喷涂点云_DENSO可达）
python run_in_denso.py --cloud "C:/Users/DYF/Desktop/点云/喷涂点云_DENSO可达/面板.ply"

# 原始点云（相机系）：轨迹必须跟着点云平移，交给脚本自动算同一份 delta
python run_in_denso.py --cloud "C:/Users/DYF/Desktop/点云/喷涂点云/面板.ply" --shift-from-cloud
```

分段手动跑：

```bash
spray_driver.exe cloud.ply traj.csv traj.ply
python traj_to_denso_plan.py --csv traj.csv --out plan.json --no-shift
DensoArmKine.exe --traj-run --traj plan.json --model VS-087 --frame tcp --traj-frame base --out r.txt

# 带碰撞（点云与轨迹必须在同一坐标系里）
DensoArmKine.exe --traj-cloud --cloud cloud.ply --traj plan.json --model VS-087 \
                 --frame tcp --traj-frame base --cloud-in-base 1 --out c.txt
```

## 2. 口径对照（三条都已核对，不是推测）

| 项 | 本仓 CSV | DensoArmKine JSON |
|---|---|---|
| 位置含义 | **TCP**（喷枪位姿点），`pos + standoff·(I,J,K)` 落在工件表面 | `Points[].Cartesian[0..2]`，TCP（`--frame tcp` 时） |
| 姿态 | `I,J,K` = 工具 +Z = 喷射方向（只有轴向，无绕轴转角） | 欧拉 `[rx,ry,rz]`，**度**，`R = Rz·Ry·Rx` |
| 绕轴滚转 | CSV 不含 | 由转换器补齐：上一点 X 投影到新切平面（最小旋转，和 planner.py 同一套） |
| 坐标系 | 点云坐标系 | `--traj-frame base` 表示已是基座系 |
| 喷枪开关 | `spray` 列：1=喷涂 / 0=空移 | `GunOpen`：true / false |

**`GunOpen=false` 的点不是"跳过点"，必须照走。** 空移段（换面、抬枪、行间跳转）
由生成器**自己插补成整条路径**，不是"只给几个拐点、中间交给控制器插补"：

- **路径**：退出当前面（沿法向）→ 抬到安全高度 → 平移 → 下降/沿法向进入下一面；
  两端本来就在件上方且直线净空足够时直连，不做无谓抬升。
- **插补**：每一段按点距（默认 4mm）铺满，**折线拐点强制保留**（不能等距重采样，
  否则会把直角切掉、路径被拉直）。
- **姿态**：从上一面喷涂末点的枪轴 slerp 到下一面首点，按累计弧长均匀分配，
  单步转角 ≤5°/点；接缝严格等于相邻喷涂点姿态（零跳变）。
- **验安全**：每条空移段插补后逐点查点云净空 + 是否穿面；**不合格就换途经**
  （反向退出 / 三个抬升高度 / 左右各三档侧绕），在可行候选里取净空最大、同净空取最短；
  绕行长度上限 = 原路径 ×2（`SPRAYTRAJ_MOVE_DETOUR` 可放宽），都不可行才退回原路径并告警。
- **安全门限**：默认取枪体安全半径（`gunClear`，20mm）；`SPRAYTRAJ_MOVE_NEED` 可抬高。
  实测开口盒子这类内腔空移的几何极限就是 20mm 出头，用 25mm 当门限会把合格路径误判成不可行。
- 插补出来的每个点随后仍进入 robot_safety 的逐点检查（IK / 连杆 / 自碰撞 / 间隙），
  所以 `--traj-run` 里空走点必须逐点走，`跳过(标记无效)` 应为 0。

实测依据（`outputs/denso_link/` 内可复算）：

- **CSV 位置是 TCP**：开口盒子 19082 点，`pos + 100·(I,J,K)` 到点云最近邻的距离中位 **0.83mm**、p90 1.77mm、99.9% 在 5mm 内。反过来 `pos` 本身距点云中位 99.99mm —— 所以 **转换时绝不能再退一个 standoff**，退两次会让 J3/J5 大面积超限（盒子限位内 36.3% → 换算前只有 20.4%）。
- **欧拉约定**：`Kinematics/EulerZyx.cs` 明确 `R = Rz·Ry·Rx`，`FromAnglesDeg` 说明单位是度。
- **平移一致性**：同一条轨迹，「在平移后点云上直接生成」与「在原云上生成 + `--cloud` 补平移」的 JSON 位置最大差 **0.0009mm**、欧拉差 **0°**。

## 3. 验证结果（当前摆放：中心 = (452.5, 0, 316.75)）

`--traj-run --frame tcp --traj-frame base`，限位内可执行 = 真机能摆到的点：

| 工件 | 点数 | 限位内可执行 | 可达率 |
|---|---|---|---|
| 面板 | 2451 | 2362 | 96.4% |
| 大平面(现场点云) | 10009 | 9021 | 90.1% |
| 槽型盘 | 1826 | 1327 | 72.7% |
| 曲型板工件 | 2490 | 1749 | 70.2% |
| 折型板(现场点云) | 1970 | 1294 | 65.7% |
| 长板宽板(中间残缺) | 12476 | 7948 | 63.7% |
| 开口盒子(带壁厚4mm) | 19082 | 6931 | 36.3% |

超限几乎都集中在 **J5**（肘腕俯仰），其次是 J3。

## 4. 关键发现：工具长度两边不一致（这是可达率的主因）

| 一侧 | VS-087 工具长度 | 来源 |
|---|---|---|
| 本仓 `robot_safety/robot.py` | **150mm**（`VS087_TOOL_LEN`） | URDF 默认值 |
| DensoArmKine | **453.08mm**（偏移 −2.21 / +6.77） | `%APPDATA%\DensoArmKine\settings.json`，实测标定值 |

同一份曲型板轨迹、同一套逆解，只换工具长度：

| 工具长度 | 限位内可执行（100 点抽样） |
|---|---|
| 150mm | **92.0%** |
| 453.08mm | **70.0%** ← 与 DensoArmKine 实测 70.2% 吻合 |

结论：可达率的差异**不是转换错**，而是两边对「枪有多长」的假设差 303mm。腕心位置 = TCP − 工具长，453mm 的枪把腕心又往后推 300mm，J5 最先顶到限位。

两条路，取决于现场真枪长度：

1. **真枪就是 453mm** → 本仓 `VS087_TOOL_LEN` 应改成 453 并重算安全层；同时 100mm 喷距意味着腕心在工件后方 553mm，工件中心要往外放（实测 dx=+300 时 73.2%，比放近好一点，但仍然不高）。
2. **真枪 150mm（本仓假设）** → 在 DensoArmKine 里把 VS-087 的 `LengthZ` 改成 150，可达率立刻回到 92% 量级；碰撞判定的「喷嘴 NEAR」也会消失（那个 5.52mm 是 453mm 枪 + 喷嘴建模再外延 18%≈81mm 造成的，实测喷嘴尖到面 65.6mm）。

## 5. 已知限制

- `spray_driver.exe` **只认 ascii PLY**。`喷涂点云` 里 `垃圾桶(模拟生成).ply`、`开口盒子(模拟生成).ply` 是 binary，要先转 ascii。
- 轨迹里同一工件的不同面是**拼接**的：开口盒子一处 TCP 突变 764mm、相邻枪轴 180° 翻转 8 处。真机需要过渡段，DensoArmKine 会把这些点记为「疑似换分支」。
- 开口盒子 36.3% 是几何必然：400×300×300 的箱内壁，150mm 枪（甚至 453mm 枪）要在箱内转向，J5 顶死。这类件必须靠变位机/外部轴，或者分站位。
- 本仓 `robot_safety/cloud.py` 依赖 scipy，托管 python 3.13 没装；`C:/Python38/python.exe`（scipy 1.10.1 + matplotlib 3.7.5）可直接跑。
