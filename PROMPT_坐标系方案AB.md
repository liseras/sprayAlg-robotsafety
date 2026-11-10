# 新窗口提示词：倾角 +Z 与手眼矩阵 —— 方案 A / 方案 B 二选一

直接整段复制粘贴。

---

```
项目：D:\Code\喷涂算法自研_robot-safety\喷涂算法自研_robot-safety
核心：cpp\SprayTrajCore\spray_traj.cpp（轨迹生成）
      cpp\SprayTrajCore\robot_safety.cpp / .h（机器人安全层，C ABI）
      wpf\SprayTrajViewer\MainWindow.xaml.cs（界面 + 调用）

【问题】
倾角逻辑把「向上」写死成点云系的 +Z，但点云是相机系，+Z 是相机光轴不是重力。
三处硬编码（已核对行号）：
  robot_safety.cpp:1308  dirOf(...)                      倾斜方向
  robot_safety.cpp:1323  tiltAxisAndDir() { V3 up(0,0,1); }   把法向往这个方向扳
  robot_safety.cpp:1339  tiltNeedDeg() { n.z ... n.z < -0.7 }  倾角大小 + 悬垂面判定
  robot_safety.cpp:1731  baseAuto 猜基座时 z(0,0,1) 当基座 Z 轴

现有手眼矩阵（T_base_cam，p_base = R·p_cam + t，已验正交 det=1）：
  文件 D:\Code\cloudepoint\点云\点云\喷涂点云_DENSO可达\槽型盘handeye_cam2base.json
  R = [[ 0.21994,  0.26892, -0.93771],
       [-0.06141,  0.96316,  0.26181],
       [ 0.97358,  0.00000,  0.22835]]
  t = [1573.33, -266.67, -120.00]
  该文件自带 Notice：「手工录入·滑条对齐，非标定级，典型 3~10mm / 0.5~2°」
  重力向上在相机系 = R 的第 3 行 = (0.9736, 0, 0.2284)，与相机 +Z 夹角 76.8°

【现状（不要推翻的既有事实）】
1. 点云载入后一次都没变换过：MainWindow.xaml.cs:454  _cloud = PlyReader.Read(path);
2. 基座位姿是单独传的：MainWindow.xaml.cs:1411 BuildBasePose()，
   CmbBaseMode.SelectedIndex != 1 时返回 null → robot_safety.cpp:1715 走 baseAuto 分支
   （按质心猜基座，日志打 baseAuto=1，「可达率仅作预判」）
3. 安全层内部已经做反向变换，不需要新写：
   robot_safety.cpp:1153  toBaseP() / toBaseR()    云系 → 基座系（IK 在这里解）
   robot_safety.cpp:1157  toCloudP() / toCloudR()  基座系 → 云系（碰撞体搬回来求距离）
4. 法向是从点重建的：spray_traj.cpp:1170 CloudNormals::build（KNN K=25 + PCA + 保边平滑）。
   所以点云搬动后法向自动跟着变，**不要"点转一次、法向再转一次"**，会转两遍。
5. 布行不受坐标系影响（这条容易搞错，务必看）：genRows 的入参是 P2，是**面局部基**
   (u,v) 坐标（spray_traj.cpp:4907-4908  P2.push_back({dot(q,u),dot(q,v)})），
   行方向靠遍历 theta 优化（:5008/:5046/:5056）。刚体变换不改平面内几何 → 行距、
   覆盖、分面全部不变。唯一可能的副作用：:4075 / :4901 用 zaxis{0,0,1} 构造局部基
   u = cross(z, n)，搬动后 u/v 朝向会变，可能让"行方向"整体差 90°。不是"歪"。

【方案 A：点云不动，给安全层传重力方向】
  A1. robot_safety.h 的 RobotSafetyParams 末尾加 double upCloud[3]（点云系重力向上，全 0 = 旧口径）
  A2. robot_safety.cpp 入口解析一次：upW = 有效 ? unit(upCloud) : V3(0,0,1)
  A3. 三处替换：
        tiltAxisAndDir(n, upW, axis, uUp)        // 替换 :1325 的 V3 up(0,0,1)
        tiltNeedDeg(n, upW, p)：n.z → dot(n, upW)；n.z < -0.7 → dot(n, upW) < -0.7
        baseAuto 的 z 轴 → upW（:1738-1744）
  A4. WPF 侧从手眼矩阵算 upCloud = R 的第 3 行，或让用户直接填 3 个数
  ⚠ ABI：robot_safety.h:94 明确写了 sizeof(RobotSafetyParams)==152，fastPose 是塞在
     尾部填充里的。加 3 个 double 会破 ABI，必须同步改 Native.cs 并核对 --abi 打印。
     用 safety_driver.exe --abi 可以打出 size/offset 对照。
  ⚠ 不修"基座在哪"这件事，baseAuto 照旧；需要可达性就必须另给 basePose16。

【方案 B：点云搬到基座系】
  B1. 在**载入处**插一次刚体变换（唯一插入点），不要放在生成轨迹之后：
        MainWindow.xaml.cs:449-454 LoadCloud() 里 _cloud = PlyReader.Read(path) 之后
        对 _cloud 每点做 p_base = R·p_cam + t
  B2. basePose16 改传显式单位阵（不能传 null：null 会触发 baseAuto 再猜一次，等于白搬）
  B3. 之后全流程照旧，robot_safety.cpp 一行不用改
  ⚠ 这样 +Z 天然正确（基座系 +Z = 重力向上，前提是机器人地面正装），
     A 方案要改的四处全部连带解决，也不碰 ABI
  ⚠ 平移 t 要给对。用文件里那份 t 试算煤气灶背板（301334 点，相机系 2280~2680 / 984~1098）：
        p_base 后质心 = (1159.6, -52.8, 2531.8)，距原点 2790mm
        VS-087 臂展 905mm → 直接全不可达
    所以要么 t 用现场摆位（align_cloud_to_robot.py:38 的 TARGET_DEFAULT =
    (0.5·reach, 0, 0.35·reach) = (452.5, 0, 316.75)），要么确认那份手眼矩阵不是这个工位的。
    只取 R、t 用摆位，是更稳的拆法。
  ⚠ 槽型盘handeye_cam2base.json 里 Notice 连写两遍：矩阵按文件原始点直接左乘；
    主程序若开「点云Y轴取反」必须先翻 diag(1,-1,1) 再乘，联调时把该开关关掉。
    界面开关是 CmbNormalFlip（MainWindow.xaml.cs:266 / :531）。

【请你做的事】
1. 先明确选 A 还是 B，给一句话理由（我倾向 B：几何上刚体变换不改任何距离/夹角，
   只换语义；且一次解决倾角方向、倾角大小、悬垂判定、基座 Z 轴四处，不破 ABI）
2. 只改选中的那一个，**不要两个都实现**
3. 改之前先核对上面列的行号（文件可能有变动）
4. 改完必须验证"不必倾斜的点结果逐字节不变"——这套代码对回归一致性很敏感，
   默认值必须等价于旧行为
5. 给出改动前后的对比：面数、行数、轨迹点数、tilt 分布、P0/P1 违规数

【测试素材】
  点云 D:\Code\cloudepoint\点云\点云\喷涂点云_DENSO可达\煤气灶背板.ply（ASCII PLY，301334 点）
  该目录下还有开口盒子/槽型盘/折型板/曲型板等回归样本
  生产参数行距 80，喷距 100（--row 80）

【不要做】
- 不要同时实现 A 和 B（会把手眼矩阵乘两次）
- 不要改 genRows 的行带逻辑（它工作在面局部基，本来就跟坐标系无关）
- 不要改残点归并容差 mergeD、平面性门 rmsCap（那是另一个已单独定位的问题）
- 不要假定手眼矩阵是标定级精度
```

---

## 附：这份提示词里我此前说错过、已在这里更正的两点

1. **布行不会因坐标系而歪。** 我前面几次说"相机系里按 XY 切带、投影到斜面会变密变疏变歪"——**错的**。`genRows` 收到的 `P2` 是面局部基 `(u,v)` 坐标（`spray_traj.cpp:4907-4908`），不是世界 XY。刚体变换不改变平面内几何，行距和覆盖率不受影响。唯一的真实副作用是 `zaxis{0,0,1}` 构造局部基时 `u/v` 朝向变化，可能让行方向整体差 90°。

2. **`RobotSafetyParams` 的 ABI 不是"可能破"，是确定会破。** `robot_safety.h:94` 自己写着 `sizeof==152`、`fastPose` 是塞在尾部填充里的。方案 A 加 3 个 `double` 必须同步改 `Native.cs`，用 `safety_driver.exe --abi` 核对。
