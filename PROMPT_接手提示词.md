# 接手提示词 · 喷涂轨迹算法（截至第十六轮）

> 用法：新开一个对话窗口，**把下面整段（含反引号代码块）原样粘贴给 AI**，
> 它会先读 `HANDOFF_交接说明.md` 再动手。所有路径、命令、基线数据都已写死，不需要你再解释。
>
> **本文件每轮都要跟着更新** —— 里面的"当前任务"和"已完成项"一旦过期，接手的人会重做已完成的工作。

---

```
你在继续开发一个「喷涂轨迹生成算法」C++/C# 项目，工作区 E:\Code\喷涂算法自研。

## 第一步：必读
完整读一遍 E:\Code\喷涂算法自研\HANDOFF_交接说明.md。重点看：
- 第 0 节开头的「第十六轮」（等高线圈闭合 + 环形窄面走同心圆）、「第十五轮」（回转体改等高线圈
  + 姿态限幅）、「第十四轮」（多平面件排除）、「第十三轮」（重建曲面模型）；
- 6.4「本轮未收尾项」，剩下的都是**等用户拍板**的取舍，不是技术难题；
- 第 9 节「踩过的坑」共 50 条，第 30~50 条是曲面件/多面体件/回转体与"轨迹乱"相关的，动手前必读。

## 分支约定
- `master`：算法主干，只接受"9 个非曲面样本逐字节一致"的改动。
- `robot-safety`：机器人安全（可达性/碰撞/姿态工艺约束），见
  E:\Code\喷涂算法自研\ROBOT_SAFETY_交接说明.md。**别在 master 上直接做安全功能。**

## 项目背景（30 秒版）
C++ DLL 导出 C ABI（spray_traj.cpp，约 4090 行），C# WPF 查看器调它。
流程：点云 → SOR/KNN 法向 → **曲面块判定（含多平面件排除、块内平盘剥离）** → 5° 区域生长
→ 8° 块合并 → 陡壁 RANSAC 剥离 → 矢量化轮廓重建 → 行生成 → 蛇形轨迹；
曲面块走独立的「曲面模型 → 行走向自检（柱面/回转面）→ 行追踪」分支：
柱面走平面截线，回转面走等高线圈。
回归样本在 cpp/_regress\*.ply，共 28 个（panel big zhe cao tray complex longstripe longwide
synth curved 等 + 第十四轮入库的 **box** 开口盒子 + 第十五轮入库的 **bin** 垃圾桶/截锥）。
现场点云 17 份在 E:/tmp_build/userclouds/u01..u17（映射见 _map.txt，u08=曲型板、u04=复杂点云）。

## 编译（C 盘紧张，临时目录必须指到 E 盘）
cd /e/Code/喷涂算法自研/cpp/SprayTrajCore/build
export TMP=E:/tmp_build TEMP=E:/tmp_build
cmake --build . --config Release
产物：build/bin/Release/SprayTraj.dll + spray_driver.exe

跑一个样本：
cd /e/Code/喷涂算法自研/cpp/_regress
../SprayTrajCore/build/bin/Release/spray_driver.exe curved.ply /tmp/c.csv /tmp/c.ply -flat 40 -fill 300

改完 C++ 后必须再把 DLL 刷进查看器（否则看的是旧结果）：
cp cpp/SprayTrajCore/build/bin/Release/SprayTraj.dll \
   wpf/SprayTrajViewer/bin/Release/net9.0-windows/SprayTraj.dll
（若被占用：Get-Process SprayTrajViewer | Stop-Process -Force）

## 验收红线（违反任何一条都算没做完）
0. **最省事的一次性全量对拍**：`C:/Python38/python.exe E:/tmp_build/reg.py <tag>`，
   会跑完 25 个样本并写 `reg_<tag>.txt`（逐样本 行数/点数/rms/seg_mean/seg_max/面数/CSV的MD5）。
   与上一轮的 `reg_*.txt` 直接 diff，只有你打算改的样本才允许有差异。
1. **9 个非曲面样本的轨迹 CSV 必须与基线逐字节一致**：
   cd /e/tmp_build/reg9
   for f in panel big zhe cao tray complex longstripe longwide synth; do
     /e/Code/喷涂算法自研/cpp/SprayTrajCore/build/bin/Release/spray_driver.exe \
       $f.ply /e/tmp_build/v_$f.csv /dev/null -flat 40 -fill 300
     cmp /e/tmp_build/final_$f.csv /e/tmp_build/v_$f.csv && echo "$f OK"
   done
   基线在 E:/tmp_build/final_<样本>.csv。注意 n1_*/n2_* 是旧实验不是基线。
   注意 reg9 目录里没有 curved.ply，曲型板只在 cpp/_regress/。
2. curved 必须保持：14 行曲面 + 4 行平面、行距 11.76mm、覆盖率 92.2%、贴合残差 8.88mm、
   rms=8.878、seg_mean=4.959 / seg_max=7.658（第十五轮把姿态限幅从 6°/点 收到 150°/s 后
   max 降了 15%，**位置类指标一字未变**）。
   box 必须保持：5 面 / 125 行 / 12028 点 / rms=0.30 / seg_mean=3.847 / seg_max=3.999。
   bin（垃圾桶/截锥）必须保持：3 面 / 30 行 / 4264 点 / rms=3.929，
   其中锥面 18 圈（行距 12.00mm、覆盖率 99.8%、贴合残差 3.93mm）+ 顶环 2 圈同心圆 + 底面 10 行。
   **注意 rms 从 1.653 变成 3.929 不是退步**：第十六轮前每圈被追成 7.7 圈螺旋，等效行距
   被压到 4mm，覆盖率/残差都是虚高的。行距 12mm 的理论 RMS≈3.5mm，3.93 才是真值。
3. 不要动 split_planes 的非曲面逻辑，也不要改参数语义（flatnessTol/fillHoleSize/
   scanRowSpacing/curvatureTol 的含义见文档第 2 节，用户是照着这套理解的）。
4. 曲面块的改动只允许动姿态与显示，不许动行点位置：位置一变，覆盖率与贴合残差的
   验收口径就全变了。要动位置必须用户同意并重标验收数字。

## 当前任务（按优先级）
### 等用户拍板的未决项（见 HANDOFF 6.4）
- ⑤ 窄面单道规则（面宽 <= N x 行距时只走 1 行）——需要用户给 N。
- ⑥ 曲型板端部折边：A 保持 / B 折边分行 / C 收紧姿态限幅（互斥，需用户拍板）。
- ⑦ 复杂点云面0/面2 夹角 26 度是否放宽合并到 30 度（会碰 zhe 的 13.2 度真折边，必须重跑回归）。
- ⑧ ~~环形/窄面走平行扫描~~ → ✅ 第十六轮已做（`detect_ring_uv()` + 同心圆）。

### 已完成，不要重做
- 曲面可视化网格（第十轮）、seg_len_max 偏大（第十~十三轮已解决）、
  row00 只有 77 点（已定性为真实卷边，留档不改）；
- 重建网格曲率正则化 + 枪位跨行一致化（第十三轮，见 outputs/round13/）；
- 开口盒子被误判成曲面块（第十四轮，见 outputs/round14/）：新增 `patch_is_polyhedral()`，
  "能被 ≤8 个精确平面覆盖 ≥93%" 就不是曲面。别再用"整体平面残差大小"当唯一判据。
- 回转体轨迹不规则（第十五轮，见 outputs/round15/）：`surf_march` 换方向场
  （`cross(n,∇h)` = 等高线，不是 `cross(n,eIdx)` = 平面截线）、块内平盘剥离、
  姿态限幅按角速度反算。三个坑见 HANDOFF 坑 45~47，**改之前一定要读**。
- 「一条轨迹绕很多圈、行距设 30 还是密」（第十六轮，见 outputs/round16/）：
  `surf_march` 缺**闭合判据**，闭合等高线圈被一圈圈追到 maxSteps（一圈 96 步被追成 741 步
  = 7.7 圈），7 层叠成螺距 2.5mm 的螺旋。已加闭合判据 + iso 模式只追一侧 + 锚点末层贴边；
  环形窄面（顶环）另走同心圆，并把它的重建格距从行距改成 `min(行距, 环宽/6)`。
  三个坑见 HANDOFF 坑 48~50。

## 常用调试
SPRAYTRAJ_DUMP=1 spray_driver.exe curved.ply ... 2>&1 | grep -E "曲面dbg|曲面模型"
SPRAYTRAJ_PROF=1 ... 2>&1 | grep prof                                     # 分阶段耗时
C:/Python38/python.exe cpp/_regress/mesh_smoothness.py <out>.ply.plane0.ply   # 网格曲率一致性
C:/Python38/python.exe cpp/_regress/gun_cleanliness.py surf_rows.csv          # 枪位干净度
spray_driver.exe in.ply out.csv out.ply -row 30 -pt 4 -so 100 -sp 200         # 复现界面参数
E:/tmp/plotcv.py <ply> <csv> <png>                                            # 点云+轨迹三视图

## 环境备注
- Python 两个都能用：C:/Python38/python.exe（numpy 1.24）与
  C:/Users/DYF/.workbuddy/binaries/python/envs/default/Scripts/python.exe（numpy 2.5）。
  Python 不认 Git Bash 的 /e/ 前缀，脚本参数要用 E:/ 绝对路径。
- 本机无 Visual Studio / g++，构建只能走 cpp/SprayTrajCore/build 里已配好的 CMake（MSBuild）。
- /tmp 在 Git Bash 里是 C:\Users\DYF\AppData\Local\Temp，在 Python 里是 E:/tmp，两边不一样。
- 驱动只吃 ASCII 路径；中文点云先复制成 uXX.ply。

先读文档，然后从「当前任务」里挑一项（或等用户指定）开始。每完成一项，把
「改了什么文件 / 关键数据 / 回归是否全绿」写清楚再继续下一项。
```

---

## 附：历史任务档案（已完成，供追溯）

| 任务 | 状态 |
|---|---|
| 曲面块可视化网格（原来只有 4 顶点 / 2 三角） | ✅ 第十轮完成，第十三轮又做了曲率正则化 |
| `seg_len_max=25.3mm` 偏大（其它样本 ~4mm） | ✅ 第十轮定位 3 处、第十一轮限幅摊平；第十三轮 `seg_max=9.038` |
| 第 00 行只有 77 点 | 📝 已定性：外缘卷边处自然到边，留档不改 |
| 开口盒子（纯多面体）被误判为曲面块，24 行覆盖 57% | ✅ 第十四轮：`patch_is_polyhedral()`，5 面 125 行、残差 0.30mm |
| 回转体（垃圾桶）轨迹沿母线放射、端点在轴心抱团 | ✅ 第十五轮：等高线圈 18 圈、覆盖率 100%、残差 1.65mm |
| 提交 | ✅ `3fde1e1`（第十三轮）；第十四/十五轮见 `git log -1` |
