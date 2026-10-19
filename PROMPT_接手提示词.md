# 接手提示词 · 喷涂轨迹算法（curved 曲面件收尾）

> 用法：新开一个对话窗口，**把下面整段（含反引号代码块）原样粘贴给 AI**，
> 它会先读 `HANDOFF_交接说明.md` 再动手。所有路径、命令、基线数据都已写死，不需要你再解释。

---

```
你在继续开发一个「喷涂轨迹生成算法」C++/C# 项目，工作区 E:\Code\喷涂算法自研。

## 第一步：必读
完整读一遍 E:\Code\喷涂算法自研\HANDOFF_交接说明.md。
重点看「6.3 第九轮：曲面件曲面路径」和「6.4 本轮未收尾项」——
那是上一轮做完的部分和明确留给你的任务。第 9 节「踩过的坑」有 29 条，
其中 23~29 条是曲面件专属的，动手前必读，能省掉几小时返工。

## 项目背景（30 秒版）
C++ DLL 导出 C ABI（spray_traj.cpp，约 2940 行），C# WPF 查看器调它。
流程：点云 → SOR/KNN 法向 → 5° 区域生长 → 8° 块合并 → 陡壁 RANSAC 剥离
→ 矢量化轮廓重建 → 行生成 → 蛇形轨迹。
10 个回归样本在 cpp/_regress\*.ply（panel big zhe cao tray complex longstripe
longwide synth curved）。

## 编译（C 盘紧张，临时目录必须指到 E 盘）
cd /e/Code/喷涂算法自研/cpp/SprayTrajCore/build
export TMP=E:/tmp_build TEMP=E:/tmp_build
cmake --build . --config Release
产物：build/bin/Release/SprayTraj.dll + spray_driver.exe

跑一个样本：
cd /e/Code/喷涂算法自研/cpp/_regress
../SprayTrajCore/build/bin/Release/spray_driver.exe curved.ply /tmp/c.csv /tmp/c.ply -flat 40 -fill 300

## 你的任务（按优先级，做完一项汇报一项）

### 任务 1：修曲面块可视化网格（最高优先）
现状：curved.ply 跑出来 `[重建] 面片0: 顶点=4 三角=2（曲面块）` —— 只有 4 个顶点。
根因：build_surface_mesh()（spray_traj.cpp 约 1930 行）调用了 build_plane_mesh()，
后者末尾的「最小外接矩形吸附」（矩形度 ≥0.85 直接替换成 4 角点）把曲面压成了平板。
轨迹本身是对的（14 行、覆盖率 92.2%），只是 WPF 里看不出曲面。

要求：把 build_surface_mesh 改成独立的高度场结构化网格，**不要**走 build_plane_mesh：
- 曲面点投到平均平面得 (u,v)；
- 网格分辨率在 ABI 预算内取最细：for K=64..8，cell=max(wU,wV)/K，
  取第一个满足 `nu*nv ≤ SPRAYTRAJ_PLANE_MAX_POLY_VERTS(384)`
  且 `2*(nu-1)*(nv-1) ≤ SPRAYTRAJ_PLANE_MAX_POLY_TRIS(768)` 的 K
  （curved 实测 K≈30、cell≈18mm，31×11 顶点、600 三角）；
- 每格取格内点云的中位高度（防单点噪声顶尖刺），空格用最多 8 轮邻域平均扩散填；
- 必须做闭运算（膨胀 1 格再腐蚀 1 格）桥接扫描线间隙，否则网格中间漏缝；
- 只对「已占用且四角都有有效高度」的格发 2 个三角；
- 顶点循环里要有 `if((int)verts.size() >= maxVerts) break;` 保护，vid[] 按 nU*nV 预分配。
验收：顶点数应到 300 量级、三角 600 量级；把 out.ply.plane0.ply 和点云叠加画图，
确认网格跟着拱面起伏而不是一块平板。加一行 `[曲面] 曲面网格: NxM 顶点=X 三角=Y` 诊断输出。

### 任务 2：查 seg_len_max=25.321mm 偏大
curved 的 OK 行是 `rms=8.878 seg_mean=6.430 seg_max=25.321`，其他样本 seg_max 都 ~4mm。
嫌疑在曲面块的 emitGap 行间插值（spray_traj.cpp 约 2649 行，d ≤ 3L 时保持 spray=1）
或面间过渡段。用 SPRAYTRAJ_DUMP=1 导出 CSV，按 spray 列分段，
定位 spray=1 但相邻点距 >8mm 的地方，判断是插值跨度太大还是过渡段漏标 spray=0。
修完 seg_max 应回落到 4mm 附近。

### 任务 3（可选）：第 00 行只有 77 点
其他曲面行 175~194 点，第 0 行贴边缘被截断。看是否需要给边缘行单独放宽 maxOff。

### 任务 4：提交
前三个任务做完后 git commit（当前 spray_traj.cpp 有未提交改动）。
commit message 写清做了什么、验收数据是什么。

## 验收红线（违反任何一条都算没做完）
1. **9 个非曲面样本的轨迹 CSV 必须与基线逐字节一致**：
   cd /e/Code/喷涂算法自研/cpp/_regress
   for f in panel big zhe cao tray complex longstripe longwide synth; do
     ../SprayTrajCore/build/bin/Release/spray_driver.exe $f.ply /tmp/v_$f.csv /dev/null -flat 40 -fill 300
     cmp /e/tmp_build/final_$f.csv /tmp/v_$f.csv && echo "$f OK"
   done
   基线在 E:/tmp_build/final_<样本>.csv。注意 n1_*/n2_* 是旧实验不是基线，别拿来对拍。
2. curved 必须保持：14 行曲面 + 4 行平面、行距 11.76mm、覆盖率 92.2%、贴合残差 8.88mm
   （允许因网格改动有小数点级波动，行数/覆盖率/残差量级不能掉）。
3. 不要动 split_planes 的非曲面逻辑，也不要改参数语义（flatnessTol/fillHoleSize/
   scanRowSpacing/curvatureTol 的含义见文档第 2 节，用户是照着这套理解的）。

## 常用调试
SPRAYTRAJ_DUMP=1 spray_driver.exe curved.ply ... 2>&1 | grep 曲面dbg   # 逐行追踪
SPRAYTRAJ_PROF=1 ... 2>&1 | grep prof                                  # 分阶段耗时
E:/tmp/plotcv.py <ply> <csv> <png>                                     # 点云+轨迹三视图（需 C:/Python38/python.exe）
uv_planeN.csv / mask_*.pgm                                              # DUMP 模式导出的重建掩膜

## 环境备注
- Python 用 C:/Python38/python.exe（已装 numpy/matplotlib）；3.13 那个没装包。
- 本机无 Visual Studio / g++，构建只能走 cpp/SprayTrajCore/build 里已配好的 CMake（MSBuild）。
- /tmp 在 Git Bash 里是 C:\Users\DYF\AppData\Local\Temp，在 Python 里是 E:/tmp，
  两边不一样，跨工具传路径时用绝对路径 E:/tmp/xxx。
- 查看器进程若占用 DLL，先 Get-Process SprayTrajViewer | Stop-Process -Force。

先读文档，然后从任务 1 开始。每完成一项，把「改了什么文件/关键数据/回归是否全绿」写清楚再继续下一项。
```

---

## 附：如果只要一件事——任务 1 的验收标准

| 项 | 现状 | 目标 |
|---|---|---|
| `[重建] 面片0` 顶点数 | 4 | **300 量级**（curved 上 31×11=341 或近似） |
| `[重建] 面片0` 三角数 | 2 | **600 量级** |
| WPF 里看到的 | 一块平板 | 贴着拱面起伏的网格 |
| 9 个非曲面样本 | 逐字节一致 | **仍逐字节一致** |
| curved 行数/覆盖率/残差 | 14 行 / 92.2% / 8.88mm | **不退化** |
