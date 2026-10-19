# 新窗口启动提示词（直接整段复制粘贴）

---

```
我要继续开发 E:\Code\喷涂算法自研 这个喷涂轨迹生成项目。请先读
E:\Code\喷涂算法自研\HANDOFF_交接说明.md 和 .workbuddy\memory\2026-10-19.md 了解现状，
然后按下面的优先级往下做。

【项目概况】
- 核心算法：cpp\SprayTrajCore\spray_traj.cpp（约 1800 行，C ABI DLL），
  流程是 体素下采样 → 统计离群滤波 → KNN PCA 法向（保边平滑）→ 5° 区域生长 →
  块级合并(8°) → 簇级合并 → 连续曲面合并 → 陡壁 RANSAC 剥离 →
  矢量化轮廓重建(网格占用/闭开运算/边界追踪/RDP/凸包规整/矩形吸附/耳切) →
  行生成(对齐矩形边) → 蛇形轨迹。
- 查看器：wpf\SprayTrajViewer（WPF + HelixToolkit 3.1.2），界面参数已全暴露。
- 测试样本：cpp\_regress\*.ply（panel/big/zhe/cao/tray/complex/longstripe/longwide/curved/synth）。
- 编译必须把临时目录指到 E 盘（C 盘空间不足）：
  cd /e/Code/喷涂算法自研/cpp/SprayTrajCore/build && TMP=E:/tmp_build TEMP=E:/tmp_build cmake --build . --config Release
- 同步 DLL 到查看器前必须先杀进程：Get-Process SprayTrajViewer | Stop-Process -Force

【参数语义，改动时不能偏离】
- 平面拟合阈值 flatnessTol(默认3mm) = 多翘的面还算同一个面
- 孔洞填补阈值 fillHoleSize(默认15mm) = 多大的洞当实心覆盖
- scanDirMode: 0自动/1沿长边/2沿短边
- 行距 12mm，点距 4mm，离面 100mm

【当前已完成，不要推倒重来】
棕色槽型盘在 -flat 30 -fill 300 下已收敛为 2 个面：
地面 4 顶点矩形 205×309（15 行）+ 近垂直侧壁 4 顶点矩形 187×56（87.7°，4 行），
行距 3.8~4.0mm 均匀，轨迹不再有梳齿外溢。
折型板 flat=3 保留 3 面；长板从 6 个"法向几乎一致"的碎面合成 2 面。

【当前待办，按优先级】
P0 查看器最终验证没做完：SprayTrajViewer.exe 最后一次启动后进程未驻留、stdout 无日志，
   原因未定位（此前有 WPS 对话框抢占前台导致 UIA 点击错位，但不是程序问题）。
   请重编 WPF 后命令行启动并重定向输出，崩溃则查 Windows 事件日志。
   验证点：棕色槽型盘应显示"1 个完整矩形地面 + 1 条近垂直侧壁"，边线是细直线不是阶梯锯齿，
   轨迹行平行矩形边、均匀铺满、不外溢。若渲染有问题，改 MainWindow.xaml.cs 的 AddPlaneVisual。
P1 曲面件（curved 曲型板）现在 flat=3/30 都被硬切成 8 片平面（94~98 行，rms 到 5.7mm），
   需要新增"按曲率带分条"的曲面路径：沿主曲率方向把曲面切成条带，每带生成 S 行，
   而不是切成互不相连的平面。这是用户"多平面要重建好"的延伸诉求。
P1 槽型盘 cao 在 flat=30 反而出 4 面（flat=3 是 3 面），阈值放大没起到合并作用，需要查。
   用 SPRAYTRAJ_DUMP=1 看 [dbg] 连续合并跳过 / [split] 侧壁面 / [split] 面块 的输出序列。
P1 面间过渡没有避障：现在是 spray=0 直线连接，会穿过工件凸台/槽壁。
   应该先沿法向抬到安全高度 → 平面内平移 → 下降接下一面起点。
P2 棕色槽型盘 flat=3 出 6 面，需导出 tray.ply.plane*.ply 逐个确认是真实结构还是过度切分
   （用户对"面太多"很敏感）。
P2 driver.cpp 不支持非 ASCII 路径（中文文件名读不了），改用 std::filesystem/_wfopen。
P2 性能：大件单件 10~30s，百万点未测。热点是 KNN、build_plane_mesh、genRows 的 5° 穷举回退。

【方法论上的坑，务必遵守】
1. 法向平滑必须保边（只平均与自身法向 ≤30° 的邻居），无差别均值会把薄侧壁拉平。
2. 先合并再按 minPts 过滤，顺序反了侧壁弧段会全丢。
3. 合并判定必须对"固定种子法向"，不能让法向随合并漂移。
4. 合并角门限 8° 是实测标定值（盘底翘曲 0.3~1.0°，折型板折角 13.2°），别随手放大到 25°。
5. 崩缩判定的少数派容忍 8%（托盘斜坡 9.7%，折型板折面 46%），放大到 15% 会把托盘侧壁吞掉。
6. 侧壁法向必须背离工件质心定向，否则喷头伸进工件。
7. 窄条面上最小外接矩形的边顺序不可靠，必须保证"行方向沿面内更长的一维"。
8. 只有网格重建成功的面才计入 plane_count 并生成轨迹，否则前端会画成包围盒。
9. 改 ABI 结构体必须同步 wpf\SprayTrajViewer\Native.cs 的 MarshalAs SizeConst。

【每次改完的验收方式】
跑全样本双阈值回归，对照 HANDOFF 第 5 节的基线表看面数/行数是否退化：
cd cpp/SprayTrajCore/build/bin/Release
for f in panel big zhe cao tray complex longstripe longwide curved synth; do
  ./spray_driver.exe /e/Code/喷涂算法自研/cpp/_regress/$f.ply /tmp/q.csv /tmp/q.ply -flat 30 -fill 300 \
    2>&1 | grep -E "^\[重建\] 面片|^OK:";
done
可视化用 cpp\_regress 下现成脚本：final_check.py（覆盖图）、plane_rows.py（逐面铺平）、
view3d.py（3D）、row_metrics.py（行规整度）。
```

---

## 补充说明（给自己看，不用复制）

- 如果新窗口的 AI 想快速确认"当前代码到底是什么状态"，让它先跑一次
  `./spray_driver.exe tray.ply /tmp/a.csv /tmp/a.ply -flat 30 -fill 300`，
  看到 `2 面 / 19 行` 就说明代码是这一版。
- 如果看到 `6 面`，说明它拿到的是旧编译产物，先重新编译。
