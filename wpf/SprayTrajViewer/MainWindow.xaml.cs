using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using Microsoft.Win32;
using HelixToolkit.Wpf;

namespace SprayTrajViewer
{
    public partial class MainWindow : Window
    {
        private List<(double x, double y, double z)> _cloud = new();
        private string? _currentFile;
        // 最近一次分割出的所有重建平面（含法向/质心/局部基/包围盒），用于"显示重建平面"
        private List<SprayTrajPlane> _planes = new();
        // 当前重建平面的视觉元素列表（精确引用，用于增删）
        private readonly List<Visual3D> _planeVisuals = new();
        // 点云与轨迹的视觉元素引用（用于独立显隐）
        private PointsVisual3D? _cloudVisual;
        private readonly List<TubeVisual3D> _trajVisuals = new();

        public MainWindow()
        {
            InitializeComponent();
            Loaded += (_, _) =>
            {
                // 窗口布局完成后再加载与取景，避免相机没对准
                try
                {
                    OnLoadedRunCli();
                }
                catch (Exception ex)
                {
                    App.Log("FATAL/Loaded", ex);
                }
            };
        }

        /// <summary>
        /// 把 3D 视图直接渲染成 PNG。
        /// 相比外部截屏：不受窗口被遮挡、焦点被抢、用户拖动相机、缩放比例影响，可作为可复现的验收证据。
        /// </summary>
        void SaveViewportShot(string path)
        {
            Viewport.UpdateLayout();
            double w = Viewport.ActualWidth, h = Viewport.ActualHeight;
            if (w < 2 || h < 2) { App.Log("SHOT/FAIL", null, $"视口尺寸异常 {w}x{h}"); return; }

            var rtb = new RenderTargetBitmap((int)w, (int)h, 96, 96, PixelFormats.Pbgra32);
            rtb.Render(Viewport);

            var enc = new PngBitmapEncoder();
            enc.Frames.Add(BitmapFrame.Create(rtb));
            using (var fs = System.IO.File.Create(path)) enc.Save(fs);
            App.Log("SHOT", null, $"{path} ({rtb.PixelWidth}x{rtb.PixelHeight})");
        }

        // ── 命令行自运行模式 ────────────────────────────────────────────────
        // 目的：让"界面渲染是否正确"的验证不依赖 UI 自动化点击（此前被 WPS 对话框抢焦点搞错过），
        // 一条命令即可复现指定样本 + 指定参数 + 指定显示项，再用截图取证。
        //   SprayTrajViewer.exe --ply <file> --flat 30 --fill 300 --auto
        //   --show-plane 1 --show-cloud 1 --show-traj 1 --log <file>
        // 无参数时行为与原来完全一致（纯交互模式）。

        static Dictionary<string, string> ParseArgs()
        {
            var d = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            string[] a = Environment.GetCommandLineArgs();
            for (int i = 1; i < a.Length; i++)
            {
                string k = a[i];
                if (!k.StartsWith("--")) { d["pos" + d.Count] = k; continue; }
                k = k.Substring(2);
                string v = (i + 1 < a.Length && !a[i + 1].StartsWith("--")) ? a[++i] : "1";
                d[k] = v;
            }
            return d;
        }

        void OnLoadedRunCli()
        {
            var o = ParseArgs();
            bool hadArgs = o.Count > 0;

            // 参数（只覆盖显式传入的项，其余保持界面默认值）
            if (o.TryGetValue("ply", out var ply)) TxtPath.Text = ply;
            if (o.TryGetValue("row", out var v)) TxtRowSpacing.Text = v;
            if (o.TryGetValue("pt", out v)) TxtPointSpacing.Text = v;
            if (o.TryGetValue("standoff", out v)) TxtStandoff.Text = v;
            if (o.TryGetValue("speed", out v)) TxtSpeed.Text = v;
            if (o.TryGetValue("turn", out v)) TxtTurnFactor.Text = v;
            if (o.TryGetValue("flat", out v)) TxtFlatness.Text = v;
            if (o.TryGetValue("curv", out v)) TxtCurvature.Text = v;
            if (o.TryGetValue("sorK", out v)) TxtSorK.Text = v;
            if (o.TryGetValue("sorZ", out v)) TxtSorZ.Text = v;
            if (o.TryGetValue("fill", out v)) TxtFillHole.Text = v;
            if (o.TryGetValue("dir", out v)) CmbScanDir.SelectedIndex = Math.Clamp(int.TryParse(v, out var di) ? di : 0, 0, 2);
            if (o.TryGetValue("flip", out v)) CmbNormalFlip.SelectedIndex = Math.Clamp(int.TryParse(v, out var fi) ? fi : 0, 0, 1);

            // 显示项（显式给值时才改，避免把用户默认勾选状态打乱）
            if (o.TryGetValue("show-cloud", out v)) ChkShowCloud.IsChecked = v != "0";
            if (o.TryGetValue("show-traj", out v)) ChkShowTraj.IsChecked = v != "0";
            if (o.TryGetValue("show-plane", out v)) ChkShowPlane.IsChecked = v != "0";

            App.Log("CLI", null, $"ply={TxtPath.Text}");
            App.Log("CLI", null,
                $"flat={TxtFlatness.Text} curv={TxtCurvature.Text} fill={TxtFillHole.Text} " +
                $"row={TxtRowSpacing.Text} pt={TxtPointSpacing.Text} standoff={TxtStandoff.Text} " +
                $"dir={CmbScanDir.SelectedIndex} flip={CmbNormalFlip.SelectedIndex}");
            App.Log("CLI", null,
                $"show: cloud={ChkShowCloud.IsChecked} traj={ChkShowTraj.IsChecked} plane={ChkShowPlane.IsChecked}");

            if (!System.IO.File.Exists(TxtPath.Text))
            {
                App.Log("CLI/WARN", null, "点云文件不存在：" + TxtPath.Text);
            }
            else
            {
                LoadCloud(TxtPath.Text);
            }

            if (hadArgs && o.ContainsKey("auto"))
            {
                // 让窗口先完成一次布局/绘制，再跑算法与渲染，避免在 Loaded 里做重活把首帧拖成黑屏
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    try
                    {
                        BtnGenerate_Click(this, new RoutedEventArgs());
                        App.Log("CLI", null, "auto-generate 完成");
                    }
                    catch (Exception ex)
                    {
                        App.Log("FATAL/AutoGenerate", ex);
                    }

                    if (o.TryGetValue("view", out var vm)) ApplyViewPreset(vm);

                    if (o.TryGetValue("shot", out var shotPath))
                    {
                        // 等一帧稳定（相机/布局/3D 绘制都完成）再截，避免抓到半成品
                        var timer = new System.Windows.Threading.DispatcherTimer
                        {
                            Interval = TimeSpan.FromMilliseconds(700)
                        };
                        timer.Tick += (_, _) =>
                        {
                            timer.Stop();
                            try { SaveViewportShot(shotPath); }
                            catch (Exception ex) { App.Log("SHOT/FAIL", ex); }
                            if (o.ContainsKey("exit"))
                            {
                                App.Log("CLI", null, "shot 完成，退出");
                                Application.Current.Shutdown(0);
                            }
                        };
                        timer.Start();
                    }
                }), System.Windows.Threading.DispatcherPriority.Background);
            }
            else if (o.TryGetValue("view", out var v2))
            {
                ApplyViewPreset(v2);
            }
        }

        void BtnBrowse_Click(object sender, RoutedEventArgs e)
        {
            var dlg = new OpenFileDialog { Filter = "PLY 点云|*.ply|所有文件|*.*" };
            if (dlg.ShowDialog() == true)
            {
                TxtPath.Text = dlg.FileName;
                LoadCloud(dlg.FileName);
            }
        }

        void LoadCloud(string path)
        {
            try
            {
                _cloud = PlyReader.Read(path);
                _currentFile = path;
                RenderCloud();
                TxtHint.Text = $"点云已加载：{_cloud.Count} 点";
                Title = $"喷涂轨迹生成器 - {System.IO.Path.GetFileName(path)}";
                App.Log("LOAD", null, $"{path} -> {_cloud.Count} 点");
            }
            catch (Exception ex)
            {
                App.Log("LOAD/FAIL", ex);
                MessageBox.Show($"点云读取失败：{ex.Message}", "错误", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        void RenderCloud()
        {
            if (_cloud.Count == 0) return;
            // 移除旧点云（保留轨迹、平面）
            if (_cloudVisual != null) { Viewport.Children.Remove(_cloudVisual); _cloudVisual = null; }

            var pts = new Point3DCollection();
            foreach (var p in _cloud.Take(200000))
                pts.Add(new Point3D(p.x, p.y, p.z));

            _cloudVisual = new PointsVisual3D
            {
                Points = pts,
                Color = Color.FromRgb(224, 64, 224),
                Size = 2.5
            };
            Viewport.Children.Add(_cloudVisual);
            Viewport.ZoomExtents(200);
            ApplyCloudVisibility();
        }

        void BtnGenerate_Click(object sender, RoutedEventArgs e)
        {
            if (_cloud.Count == 0)
            {
                MessageBox.Show("请先加载点云文件", "提示", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            SprayTrajParams p = new()
            {
                scanRowSpacing = ParseDouble(TxtRowSpacing.Text, 12),
                pointSpacing   = ParseDouble(TxtPointSpacing.Text, 4),
                standoff       = ParseDouble(TxtStandoff.Text, 100),
                speed          = ParseDouble(TxtSpeed.Text, 200),
                turnRadiusFactor = ParseDouble(TxtTurnFactor.Text, 0.5),
                flatnessTol    = ParseDouble(TxtFlatness.Text, 3),
                curvatureTol   = ParseDouble(TxtCurvature.Text, 5),
                sorK           = (int)ParseDouble(TxtSorK.Text, 20),
                sorZ           = ParseDouble(TxtSorZ.Text, 2.5),
                normalFlip     = CmbNormalFlip.SelectedIndex == 1 ? 1 : 0,
                fillHoleSize   = ParseDouble(TxtFillHole.Text, 15),
                scanDirMode    = CmbScanDir.SelectedIndex,
            };

            // 组织 xyz 数组
            var xyz = new double[_cloud.Count * 3];
            for (int i = 0; i < _cloud.Count; i++)
            {
                xyz[i * 3]     = _cloud[i].x;
                xyz[i * 3 + 1] = _cloud[i].y;
                xyz[i * 3 + 2] = _cloud[i].z;
            }

            SprayTrajResult res;
            try
            {
                res = Native.spray_traj_generate(xyz, _cloud.Count, ref p);
            }
            catch (DllNotFoundException)
            {
                MessageBox.Show("未找到 SprayTraj.dll，请先编译 C++ 核心库并复制到程序目录。", "错误", MessageBoxButton.OK, MessageBoxImage.Error);
                return;
            }
            catch (Exception ex)
            {
                MessageBox.Show($"调用失败：{ex.Message}", "错误", MessageBoxButton.OK, MessageBoxImage.Error);
                return;
            }

            if (res.error != IntPtr.Zero)
            {
                string err = Marshal.PtrToStringAnsi(res.error) ?? "未知错误";
                TxtResult.Text = $"失败：{err}";
                Native.spray_traj_free(ref res);
                return;
            }

            // 读取轨迹点
            var traj = new List<SprayTrajPoint>();
            for (int i = 0; i < res.count; i++)
            {
                IntPtr ptr = IntPtr.Add(res.points, i * Marshal.SizeOf<SprayTrajPoint>());
                traj.Add(Marshal.PtrToStructure<SprayTrajPoint>(ptr));
            }

            RenderTrajectory(traj);

            // 保存所有分割平面参数，供"显示重建平面"使用
            _planes.Clear();
            for (int i = 0; i < res.plane_count && i < 8; i++)
                _planes.Add(res.planes[i]);
            // 兼容：若 plane_count 未填充（旧 DLL），回退用主平面
            if (_planes.Count == 0)
            {
                _planes.Add(new SprayTrajPlane
                {
                    nx = res.normal_x, ny = res.normal_y, nz = res.normal_z,
                    cx = res.plane_cx, cy = res.plane_cy, cz = res.plane_cz,
                    count = res.inlier_points,
                    minU = -1, maxU = 1, minV = -1, maxV = 1
                });
            }

            App.Log("GEN", null,
                $"plane_count={res.plane_count} rows={res.rows} points={res.count} " +
                $"rms={res.plane_rms:F3} inlier={res.inlier_points} " +
                $"seg={res.seg_len_mean:F3}/{res.seg_len_max:F3}");
            LogPlaneSummary();

            RenderPlane();

            TxtResult.Text =
                $"✓ 生成成功\n" +
                $"  分割面数：{Math.Max(res.plane_count, 1)}\n" +
                $"  行数：{res.rows}\n" +
                $"  轨迹点：{res.count}\n" +
                $"  平面残差 RMS：{res.plane_rms:F3} mm\n" +
                $"  主面法向：({res.normal_x:F3}, {res.normal_y:F3}, {res.normal_z:F3})\n" +
                $"  内点：{res.inlier_points}（占去噪后 {(double)res.inlier_points / Math.Max(1, res.input_points - res.removed_points) * 100:F1}%）\n" +
                $"  点距：均值 {res.seg_len_mean:F3} / 最大 {res.seg_len_max:F3} mm\n" +
                $"  输入点：{res.input_points}，去噪剔除：{res.removed_points}\n" +
                $"  红色=喷涂行，灰色=空移（面间抬枪绕行）";

            TxtHint.Text = $"轨迹已生成：{res.count} 点，{res.rows} 行";

            Native.spray_traj_free(ref res);

            // 同时写出 CSV 供外部使用
            try { ExportCsv(traj); } catch { }
        }

        void RenderTrajectory(List<SprayTrajPoint> traj)
        {
            if (traj.Count == 0) return;
            // 移除旧轨迹（保留点云与平面）
            foreach (var t in _trajVisuals) Viewport.Children.Remove(t);
            _trajVisuals.Clear();

            // 按"段"分类，而不是按"点"分类：
            //   只有【两端都是 spray=1】的段才算喷涂段（红色）；
            //   只要有一端是 spray=0，这一段就是空移（灰色）。
            // 为什么：关枪/开枪发生在点上，枪在 A→B 这一段上的状态是"关"。
            // 早先按点分类，会把"抬枪上升"和"下降到首喷点"这两段染成红色，看着像喷着飞过去。
            const double breakDist = 15.0;   // 相邻点距约 4mm，超过 15mm 视为断开（不应出现）

            var runCls = new List<int>();                 // 1=喷涂 0=空移
            var runPts = new List<List<Point3D>>();
            List<Point3D> acc = null;
            int accCls = -1;

            for (int i = 1; i < traj.Count; i++)
            {
                var a = new Point3D(traj[i - 1].x, traj[i - 1].y, traj[i - 1].z);
                var b = new Point3D(traj[i].x, traj[i].y, traj[i].z);
                double d = (b - a).Length;
                if (d < 1e-6) continue;                    // 重复点

                int cls = (traj[i].spray == 1 && traj[i - 1].spray == 1) ? 1 : 0;
                bool jump = d > breakDist;

                if (acc == null || cls != accCls || jump)
                {
                    if (acc != null && acc.Count >= 2) { runCls.Add(accCls); runPts.Add(acc); }
                    acc = new List<Point3D> { a, b };
                    accCls = cls;
                }
                else acc.Add(b);
            }
            if (acc != null && acc.Count >= 2) { runCls.Add(accCls); runPts.Add(acc); }

            for (int i = 0; i < runPts.Count; i++)
            {
                bool spraying = runCls[i] == 1;
                _trajVisuals.Add(new TubeVisual3D
                {
                    Path = new Point3DCollection(runPts[i]),
                    Diameter = spraying ? 2.5 : 1.1,
                    Fill = spraying
                        ? new SolidColorBrush(Color.FromRgb(255, 48, 32))
                        : new SolidColorBrush(Color.FromRgb(150, 158, 170)),
                    ThetaDiv = spraying ? 6 : 4
                });
            }
            foreach (var t in _trajVisuals) Viewport.Children.Add(t);

            Viewport.ZoomExtents(0);   // 立即取景到全部内容
            ApplyTrajVisibility();

            App.Log("TRAJ", null,
                $"段数={runPts.Count}（喷涂 {runCls.Count(c => c == 1)} / 空移 {runCls.Count(c => c == 0)}）");
        }

        // 记录 C++ 侧实际交过来的面片几何：用于区分"算法分错了"还是"查看器画错了"
        void LogPlaneSummary()
        {
            for (int i = 0; i < _planes.Count; i++)
            {
                var p = _planes[i];
                double nx = p.nx, ny = p.ny, nz = p.nz;
                double angZ = Math.Acos(Math.Clamp(Math.Abs(nz) / Math.Max(1e-9, Math.Sqrt(nx * nx + ny * ny + nz * nz)), 0, 1)) * 180.0 / Math.PI;

                // 顶点包围盒尺寸（以质心为原点的局部坐标，直接给世界坐标尺寸更有用）
                double minX = double.MaxValue, maxX = double.MinValue;
                double minY = double.MaxValue, maxY = double.MinValue;
                double minZ = double.MaxValue, maxZ = double.MinValue;
                int nv = Math.Min(p.polyVertCount, p.polyVerts?.Length / 3 ?? 0);
                for (int k = 0; k < nv; k++)
                {
                    double x = p.polyVerts![k * 3], y = p.polyVerts[k * 3 + 1], z = p.polyVerts[k * 3 + 2];
                    if (x < minX) minX = x; if (x > maxX) maxX = x;
                    if (y < minY) minY = y; if (y > maxY) maxY = y;
                    if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
                }
                string size = nv > 0
                    ? $"尺寸={maxX - minX:F1}x{maxY - minY:F1}x{maxZ - minZ:F1}"
                    : "尺寸=<无顶点>";

                App.Log("PLANE", null,
                    $"#{i} 顶点={p.polyVertCount} 三角={p.polyTriCount} {size} " +
                    $"内点={p.count} 与Z夹角={angZ:F1}° 质心=({p.cx:F1},{p.cy:F1},{p.cz:F1}) " +
                    $"法向=({nx:F3},{ny:F3},{nz:F3})");
            }
        }

        // 渲染重建平面（半透明网格面片）
        void RenderPlane()
        {
            // 精确移除旧平面（只删平面的 ModelVisual3D，不动点云/轨迹）
            foreach (var v in _planeVisuals) Viewport.Children.Remove(v);
            _planeVisuals.Clear();

            if (!ChkShowPlane.IsChecked.GetValueOrDefault()) return;
            foreach (var p in _planes)
                AddPlaneVisual(p);
        }

        void AddPlaneVisual(SprayTrajPlane p)
        {
            var n = new Vector3D(p.nx, p.ny, p.nz);
            if (n.Length > 1e-9) n.Normalize();
            var c = new Point3D(p.cx, p.cy, p.cz);

            Vector3D u = new Vector3D(p.ux, p.uy, p.uz);
            Vector3D v = new Vector3D(p.vx, p.vy, p.vz);
            if (u.Length < 1e-6 || v.Length < 1e-6)
            {
                u = Vector3D.CrossProduct(new Vector3D(0, 0, 1), n);
                if (u.Length < 1e-6) u = new Vector3D(1, 0, 0);
                u.Normalize();
                v = Vector3D.CrossProduct(n, u);
                v.Normalize();
            }
            else { u.Normalize(); v.Normalize(); }

            // 直接用核心库重建好的三角网格渲染：
            // 轮廓已矢量化（直边/矩形吸附），不再是逐 cell 小方片拼出来的锯齿拼块。
            var mesh = new MeshGeometry3D();
            int nv = p.polyVertCount, nt = p.polyTriCount;
            bool usedMesh = p.polyVerts != null && p.polyTris != null && nv >= 3 && nt >= 1;
            if (usedMesh)
            {
                for (int i = 0; i < nv && i * 3 + 2 < p.polyVerts.Length; i++)
                    mesh.Positions.Add(new Point3D(p.polyVerts[i * 3], p.polyVerts[i * 3 + 1], p.polyVerts[i * 3 + 2]));
                for (int t = 0; t < nt && t * 3 + 2 < p.polyTris.Length; t++)
                {
                    int a = p.polyTris[t * 3], b = p.polyTris[t * 3 + 1], cc = p.polyTris[t * 3 + 2];
                    if (a < 0 || b < 0 || cc < 0 || a >= nv || b >= nv || cc >= nv) continue;
                    mesh.TriangleIndices.Add(a);
                    mesh.TriangleIndices.Add(b);
                    mesh.TriangleIndices.Add(cc);
                }
            }
            else
            {
                // 兜底：无网格数据时画包围盒矩形
                double minU = p.minU, maxU = p.maxU, minV = p.minV, maxV = p.maxV;
                if (maxU - minU < 1e-6 || maxV - minV < 1e-6) return;
                mesh.Positions.Add(c + u * minU + v * minV);
                mesh.Positions.Add(c + u * maxU + v * minV);
                mesh.Positions.Add(c + u * maxU + v * maxV);
                mesh.Positions.Add(c + u * minU + v * maxV);
                mesh.TriangleIndices.Add(0); mesh.TriangleIndices.Add(1); mesh.TriangleIndices.Add(2);
                mesh.TriangleIndices.Add(0); mesh.TriangleIndices.Add(2); mesh.TriangleIndices.Add(3);
            }
            if (mesh.Positions.Count == 0 || mesh.TriangleIndices.Count == 0) return;

            var mat = new DiffuseMaterial(new SolidColorBrush(Color.FromArgb(78, 64, 156, 255)));
            var model = new GeometryModel3D(mesh, mat) { BackMaterial = mat };
            var visual = new ModelVisual3D { Content = model };
            _planeVisuals.Add(visual);
            Viewport.Children.Add(visual);

            // 轮廓线：三角网格中"只被一个三角形使用"的边就是边界（含孔洞边界），
            // 描一圈深色线让重建面的形状一目了然，不会糊成一团半透明蓝色。
            var edgeCount = new Dictionary<(int, int), int>();
            for (int t = 0; t + 2 < mesh.TriangleIndices.Count; t += 3)
            {
                int a = mesh.TriangleIndices[t], b = mesh.TriangleIndices[t + 1], cc = mesh.TriangleIndices[t + 2];
                foreach (var (i, j) in new[] { (a, b), (b, cc), (cc, a) })
                {
                    var key = i < j ? (i, j) : (j, i);
                    edgeCount[key] = edgeCount.TryGetValue(key, out var cnt) ? cnt + 1 : 1;
                }
            }
            var pts = new Point3DCollection();
            foreach (var kv in edgeCount)
            {
                if (kv.Value != 1) continue;
                pts.Add(mesh.Positions[kv.Key.Item1]);
                pts.Add(mesh.Positions[kv.Key.Item2]);
            }
            if (pts.Count > 0)
            {
                var lines = new LinesVisual3D
                {
                    Points = pts,
                    Color = Color.FromRgb(32, 82, 150),
                    Thickness = 1.6
                };
                _planeVisuals.Add(lines);
                Viewport.Children.Add(lines);
            }

            App.Log("RENDER", null,
                $"{(usedMesh ? "网格" : "包围盒兜底")} 顶点={mesh.Positions.Count} " +
                $"三角={mesh.TriangleIndices.Count / 3} 边界边={pts.Count / 2} " +
                $"尺寸={FmtLen(mesh)}");
        }

        void ExportCsv(List<SprayTrajPoint> traj)
        {
            var path = System.IO.Path.Combine(
                System.IO.Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location)!,
                "trajectory.csv");
            var sb = new System.Text.StringBuilder("X,Y,Z,I,J,K,speed,spray\n");
            foreach (var p in traj)
                sb.Append($"{p.x:F3},{p.y:F3},{p.z:F3},{p.i:F4},{p.j:F4},{p.k:F4},{p.speed:F1},{p.spray}\n");
            System.IO.File.WriteAllText(path, sb.ToString());
            TxtHint.Text += $"\nCSV 已导出：{path}";
        }

        static double ParseDouble(string s, double def)
        {
            return double.TryParse(s, out var v) ? v : def;
        }

        static string FmtLen(MeshGeometry3D m)
        {
            if (m.Positions.Count == 0) return "-";
            double minX = double.MaxValue, maxX = double.MinValue;
            double minY = double.MaxValue, maxY = double.MinValue;
            double minZ = double.MaxValue, maxZ = double.MinValue;
            foreach (var p in m.Positions)
            {
                if (p.X < minX) minX = p.X; if (p.X > maxX) maxX = p.X;
                if (p.Y < minY) minY = p.Y; if (p.Y > maxY) maxY = p.Y;
                if (p.Z < minZ) minZ = p.Z; if (p.Z > maxZ) maxZ = p.Z;
            }
            return $"{maxX - minX:F1}x{maxY - minY:F1}x{maxZ - minZ:F1}";
        }

        /// <summary>
        /// 场景包围盒：直接对点云、面片顶点、轨迹点求并集。
        /// 不用 HelixToolkit 的 Viewport3DHelper.GetBounds（3.1.2 没有该 API），也不依赖视觉树遍历，
        /// 便于在不同 Helix 版本下保持行为一致。
        /// </summary>
        Rect3D? ComputeSceneBounds()
        {
            double minX = double.MaxValue, maxX = double.MinValue;
            double minY = double.MaxValue, maxY = double.MinValue;
            double minZ = double.MaxValue, maxZ = double.MinValue;
            bool any = false;

            void Acc(double x, double y, double z)
            {
                if (double.IsNaN(x) || double.IsNaN(y) || double.IsNaN(z)) return;
                if (x < minX) minX = x; if (x > maxX) maxX = x;
                if (y < minY) minY = y; if (y > maxY) maxY = y;
                if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
                any = true;
            }

            foreach (var p in _cloud) Acc(p.x, p.y, p.z);
            foreach (var pl in _planes)
            {
                int nv = Math.Min(pl.polyVertCount, pl.polyVerts?.Length / 3 ?? 0);
                for (int k = 0; k < nv; k++)
                    Acc(pl.polyVerts![k * 3], pl.polyVerts[k * 3 + 1], pl.polyVerts[k * 3 + 2]);
            }
            foreach (var t in _trajVisuals)
                foreach (var p in t.Path) Acc(p.X, p.Y, p.Z);

            if (!any) return null;
            return new Rect3D(minX, minY, minZ, maxX - minX, maxY - minY, maxZ - minZ);
        }

        /// <summary>
        /// 相机预设：让界面验证可复现（同一样本每次都用同一视角截图）。
        /// 轴向视图（top/front/side）用正交相机：轨迹悬在离面 standoff 处，
        /// 透视相机下会被按距离放大（615/(615-100)≈1.19 倍），看起来像"轨迹溢出重建面"——
        /// 那是投影失真，不是算法外溢。正交投影没有这个问题，可用于测量/验收。
        /// iso 仍用透视，符合人眼观感。
        /// </summary>
        void ApplyViewPreset(string mode)
        {
            if (string.IsNullOrWhiteSpace(mode)) return;

            var b = ComputeSceneBounds();
            if (b is not Rect3D r || r.IsEmpty) { App.Log("VIEW/WARN", null, "场景为空，跳过视角预设"); return; }

            var center = new Point3D(r.X + r.SizeX / 2, r.Y + r.SizeY / 2, r.Z + r.SizeZ / 2);
            double span = Math.Max(Math.Max(r.SizeX, r.SizeY), r.SizeZ);
            if (span <= 1e-6) return;

            string m = mode.ToLowerInvariant();
            Vector3D dir = m switch
            {
                "top"    => new Vector3D(0, 0, 1),
                "bottom" => new Vector3D(0, 0, -1),
                "front"  => new Vector3D(0, -1, 0.0001),
                "side"   => new Vector3D(1, 0, 0.0001),
                _        => new Vector3D(-0.55, -0.72, 0.42),   // iso
            };
            dir.Normalize();
            Vector3D up = Math.Abs(dir.Z) > 0.99 ? new Vector3D(0, 1, 0) : new Vector3D(0, 0, 1);
            bool ortho = m is "top" or "bottom" or "front" or "side";

            if (ortho)
            {
                // 正交：相机推得足够远（正交投影与距离无关），宽度略大于场景最大跨度
                double back = Math.Max(span * 3.0, 1000.0);
                var cam = new OrthographicCamera
                {
                    Position = center - dir * back,
                    LookDirection = dir * back,
                    UpDirection = up,
                    Width = span * 1.30,
                    NearPlaneDistance = 1.0,
                    FarPlaneDistance = back * 4.0,
                };
                Viewport.Camera = cam;
                App.Log("VIEW", null,
                    $"正交 {m} center=({center.X:F1},{center.Y:F1},{center.Z:F1}) span={span:F1} width={span * 1.30:F1}");
            }
            else
            {
                double halfFov = 45.0;
                double dist = span / (2.0 * Math.Tan(halfFov * Math.PI / 360.0)) * 1.25;
                var cam = new PerspectiveCamera
                {
                    Position = center - dir * dist,
                    LookDirection = dir * dist,
                    UpDirection = up,
                    FieldOfView = halfFov,
                    NearPlaneDistance = Math.Max(0.5, dist * 0.001),
                    FarPlaneDistance = dist * 20.0,
                };
                Viewport.Camera = cam;
                App.Log("VIEW", null,
                    $"透视 {m} center=({center.X:F1},{center.Y:F1},{center.Z:F1}) span={span:F1} dist={dist:F1}");
            }
        }

        void ChkShowPlane_Changed(object sender, RoutedEventArgs e)
        {
            RenderPlane();
        }

        void ChkShowCloud_Changed(object sender, RoutedEventArgs e)
        {
            ApplyCloudVisibility();
        }

        void ChkShowTraj_Changed(object sender, RoutedEventArgs e)
        {
            ApplyTrajVisibility();
        }

        // Visual3D 无 Visibility 属性，通过增删 Children 实现显隐
        void ApplyCloudVisibility()
        {
            if (_cloudVisual == null) return;
            bool show = ChkShowCloud.IsChecked.GetValueOrDefault(true);
            bool inTree = Viewport.Children.Contains(_cloudVisual);
            if (show && !inTree) Viewport.Children.Add(_cloudVisual);
            else if (!show && inTree) Viewport.Children.Remove(_cloudVisual);
        }

        void ApplyTrajVisibility()
        {
            bool show = ChkShowTraj.IsChecked.GetValueOrDefault(true);
            foreach (var t in _trajVisuals)
            {
                bool inTree = Viewport.Children.Contains(t);
                if (show && !inTree) Viewport.Children.Add(t);
                else if (!show && inTree) Viewport.Children.Remove(t);
            }
        }
    }
}
