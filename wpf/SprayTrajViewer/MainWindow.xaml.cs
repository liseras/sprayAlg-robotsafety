using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
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

        // 可喷面角度窗口，必须与 C++ 的 SPRAYTRAJ_NORMAL_WIN_DEG 保持一致：
        // 法向与 +Z 夹角 ≤ 此值才列入清单。130° = 朝上、侧向、最多向下倾 40° 都算可喷，
        // 开口盒底板外侧（180°）被排除。
        private const double NormalWinDeg = 130.0;
        // 当前重建平面的视觉元素列表（精确引用，用于增删）
        private readonly List<Visual3D> _planeVisuals = new();
        // 点云与轨迹的视觉元素引用（用于独立显隐）
        private PointsVisual3D? _cloudVisual;
        private readonly List<TubeVisual3D> _trajVisuals = new();
        // 命令行 --keep-cloud：生成后不自动收起点云，用于"点云抽样是否正确"的验收截图
        private bool _keepCloudOnGenerate;

        public MainWindow()
        {
            InitializeComponent();
            // ABI 自检：与 C++ `spray_safety.exe --abi` 打印的 sizeof 对齐
            // （Plane=18664 Point=64 Result=298744 Params=128；SafetyParams=152 SafetyPoint=216 SafetyResult=304）
            try
            {
                App.Log("ABI", null,
                    $"Plane={Marshal.SizeOf<SprayTrajPlane>()} " +
                    $"Point={Marshal.SizeOf<SprayTrajPoint>()} " +
                    $"Result={Marshal.SizeOf<SprayTrajResult>()} " +
                    $"Params={Marshal.SizeOf<SprayTrajParams>()} | " +
                    $"SafetyParams={Marshal.SizeOf<RobotSafetyParams>()} " +
                    $"SafetyPoint={Marshal.SizeOf<RobotSafetyPoint>()} " +
                    $"SafetyResult={Marshal.SizeOf<RobotSafetyResult>()} " +
                    $"safety@{Marshal.OffsetOf<RobotSafetyParams>(nameof(RobotSafetyParams.safety))} " +
                    $"hardClearance@{Marshal.OffsetOf<RobotSafetyParams>(nameof(RobotSafetyParams.hardClearance))}");
            }
            catch (Exception ex) { App.Log("ABI/FAIL", ex); }
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

        /// <summary>
        /// 整窗截图：包含左侧参数面板与右上角安全结论卡。
        /// Viewport 截图只含 3D 内容，验收「卡片显示了什么」必须用这个。
        /// </summary>
        void SaveWindowShot(string path)
        {
            UpdateLayout();
            double w = ActualWidth, h = ActualHeight;
            if (w < 2 || h < 2) { App.Log("SHOT-WIN/FAIL", null, $"窗口尺寸异常 {w}x{h}"); return; }

            var rtb = new RenderTargetBitmap((int)w, (int)h, 96, 96, PixelFormats.Pbgra32);
            rtb.Render(this);

            var enc = new PngBitmapEncoder();
            enc.Frames.Add(BitmapFrame.Create(rtb));
            using (var fs = System.IO.File.Create(path)) enc.Save(fs);
            App.Log("SHOT-WIN", null, $"{path} ({rtb.PixelWidth}x{rtb.PixelHeight})");
        }

        // ── 命令行自运行模式 ────────────────────────────────────────────────
        // 目的：让"界面渲染是否正确"的验证不依赖 UI 自动化点击（此前被 WPS 对话框抢焦点搞错过），
        // 一条命令即可复现指定样本 + 指定参数 + 指定显示项，再用截图取证。
        //   SprayTrajViewer.exe --ply <file> --flat 30 --fill 300 --auto
        //   --show-plane 1 --show-cloud 1 --show-traj 1 --log <file>
        //   --view top|bottom|front|side|iso|face   （face=沿面0法向正视，曲面件验收用）
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
            _keepCloudOnGenerate = o.ContainsKey("keep-cloud");

            // 参数（只覆盖显式传入的项，其余保持界面默认值）
            if (o.TryGetValue("ply", out var ply)) TxtPath.Text = ply;
            if (o.TryGetValue("row", out var v)) TxtRowSpacing.Text = v;
            if (o.TryGetValue("pt", out v)) TxtPointSpacing.Text = v;
            if (o.TryGetValue("standoff", out v)) TxtStandoff.Text = v;
            if (o.TryGetValue("gc", out v)) TxtGunClear.Text = v;
            if (o.TryGetValue("ee", out v)) TxtEdgeExtend.Text = v;
            if (o.TryGetValue("ea", out v)) TxtEdgeExtendAcross.Text = v;
            if (o.TryGetValue("st", out v)) TxtShellThick.Text = v;
            if (o.TryGetValue("speed", out v)) TxtSpeed.Text = v;
            if (o.TryGetValue("turn", out v)) TxtTurnFactor.Text = v;
            if (o.TryGetValue("flat", out v)) TxtFlatness.Text = v;
            if (o.TryGetValue("curv", out v)) TxtCurvature.Text = v;
            if (o.TryGetValue("sorK", out v)) TxtSorK.Text = v;
            if (o.TryGetValue("sorZ", out v)) TxtSorZ.Text = v;
            if (o.TryGetValue("fill", out v)) TxtFillHole.Text = v;
            if (o.TryGetValue("dir", out v)) CmbScanDir.SelectedIndex = Math.Clamp(int.TryParse(v, out var di) ? di : 0, 0, 2);
            if (o.TryGetValue("flip", out v)) CmbNormalFlip.SelectedIndex = Math.Clamp(int.TryParse(v, out var fi) ? fi : 0, 0, 1);
            // --tab 1：启动就停在「机器人安全」页（UI 探针截图用）
            if (o.TryGetValue("tab", out v)) Tabs.SelectedIndex = Math.Clamp(int.TryParse(v, out var ti) ? ti : 0, 0, 1);
            // --azmode 1：切回「只用向上」的旧口径（与 spray_safety -azmode 1 等价）
            if (o.TryGetValue("azmode", out v)) CmbAzMode.SelectedIndex = (int.TryParse(v, out var azm) && azm == 1) ? 1 : 0;
            // 注意：--safety 是"跑安全评估"的开关（无值），不能拿它覆盖告警阈值，
            // 否则 ParseArgs 会把无值开关当成 "1"，把 20mm 告警阈值误改成 1mm。
            // 要覆盖阈值请用 --safety-gap 20。
            if (o.TryGetValue("safety-gap", out v)) TxtSafetyGap.Text = v;
            if (o.TryGetValue("hardclr", out v)) TxtHardClearance.Text = v;

            // 显示项（显式给值时才改，避免把用户默认勾选状态打乱）
            if (o.TryGetValue("show-cloud", out v)) ChkShowCloud.IsChecked = v != "0";
            if (o.TryGetValue("show-traj", out v)) ChkShowTraj.IsChecked = v != "0";
            if (o.TryGetValue("show-plane", out v)) ChkShowPlane.IsChecked = v != "0";

            App.Log("CLI", null, $"ply={TxtPath.Text}");
            App.Log("CLI", null,
                $"flat={TxtFlatness.Text} curv={TxtCurvature.Text} fill={TxtFillHole.Text} " +
                $"row={TxtRowSpacing.Text} pt={TxtPointSpacing.Text} standoff={TxtStandoff.Text} " +
                $"gc={TxtGunClear.Text} ee={TxtEdgeExtend.Text} ea={TxtEdgeExtendAcross.Text} st={TxtShellThick.Text} " +
                $"dir={CmbScanDir.SelectedIndex} flip={CmbNormalFlip.SelectedIndex}");
            App.Log("CLI", null,
                $"show: cloud={ChkShowCloud.IsChecked} traj={ChkShowTraj.IsChecked} plane={ChkShowPlane.IsChecked}");
            App.Log("CLI", null, $"safetyGap={TxtSafetyGap.Text} hardClr={TxtHardClearance.Text}");

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
                // CLI 模式由这里控制时序（生成 → 可选安全评估 → 视角 → 截图），
                // 关掉「生成后自动评估」，避免它和截图抢时序。
                _suppressAutoSafety = true;
                // 让窗口先完成一次布局/绘制，再跑算法与渲染，避免在 Loaded 里做重活把首帧拖成黑屏
                Dispatcher.BeginInvoke(new Action(async () =>
                {
                    try
                    {
                        BtnGenerate_Click(this, new RoutedEventArgs());
                        App.Log("CLI", null, "auto-generate 完成");

                        // --safety：跑一次安全评估，结论会打到日志，便于无人值守验收
                        if (o.ContainsKey("safety") && _lastTraj.Count > 0)
                        {
                            await RunSafetyCheckAsync(true);
                            App.Log("CLI", null, "auto-safety 完成");

                            // --export-fixed：把倾斜修正后的轨迹落盘（等同界面上的导出按钮）
                            if (o.TryGetValue("export-fixed", out var fx) && _safetyPts.Length > 0)
                            {
                                WriteFixedCsv(fx);
                                App.Log("CLI", null, $"export-fixed -> {fx}（{_safetyPts.Length} 点）");
                            }
                        }
                    }
                    catch (Exception ex)
                    {
                        App.Log("FATAL/AutoGenerate", ex);
                    }

                    if (o.TryGetValue("view", out var vm)) ApplyViewPreset(vm);

                    bool wantShot = o.ContainsKey("shot") || o.ContainsKey("shot-win");
                    if (wantShot)
                    {
                        // 等一帧稳定（相机/布局/3D 绘制都完成）再截，避免抓到半成品
                        var timer = new System.Windows.Threading.DispatcherTimer
                        {
                            Interval = TimeSpan.FromMilliseconds(700)
                        };
                        timer.Tick += (_, _) =>
                        {
                            timer.Stop();
                            // shot = 只截 3D 视图（看几何）；shot-win = 整窗（看界面与结论卡）
                            if (o.TryGetValue("shot", out var shotPath))
                                try { SaveViewportShot(shotPath); } catch (Exception ex) { App.Log("SHOT/FAIL", ex); }
                            if (o.TryGetValue("shot-win", out var winPath))
                                try { SaveWindowShot(winPath); } catch (Exception ex) { App.Log("SHOT-WIN/FAIL", ex); }
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

        /// <summary>
        /// 浏览对话框的起始目录。优先级：
        ///   1) 当前路径框里文件所在目录（用户此刻的意图，最贴近他要找的地方）；
        ///   2) 上次成功打开过的目录（换了一批文件、路径框已清空时用）；
        ///   3) 文档目录（首次运行兜底）。
        /// 任一候选目录不存在就跳到下一档，最后再退回文档，避免对话框报错或乱跳。
        /// </summary>
        string PickStartDir()
        {
            string docs = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);

            string? inBox = null;
            try
            {
                string t = TxtPath.Text?.Trim();
                if (!string.IsNullOrEmpty(t))
                    inBox = System.IO.Path.GetDirectoryName(System.IO.Path.GetFullPath(t));
            }
            catch { }
            if (!string.IsNullOrEmpty(inBox) && System.IO.Directory.Exists(inBox)) return inBox;

            string? remembered = UiState.Get("lastCloudDir");
            if (!string.IsNullOrEmpty(remembered))
            {
                // 用户可能把文件夹改名/移到别处，逐级上退找一个还存在的
                var d = new System.IO.DirectoryInfo(remembered);
                while (d != null && !d.Exists) d = d.Parent;
                if (d != null) return d.FullName;
            }

            return System.IO.Directory.Exists(docs) ? docs : System.IO.Path.GetTempPath();
        }

        void BtnBrowse_Click(object sender, RoutedEventArgs e)
        {
            var dlg = new OpenFileDialog
            {
                Filter = "PLY 点云|*.ply|所有文件|*.*",
                // 必须显式给 InitialDirectory：.NET Core 下的 WPF 对话框不会自己记忆，
                // 不给就每次都从"文档"开始（这正是此前"每次点开都到这里"的原因）。
                InitialDirectory = PickStartDir(),
                // 不要设 true：true 会把进程当前目录恢复原值并绕开系统记忆，反而更死。
                RestoreDirectory = false,
            };
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
                var swLoad = System.Diagnostics.Stopwatch.StartNew();
                _cloud = PlyReader.Read(path);
                swLoad.Stop();
                _currentFile = path;
                // 记住这次文件的目录：下次点"浏览"直接从这里开始（含命令行 --ply 打开的情况）
                try
                {
                    var dir = System.IO.Path.GetDirectoryName(System.IO.Path.GetFullPath(path));
                    if (!string.IsNullOrEmpty(dir)) UiState.Set("lastCloudDir", dir);
                }
                catch { }
                RenderCloud();
                TxtHint.Text = $"点云已加载：{_cloud.Count} 点（读取 {swLoad.ElapsedMilliseconds} ms）";
                Title = $"喷涂轨迹生成器 - {System.IO.Path.GetFileName(path)}";
                App.Log("LOAD", null, $"{path} -> {_cloud.Count} 点，读取耗时 {swLoad.ElapsedMilliseconds} ms");
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

            // 大点云必须**等距抽样**，不能 Take(N)。
            // 现场点云是逐扫描线存盘的，Take(200000) 取到的是文件开头那一段连续区域，
            // 601k 点的工件只画前 200k，看起来就像"点云只盖住一条斜带"，
            // 而且这团点是不透明地压在轨迹上面 → 用户以为"轨迹只喷了左半边"。
            // 抽样密度上限按"视觉够用"定：20 万点已能把轮廓看清楚。
            const int MaxDraw = 200000;
            int step = Math.Max(1, (_cloud.Count + MaxDraw - 1) / MaxDraw);
            var pts = new Point3DCollection(_cloud.Count / step + 8);
            for (int i = 0; i < _cloud.Count; i += step)
            {
                var p = _cloud[i];
                pts.Add(new Point3D(p.x, p.y, p.z));
            }

            _cloudVisual = new PointsVisual3D
            {
                Points = pts,
                Color = Color.FromRgb(224, 64, 224),
                // 点越密，点径越小，否则整片糊成色块把轨迹/面片全遮住
                Size = pts.Count > 120000 ? 1.2 : (pts.Count > 40000 ? 1.8 : 2.5)
            };
            Viewport.Children.Add(_cloudVisual);
            Viewport.ZoomExtents(200);
            ApplyCloudVisibility();
            App.Log("CLOUD", null,
                $"绘制 {pts.Count}/{_cloud.Count} 点（等距抽样 step={step}，点径={_cloudVisual.Size}）");
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
                flatnessTol    = ParseDouble(TxtFlatness.Text, 40),
                curvatureTol   = ParseDouble(TxtCurvature.Text, 5),
                sorK           = (int)ParseDouble(TxtSorK.Text, 20),
                sorZ           = ParseDouble(TxtSorZ.Text, 2.5),
                normalFlip     = CmbNormalFlip.SelectedIndex == 1 ? 1 : 0,
                fillHoleSize   = ParseDouble(TxtFillHole.Text, 300),
                scanDirMode    = CmbScanDir.SelectedIndex,
                gunClear       = ParseDouble(TxtGunClear.Text, 20),
                edgeExtend     = ParseDouble(TxtEdgeExtend.Text, -1),
                edgeExtendAcross = ParseDouble(TxtEdgeExtendAcross.Text, -1),
                shellThick     = ParseDouble(TxtShellThick.Text, 0),
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
            var swGen = System.Diagnostics.Stopwatch.StartNew();
            try
            {
                res = Native.spray_traj_generate(xyz, _cloud.Count, ref p);
                swGen.Stop();
                App.Log("GEN", null, $"算法耗时 {swGen.ElapsedMilliseconds} ms（{_cloud.Count} 点）");
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
            // 留一份给安全评估用：res.points 马上就要被 spray_traj_free 释放
            _lastTraj = traj;

            // 保存所有分割平面参数，供"显示重建平面"使用
            // 上限取数组槽位数，不能写死 8：薄壁件拆层后面数会超过 8（10 面很常见），
            // 写死会让后几个面既不渲染也不进清单，看上去像"内壁丢了"。
            _planes.Clear();
            for (int i = 0; i < res.plane_count && i < res.planes.Length; i++)
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
            FillFaceTable();

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
                $"  红色=喷涂行，灰色=空移（面间抬枪绕行）\n" +
                $"  蓝色面=外壁，橙色面=内壁";

            TxtHint.Text = $"轨迹已生成：{res.count} 点，{res.rows} 行";

            // 生成完自动收起原始点云：几十万点的洋红点云是不透明的，会整片压在
            // 轨迹和重建面上面，看上去像"轨迹只喷了一半"。点云随时可以再勾回来对比。
            if (ChkShowCloud.IsChecked == true && !_keepCloudOnGenerate)
            {
                ChkShowCloud.IsChecked = false;
                App.Log("CLOUD", null, "生成后自动隐藏原始点云（可在显示控制里重新勾选）");
            }

            Native.spray_traj_free(ref res);

            // 同时写出 CSV 供外部使用
            try { ExportCsv(traj); } catch { }

            // 生成完顺手做一次安全评估：轨迹能生成 ≠ 机器人做得到，
            // 少了这一步，用户看到的漂亮轨迹可能是不可达或撞工件的。
            if (ChkAutoSafety.IsChecked == true && !_suppressAutoSafety)
                _ = RunSafetyCheckAsync(true);
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
                    ThetaDiv = spraying ? 12 : 8
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
                    $"法向=({nx:F3},{ny:F3},{nz:F3}) 真实=({p.n0x:F3},{p.n0y:F3},{p.n0z:F3}) " +
                    $"分类={p.side} 内壁={(p.isInner > 0 ? "是" : (p.isInner < 0 ? "未分" : "否"))} 依据={p.sideSrc}");
            }
        }

        // ---- 面级法向清单（内壁/外壁）----
        // 数据全部来自 C++ 侧：nx,ny,nz 是喷距法向（轨迹沿它抬升），
        // n0 是射线自由空间判出的真实可喷侧法向。两者点积为负 = 该面喷距方向反了。
        private sealed class FaceRow
        {
            public int Face { get; set; }
            public string Normal { get; set; } = "";
            public string AngZ { get; set; } = "";
            public string Tilt { get; set; } = "";
            public string Kind { get; set; } = "";
            public string Inner { get; set; } = "";
            public string Src { get; set; } = "";
            public string Free { get; set; } = "";
            public string Layer { get; set; } = "";
            public string Note { get; set; } = "";
        }

        void FillFaceTable()
        {
            var rows = new List<FaceRow>();
            int skipped = 0, nInner = 0, nOuter = 0, nSingle = 0;
            for (int i = 0; i < _planes.Count; i++)
            {
                var p = _planes[i];
                // 口径与 C++ 一致：与 +Z 夹角 ≤ NormalWinDeg 的面全部列出（内壁、外壁一起），
                // 完全朝下的面（开口盒底板外侧，180°）不列入。
                if (p.angZ > NormalWinDeg) { skipped++; continue; }
                // 单层壳（数据里这一面只有一层点）：单独计数，只用来在统计栏里提示
                // "这一面没有第二层点"。但**不影响内/外壁计数** —— 射线已经给出了
                // 可喷侧，是内壁就是内壁；把单层壳排除在统计外，反而会出现
                // "内壁 0 / 外壁 4"这种更离谱的结果。
                bool dual = p.sideSrc == 5;
                bool single = p.layerCnt == 1 && !dual;
                if (single) nSingle++;
                if (p.isInner > 0) nInner++;
                else if (p.isInner == 0) nOuter++;
                bool hasDiag = !(p.n0x == 0 && p.n0y == 0 && p.n0z == 0);
                double flip = p.nx * p.n0x + p.ny * p.n0y + p.nz * p.n0z;
                string kind = p.side > 0 ? "正向" : (p.side < 0 ? "反向" : "竖直");
                string src = p.sideSrc == 1 ? "射线" : (p.sideSrc == 3 ? "强制"
                           : (p.sideSrc == 4 ? "薄壁" : (p.sideSrc == 5 ? "双面" : "兜底")));
                string fp = p.freePlus < 0 ? "∞" : p.freePlus.ToString("F0");
                string fm = p.freeMinus < 0 ? "∞" : p.freeMinus.ToString("F0");

                string note;
                if (!hasDiag) note = "无诊断数据（旧 DLL）";
                else if (dual) note = $"薄板双面：按壁厚 {p.layerThick:F1}mm 补出的另一面";
                else if (single) note = "单层壳：这一面在点云里只有一层点，没有壁厚";
                else if (p.sideSrc == 4) note = p.isInner > 0
                    ? "薄壁拆层：内壁（法向朝腔体）"
                    : "薄壁拆层：外壁（法向朝外）";
                else if (p.sideSrc == 2) note = "两侧都自由（单层壳），待人工确认";
                else if (flip < 0) note = "内外反：喷距方向与真实可喷侧相反";
                else note = "喷距方向一致";

                rows.Add(new FaceRow
                {
                    Face = i,
                    Normal = $"({p.n0x:F3}, {p.n0y:F3}, {p.n0z:F3})",
                    AngZ = $"{p.angZ:F1}",
                    Tilt = $"{p.tiltDeg:F1}",
                    Kind = kind,
                    Inner = p.isInner > 0 ? "是" : (p.isInner < 0 ? "未分" : "否"),
                    Src = src,
                    Free = $"{fp} / {fm}",
                    // 层/跨度：单层 = 点云在这面只有一层点（跨度≈噪声/面起伏，没有壁厚）；
                    // 2层 = 沿法向有两片点，跨度就是壁厚。未判 = 曲面块或旧 DLL。
                    Layer = dual ? $"双面 {p.layerThick:F1}"
                          : (p.layerCnt == 2 ? $"2层 {p.layerThick:F1}"
                          : (p.layerCnt == 1 ? $"单层 {p.layerThick:F1}" : "未判")),
                    Note = note
                });
            }
            GridFaces.ItemsSource = rows;
            string stat;
            if (rows.Count > 0 && nSingle == rows.Count)
            {
                // 一面都没补：这份点云的每个面只有一层点，没有壁厚，也就没有"外壁"这种面。
                // 直接说清楚，比"内壁 5 / 外壁 0"更容易看懂。
                stat = $"面级法向清单：窗口内 {rows.Count} 面 —— 全部是单层壳" +
                       $"（每个面只有一层点、没有壁厚，点云里不含第二层；" +
                       $"要按板的两面出轨迹，把「薄板双面壁厚」填上板厚）；" +
                       $"分割共 {_planes.Count} 面";
            }
            else
            {
                stat = $"面级法向清单：窗口内 {rows.Count} 面（内壁 {nInner} / 外壁 {nOuter}）" +
                       $"，分割共 {_planes.Count} 面";
                if (nSingle > 0) stat += $"；其中 {nSingle} 面数据只有一层点（该面的另一侧未扫描到）";
            }
            if (skipped > 0) stat += $"，滤除 {skipped} 面（与 +Z 夹角 >{NormalWinDeg:F0}°，朝下）";
            if (TxtFaceStat != null) TxtFaceStat.Text = stat;
            App.Log("FACES", null, stat);
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

            // 内壁/外壁用不同颜色：薄壁件拆层后内外两面几乎重合，
            // 同色半透明叠在一起看不出是两块，容易被当成"只画了外壁"。
            Color fill = p.isInner > 0 ? Color.FromArgb(88, 255, 148, 52)
                       : p.isInner < 0 ? Color.FromArgb(88, 170, 170, 170)
                                       : Color.FromArgb(78, 64, 156, 255);
            var mat = new DiffuseMaterial(new SolidColorBrush(fill));
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
                    Color = p.isInner > 0 ? Color.FromRgb(176, 88, 16) : Color.FromRgb(32, 82, 150),
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

        // ══════════════════════════════════════════════════════════════════
        //  机器人执行安全
        //
        //  职责分离：spray_traj 只管「喷得对不对」，这一段只管「机器人能不能
        //  安全做到」——可达、关节限位、枪体/连杆与工件碰撞、连杆自碰撞、喷角。
        //
        //  生成结果一个字不改：先按原样画轨迹，再把安全层的修正轨迹用青色叠
        //  上去。哪个点被倾斜了、倾斜后间隙变大还是变小，界面上直接看得见。
        // ══════════════════════════════════════════════════════════════════

        private List<SprayTrajPoint> _lastTraj = new();
        private RobotSafetyPoint[] _safetyPts = Array.Empty<RobotSafetyPoint>();
        private readonly List<Visual3D> _riskVisuals = new();
        private LinesVisual3D? _fixedVisual;
        private bool _suppressAutoSafety;

        // 严重度分组，与 C++ 侧 RS_VIOL_* 的 P0/P1/P2 划分一致
        const uint MaskP0 = RobotSafetyRule.IkFail | RobotSafetyRule.LinkHit | RobotSafetyRule.ToolHit
                          | RobotSafetyRule.SelfCollision;
        const uint MaskP1 = RobotSafetyRule.Clearance | RobotSafetyRule.JointLimit | RobotSafetyRule.SprayAngle;
        const uint MaskP2 = RobotSafetyRule.NearSingular | RobotSafetyRule.TiltJump;

        async void BtnSafety_Click(object sender, RoutedEventArgs e)
        {
            try { await RunSafetyCheckAsync(false); }
            catch (Exception ex)
            {
                App.Log("SAFETY/FAIL", ex);
                MessageBox.Show($"安全评估失败：{ex.Message}", "错误", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        void CmbBaseMode_Changed(object sender, System.Windows.Controls.SelectionChangedEventArgs e)
        {
            if (PanelBaseManual == null) return;
            PanelBaseManual.Visibility = CmbBaseMode.SelectedIndex == 1
                ? Visibility.Visible : Visibility.Collapsed;
        }

        /// <summary>
        /// 手输基座位姿 → T_base_in_cloud（4x4 行主序）。
        /// 约定：基座 +Z 竖直向上；+X 默认指向工件质心，再用「朝向偏转」绕 Z 微调。
        /// 自动模式返回 null，交给核心库按工件质心假设（报告里会明确标注是猜的）。
        /// </summary>
        double[]? BuildBasePose()
        {
            if (CmbBaseMode.SelectedIndex != 1) return null;

            double x = ParseDouble(TxtBaseX.Text, 0);
            double y = ParseDouble(TxtBaseY.Text, 0);
            double z = ParseDouble(TxtBaseZ.Text, 0);
            double yaw = ParseDouble(TxtBaseYaw.Text, 0) * Math.PI / 180.0;

            // 质心（水平投影）；大点云等距抽样，够用且不卡界面
            double cx = 0, cy = 0;
            int step = Math.Max(1, _cloud.Count / 4000);
            int cnt = 0;
            for (int i = 0; i < _cloud.Count; i += step) { cx += _cloud[i].x; cy += _cloud[i].y; cnt++; }
            if (cnt > 0) { cx /= cnt; cy /= cnt; }

            double dx = cx - x, dy = cy - y;
            double len = Math.Sqrt(dx * dx + dy * dy);
            if (len < 1e-6) { dx = 1; dy = 0; } else { dx /= len; dy /= len; }

            double c = Math.Cos(yaw), s = Math.Sin(yaw);
            double ax = c * dx - s * dy;      // 基座 +X（绕 Z 转 yaw 后）
            double ay = s * dx + c * dy;
            // +Y = +Z × +X = (-ay, ax, 0)
            return new double[]
            {
                ax, -ay, 0, x,
                ay,  ax, 0, y,
                 0,   0, 1, z,
                 0,   0, 0, 1,
            };
        }

        async Task RunSafetyCheckAsync(bool auto)
        {
            if (_lastTraj.Count == 0)
            {
                if (!auto)
                    MessageBox.Show("请先生成轨迹，再做安全评估。", "提示",
                        MessageBoxButton.OK, MessageBoxImage.Information);
                return;
            }
            if (_cloud.Count == 0)
            {
                if (!auto) MessageBox.Show("请先加载点云。", "提示", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            var p = Native.robot_safety_default_params();
            // 喷距必须与生成时一致 —— 倾斜是绕面点的球面旋转，喷距错了球心就错，
            // 后面「可达 / 碰撞」的结论全部不成立，所以这里强制取界面上的生成值。
            p.standoff = ParseDouble(TxtStandoff.Text, p.standoff);
            p.safety = ParseDouble(TxtSafetyGap.Text, p.safety);
            // 绝对安全距离/回退硬下限：只有低于它才回退，间隙略小但高于它时保留平滑。
            p.hardClearance = Math.Max(0, ParseDouble(TxtHardClearance.Text, p.hardClearance));
            p.toolLen = ParseDouble(TxtToolLen.Text, p.toolLen);
            p.tiltMaxDeg = ParseDouble(TxtTiltMax.Text, p.tiltMaxDeg);
            p.maxSprayAngleDeg = ParseDouble(TxtMaxSprayAngle.Text, p.maxSprayAngleDeg);
            p.jointMarginMin = ParseDouble(TxtJointMargin.Text, p.jointMarginMin);
            // 并行线程数：默认取"逻辑核数 - 1"，留一个核给界面线程，否则评估时窗口卡顿。
            // 单线程（=1）结果与命令行逐字节一致，用于要和存档/对拍对齐的场合。
            int defThreads = Math.Max(1, Environment.ProcessorCount - 1);
            p.threads = Math.Clamp((int)ParseDouble(TxtThreads.Text, defThreads), 1, 64);
            // 点云体素：0=全量精判。>0 只用于在线预判，薄壁件会偏乐观（见 README 的量化数据）。
            p.voxelMm = Math.Max(0, ParseDouble(TxtVoxel.Text, 0));
            // 进刀侧：0=按机器人站位逐行择优（默认），1=只用「向上」（旧口径，用于和历史存档对拍）
            p.azMode = CmbAzMode.SelectedIndex == 1 ? 1 : 0;

            double[]? basePose = BuildBasePose();
            var pts = _lastTraj.ToArray();
            var cloudXyz = new double[_cloud.Count * 3];
            for (int i = 0; i < _cloud.Count; i++)
            {
                cloudXyz[i * 3] = _cloud[i].x;
                cloudXyz[i * 3 + 1] = _cloud[i].y;
                cloudXyz[i * 3 + 2] = _cloud[i].z;
            }

            BtnSafety.IsEnabled = false;
            TxtHint.Text = $"机器人安全评估中…（轨迹 {pts.Length} 点 / 点云 {_cloud.Count} 点）";
            var sw = System.Diagnostics.Stopwatch.StartNew();

            RobotSafetyPoint[] got = Array.Empty<RobotSafetyPoint>();
            RobotSafetyResult st = default;
            RobotSafetyResult sum = default;   // free 会把 st.count 清 0，标量先留一份副本
            string? err = null;
            try
            {
                st = await Task.Run(() =>
                    Native.robot_safety_check(pts, pts.Length, cloudXyz, _cloud.Count, basePose, ref p));
                sw.Stop();

                if (st.error != IntPtr.Zero)
                {
                    err = Marshal.PtrToStringAnsi(st.error) ?? "未知错误";
                }
                else
                {
                    sum = st;
                    got = new RobotSafetyPoint[st.count];
                    int sz = Marshal.SizeOf<RobotSafetyPoint>();
                    for (int i = 0; i < st.count; i++)
                        got[i] = Marshal.PtrToStructure<RobotSafetyPoint>(IntPtr.Add(st.points, i * sz));
                }
            }
            catch (EntryPointNotFoundException)
            {
                err = "当前 SprayTraj.dll 不含安全层导出（robot_safety_check）。请重新编译 C++ 核心库并复制到程序目录。";
            }
            finally
            {
                try { Native.robot_safety_free(ref st); } catch { }
                BtnSafety.IsEnabled = true;
            }

            if (err != null)
            {
                TxtHint.Text = "安全评估失败：" + err;
                ShowSafetyCard("ERROR", Color.FromRgb(248, 113, 113), "", err);
                App.Log("SAFETY/FAIL", null, err);
                return;
            }

            // 按最高严重度归类（一个点可能同时命中多条规则）
            int nP0 = 0, nP1 = 0, nP2 = 0;
            foreach (var q in got)
            {
                if ((q.violation & MaskP0) != 0) nP0++;
                else if ((q.violation & MaskP1) != 0) nP1++;
                else if ((q.violation & MaskP2) != 0) nP2++;
            }

            string verdict = nP0 > 0 ? "REJECT" : (nP1 > 0 ? "PASS WITH WARNINGS" : "PASS");
            Color vColor = nP0 > 0 ? Color.FromRgb(248, 113, 113)
                         : nP1 > 0 ? Color.FromRgb(251, 191, 36)
                                   : Color.FromRgb(74, 222, 128);

            int[] b = sum.byRule ?? new int[10];
            string Rule(int bit) => b.Length > bit ? b[bit].ToString() : "0";

            // 最危险的连杆段与枪体段：现场最想知道的不是"最小距离是多少"，
            // 而是"哪一段撞了"——枪体撞改喷距/角度，连杆撞要改站位或整条姿态。
            int linkIdx = -1, toolIdx = -1;
            double linkMin = double.PositiveInfinity, toolMin = double.PositiveInfinity;
            for (int i = 0; i < got.Length; i++)
            {
                if (got[i].clearanceLinkMm < linkMin) { linkMin = got[i].clearanceLinkMm; linkIdx = i; }
                if (got[i].clearanceToolMm < toolMin) { toolMin = got[i].clearanceToolMm; toolIdx = i; }
            }
            string LinkDesc() => linkIdx >= 0 && double.IsFinite(linkMin)
                ? $"{SegKind.Name(got[linkIdx].linkKind, got[linkIdx].linkJ0, got[linkIdx].linkJ1)} {linkMin:F2}mm"
                : "—";
            string ToolDesc() => toolIdx >= 0 && double.IsFinite(toolMin)
                ? $"{SegKind.Name(got[toolIdx].worstKind, got[toolIdx].worstJ0, got[toolIdx].worstJ1)} {toolMin:F2}mm"
                : "—";

            var sb = new System.Text.StringBuilder();
            sb.AppendLine($"可达 {sum.reachable}/{sum.count}" +
                          $"（{(sum.count > 0 ? 100.0 * sum.reachable / sum.count : 0):F1}%）");
            sb.AppendLine($"最小表面距离 {sum.minClearanceMm:F2} mm" +
                          (sum.minClearanceMm < 0 ? "  ← 穿透" : ""));
            sb.AppendLine($"　枪体 {ToolDesc()}　连杆 {LinkDesc()}");
            sb.AppendLine($"喷角 均值 {sum.tiltMeanDeg:F1}° · p95 {sum.tiltP95Deg:F1}° · 最大 {sum.tiltMaxDeg:F1}°");
            // 进刀侧：说明"姿势是谁定的"——重力（向上）还是机器人站位。
            // 只报**实际生效**的行数：改判了但被逐点回退掉的不算（否则数字会虚高）。
            sb.AppendLine(p.azMode == 1
                ? $"进刀侧：全部向上（旧口径，{sum.azRows} 行）"
                : $"进刀侧：{sum.azRows} 行中 {sum.azStanceRows} 行按站位侧（其余沿法向上）");
            sb.AppendLine($"喷距守恒偏差 {sum.standoffMaxErrMm:E2} mm");
            sb.AppendLine($"回退策略：绝对安全距离 {p.hardClearance:F1} mm（硬约束通过即保留平滑）");
            if (sum.poseReverted > 0)
                sb.AppendLine($"平滑回退 {sum.poseReverted} 点（仅低于硬下限/碰撞/自碰撞/关节超限）");
            else
                sb.AppendLine("平滑回退 0 点（全部硬约束通过）");
            sb.AppendLine($"违规：不可达 {Rule(0)} · 连杆撞 {Rule(8)} · 枪体撞 {Rule(9)} · 自碰撞 {Rule(3)}");
            sb.AppendLine($"　　　间隙 {Rule(2)} · 关节 {Rule(4)} · 喷角 {Rule(5)} · 奇异 {Rule(6)} · 跳变 {Rule(7)}");
            sb.AppendLine(sum.baseAuto != 0
                ? $"基座：自动假设 ({sum.basePose[3]:F0}, {sum.basePose[7]:F0}, {sum.basePose[11]:F0}) —— 未标定，可达率仅作预判"
                : $"基座：手动 ({sum.basePose[3]:F0}, {sum.basePose[7]:F0}, {sum.basePose[11]:F0})");

            ShowSafetyCard(verdict, vColor, $"耗时 {sw.ElapsedMilliseconds} ms", sb.ToString().TrimEnd());

            _safetyPts = got;
            RenderSafetyOverlay(got, sum.tiltedPoints > 0);

            TxtHint.Text = $"安全评估：{verdict} · 可达 {sum.reachable}/{sum.count} · 最小间隙 {sum.minClearanceMm:F1}mm";
            App.Log("SAFETY", null,
                $"verdict={verdict} p0={nP0} p1={nP1} p2={nP2} count={sum.count} reachable={sum.reachable} " +
                $"minClr={sum.minClearanceMm:F2} hardClr={p.hardClearance:F2} " +
                $"tool={sum.minClearanceToolMm:F2} link={sum.minClearanceLinkMm:F2} " +
                $"linkSeg={LinkDesc()} " +
                $"tilt(mean/p95/max)={sum.tiltMeanDeg:F1}/{sum.tiltP95Deg:F1}/{sum.tiltMaxDeg:F1} " +
                $"tilted={sum.tiltedPoints} reverted={sum.poseReverted} standoffErr={sum.standoffMaxErrMm:E2} " +
                $"baseAuto={sum.baseAuto} cloudPts={sum.cloudPoints} l0={sum.l0Issues} " +
                $"threads={p.threads} voxel={p.voxelMm:F1} elapsed={sw.ElapsedMilliseconds}ms");
        }

        void ShowSafetyCard(string verdict, Color color, string cost, string body)
        {
            CardSafety.Visibility = Visibility.Visible;
            TxtSafetyVerdict.Text = verdict;
            TxtSafetyVerdict.Foreground = new SolidColorBrush(color);
            TxtSafetyCost.Text = cost;
            TxtSafetyBody.Text = body;
            // 「机器人安全」页签里放一份同样的文字：3D 视图被相机拖动时卡片仍可见，
            // 但用户也可能就待在参数页，不该为了看结论切回去。
            TxtSafetyInline.Text = string.IsNullOrEmpty(cost)
                ? $"{verdict}\n{body}"
                : $"{verdict}（{cost}）\n{body}";
        }

        /// <summary>
        /// 把风险点分三层画，否则会糊成一片看不出重点：
        ///   红（大）   P0 碰撞 / 不可达 / 自碰撞 —— 这里必须停线
        ///   橙（中）   间隙不足 / 关节余量不足 —— 局部风险，位置有意义
        ///   黄（小）   喷角 / 奇异 / 倾角跳变 —— 趋势性指标，成百上千个点，
        ///              抽稀到 ≤400 个只为"告诉你在哪一段"，具体数量看结论卡
        /// </summary>
        void RenderSafetyOverlay(RobotSafetyPoint[] arr, bool showFixed)
        {
            foreach (var v in _riskVisuals) Viewport.Children.Remove(v);
            _riskVisuals.Clear();
            if (_fixedVisual != null)
            {
                Viewport.Children.Remove(_fixedVisual);
                _fixedVisual = null;
            }
            if (arr.Length == 0) return;

            const uint MaskLocal = RobotSafetyRule.Clearance | RobotSafetyRule.JointLimit;

            var red = new Point3DCollection();
            var orange = new Point3DCollection();
            var yellowAll = new List<Point3D>();

            foreach (var q in arr)
            {
                if (q.violation == 0) continue;
                var pt = new Point3D(q.x, q.y, q.z);
                if ((q.violation & MaskP0) != 0) red.Add(pt);
                else if ((q.violation & MaskLocal) != 0) orange.Add(pt);
                else yellowAll.Add(pt);
            }

            void AddRisk(System.Windows.Media.Media3D.Point3DCollection c, Color col, double size)
            {
                if (c.Count == 0) return;
                var v = new PointsVisual3D { Points = c, Color = col, Size = size };
                _riskVisuals.Add(v);
                Viewport.Children.Add(v);
            }
            AddRisk(red, Color.FromRgb(255, 48, 32), 13);
            AddRisk(orange, Color.FromRgb(251, 146, 60), 8);

            // 趋势性风险点多到几千个，全画出来只会盖住工件，抽稀保留分布位置
            int nYellowDraw = 0;
            if (yellowAll.Count > 0)
            {
                int step = Math.Max(1, yellowAll.Count / 400);
                var y = new Point3DCollection(yellowAll.Count / step + 2);
                for (int i = 0; i < yellowAll.Count; i += step) { y.Add(yellowAll[i]); nYellowDraw++; }
                AddRisk(y, Color.FromRgb(250, 204, 21), 4);
            }

            // 修正后轨迹：是否连线看**原轨迹**相邻点距（>15mm 即跨行/空移，不连），
            // 不能用修正后的距离——倾斜会把点整体挪开，按修正后距离判断会把每行都判成断开。
            if (showFixed && _lastTraj.Count == arr.Length)
            {
                var pts = new Point3DCollection(arr.Length * 2);
                for (int i = 1; i < arr.Length; i++)
                {
                    var o0 = _lastTraj[i - 1];
                    var o1 = _lastTraj[i];
                    double dx = o1.x - o0.x, dy = o1.y - o0.y, dz = o1.z - o0.z;
                    if (Math.Sqrt(dx * dx + dy * dy + dz * dz) > 15.0) continue;
                    pts.Add(new Point3D(arr[i - 1].x, arr[i - 1].y, arr[i - 1].z));
                    pts.Add(new Point3D(arr[i].x, arr[i].y, arr[i].z));
                }
                if (pts.Count > 0)
                {
                    _fixedVisual = new LinesVisual3D
                    {
                        Points = pts,
                        Color = Color.FromArgb(205, 34, 211, 238),
                        Thickness = 1.1
                    };
                    Viewport.Children.Add(_fixedVisual);
                }
            }

            ApplyRiskVisibility();
            ApplyFixedVisibility();
            App.Log("SAFETY/DRAW", null,
                $"风险点 红={red.Count} 橙={orange.Count} 黄={nYellowDraw}/{yellowAll.Count}（抽稀）；" +
                $"修正轨迹线段={(_fixedVisual?.Points.Count ?? 0) / 2}");
        }

        void ChkShowRisk_Changed(object sender, RoutedEventArgs e) => ApplyRiskVisibility();
        void ChkShowFixed_Changed(object sender, RoutedEventArgs e) => ApplyFixedVisibility();

        void ApplyRiskVisibility()
        {
            if (Viewport == null) return;
            bool show = ChkShowRisk.IsChecked.GetValueOrDefault(true);
            foreach (var v in _riskVisuals)
            {
                bool inTree = Viewport.Children.Contains(v);
                if (show && !inTree) Viewport.Children.Add(v);
                else if (!show && inTree) Viewport.Children.Remove(v);
            }
        }

        void ApplyFixedVisibility()
        {
            if (Viewport == null || _fixedVisual == null) return;
            bool show = ChkShowFixed.IsChecked.GetValueOrDefault(true);
            bool inTree = Viewport.Children.Contains(_fixedVisual);
            if (show && !inTree) Viewport.Children.Add(_fixedVisual);
            else if (!show && inTree) Viewport.Children.Remove(_fixedVisual);
        }

        /// <summary>
        /// 导出修正后的轨迹。列序与命令行版 spray_safety.exe 完全一致，
        /// 下游（机器人程序 / 离线校验脚本）不区分数据来源。
        /// </summary>
        void BtnExportFixed_Click(object sender, RoutedEventArgs e)
        {
            if (_safetyPts.Length == 0)
            {
                MessageBox.Show("请先做一次安全评估。", "提示", MessageBoxButton.OK, MessageBoxImage.Information);
                return;
            }

            var dlg = new SaveFileDialog
            {
                Filter = "CSV 轨迹|*.csv|所有文件|*.*",
                FileName = "trajectory_safe.csv",
            };
            if (dlg.ShowDialog() != true) return;

            try
            {
                WriteFixedCsv(dlg.FileName);
                TxtHint.Text = $"修正后轨迹已导出：{dlg.FileName}（{_safetyPts.Length} 点）";
                App.Log("SAFETY/EXPORT", null, $"{dlg.FileName} {_safetyPts.Length} 点");
            }
            catch (Exception ex)
            {
                MessageBox.Show($"导出失败：{ex.Message}", "错误", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        /// <summary>写修正后轨迹 CSV（无对话框版本，供命令行验收复用）。</summary>
        void WriteFixedCsv(string path)
        {
            var sb = new System.Text.StringBuilder();
            sb.Append("X,Y,Z,I,J,K,speed,spray,tilt_deg,clearance_mm,self_gap_mm,margin_deg,")
              .Append("standoff_err_mm,lifted_mm,cond,reverted,J1_deg,J2_deg,J3_deg,J4_deg,J5_deg,J6_deg,violation\n");
            foreach (var q in _safetyPts)
            {
                sb.Append($"{q.x:F3},{q.y:F3},{q.z:F3},{q.i:F4},{q.j:F4},{q.k:F4},{q.speed:F1},{q.spray},");
                sb.Append($"{q.tiltDeg:F2},{q.clearanceMm:F2},{q.selfGapMm:F2},{q.marginDeg:F2},");
                sb.Append($"{q.standoffErrMm:E6},{q.liftedMm:F1},{q.cond:F1},{q.reverted}");
                var th = q.theta ?? new double[6];
                for (int k = 0; k < 6; k++)
                    sb.Append($",{(k < th.Length ? th[k] : 0.0):F3}");
                sb.Append($",{RuleText(q.violation)}\n");
            }
            System.IO.File.WriteAllText(path, sb.ToString(), new System.Text.UTF8Encoding(true));
        }

        static string RuleText(uint v)
        {
            if (v == 0) return "";
            var parts = new List<string>(4);
            if ((v & RobotSafetyRule.IkFail) != 0) parts.Add("IK_FAIL");
            if ((v & RobotSafetyRule.LinkHit) != 0) parts.Add("LINK_HIT");
            if ((v & RobotSafetyRule.ToolHit) != 0) parts.Add("TOOL_HIT");
            if ((v & RobotSafetyRule.Clearance) != 0) parts.Add("CLEARANCE");
            if ((v & RobotSafetyRule.SelfCollision) != 0) parts.Add("SELF_COLLISION");
            if ((v & RobotSafetyRule.JointLimit) != 0) parts.Add("JOINT_LIMIT");
            if ((v & RobotSafetyRule.SprayAngle) != 0) parts.Add("SPRAY_ANGLE");
            if ((v & RobotSafetyRule.NearSingular) != 0) parts.Add("NEAR_SINGULAR");
            if ((v & RobotSafetyRule.TiltJump) != 0) parts.Add("TILT_JUMP");
            return string.Join("|", parts);
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
        /// face（第十轮新增）：沿**面0自身法向**正视。曲面件的工件坐标系里没有"正对"方向，
        /// top/front/side 一律是斜视，验收"行是否铺满面、行距是否均匀"必须用 face。
        /// </summary>
        void ApplyViewPreset(string mode)
        {
            if (string.IsNullOrWhiteSpace(mode)) return;

            // face：沿面0法向正视（正交），取景只按面0的顶点范围，避免被点云/安全高度拉远
            if (string.Equals(mode, "face", StringComparison.OrdinalIgnoreCase))
            {
                if (_planes.Count == 0) { App.Log("VIEW/WARN", null, "face: 没有面片，跳过"); return; }
                var p0 = _planes[0];
                var nf = new Vector3D(p0.nx, p0.ny, p0.nz);
                if (nf.Length < 1e-9) { App.Log("VIEW/WARN", null, "face: 面0法向退化，跳过"); return; }
                nf.Normalize();

                Vector3D u3 = new Vector3D(p0.ux, p0.uy, p0.uz);
                Vector3D v3 = new Vector3D(p0.vx, p0.vy, p0.vz);
                if (u3.Length < 1e-6) u3 = Vector3D.CrossProduct(new Vector3D(0, 0, 1), nf);
                if (u3.Length < 1e-6) u3 = new Vector3D(1, 0, 0);
                u3.Normalize();
                if (v3.Length < 1e-6) v3 = Vector3D.CrossProduct(nf, u3);
                v3.Normalize();

                int nv0 = Math.Min(p0.polyVertCount, p0.polyVerts?.Length / 3 ?? 0);
                if (nv0 < 3) { App.Log("VIEW/WARN", null, "face: 面0顶点不足，跳过"); return; }
                var c0 = new Point3D(p0.cx, p0.cy, p0.cz);
                double minU = double.MaxValue, maxU = double.MinValue, minV = double.MaxValue, maxV = double.MinValue;
                for (int k = 0; k < nv0; k++)
                {
                    var q = new Vector3D(p0.polyVerts![k * 3] - c0.X, p0.polyVerts[k * 3 + 1] - c0.Y, p0.polyVerts[k * 3 + 2] - c0.Z);
                    double uu = Vector3D.DotProduct(q, u3), vv = Vector3D.DotProduct(q, v3);
                    if (uu < minU) minU = uu; if (uu > maxU) maxU = uu;
                    if (vv < minV) minV = vv; if (vv > maxV) maxV = vv;
                }
                double span0 = Math.Max(maxU - minU, maxV - minV);
                if (span0 <= 1e-6) { App.Log("VIEW/WARN", null, "face: 面0跨度为 0，跳过"); return; }

                // 相机站在面的 +n 一侧（枪位就在这一侧：轨迹 = 面点 + standoff·n）
                double back0 = Math.Max(span0 * 3.0, 1000.0);
                Viewport.Camera = new OrthographicCamera
                {
                    Position = c0 + nf * back0,
                    LookDirection = -nf * back0,
                    UpDirection = v3,
                    Width = span0 * 1.15,
                    NearPlaneDistance = 1.0,
                    FarPlaneDistance = back0 * 4.0,
                };
                App.Log("VIEW", null,
                    $"正交 face n=({nf.X:F3},{nf.Y:F3},{nf.Z:F3}) center=({c0.X:F1},{c0.Y:F1},{c0.Z:F1}) span={span0:F1} width={span0 * 1.15:F1}");
                return;
            }

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
