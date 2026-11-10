using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Input;
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
        // 喷涂段（spray=1，枪开）：红色
        private readonly List<TubeVisual3D> _trajVisuals = new();
        // 空移段（spray=0，枪关）：亮蓝色。单独留引用是必须的 ——
        // 空移段是"过渡动作"，和有效喷涂混在一起看，根本分不清机器人到底在哪儿喷，
        // 也不方便只看空移路径判断它绕得对不对。
        private readonly List<TubeVisual3D> _moveVisuals = new();
        // 配色统一放这里：喷涂红 / 空移亮蓝。
        // 空移不能用灰色 —— 点云本身是灰的，灰管子在深色背景上跟点云糊成一片，
        // 等于没标。亮蓝和点云灰、修正轨迹青（#22D3EE）、风险点红橙黄都不撞。
        static readonly Color SprayColor = Color.FromRgb(255, 48, 32);
        static readonly Color MoveColor = Color.FromRgb(64, 176, 240);
        // 命令行 --keep-cloud：生成后不自动收起点云，用于"点云抽样是否正确"的验收截图
        private bool _keepCloudOnGenerate;
        // 自动导出目录的解析结果：点云目录不可写（安全软件拦桌面、只读盘）时回退到程序目录，
        // 解析一次就缓存，避免每次生成都写探测文件；换点云时清空。
        private string? _exportDirResolved;
        // 命令行 --tiltk：倾角系数（δ=90° 侧壁时给多少度）。NaN = 用库默认（60）。
        private double _tiltCoefDeg = double.NaN;
        // 上一次评估**实际用上**的倾角系数与上限。写进 plan.json 的 Meta 里：
        // 同一件用不同系数导出的 json，光看 Points/GunOpen 完全分不出来，而枪的姿态差几十度。
        private double _effTiltK = double.NaN;
        private double _effTiltMaxDeg = double.NaN;

        public MainWindow()
        {
            InitializeComponent();
            InstallSafeCopyCommand();
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
        /// 接管 3D 视口的 Ctrl+C 复制命令（替换 HelixToolkit 自带的 Copy 绑定）。
        /// 起因（2026-10-28 闪退）：HelixToolkit 3.1.2 的 Copy 直接调 Clipboard.SetImage，
        /// .NET 9 已移除 BinaryFormatter，WPF 剪贴板在这条路径上会抛 PlatformNotSupportedException；
        /// 原绑定没有 try/catch，异常冒到 Dispatcher 后进程直接退出。
        /// 这里改成自家实现：渲染失败、剪贴板失败都只提示，绝不终止进程；
        /// 剪贴板不可用时退化为把视图存成 PNG。
        /// </summary>
        void InstallSafeCopyCommand()
        {
            CommandBinding? old = Viewport.CommandBindings.OfType<CommandBinding>()
                .FirstOrDefault(b => b.Command == ApplicationCommands.Copy);
            if (old != null) Viewport.CommandBindings.Remove(old);
            Viewport.CommandBindings.Add(new CommandBinding(ApplicationCommands.Copy, SafeCopyViewportCommand));
        }

        void SafeCopyViewportCommand(object sender, ExecutedRoutedEventArgs e)
        {
            e.Handled = true;
            var v = Viewport.Viewport;
            if (v.ActualWidth < 2 || v.ActualHeight < 2)
            {
                TxtHint.Text = "3D 视图还没布局好，暂时无法复制";
                return;
            }

            BitmapSource bmp;
            try
            {
                // 与 HelixToolkit 默认行为一致：2 倍超采样、白底
                bmp = v.RenderBitmap(v.ActualWidth * 2, v.ActualHeight * 2, Brushes.White, 2);
            }
            catch (Exception ex)
            {
                App.Log("COPY/RENDER-FAIL", ex);
                TxtHint.Text = "复制失败：渲染 3D 视图出错（详见日志）";
                return;
            }

            // 剪贴板经常被别的进程短暂占用（OpenClipboard 返回 CLIPBRD_E_CANT_OPEN）。
            // .NET 9 的 WPF 在这种状态下会继续走 OLE 持久化、抛 BinaryFormatter 已移除的异常
            // （2026-10-28 闪退就是这么来的）。所以先重试几次，仍失败再落盘 PNG，绝不外抛。
            Exception? clipErr = null;
            for (int attempt = 0; attempt < 5; attempt++)
            {
                try
                {
                    Clipboard.SetImage(bmp);
                    clipErr = null;
                    break;
                }
                catch (Exception ex)
                {
                    clipErr = ex;
                    System.Threading.Thread.Sleep(120);
                }
            }

            if (clipErr == null)
            {
                TxtHint.Text = "已复制 3D 视图到剪贴板";
                App.Log("COPY/OK", null, $"{bmp.PixelWidth}x{bmp.PixelHeight}");
                return;
            }

            App.Log("COPY/CLIPBOARD-FAIL", clipErr);
            // 兜底：先 TEMP，写不进去再退到桌面，保证图能落下来
            string[] dirs =
            {
                System.IO.Path.GetTempPath(),
                Environment.GetFolderPath(Environment.SpecialFolder.DesktopDirectory),
            };
            string name = "SprayTrajViewer_copy_" + DateTime.Now.ToString("HHmmss") + ".png";
            foreach (string dir in dirs)
            {
                try
                {
                    string path = System.IO.Path.Combine(dir, name);
                    var enc = new PngBitmapEncoder();
                    enc.Frames.Add(BitmapFrame.Create(bmp));
                    using (var fs = System.IO.File.Create(path)) enc.Save(fs);
                    TxtHint.Text = "剪贴板被占用，已把视图存成 PNG：" + path;
                    App.Log("COPY/FALLBACK", null, path);
                    return;
                }
                catch (Exception ex2)
                {
                    App.Log("COPY/FALLBACK-FAIL", ex2);
                }
            }

            TxtHint.Text = "复制失败（剪贴板被占用，PNG 也写不出去）：" + clipErr.Message;
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
            // --reach 0：不要求可达，只查几何（等价于取消勾选「要求可达」，与 spray_safety -reach 0 一致）。
            // 不传该参数 = 保持界面默认（勾选），所以历史脚本的验收结果不会变。
            if (o.TryGetValue("reach", out v)) ChkRequireReach.IsChecked = !(int.TryParse(v, out var rq) && rq == 0);
            // --fast 0|1：姿态求解策略（等价于勾选/取消「单次逆解」，与 spray_safety -fast 一致）。
            // 不传 = 保持界面默认（勾选 = 单次逆解）。
            if (o.TryGetValue("fast", out v)) ChkFastPose.IsChecked = !(int.TryParse(v, out var fp) && fp == 0);
            // --tiltmode 0|1：1 = 倾角逐点连续（tilt=60°×δ/90，无面级常量），与 spray_safety -tiltmode 一致
            if (o.TryGetValue("tiltmode", out v)) ChkTiltPerPoint.IsChecked = int.TryParse(v, out var tm) && tm != 0;
            // --tiltk 60：倾角系数（δ=90° 侧壁时给多少度）。不传 = 库默认 60，等价 spray_safety -tilt。
            // 做"系数扫一遍"的 A/B 时用它，不必改代码重编。
            if (o.TryGetValue("tiltk", out v) &&
                double.TryParse(v, System.Globalization.NumberStyles.Float,
                                System.Globalization.CultureInfo.InvariantCulture, out var tk) && tk > 0)
                _tiltCoefDeg = tk;
            // --auto-export 0|1：自动导出（点云同名 CSV + 模拟平台 JSON）
            if (o.TryGetValue("auto-export", out v)) ChkAutoExport.IsChecked = !(int.TryParse(v, out var ae) && ae == 0);
            // 注意：--safety 是"跑安全评估"的开关（无值），不能拿它覆盖告警阈值，
            // 否则 ParseArgs 会把无值开关当成 "1"，把 20mm 告警阈值误改成 1mm。
            // 要覆盖阈值请用 --safety-gap 20。
            if (o.TryGetValue("safety-gap", out v)) TxtSafetyGap.Text = v;
            if (o.TryGetValue("hardclr", out v)) TxtHardClearance.Text = v;

            // --handeye <矩阵文件>：加载点云后按手眼标定搬到基座系再挪到可达位。
            // 无人值守验收要用它——否则命令行跑出来的结论全是相机系口径的，
            // 和界面上手点出来的对不上。给值即自动勾选（省掉 --handeye 1 这种多余写法）。
            if (o.TryGetValue("handeye", out v) && v != "0")
            {
                if (v != "1") TxtHandEye.Text = v;
                ChkHandEye.IsChecked = true;
            }
            // --facex 450 / --floorz -200：可达位摆位（仅在启用 --handeye 时生效）
            // --floorz 也接受字面量 auto：按工件自身高度让中心落到 0.35×reach。
            if (o.TryGetValue("facex", out v) &&
                double.TryParse(v, System.Globalization.NumberStyles.Float,
                                System.Globalization.CultureInfo.InvariantCulture, out var fxv))
                TxtFaceX.Text = fxv.ToString(System.Globalization.CultureInfo.InvariantCulture);
            if (o.TryGetValue("floorz", out v))
            {
                TxtFloorZ.Text = v.Trim().Equals("auto", StringComparison.OrdinalIgnoreCase)
                    ? "auto"
                    : (double.TryParse(v, System.Globalization.NumberStyles.Float,
                                       System.Globalization.CultureInfo.InvariantCulture, out var fzv)
                        ? fzv.ToString(System.Globalization.CultureInfo.InvariantCulture)
                        : TxtFloorZ.Text);
            }
            // --no-placement / --placement 0：跳过可达位平移，直接用手眼算出的基座坐标。
            // 默认不跳过（保持历史口径），要"看真实位置"的验收必须显式打开。
            if (o.TryGetValue("no-placement", out v))
                ChkSkipPlacement.IsChecked = v != "0";
            else if (o.TryGetValue("placement", out v) && int.TryParse(v, out var pl))
                ChkSkipPlacement.IsChecked = pl == 0;
            // --base-frame 1：点云已在基座系（恒等）。手眼变换会自动置此项，
            // 单独给它是为了处理「点云早就被 python 脚本变换过」的存量文件。
            if (o.TryGetValue("base-frame", out v) && int.TryParse(v, out var bf) && bf != 0)
                CmbBaseMode.SelectedIndex = 2;
            // --up +Z|-Y|+Y|-Z|custom：重力向上。手眼变换后点云在基座系，必须是 +Z；
            // 存量相机系点云才需要指定别的方向。custom 时读 --upxyz x y z。
            if (o.TryGetValue("up", out v))
            {
                switch (v.ToUpperInvariant())
                {
                    case "+Z": case "Z": CmbUpDir.SelectedIndex = 0; break;
                    case "-Y": CmbUpDir.SelectedIndex = 1; break;
                    case "+Y": case "Y": CmbUpDir.SelectedIndex = 2; break;
                    case "-Z": CmbUpDir.SelectedIndex = 3; break;
                    case "CUSTOM":
                        CmbUpDir.SelectedIndex = 4;
                        if (o.TryGetValue("upxyz", out var uxyz))
                        {
                            var t = uxyz.Split(new[] { ',', ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            if (t.Length >= 3)
                            {
                                TxtUpX.Text = t[0]; TxtUpY.Text = t[1]; TxtUpZ.Text = t[2];
                            }
                        }
                        break;
                }
            }

            // --sprayside 已废弃：可喷侧在相机系恒为 −Z，核心库默认按此处理，不再可配。
            if (o.TryGetValue("sprayside", out v))
                App.Log("GEN", null, $"忽略 --sprayside {v}：可喷侧恒为 −Z（朝相机侧）");

            // 显示项（显式给值时才改，避免把用户默认勾选状态打乱）
            if (o.TryGetValue("show-cloud", out v)) ChkShowCloud.IsChecked = v != "0";
            if (o.TryGetValue("show-traj", out v)) ChkShowTraj.IsChecked = v != "0";
            if (o.TryGetValue("show-plane", out v)) ChkShowPlane.IsChecked = v != "0";
            // --show-move 0：只留喷涂段（红），空移段（蓝）藏起来 —— 查行距/覆盖时用。
            // 修正层也会跟着藏（青/紫成对），保证两层口径一致。
            if (o.TryGetValue("show-move", out v)) ChkShowMove.IsChecked = v != "0";
            // --show-fixed 0|1：修正后轨迹（喷涂青 / 空移紫）
            if (o.TryGetValue("show-fixed", out v)) ChkShowFixed.IsChecked = v != "0";
            // --show-risk 0：藏掉风险点，只看轨迹本身（截图/讲解时用）
            if (o.TryGetValue("show-risk", out v)) ChkShowRisk.IsChecked = v != "0";

            App.Log("CLI", null, $"ply={TxtPath.Text}");
            App.Log("CLI", null,
                $"flat={TxtFlatness.Text} curv={TxtCurvature.Text} fill={TxtFillHole.Text} " +
                $"row={TxtRowSpacing.Text} pt={TxtPointSpacing.Text} standoff={TxtStandoff.Text} " +
                $"gc={TxtGunClear.Text} ee={TxtEdgeExtend.Text} ea={TxtEdgeExtendAcross.Text} st={TxtShellThick.Text} " +
                $"dir={CmbScanDir.SelectedIndex} flip={CmbNormalFlip.SelectedIndex}");
            App.Log("CLI", null,
                $"show: cloud={ChkShowCloud.IsChecked} traj={ChkShowTraj.IsChecked} plane={ChkShowPlane.IsChecked}");
            App.Log("CLI", null, $"safetyGap={TxtSafetyGap.Text} hardClr={TxtHardClearance.Text} requireReach={ChkRequireReach.IsChecked}");
            App.Log("CLI", null, $"fastPose={ChkFastPose.IsChecked} autoExport={ChkAutoExport.IsChecked}");

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

                    // 这个标志只该抑制**本次启动**的那一次自动评估（避免和截图抢时序）。
                    // 以前置 true 之后永不复位，于是"从命令行打开过一次"的窗口里，手动点生成
                    // 再也不会跑安全评估 —— 而界面上「生成后自动评估」还是勾着的：
                    // 结果是自动导出把**未修正**的 plan.json 写进 traj_out，还顺手删掉了上一次的
                    // _safe.csv，用户看到的修正层自然一直没有（16:18 垃圾桶就这么丢了一次）。
                    _suppressAutoSafety = false;

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

        void BtnHandEyeBrowse_Click(object sender, RoutedEventArgs e)
        {
            var dlg = new OpenFileDialog
            {
                Title = "选择手眼标定矩阵（相机系 → 基座系）",
                Filter = "标定矩阵 (*.txt;*.csv;*.json)|*.txt;*.csv;*.json|所有文件 (*.*)|*.*",
                CheckFileExists = true,
                RestoreDirectory = false,
            };
            if (dlg.ShowDialog() == true)
            {
                TxtHandEye.Text = dlg.FileName;
                // 选完就立刻生效：让人马上看见"矩阵到底有没有用上"，
                // 而不是等到下一次加载点云才发现勾选框没勾。
                ApplyHandEyeIfNeeded();
            }
        }

        /// <summary>
        /// 勾选框变化：只重算「基座位姿 + 重力向上」这两项口径，点云不动。
        ///
        /// 这里不需要"退回原始点云"那一步 —— 新方案下点云从来没被改过，
        /// 勾/取消勾影响的是安全层与导出用的是哪个 T_base_in_cloud，
        /// 以及倾角按哪个方向当"上"。
        /// </summary>
        void ChkHandEye_Changed(object sender, RoutedEventArgs e)
        {
            if (_applyingHandEye) return;

            if (ChkHandEye.IsChecked != true)
            {
                _handEyeM = null;
                _handEyeBasePose = null;
                _handEyeExportPose = null;
                CmbBaseMode.SelectedIndex = 0;   // 回自动假设
                CmbUpDir.SelectedIndex = 0;      // 回 +Z
                UpdatePlacementInputsEnabled();  // 手眼关了 ⇒ faceX/floorZ 本就不起作用，统一灰掉
                InvalidateTrajectory("手眼标定已关闭");
                // 点云/平面要退回相机系显示 —— 少了这一句，界面还停在基座系的画面上，
                // 而 _handEyeExportPose 已经是 null，后续任何重画都会跳变一次。
                RefreshDisplayFrame("手眼标定已关闭");
                TxtHint.Text = "手眼标定已关闭：安全层改回按界面「基座位姿」推算；视图回到相机系";
                App.Log("HANDEYE/OFF", null, "关闭手眼，点云始终是相机系原始点云（未被改过）");
                return;
            }

            // 勾上但还没有点云：口径已经该切了（faceX/floorZ 变可用），只是没有矩阵可算
            UpdatePlacementInputsEnabled();
            if (_cloud.Count == 0) return;
            ApplyHandEyeIfNeeded();
        }

        /// <summary>
        /// 「正面到基座」/「最低点高度」改完（失焦）→ 重算可达位平移并重画视图。
        ///
        /// 这两个数只影响**平移**（把工件整体挪到机器人舒适区），不改变姿态，
        /// 所以只需重跑 ComputeBasePoseToReachable 再重画，不必重算点云与算法。
        /// 用 LostFocus 而不是 TextChanged：逐字符触发会让工件在输入 "450" 的过程中
        /// 先跳到 4mm、再跳到 45mm，视图一直在闪。
        ///
        /// 手眼没启用时这两个输入框在数学上不起作用（没有 T_cam_in_base 可叠加），
        /// 直接返回，不去打扰用户。
        /// </summary>
        void TxtPlacement_LostFocus(object sender, RoutedEventArgs e)
        {
            if (_applyingHandEye) return;
            // 手眼未生效 → 可达位平移无意义，静默跳过
            if (ChkHandEye.IsChecked != true) return;
            if (_handEyeM == null) return;
            // 勾了「跳过可达位平移」→ 这两个输入框不参与计算，改了也没用，不去重算
            if (ChkSkipPlacement.IsChecked == true) return;

            // 记下改前的值，避免"点了下框又离开"也触发一次全量重画
            string sig = $"{TxtFaceX.Text?.Trim()}|{TxtFloorZ.Text?.Trim()}";
            if (sig == _placementSig) return;
            _placementSig = sig;

            // 平移变了 ⇒ 这些坐标系下的结论全部作废，交给 ApplyHandEyeIfNeeded 统一处理
            // （它会重算矩阵 → 切口径 → InvalidateTrajectory → RefreshDisplayFrame）
            ApplyHandEyeIfNeeded();
        }

        /// <summary>上一次「正面到基座 / 最低点高度」的文本快照，用来判断是否真的改了。</summary>
        private string _placementSig = "";

        /// <summary>
        /// 「跳过可达位平移」勾选变化：重算平移 + 重画视图。
        ///
        /// 为什么必须重画：这个开关只改 dx/dy/dz，不改点云本身。不重画的话
        /// 屏幕上工件一点不动，人会以为"勾了没反应"（矩阵字段已经变了，没人消费）。
        /// 旧轨迹也要作废 —— 它是按旧摆放算出来的，留着会与新摆放差好几百毫米。
        /// </summary>
        void ChkSkipPlacement_Changed(object sender, RoutedEventArgs e)
        {
            if (_applyingHandEye) return;

            // 输入框可用态跟着走：勾上后 faceX / floorZ 不参与计算，灰掉才不会误导
            UpdatePlacementInputsEnabled();

            if (ChkHandEye.IsChecked != true || _handEyeM == null) return;

            InvalidateTrajectory("可达位平移口径已变更");
            ApplyHandEyeIfNeeded();
        }

        /// <summary>
        /// 「正面到基座 / 最低点高度」的可用态：跳过平移时它们不参与计算，灰掉。
        /// 只在手眼生效时才有意义 —— 手眼没开时这两个框本来就不起作用。
        /// </summary>
        void UpdatePlacementInputsEnabled()
        {
            if (TxtFaceX == null || TxtFloorZ == null || ChkSkipPlacement == null) return;
            bool on = ChkHandEye.IsChecked == true && ChkSkipPlacement.IsChecked != true;
            TxtFaceX.IsEnabled = on;
            TxtFloorZ.IsEnabled = on;
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
                _exportDirResolved = null;
                // 原样留一份：手眼变换每次从它重算，保证反复触发不会叠加。
                // 必须是**拷贝**而不是引用赋值 —— 同一份 List 被 _cloud 原地改写时
                // （比如以后加去噪/平滑）会把"原始点云"一起污染掉，那时就再也退不回去了。
                _rawCloud = new List<(double, double, double)>(_cloud);
                // 记住这次文件的目录：下次点"浏览"直接从这里开始（含命令行 --ply 打开的情况）
                try
                {
                    var dir = System.IO.Path.GetDirectoryName(System.IO.Path.GetFullPath(path));
                    if (!string.IsNullOrEmpty(dir)) UiState.Set("lastCloudDir", dir);
                }
                catch { }
                RenderCloud();
                TxtHint.Text = $"点云已加载：{_cloud.Count} 点（读取 {swLoad.ElapsedMilliseconds} ms）";
                // 坐标系标签在点云加载时就同步一次，免得首次生成前面板一直停在初始文案
                if (TxtFrameTag != null) TxtFrameTag.Text = $"视图坐标系：{DisplayFrameName()}";
                ShowOutDir();
                Title = $"喷涂轨迹生成器 - {System.IO.Path.GetFileName(path)}";
                App.Log("LOAD", null, $"{path} -> {_cloud.Count} 点，读取耗时 {swLoad.ElapsedMilliseconds} ms");

                // 新方案下这里**不改点云**，只是把「基座位姿 + 重力向上」这两项口径
                // 与手眼矩阵对齐。点云原样留在相机系，生成层口径与不勾手眼时完全一致
                // （实测同一份点云：勾/不勾都产出 12571 点；旧方案预变换会变成 14025 点）。
                ApplyHandEyeIfNeeded();
            }
            catch (Exception ex)
            {
                App.Log("LOAD/FAIL", ex);
                MessageBox.Show($"点云读取失败：{ex.Message}", "错误", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        // ==================== 手眼标定：只算 basePose，不动点云 ====================
        //
        // 为什么**不能**把点云预先变换到基座系（2026-11-10 修正）：
        //   生成层的可喷窗口是「面法向与 +Z 夹角 ≤ 130°」，这里的 +Z 是**点云自己的坐标轴**。
        //   相机正对着拍时，点云 +Z 就是光轴，也就是"拍摄方向"，主面 angZ ≈ 4°（正对相机）。
        //   一旦把手眼变换提前做掉，点云的 +Z 变成了"重力向上"，同一个主面的 angZ 从 4.3°
        //   跳到 86.1°，算法眼里的"正对相机的板"变成了"竖直侧壁"——面清单、布行方向、
        //   轨迹长度全变（实测 12571 点 → 14025 点），而且竖直面会掉进 130° 窗口的盲区
        //   （正反对称，翻不翻都放行），薄壁误切造出来的反向面因此拦不住。
        //
        // 正确分工：
        //   · 生成层：吃**相机系原始点云**，+Z 就是拍摄方向，语义与原设计一致；
        //   · 安全层：拿 basePose（点云→基座）和 upCloud（点云里哪边是重力向上）自己换算；
        //   · 导出：轨迹本来在点云系，写 plan.json / CSV 前用 basePose 搬到基座系。
        // 这样"要不要给机器人执行"由安全层和导出环节负责，生成层永远只关心"哪个面朝着相机"。

        /// <summary>原始点云（相机系）。生成层与安全层都用它，点云本身永不被改写。</summary>
        private List<(double x, double y, double z)>? _rawCloud;

        /// <summary>防重入：勾选框引发重载 → 重载又改勾选框，会递归。</summary>
        private bool _applyingHandEye;

        /// <summary>当前生效的手眼矩阵 T_cam_in_base（p_base = R·p_cam + t，4×4 行主序）；null = 未启用。</summary>
        private double[]? _handEyeM;

        /// <summary>
        /// T_base_in_cloud（= 手眼矩阵的逆 + 可达位平移），直接喂安全层。
        /// 注意方向：C++ 的 toBaseP(p) = Rᵀ·(p − t)，与手眼文件**互为逆**。
        /// null = 未启用手眼（那时走界面上的「基座位姿」选项）。
        /// </summary>
        private double[]? _handEyeBasePose;

        /// <summary>
        /// T_cam_in_base（= 手眼矩阵 + 可达位平移），导出时把轨迹从点云系搬到基座系。
        /// 与 _handEyeBasePose 互为逆矩阵。
        /// </summary>
        private double[]? _handEyeExportPose;

        /// <summary>
        /// 手眼配置总入口（只做步骤调度）。
        ///
        /// 注意它**不碰 _cloud**。它做三件事：读矩阵、把它记成 _handEyeM / _handEyeBasePose、
        /// 把界面的「基座位姿」与「重力向上」设成与矩阵配套的口径。
        /// </summary>
        void ApplyHandEyeIfNeeded()
        {
            if (ChkHandEye.IsChecked != true)
            {
                _handEyeM = null;
                _handEyeBasePose = null;
                _handEyeExportPose = null;
                return;
            }
            if (_applyingHandEye) return;

            try
            {
                _applyingHandEye = true;

                // ---------- 步骤 1：读手眼矩阵 p_base = R·p_cam + t ----------
                string? err;
                double[]? M = ReadHandEyeMatrix(TxtHandEye.Text, out err);
                if (M == null)
                {
                    _handEyeM = null;
                    _handEyeBasePose = null;
                    _handEyeExportPose = null;
                    TxtHint.Text = "手眼标定未生效：" + err;
                    App.Log("HANDEYE/SKIP", null, err ?? "");
                    return;
                }
                _handEyeM = M;

                // ---------- 步骤 2：重力向上在点云系里指向哪（R 第 3 行）----------
                LogUpVectorInCloud(M);

                // ---------- 步骤 3：叠加可达位平移，得到安全层与导出各自的矩阵 ----------
                // 两个矩阵互为逆：安全层要 T_base_in_cloud，导出要 T_cam_in_base。
                // 只改矩阵、点云一个点都不动 —— 生成层拿到的永远是相机系原始点云。
                ComputeBasePoseToReachable(M, out var bp, out var ep);
                _handEyeBasePose = bp;
                _handEyeExportPose = ep;
                // 记下当前摆放参数，供 TxtPlacement_LostFocus 判断"到底改没改"。
                // 放在这里而不是构造函数：用户换了手眼矩阵文件但没改摆放参数时，
                // 不该被当成"摆放改了"再重画一次。
                _placementSig = $"{TxtFaceX.Text?.Trim()}|{TxtFloorZ.Text?.Trim()}";

                // ---------- 步骤 4：把界面口径设成与该矩阵配套 ----------
                SwitchToHandEyeConventions(M);
                // 勾了「跳过可达位平移」时 faceX / floorZ 不参与计算，灰掉避免误导
                UpdatePlacementInputsEnabled();

                // ---------- 步骤 5：按新坐标系重画视图 ----------
                // 必须先重画再汇报：这条路径上头一次生效时视图还是相机系的，
                // 不重画的话"勾了手眼"看上去毫无变化（矩阵只写进了字段，没人消费）。
                RefreshDisplayFrame("手眼标定生效");

                // ---------- 步骤 6：汇报 ----------
                TxtHint.Text = $"手眼标定已生效：视图已切到{DisplayFrameName()}，点云内部保持相机系；" +
                               $"重力向上 = 点云 ({M[8]:F2}, {M[9]:F2}, {M[10]:F2})";
                App.Log("HANDEYE", null,
                    "点云保持相机系，生成层口径不变；安全层用 basePose + upCloud 换算到基座系");
            }
            catch (Exception ex)
            {
                _handEyeM = null;
                _handEyeBasePose = null;
                _handEyeExportPose = null;
                App.Log("HANDEYE/FAIL", ex);
                MessageBox.Show($"手眼标定读取失败：{ex.Message}", "错误", MessageBoxButton.OK, MessageBoxImage.Error);
            }
            finally
            {
                _applyingHandEye = false;
            }
        }


        /// <summary>
        /// 步骤 1：读手眼标定矩阵，返回 4×4 行主序（16 个 double）。
        ///
        /// 约定 p_base = R·p_cam + t，文件可以是 4×4 或 3×4；<c>#</c> 之后当注释。
        /// 校验只告警不拦：现场标定本身就有几度误差，硬拦下来反而没法用；
        /// 但偏离正交到 1e-2 以上基本就是抄错/行列颠倒，必须让人看见。
        /// </summary>
        double[]? ReadHandEyeMatrix(string spec, out string? err)
        {
            err = null;
            string path = (spec ?? "").Trim();
            if (path.Length == 0) { err = "未填写手眼矩阵文件"; return null; }

            // 相对路径优先按程序目录找，其次按当前目录 —— 双击启动和调试启动的 cwd 不同，
            // 只认一个会让人莫名其妙地"文件明明在却读不到"。
            if (!System.IO.Path.IsPathRooted(path))
            {
                string baseDir = AppDomain.CurrentDomain.BaseDirectory;
                string p1 = System.IO.Path.Combine(baseDir, path);
                if (System.IO.File.Exists(p1)) path = p1;
            }
            if (!System.IO.File.Exists(path)) { err = $"手眼矩阵文件不存在：{spec}"; return null; }

            string raw = System.IO.File.ReadAllText(path);
            raw = string.Join("\n", raw.Split('\n').Select(ln =>
            {
                int k = ln.IndexOf('#');
                return k >= 0 ? ln.Substring(0, k) : ln;
            }));

            var nums = new List<double>();
            foreach (System.Text.RegularExpressions.Match m in
                     System.Text.RegularExpressions.Regex.Matches(raw, @"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"))
            {
                if (double.TryParse(m.Value, System.Globalization.NumberStyles.Float,
                                    System.Globalization.CultureInfo.InvariantCulture, out double v))
                    nums.Add(v);
            }
            if (nums.Count < 12)
            {
                err = $"手眼矩阵至少要 12 个数（3×4），文件里只解析到 {nums.Count} 个";
                return null;
            }

            var M = new double[16];
            M[15] = 1.0;
            if (nums.Count >= 16)
                for (int i = 0; i < 16; i++) M[i] = nums[i];
            else
                for (int i = 0; i < 12; i++) M[i] = nums[i];

            // 正交性/行列式自检：det 不接近 +1 说明不是纯旋转（含缩放或镜像），
            // 那种矩阵把点云搬过去会变形，比例、尺寸全错。
            double o = 0;
            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 3; c++)
                {
                    double dot = 0;
                    for (int k = 0; k < 3; k++) dot += M[r * 4 + k] * M[c * 4 + k];
                    o = Math.Max(o, Math.Abs(dot - (r == c ? 1.0 : 0.0)));
                }
            double det =
                M[0] * (M[5] * M[10] - M[6] * M[9]) -
                M[1] * (M[4] * M[10] - M[6] * M[8]) +
                M[2] * (M[4] * M[9] - M[5] * M[8]);

            App.Log("HANDEYE/M", null,
                $"{path} |RR^T-I|max={o:G4} det={det:F6} t=({M[3]:F2},{M[7]:F2},{M[11]:F2})");
            if (o > 1e-2)
            {
                MessageBox.Show(
                    $"手眼矩阵的旋转部分不像纯旋转（|RRᵀ−I| = {o:F4g}）。\n\n" +
                    "常见原因：行列抄颠倒、把 3×4 当成 4×3 读、或者混进了缩放。\n" +
                    "点云搬过去会变形，尺寸和倾角都会偏。请先核对矩阵再继续。",
                    "手眼矩阵校验", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
            if (Math.Abs(det - 1.0) > 1e-2)
            {
                App.Log("HANDEYE/DET", null, $"det={det:F4f} ≠ +1");
                MessageBox.Show(
                    $"手眼矩阵 det = {det:F4f}，不等于 +1。\n\n" +
                    "det = −1 表示含镜像（手性反了），det 明显偏离 1 表示含缩放。\n" +
                    "这两种情况点云都会变形，建议先修正标定再继续。",
                    "手眼矩阵校验", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
            return M;
        }

        /// <summary>
        /// 步骤 2：算出「重力向上在点云系里指向哪」，记进日志供人工核对。
        ///
        /// 已知 p_base = R·p_cam + t，则基座 +Z（重力向上）在点云系里 = **Rᵀ·ẑ**，
        /// 写成行向量正好是 **R 的第 3 行** (M[8], M[9], M[10])。
        /// 不是第 3 列 —— 第 3 列是"相机光轴在基座系里的朝向"，两者在本机实测差 149.4°，
        /// 拿反了整件工件的倾角口径会全部反过来，而且位置、可达率都正常，只表现为
        /// "明明立着挂的板，算法当成平放"。
        ///
        /// 本标定文件算出来是 (0,1,0)：点云 y 轴竖直向上，与「工件立着挂、相机正对拍」一致。
        /// </summary>
        void LogUpVectorInCloud(double[] M)
        {
            double ux = M[8], uy = M[9], uz = M[10];
            double L = Math.Sqrt(ux * ux + uy * uy + uz * uz);
            if (L < 1e-9) return;
            App.Log("HANDEYE/UP", null,
                $"重力向上(点云系) = R 第3行 = ({ux / L:F4}, {uy / L:F4}, {uz / L:F4})");
        }

        /// <summary>
        /// 步骤 3：算可达位平移，并同时产出**两个方向相反**的矩阵。
        ///
        /// ⚠ 这里有个必须讲清楚的坑：手眼文件与安全层用的不是同一个方向的矩阵。
        ///
        ///   手眼文件（标定软件给的）：p_base = R·p_cam + t        ← T_cam_in_base
        ///   安全层 basePose（C++ robot_safety.cpp:1160）：
        ///       toBaseP(p) = Rᵀ·(p − t)                          ← T_base_in_cloud
        /// 两者互为逆。把文件里的矩阵原样喂进安全层，等于传了个逆矩阵，
        /// 表现是 **reachable=0、全部不可达、最小间隙 −6.2mm**（本机实测），
        /// 但 baseAuto=0 说明"矩阵确实传进去了"—— 很容易误判成"点云不适配这台机器人"。
        ///
        /// 换算（要满足 toBaseP(p) = R_h·p + t_h + d）：
        ///   R_base = R_hᵀ
        ///   t_base = −R_hᵀ·(t_h + d)
        ///
        /// 平移口径（与 align_cloud_to_robot.py 一致，只是这里算而不改点云）：
        ///   · 正面（min X 那一侧）贴到 faceX —— 决定机器人要伸多远；
        ///   · Y 居中到 0 —— 别让工件偏向一侧，平白吃掉一侧的臂展；
        ///   · 最低点落到 floorZ —— 决定工件高低，立式件靠它把中心顶到舒适高度。
        ///
        /// 为什么不是"质心对准舒适区"：质心对中了，正面还得看它有多厚。
        /// 臂展约束卡的是**正面**到基座的距离（腕心 = 正面 − 喷距），
        /// 质心对中反而在厚件上把正面推远。所以这里显式给正面位置。
        /// </summary>
        /// <param name="M">手眼矩阵 T_cam_in_base（p_base = R·p_cam + t）</param>
        /// <param name="basePoseOut">T_base_in_cloud，喂安全层</param>
        /// <param name="exportOut">T_cam_in_base（含可达位平移），导出时把轨迹搬到基座系</param>
        void ComputeBasePoseToReachable(double[] M, out double[] basePoseOut, out double[] exportOut)
        {
            // ---------- 步骤 1：算手眼变换后的基座系包围盒 ----------
            // 正面/最低点是**基座系**里说的，不能拿相机系包围盒直接减 faceX
            // （两套轴的方向都不同）。
            var box = ComputeBaseBox(M);

            // ---------- 步骤 2：算可达位平移 ----------
            var shift = ComputePlacementShift(box);

            // ---------- 步骤 3：组装导出用矩阵（T_cam_in_base + 平移） ----------
            exportOut = BuildExportPose(M, shift);

            // ---------- 步骤 4：组装安全层用矩阵（转置 + 反解平移） ----------
            basePoseOut = BuildBasePose(M, exportOut);

            // ---------- 步骤 5：互逆自检 ----------
            VerifyPosePair(basePoseOut, exportOut);
        }

        /// <summary>VS-087 臂展（mm）。舒适区口径与 shift_cloud_to_denso.py / planner.py::auto_base 一致。</summary>
        private const double RobotReachMm = 905.0;

        /// <summary>
        /// 步骤 1：手眼变换后的基座系包围盒。
        ///
        /// 空点云时各分量返回 ±Infinity/∓Infinity，由 <see cref="ComputePlacementShift"/>
        /// 判 IsFinite 后按「不平移」处理 —— 与旧行为一致，不会因为除零或 NaN 把矩阵搞坏。
        /// </summary>
        (double loX, double hiX, double loY, double hiY, double loZ, double hiZ) ComputeBaseBox(double[] M)
        {
            double loX = double.PositiveInfinity, hiX = double.NegativeInfinity;
            double loY = double.PositiveInfinity, hiY = double.NegativeInfinity;
            double loZ = double.PositiveInfinity, hiZ = double.NegativeInfinity;
            var src = _rawCloud ?? _cloud;
            for (int i = 0; i < src.Count; i++)
            {
                var p = src[i];
                double x = M[0] * p.x + M[1] * p.y + M[2] * p.z + M[3];
                double y = M[4] * p.x + M[5] * p.y + M[6] * p.z + M[7];
                double z = M[8] * p.x + M[9] * p.y + M[10] * p.z + M[11];
                if (x < loX) loX = x; if (x > hiX) hiX = x;
                if (y < loY) loY = y; if (y > hiY) hiY = y;
                if (z < loZ) loZ = z; if (z > hiZ) hiZ = z;
            }
            return (loX, hiX, loY, hiY, loZ, hiZ);
        }

        /// <summary>
        /// 步骤 2：算「可达位平移」这一层的位移量 (dx, dy, dz)。
        ///
        /// ⚠ 这一层与手眼换系是**两件事**，必须分开理解：
        ///   · 手眼矩阵 p_base = R·p_cam + t 做完，点云已经在基座系了，坐标是被测出来的；
        ///   · 本方法再做一次整体平移，把工件挪到机器人「舒适区」—— 这个位置是**规划出来的**，
        ///     不是测出来的，所以它当然会覆盖掉手眼给出的真实位置。
        ///
        /// 勾了「跳过可达位平移」时直接返回全 0：手眼结果就是终点，位置完全可追溯。
        /// 没有这个开关，矮件会被 floorZ 的绝对值压得很低——实测 476 高的槽型盘悬挂
        /// 套用给 860 背板配的 floorZ=−200，中心只剩 +38mm，看起来就像掉到地上。
        /// </summary>
        (double dx, double dy, double dz) ComputePlacementShift(
            (double loX, double hiX, double loY, double hiY, double loZ, double hiZ) box)
        {
            if (!double.IsFinite(box.loX))
            {
                App.Log("HANDEYE/BOX", null, "点云为空，可达位平移取 0");
                return (0, 0, 0);
            }

            // ---------- 分支 A：跳过平移，手眼算到哪就是哪 ----------
            if (ChkSkipPlacement.IsChecked == true)
            {
                App.Log("HANDEYE/BOX", null,
                    $"跳过可达位平移：基座系包围盒 X[{box.loX:F0},{box.hiX:F0}] " +
                    $"Y[{box.loY:F0},{box.hiY:F0}] Z[{box.loZ:F0},{box.hiZ:F0}] " +
                    $"尺寸 {box.hiX - box.loX:F0}×{box.hiY - box.loY:F0}×{box.hiZ - box.loZ:F0} " +
                    "位移(0,0,0)（纯手眼结果）");
                return (0, 0, 0);
            }

            // ---------- 分支 B：按界面参数挪到可达位 ----------
            double floorZ = ResolveFloorZ(box);
            double dx = ParseDouble(TxtFaceX.Text, 450) - box.loX;   // 正面贴到 faceX
            double dy = -0.5 * (box.loY + box.hiY);                  // Y 居中
            double dz = floorZ - box.loZ;                            // 最低点落到 floorZ

            App.Log("HANDEYE/BOX", null,
                $"基座系包围盒 X[{box.loX + dx:F0},{box.hiX + dx:F0}] " +
                $"Y[{box.loY + dy:F0},{box.hiY + dy:F0}] " +
                $"Z[{box.loZ + dz:F0},{box.hiZ + dz:F0}] " +
                $"尺寸 {box.hiX - box.loX:F0}×{box.hiY - box.loY:F0}×{box.hiZ - box.loZ:F0} " +
                $"位移({dx:F0},{dy:F0},{dz:F0})");
            return (dx, dy, dz);
        }

        /// <summary>
        /// 「最低点高度」的取值解析：支持 <c>auto</c>（或留空）表示按工件自身高度自动算。
        ///
        /// auto 的口径 = 让工件**中心**落在舒适高度 0.35×reach（VS-087 ≈ 316.75mm）。
        /// 为什么不再写死 −200：那是给 860 高背板配的绝对值，套到任意高度的工件上都会偏；
        /// auto 只关心「中心落在哪」，换工件高度不用重新手算，也不会出现矮件被压到地上。
        /// </summary>
        double ResolveFloorZ((double loX, double hiX, double loY, double hiY, double loZ, double hiZ) box)
        {
            string raw = (TxtFloorZ.Text ?? "").Trim();
            if (raw.Length == 0 || raw.Equals("auto", StringComparison.OrdinalIgnoreCase))
            {
                double comfort = 0.35 * RobotReachMm;
                double h = box.hiZ - box.loZ;
                double floorZ = comfort - 0.5 * h;
                App.Log("HANDEYE/BOX", null,
                    $"最低点高度=auto：高度 {h:F1} → floorZ={floorZ:F1}（中心落到 {comfort:F1}）");
                return floorZ;
            }
            return ParseDouble(raw, -200);
        }

        /// <summary>步骤 3：导出用矩阵 = 手眼 R/t 叠加可达位平移。</summary>
        double[] BuildExportPose(double[] M, (double dx, double dy, double dz) shift)
        {
            double tx = M[3] + shift.dx, ty = M[7] + shift.dy, tz = M[11] + shift.dz;
            return new double[]
            {
                M[0], M[1], M[2],  tx,
                M[4], M[5], M[6],  ty,
                M[8], M[9], M[10], tz,
                0,     0,     0,     1,
            };
        }

        /// <summary>步骤 4：安全层用矩阵 —— 转置 + 反解平移，t_base = −Rᵀ·(t_h + d)，逐行展开（Rᵀ 的第 r 行 = R 的第 r 列）。</summary>
        double[] BuildBasePose(double[] M, double[] exportPose)
        {
            double tx = exportPose[3], ty = exportPose[7], tz = exportPose[11];
            double bx = -(M[0] * tx + M[4] * ty + M[8] * tz);
            double by = -(M[1] * tx + M[5] * ty + M[9] * tz);
            double bz = -(M[2] * tx + M[6] * ty + M[10] * tz);
            return new double[]
            {
                M[0], M[4], M[8],  bx,
                M[1], M[5], M[9],  by,
                M[2], M[6], M[10], bz,
                0,     0,     0,     1,
            };
        }

        /// <summary>
        /// 步骤 5：自检 —— 两个矩阵必须互逆（basePose · exportPose ≈ I）。
        /// 不核对这一条，方向搞反只会表现为「全部不可达」，看不出是矩阵方向错
        /// （实测就是这么被误导过一次：reachable=0、minClr=−6.2mm，看着像点云不适配）。
        ///
        /// 阈值取 1e-2：本标定文件是屏幕抄录、两位小数，|RRᵀ−I|max 本身就有 1.6e-3，
        /// 该残差会原样出现在这里。方向反会到 O(1)，与 1e-3 差三个数量级，不会混淆。
        /// </summary>
        void VerifyPosePair(double[] basePoseOut, double[] exportOut)
        {
            double errMax = 0;
            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 4; c++)
                {
                    double s = 0;
                    for (int k = 0; k < 4; k++) s += basePoseOut[r * 4 + k] * exportOut[k * 4 + c];
                    double want = (r == c) ? 1.0 : 0.0;
                    errMax = Math.Max(errMax, Math.Abs(s - want));
                }
            App.Log("HANDEYE/CHK", null,
                $"basePose·exportPose 与单位阵最大偏差 = {errMax:E2}" +
                (errMax > 1e-2 ? "  ⚠ 方向/平移不自洽" : " ✓（残差来自标定矩阵非严格正交）"));
        }

        /// <summary>
        /// 步骤 4：把「基座位姿」与「重力向上」设成与手眼矩阵配套的口径。
        ///
        /// 要改这两项，是因为它们**描述的是同一件事的两面**：
        ///   · 基座位姿 = 点云 → 基座的变换（手眼矩阵 + 可达位平移），界面无法用「手输 x/y/z + 偏转」
        ///     表达一个任意旋转，所以这里只把平移分量填进去供人看，真正的矩阵走 _handEyeM；
        ///   · 重力向上 = 该矩阵 R 的第 3 行（本例 +Y），点云里竖直向上的那根轴。
        /// 只改一个不改另一个，是最难看出来的错：位置全对、可达率也正常，
        /// 只有倾角相对一个歪掉 90° 的"向上"算 —— 表现为"明明立着挂的板当成平放"。
        /// </summary>
        void SwitchToHandEyeConventions(double[] M)
        {
            // 基座位姿：切到「手动输入」并把原点喂进去。填的是 _handEyeBasePose 的平移列，
            // 也就是真正喂给安全层的那个数 —— 不能填 M 的 t，那是反的（见该字段的注释）。
            // 旋转无法用「原点 + 朝向偏转」表达，所以这一项只作核对用，真值走 _handEyeBasePose。
            var B = _handEyeBasePose ?? M;
            CmbBaseMode.SelectedIndex = 1;
            TxtBaseX.Text = B[3].ToString("F1", System.Globalization.CultureInfo.InvariantCulture);
            TxtBaseY.Text = B[7].ToString("F1", System.Globalization.CultureInfo.InvariantCulture);
            TxtBaseZ.Text = B[11].ToString("F1", System.Globalization.CultureInfo.InvariantCulture);
            TxtBaseYaw.Text = "0";

            // 重力向上：R 第 3 行。已是 ±X/±Y/±Z 的整数轴时选预设项，否则落到「自定义」。
            double ux = M[8], uy = M[9], uz = M[10];
            double L = Math.Sqrt(ux * ux + uy * uy + uz * uz);
            if (L > 1e-9) { ux /= L; uy /= L; uz /= L; }
            int preset = -1;
            if (Math.Abs(ux) < 0.05 && Math.Abs(uy) < 0.05 && uz > 0.95) preset = 0;        // +Z
            else if (Math.Abs(ux) < 0.05 && Math.Abs(uy) < 0.05 && uz < -0.95) preset = 3;   // −Z
            else if (Math.Abs(ux) < 0.05 && Math.Abs(uz) < 0.05 && uy < -0.95) preset = 1;   // −Y
            else if (Math.Abs(ux) < 0.05 && Math.Abs(uz) < 0.05 && uy > 0.95) preset = 2;    // +Y
            if (preset >= 0)
            {
                CmbUpDir.SelectedIndex = preset;
            }
            else
            {
                CmbUpDir.SelectedIndex = 4;   // 自定义三分量
                TxtUpX.Text = ux.ToString("F4", System.Globalization.CultureInfo.InvariantCulture);
                TxtUpY.Text = uy.ToString("F4", System.Globalization.CultureInfo.InvariantCulture);
                TxtUpZ.Text = uz.ToString("F4", System.Globalization.CultureInfo.InvariantCulture);
            }

            // 旧的轨迹是按旧口径算的，留着会和新口径混在一起
            InvalidateTrajectory("手眼标定口径已变更");
        }

        /// <summary>
        /// 坐标系变了，旧的轨迹与评估结果就作废了。
        ///
        /// 必须主动清掉，不能留着：轨迹是用**相机系**点云算出来的，
        /// 点云搬到基座系后两者相差一个刚体变换（这里约 1.5 米平移 + 90° 级旋转）。
        /// 留在屏幕上会看到"轨迹浮在工件外面一大截"，而结论卡里还是上一次的 PASS ——
        /// 最容易误导人的状态：看上去有结果，其实对不上。
        /// </summary>
        void InvalidateTrajectory(string reason)
        {
            if (_lastTraj.Count == 0 && _safetyPts.Length == 0) return;

            _lastTraj.Clear();
            _safetyPts = Array.Empty<RobotSafetyPoint>();

            foreach (var t in _trajVisuals) Viewport.Children.Remove(t);
            foreach (var t in _moveVisuals) Viewport.Children.Remove(t);
            _trajVisuals.Clear();
            _moveVisuals.Clear();

            // 风险点与修正层也是按**上一次的坐标系**烘死进 Viewport 的，
            // 只清轨迹会留下一团红点/绿线停在旧位置 —— 换坐标系时尤其明显
            // （点云已经在基座系，风险点还堆在相机系的坐标上）。
            foreach (var v in _riskVisuals) Viewport.Children.Remove(v);
            _riskVisuals.Clear();
            if (_fixedSprayVisual != null) { Viewport.Children.Remove(_fixedSprayVisual); _fixedSprayVisual = null; }
            if (_fixedMoveVisual != null) { Viewport.Children.Remove(_fixedMoveVisual); _fixedMoveVisual = null; }
            if (_fixedDroppedVisual != null) { Viewport.Children.Remove(_fixedDroppedVisual); _fixedDroppedVisual = null; }

            TxtResult.Text = "";
            TxtSafetyInline.Text = "尚未评估（" + reason + "，请重新生成轨迹）";
            App.Log("INVALIDATE", null, reason + " → 已清除旧轨迹与评估结果");
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
                // 显示坐标系：手眼生效时搬到基座系，否则原样（见「显示坐标系」段）
                var d = DispP(p.x, p.y, p.z);
                pts.Add(new Point3D(d.x, d.y, d.z));
            }

            _cloudVisual = new PointsVisual3D
            {
                Points = pts,
                Color = Color.FromRgb(224, 64, 224),
                // 点越密，点径越小，否则整片糊成色块把轨迹/面片全遮住
                Size = pts.Count > 120000 ? 1.2 : (pts.Count > 40000 ? 1.8 : 2.5)
            };
            Viewport.Children.Add(_cloudVisual);
            // 每次重画都重设取景：这个函数只在"点云换了"或"坐标系换了"时被调用，
            // 两种情况下旧镜头都不再适用（换坐标系时工件整体平移约 1.5 米并转 90° 级，
            // 保住旧镜头会让人看到一片空白，以为程序挂了）。
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
            // 可喷侧方向不需要选：相机系里 +Z 是光轴、指向景物，相机只看得到朝向自己那一面
            // （n_z < 0），枪又和相机同在一个腕上，所以可喷侧恒定是 −Z。核心库默认就按这个
            // 处理（入口把 −Z 旋到 +Z 算，出口再转回来），这里不再设 SPRAYTRAJ_UP。
            // 注意与上面的「重力向上」区分：那个只喂安全层（倾角/碰撞），与轨迹定向无关。
            App.Log("GEN", null, "可喷侧方向 up=-Z（相机系恒定）");
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
                $"  红色=喷涂行，蓝色=空移（面间抬枪绕行，枪关）\n" +
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

            // 自动导出：生成这一刻先落一份**未修正**的轨迹 + plan.json，
            // 安全评估跑完后再覆盖成修正版（见 RunSafetyCheckAsync 末尾的 AutoExport(true)）。
            // 这样即使用户关掉「自动评估」，DensoArmKine 也有东西可跑。
            //
            // 注意：这里**不再**额外往 exe 目录写一份 trajectory.csv。
            // 那一路是早期没有 traj_out 时的临时出口，两个问题：
            //   1) 提示条写"CSV 已导出：…\bin\Release\net9.0-windows\trajectory.csv"，
            //      看着像"只导了 CSV、没有 JSON"，而真正的 json 在同名目录里；
            //   2) 所有工件都覆盖同一个文件，事后分不清谁是谁。
            AutoExport(false);

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
            foreach (var t in _moveVisuals) Viewport.Children.Remove(t);
            _trajVisuals.Clear();
            _moveVisuals.Clear();

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
                // 显示坐标系：轨迹点也要搬（见「显示坐标系」段）。搬的是**拷贝**，
                // traj 本身（= _lastTraj）不动 —— 生成层与导出用的仍是相机系原值。
                var pa = DispP(traj[i - 1].x, traj[i - 1].y, traj[i - 1].z);
                var pb = DispP(traj[i].x, traj[i].y, traj[i].z);
                var a = new Point3D(pa.x, pa.y, pa.z);
                var b = new Point3D(pb.x, pb.y, pb.z);
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
                // 空移段也画出来（不能跳点），但配色/粗细和喷涂段分开：
                //   喷涂：红、粗 2.5 —— 有效动作
                //   空移：亮蓝、细 1.7 —— 过渡动作（抬枪/跨面/下枪），比原来的灰更容易看见
                var tube = new TubeVisual3D
                {
                    Path = new Point3DCollection(runPts[i]),
                    Diameter = spraying ? 2.5 : 1.7,
                    Fill = new SolidColorBrush(spraying ? SprayColor : MoveColor),
                    ThetaDiv = spraying ? 12 : 10
                };
                if (spraying) _trajVisuals.Add(tube); else _moveVisuals.Add(tube);
            }
            foreach (var t in _trajVisuals) Viewport.Children.Add(t);
            foreach (var t in _moveVisuals) Viewport.Children.Add(t);

            Viewport.ZoomExtents(0);   // 立即取景到全部内容
            ApplyTrajVisibility();

            App.Log("TRAJ", null,
                $"段数={runPts.Count}（喷涂 {runCls.Count(c => c == 1)} / 空移 {runCls.Count(c => c == 0)}，" +
                $"空移点 {traj.Count(p => p.spray != 1)}/{traj.Count}）");
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
                           : (p.sideSrc == 4 ? "薄壁" : (p.sideSrc == 5 ? "双面"
                           : (p.sideSrc == 6 ? "同向" : "兜底"))));
                string fp = p.freePlus < 0 ? "∞" : p.freePlus.ToString("F0");
                string fm = p.freeMinus < 0 ? "∞" : p.freeMinus.ToString("F0");

                string note;
                if (!hasDiag) note = "无诊断数据（旧 DLL）";
                else if (dual) note = $"薄板双面：按壁厚 {p.layerThick:F1}mm 补出的另一面";
                else if (single) note = "单层壳：这一面在点云里只有一层点，没有壁厚";
                else if (p.sideSrc == 6) note = "拆层同向：两片朝向一致，是台阶/两层平台，不是壁的两侧（各自布行）";
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

        void AddPlaneVisual(SprayTrajPlane pIn)
        {
            // 显示坐标系：面片整体搬到基座系（法向/基向量/顶点一起转，见「显示坐标系」段）。
            // 不转的话就会出现"轨迹在基座系、面片还在相机系"——面片会飘到工件外面去。
            var p = DispPlane(pIn);

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

        /// <summary>
        /// 写原始轨迹 CSV（8 列，与本仓 spray_driver 同口径）。
        /// 传进来的 <paramref name="traj"/> 必须已经在**目标坐标系**里了（见 ExportPose）。
        /// </summary>
        static void WriteTrajCsv(string path, List<SprayTrajPoint> traj)
        {
            var sb = new System.Text.StringBuilder("X,Y,Z,I,J,K,speed,spray\n");
            foreach (var p in traj)
                sb.Append($"{p.x:F3},{p.y:F3},{p.z:F3},{p.i:F4},{p.j:F4},{p.k:F4},{p.speed:F1},{p.spray}\n");
            System.IO.File.WriteAllText(path, sb.ToString());
        }

        // ══════════════════════════════════════════════════════════════════
        //  导出前的坐标系搬运
        //
        //  为什么要单独一步：生成层吃的是**相机系原始点云**，它吐出来的轨迹自然也在相机系。
        //  但 DensoArmKine / 示教器要的是**基座系**位姿，两者相差一个刚体变换。
        //  所以「不改点云」这个决定并不等于"导出也不用管坐标系" ——
        //  恰恰相反，坐标系这件事被**整体推迟**到了这一步：点云和轨迹全程留在相机系，
        //  只有落盘那一刻才搬到基座系。
        //
        //  只在于手眼标定生效时才搬：没有标定时根本不知道基座在哪，
        //  硬用一个"按质心猜"的位姿去搬，会让原来能对上的 CSV 突然对不上。
        // ══════════════════════════════════════════════════════════════════

        /// <summary>
        /// 导出/画图要用的 T_cam_in_base（把点云系轨迹搬到基座系）；null = 不搬（保持相机系）。
        /// </summary>
        double[]? ExportPose()
        {
            if (_handEyeExportPose == null) return null;
            return _handEyeExportPose;
        }

        /// <summary>
        /// 把一条轨迹从点云系搬到基座系：位置用 p' = R·p + t，朝向只用 R 转（不平移）。
        ///
        /// IJK 是**方向**，不是点 —— 给它加上平移分量会让枪轴凭空歪掉，
        /// 而且位置完全正确、只有姿态错，是最难发现的一类错。
        /// </summary>
        static List<SprayTrajPoint> ToBaseFrame(List<SprayTrajPoint> src, double[] M)
        {
            var dst = new List<SprayTrajPoint>(src.Count);
            foreach (var q in src)
            {
                var r = q;
                r.x = M[0] * q.x + M[1] * q.y + M[2] * q.z + M[3];
                r.y = M[4] * q.x + M[5] * q.y + M[6] * q.z + M[7];
                r.z = M[8] * q.x + M[9] * q.y + M[10] * q.z + M[11];
                r.i = M[0] * q.i + M[1] * q.j + M[2] * q.k;
                r.j = M[4] * q.i + M[5] * q.j + M[6] * q.k;
                r.k = M[8] * q.i + M[9] * q.j + M[10] * q.k;
                dst.Add(r);
            }
            return dst;
        }

        /// <summary>同 ToBaseFrame，但作用于安全层输出（字段多，只改位置与朝向）。</summary>
        static RobotSafetyPoint[] ToBaseFrame(RobotSafetyPoint[] src, double[] M)
        {
            var dst = new RobotSafetyPoint[src.Length];
            for (int n = 0; n < src.Length; n++)
            {
                var q = src[n];
                double x = q.x, y = q.y, z = q.z, i = q.i, j = q.j, k = q.k;
                q.x = M[0] * x + M[1] * y + M[2] * z + M[3];
                q.y = M[4] * x + M[5] * y + M[6] * z + M[7];
                q.z = M[8] * x + M[9] * y + M[10] * z + M[11];
                q.i = M[0] * i + M[1] * j + M[2] * k;
                q.j = M[4] * i + M[5] * j + M[6] * k;
                q.k = M[8] * i + M[9] * j + M[10] * k;
                dst[n] = q;
            }
            return dst;
        }

        /// <summary>导出用的算法轨迹：手眼生效时搬到基座系，否则原样返回（不拷贝）。</summary>
        List<SprayTrajPoint> ExportTraj()
        {
            var M = ExportPose();
            return M == null ? _lastTraj : ToBaseFrame(_lastTraj, M);
        }

        /// <summary>导出用的修正后轨迹：手眼生效时搬到基座系，否则原样返回。</summary>
        RobotSafetyPoint[] ExportSafetyPts()
        {
            var M = ExportPose();
            return M == null ? _safetyPts : ToBaseFrame(_safetyPts, M);
        }

        // ══════════════════════════════════════════════════════════════════
        //  显示坐标系（3D 视图）—— 与内部数据流是**两件事**
        //
        //  内部数据流（不能动，动了算法口径就变）：
        //    点云 → 相机系；生成层吃相机系（可喷窗口「法向与 +Z 夹角 ≤130°」里的 +Z
        //    是**点云的拍摄方向**，一换到基座系，正对相机的板就被当成竖直侧壁，
        //    面清单/布行方向/轨迹长度全变，实测 12571 → 14025 点）。
        //
        //  显示层（本段负责）：
        //    用户要看的是「工件在机器人工作空间里的什么样、轨迹离它多远、会不会撞」，
        //    这件事**必须**在基座系里看。相机系的点云配基座系的轨迹，看起来就是
        //    "轨迹浮在工件外一大截"，而且没有任何参照物能让人发现是坐标系错了。
        //
        //  所以：手眼标定生效时，视图里**四类几何全部搬到基座系**——原始点云、算法轨迹、
        //  重建平面、修正轨迹与风险点。它们用同一个矩阵，所以彼此仍然严丝合缝，
        //  "原轨迹 vs 修正轨迹差多少"这个比较依然成立（刚体变换不改相对关系）。
        //
        //  手眼没启用时 DisplayPose()=null，全部原样输出，与旧行为逐像素一致。
        // ══════════════════════════════════════════════════════════════════

        /// <summary>显示要用的 T_cam_in_base；null = 保持点云原坐标系（未启用手眼）。</summary>
        double[]? DisplayPose() => ExportPose();

        /// <summary>当前视图显示在哪个坐标系里（给界面/日志用）。</summary>
        string DisplayFrameName() => DisplayPose() == null ? "点云原坐标系（相机系）" : "机器人基座系";

        /// <summary>把一个点搬到显示坐标系。M=null 时原样返回，不做任何浮点运算。</summary>
        (double x, double y, double z) DispP(double x, double y, double z)
        {
            var M = DisplayPose();
            if (M == null) return (x, y, z);
            return (M[0] * x + M[1] * y + M[2] * z + M[3],
                    M[4] * x + M[5] * y + M[6] * z + M[7],
                    M[8] * x + M[9] * y + M[10] * z + M[11]);
        }

        /// <summary>把一个方向（法向/枪轴/平面基向量）搬到显示坐标系：**只转不平移**。</summary>
        (double x, double y, double z) DispV(double x, double y, double z)
        {
            var M = DisplayPose();
            if (M == null) return (x, y, z);
            return (M[0] * x + M[1] * y + M[2] * z,
                    M[4] * x + M[5] * y + M[6] * z,
                    M[8] * x + M[9] * y + M[10] * z);
        }

        /// <summary>
        /// 把面片搬到显示坐标系。结构体里有 3 个大数组（顶点/三角形），
        /// 只复制**实际用到的**顶点段（polyVertCount），不做整块 1152 拷贝。
        /// </summary>
        SprayTrajPlane DispPlane(SprayTrajPlane p)
        {
            var M = DisplayPose();
            if (M == null) return p;

            var q = p;   // 先整体复制（数组字段随后替换）

            var c = DispP(p.cx, p.cy, p.cz);
            q.cx = c.x; q.cy = c.y; q.cz = c.z;

            // 法向：n 是"喷距法向"，n0 是"真实可喷侧法向"，两个都要转，
            // 否则面清单里「内外反」的判定会拿两套坐标系的向量做点积。
            var n = DispV(p.nx, p.ny, p.nz);
            q.nx = n.x; q.ny = n.y; q.nz = n.z;
            var n0 = DispV(p.n0x, p.n0y, p.n0z);
            q.n0x = n0.x; q.n0y = n0.y; q.n0z = n0.z;

            // 面内基向量：刚体变换后仍正交归一，但指向会变，必须一起转
            var u = DispV(p.ux, p.uy, p.uz);
            q.ux = u.x; q.uy = u.y; q.uz = u.z;
            var v = DispV(p.vx, p.vy, p.vz);
            q.vx = v.x; q.vy = v.y; q.vz = v.z;

            // 顶点：只转实际顶点的前 polyVertCount 个，剩下的保持原样（没人读）
            int nv = Math.Min(p.polyVertCount, p.polyVerts?.Length / 3 ?? 0);
            if (nv > 0)
            {
                var verts = new double[p.polyVerts!.Length];
                Array.Copy(p.polyVerts, verts, p.polyVerts.Length);
                for (int i = 0; i < nv; i++)
                {
                    var t = DispP(p.polyVerts[i * 3], p.polyVerts[i * 3 + 1], p.polyVerts[i * 3 + 2]);
                    verts[i * 3] = t.x; verts[i * 3 + 1] = t.y; verts[i * 3 + 2] = t.z;
                }
                q.polyVerts = verts;
            }
            return q;
        }

        // ══════════════════════════════════════════════════════════════════
        //  点云导出：手眼变换后的基座系点云
        //
        //  为什么必须导出：CSV 里全是**基座系**坐标，但如果只给两份 CSV，
        //  模拟平台/示教器那边没有工件可看，只能凭坐标想象"轨迹在工件前面
        //  100mm 处"长什么样；出问题时也没法把轨迹和工件放在一起比对。
        //  旧版导入说明写的是"点云 = 与本目录同名的 .ply"，可那个 .ply 是
        //  **相机系**的原始文件 —— 拿去和基座系 CSV 一起看，两者差一个
        //  1.5 米平移 + 90° 级旋转，看起来就是"轨迹飞在工件外面"。
        //
        //  口径：和 ExportTraj/ExportSafetyPts 用**同一个** ExportPose()，
        //  所以落盘的三份几何（点云/算法轨迹/修正轨迹）严格同系。
        //  手眼未启用时不写这个文件 —— 那时点云本来就在相机系，
        //  再复制一份改名只会让人分不清哪份才是原文件。
        // ══════════════════════════════════════════════════════════════════

        /// <summary>
        /// 写出「手眼变换后的点云」为 ASCII PLY（基座系）。手眼未生效时返回 null。
        ///
        /// 格式与 robot_safety/cloud.py 的 write_ply_xyz 一致（ascii 1.0 + 三位小数），
        /// 这样 Python 侧、C++ 侧、模拟平台三方都能直接读，不需要为它写特例。
        /// </summary>
        string? ExportCloudPly(string dir, string bas)
        {
            var M = ExportPose();
            if (M == null) return null;              // 相机系：不写，原 .ply 就是它
            if (_rawCloud == null || _rawCloud.Count == 0) return null;

            string path = System.IO.Path.Combine(dir, bas + "_base.ply");
            var ci = System.Globalization.CultureInfo.InvariantCulture;

            // 百万点用 StringBuilder 拼会吃掉几十 MB，直接流式写。
            // 缓冲 1MB：60 万点约 30MB，逐行 Write 会被格式化开销拖到几秒。
            using (var fs = new System.IO.FileStream(path, System.IO.FileMode.Create,
                                                     System.IO.FileAccess.Write, System.IO.FileShare.None,
                                                     1 << 20))
            using (var sw = new System.IO.StreamWriter(fs, new System.Text.UTF8Encoding(false)))
            {
                sw.NewLine = "\n";
                sw.Write("ply\n");
                sw.Write("format ascii 1.0\n");
                // comment 里写明坐标系与来源：纯数字的 ply 事后没法自证是哪一步的产物
                sw.Write("comment SprayTrajViewer export: hand-eye transformed to robot base frame\n");
                sw.Write("comment T_cam_in_base applied; same matrix as *_traj.csv / *_safe.csv\n");
                sw.Write($"element vertex {_rawCloud.Count}\n");
                sw.Write("property float x\nproperty float y\nproperty float z\n");
                sw.Write("end_header\n");

                for (int i = 0; i < _rawCloud.Count; i++)
                {
                    var p = _rawCloud[i];
                    double x = M[0] * p.x + M[1] * p.y + M[2] * p.z + M[3];
                    double y = M[4] * p.x + M[5] * p.y + M[6] * p.z + M[7];
                    double z = M[8] * p.x + M[9] * p.y + M[10] * p.z + M[11];
                    sw.Write(x.ToString("F3", ci)); sw.Write(' ');
                    sw.Write(y.ToString("F3", ci)); sw.Write(' ');
                    sw.Write(z.ToString("F3", ci)); sw.Write('\n');
                }
            }

            App.Log("PLYOUT", null,
                $"{System.IO.Path.GetFileName(path)}  {_rawCloud.Count} 点（基座系，含可达位平移）");
            return path;
        }

        /// <summary>
        /// 显示坐标系变了之后，把视图里所有几何按新口径重画一遍。
        ///
        /// 为什么必须有这一步：DispP/DispV 是在**渲染的那一刻**才取的矩阵，
        /// 已经加进 Viewport 的 Visual3D 存的是烘死的 Point3D。只改 _handEyeExportPose
        /// 而不重画，屏幕上一个像素都不会动 —— 用户会以为"勾了手眼根本没反应"。
        ///
        /// 只重画**数据还在内存里**的那几层：
        ///   · 点云（_cloud）：永远有；
        ///   · 平面（_planes）：生成过才有；
        ///   · 轨迹（_lastTraj）/ 修正层（_safetyPts）：只有生成过才有。
        /// 轨迹不在这里重算 —— 换了坐标系说明口径变了，旧轨迹本来就该作废
        /// （调用方 SwitchToHandEyeConventions 会先 InvalidateTrajectory）。
        /// 这里只负责把"还该显示的东西"搬到新坐标系，不制造结果。
        /// </summary>
        void RefreshDisplayFrame(string reason)
        {
            if (Viewport == null) return;   // XAML 解析阶段可能还没建好

            // ---------- 步骤 1：点云 ----------
            // _cloud 非空才画。这里**不传 keepView** —— 换坐标系意味着整个工件
            // 平移了约 1.5 米、还带 90° 级旋转，保住旧镜头只会让人看到一片空白
            // （相机还对着原来的坐标区间）。ZoomExtents 一次最省事。
            if (_cloud.Count > 0) RenderCloud();

            // ---------- 步骤 2：重建平面 ----------
            if (_planes.Count > 0) RenderPlane();

            // ---------- 步骤 3：算法轨迹 ----------
            // 换了坐标系后旧轨迹通常已被作废（_lastTraj 为空），那就不用画；
            // 如果还在（比如只是改了 faceX 又没触发作废），照样按新系重画。
            if (_lastTraj.Count > 1) RenderTrajectory(_lastTraj);

            // ---------- 步骤 4：修正轨迹 + 风险点 ----------
            if (_safetyPts.Length > 0)
                RenderSafetyOverlay(_safetyPts, ChkShowFixed.IsChecked != false);

            // ---------- 步骤 5：取景 ----------
            // 换坐标系后工件已经挪到别处，ZoomExtents 给的默认方向往往是"从下往上看"，
            // 落地看工件看不出所以然。有面片就正视面0（就是喷的那一面），
            // 没有面片才退回 ZoomExtents 的通用取景。
            if (_planes.Count > 0) ApplyViewPreset("face");

            // ---------- 步骤 6：坐标系标签 ----------
            // 常驻显示，避免"看着像算法错了其实是坐标系没切"这类误判。
            TxtFrameTag.Text = $"视图坐标系：{DisplayFrameName()}";
            BdFrameTag.Background = new SolidColorBrush(DisplayPose() == null
                ? Color.FromArgb(0xCC, 0x37, 0x41, 0x51)      // 相机系：灰蓝
                : Color.FromArgb(0xCC, 0x06, 0x5F, 0x46));    // 基座系：墨绿，一眼区分

            App.Log("DISPFRAME", null, $"显示坐标系 → {DisplayFrameName()}（{reason}）");
        }

        // ══════════════════════════════════════════════════════════════════
        //  自动导出：<点云目录>\traj_out\<点云名>_traj.csv / _safe.csv / .plan.json
        //
        //  为什么要自动：修正后的轨迹以前只在内存里（要手点导出、还要自己选目录），
        //  换个窗口就没了 —— 这不是"用户忘点"，是流程设计缺一步。
        //
        //  为什么随点云走：点云是用户自己挑的、天然知道放在哪个工程目录里；
        //  写死在 exe 目录会把不同工件的结果混在一个 bin 文件夹里，事后分不清谁是谁。
        // ══════════════════════════════════════════════════════════════════

        /// <summary>
        /// 自动导出目录 = 点云所在目录下的 traj_out。
        /// 点云目录写不进去时（典型原因：安全软件的"文档保护"拦桌面/文档，点云在只读盘）
        /// 自动退到程序目录下的 traj_out —— 先真实写一个探针文件判断，不能等 AutoExport 抛异常。
        /// 结果缓存，换点云时由 LoadCloud 清空。
        /// </summary>
        string AutoExportDir()
        {
            if (_exportDirResolved is not null) return _exportDirResolved;

            // 1) 界面上手动指定的目录最优先（选一次就记住，换点云/重启都沿用）
            string? custom = UiState.Get("exportDir");
            if (!string.IsNullOrWhiteSpace(custom) && CanWriteTo(custom))
            {
                _exportDirResolved = custom;
                return custom;
            }

            // 2) 默认放程序目录\traj_out，**不再跟随点云目录**：
            //    桌面/文档经常被安全软件的"文档保护"拦截，写那里会静默失败
            //    （2026-10-28 用户就是因此以为"轨迹没生成"）。程序目录在 E: 盘，
            //    未签名 exe 也能写，位置固定、好找。想跟点云走就点「更改…」自己选。
            string fallback = System.IO.Path.Combine(AppContext.BaseDirectory, "traj_out");
            _exportDirResolved = fallback;
            return fallback;
        }

        /// <summary>真实写一个小文件判断目录可写：ACL 和安全软件拦截都能测出来。</summary>
        static bool CanWriteTo(string dir)
        {
            try
            {
                System.IO.Directory.CreateDirectory(dir);
                string probe = System.IO.Path.Combine(dir, ".write_probe_" + Guid.NewGuid().ToString("N") + ".tmp");
                System.IO.File.WriteAllText(probe, "");
                System.IO.File.Delete(probe);
                return true;
            }
            catch
            {
                return false;
            }
        }

        /// <summary>
        /// 把这次的导出目录写到界面上（勾选框下面那行）。
        /// 输出位置是按**点云所在目录**算的：桌面\点云\喷涂点云 和 喷涂点云_DENSO可达 里各有一份
        /// 同名的 垃圾桶(模拟生成).ply，生成完就在两处各落一个 traj_out —— 只看文件名永远说不清
        /// 刚才那一份在哪。所以这行不是"锦上添花"，是这条流程唯一能自证位置的东西。
        /// </summary>
        void ShowOutDir()
        {
            if (TxtOutDir == null) return;
            if (string.IsNullOrEmpty(_currentFile))
            {
                TxtOutDir.Text = "导出目录：加载点云后显示";
                TxtOutDir.ToolTip = null;
                return;
            }
            string dir = AutoExportDir();
            bool custom = !string.IsNullOrWhiteSpace(UiState.Get("exportDir"));
            TxtOutDir.Text = "导出目录：" + dir + (custom ? "（手动指定）" : "（默认：程序目录）");
            // _base.ply 只在手眼生效时才写，所以按当前状态决定要不要列出来
            bool hasBase = ExportPose() != null;
            TxtOutDir.ToolTip = $"{dir}\n\n本次会写：{CloudBaseName()}_traj.csv、"
                              + $"{CloudBaseName()}_safe.csv、{CloudBaseName()}.plan.json"
                              + (hasBase ? $"、{CloudBaseName()}_base.ply（基座系点云）" : "");
            App.Log("OUTDIR", null, dir);
        }

        void BtnOpenOutDir_Click(object sender, RoutedEventArgs e)
        {
            if (string.IsNullOrEmpty(_currentFile))
            {
                TxtHint.Text = "先加载点云，才知道要导出到哪个目录";
                return;
            }
            try
            {
                string dir = AutoExportDir();
                System.IO.Directory.CreateDirectory(dir);   // 还没生成过也让它能打开（空目录）
                System.Diagnostics.Process.Start(new System.Diagnostics.ProcessStartInfo
                {
                    FileName = "explorer.exe",
                    Arguments = "\"" + dir + "\"",
                    UseShellExecute = true
                });
                App.Log("OUTDIR", null, "打开 " + dir);
            }
            catch (Exception ex)
            {
                App.Log("OUTDIR/FAIL", ex);
                TxtHint.Text = "打开目录失败：" + ex.Message;
            }
        }

        /// <summary>
        /// 手动指定自动导出目录（安全软件拦桌面时改到可写盘，或任何你想放的工程目录）。
        /// 选择存 UiState，换点云/重启都沿用；改回默认就选程序目录下的 traj_out。
        /// </summary>
        void BtnSetOutDir_Click(object sender, RoutedEventArgs e)
        {
            var dlg = new OpenFolderDialog
            {
                Title = "选择自动导出目录",
                InitialDirectory = System.IO.Directory.Exists(AutoExportDir())
                    ? AutoExportDir()
                    : AppContext.BaseDirectory,
            };
            if (dlg.ShowDialog() != true) return;

            UiState.Set("exportDir", dlg.FolderName);
            _exportDirResolved = null;
            ShowOutDir();
            TxtHint.Text = "导出目录已改为：" + dlg.FolderName;
            App.Log("OUTDIR", null, "用户指定 " + dlg.FolderName);
        }

        /// <summary>导出用的基名 = 点云文件名（去扩展名）。</summary>
        string CloudBaseName()
        {
            string n = System.IO.Path.GetFileNameWithoutExtension(_currentFile ?? "");
            return string.IsNullOrWhiteSpace(n) ? "trajectory" : n;
        }

        /// <summary>
        /// 自动导出。<paramref name="withSafety"/>=true 时额外写修正后轨迹，
        /// 并把 plan.json 换成**修正后**的位姿（模拟平台要跑的就是这一份）。
        /// </summary>
        void AutoExport(bool withSafety)
        {
            if (ChkAutoExport.IsChecked != true) return;
            if (string.IsNullOrEmpty(_currentFile) || _lastTraj.Count == 0) return;

            try
            {
                string dir = AutoExportDir();
                System.IO.Directory.CreateDirectory(dir);
                string bas = CloudBaseName();
                var wrote = new List<string>();

                string fTraj = System.IO.Path.Combine(dir, bas + "_traj.csv");
                WriteTrajCsv(fTraj, ExportTraj());
                wrote.Add(System.IO.Path.GetFileName(fTraj));

                // 基座系点云：手眼生效时必须有，否则两份基座系 CSV 旁边没有工件可对照
                // （导入说明里"点云 = 同名的 .ply"指的是相机系原文件，坐标系对不上）。
                // 手眼未启用时返回 null，不写 —— 那时原 .ply 本身就是正确的参照。
                string? fPly = ExportCloudPly(dir, bas);
                if (fPly != null) wrote.Add(System.IO.Path.GetFileName(fPly));
                else
                {
                    // 手眼关了但上一次的 _base.ply 还留着 → 那份是旧矩阵的产物，
                    // 留在同一目录里会被误认为本次结果（和 _safe.csv 同一个坑）。
                    string stale = System.IO.Path.Combine(dir, bas + "_base.ply");
                    if (System.IO.File.Exists(stale))
                    {
                        try { System.IO.File.Delete(stale); } catch { }
                    }
                }

                // 重新生成轨迹时，上一次的修正结果已经作废（它对应的是另一条轨迹）。
                // 不清掉会留下"新 _traj.csv + 旧 _safe.csv"这种看不出错、但会把旧姿态
                // 送进模拟平台的组合。
                string fSafeOld = System.IO.Path.Combine(dir, bas + "_safe.csv");
                if (System.IO.File.Exists(fSafeOld))
                {
                    try { System.IO.File.Delete(fSafeOld); } catch { }
                }

                // ---------- 两份 plan.json：算法轨迹 / 修正轨迹，并行存在 ----------
                // 以前是同一个文件名靠 useSafe 覆盖，于是"生成 → 评估"跑完后只剩修正版，
                // 算法原始轨迹那份 json 就没了；想对比"安全层到底改了什么"只能翻 CSV。
                // 现在固定拆成两个名字，两份都能直接导进 DensoArmKine 对比着跑。
                string fJsonRaw = System.IO.Path.Combine(dir, bas + "_traj.plan.json");
                WriteDensoPlanJson(fJsonRaw, bas + "_traj", false);
                wrote.Add(System.IO.Path.GetFileName(fJsonRaw));

                bool useSafe = withSafety && _safetyPts.Length > 0;
                if (useSafe)
                {
                    string fSafe = System.IO.Path.Combine(dir, bas + "_safe.csv");
                    WriteFixedCsv(fSafe);
                    wrote.Add(System.IO.Path.GetFileName(fSafe));

                    string fJsonSafe = System.IO.Path.Combine(dir, bas + "_safe.plan.json");
                    WriteDensoPlanJson(fJsonSafe, bas + "_safe", true);
                    wrote.Add(System.IO.Path.GetFileName(fJsonSafe));
                }
                else
                {
                    // 没有修正结果时删掉上一轮的 _safe.plan.json，否则 DensoArmKine 里
                    // 会看到一份"属于另一条轨迹"的修正版，和旧 _safe.csv 同一个坑。
                    string stale = System.IO.Path.Combine(dir, bas + "_safe.plan.json");
                    if (System.IO.File.Exists(stale)) { try { System.IO.File.Delete(stale); } catch { } }
                }

                EnsureExportReadme(dir);

                string msg = $"自动导出 → {dir}（{string.Join("、", wrote)}）"
                           + (useSafe ? " · 修正后轨迹" : " · 未修正（算法轨迹）");
                TxtHint.Text = TxtHint.Text.Split('\n')[0] + "\n" + msg;
                if (TxtOutDir != null) TxtOutDir.Text = "导出目录：" + dir + "（已写入）";
                App.Log("AUTOEXPORT", null, msg);
                if (!useSafe && ChkAutoSafety.IsChecked == true)
                    App.Log("AUTOEXPORT", null,
                        "警告：勾了「自动评估」但本次没有可用的修正结果，json 写的是算法轨迹（Meta.Source=raw）");
            }
            catch (Exception ex)
            {
                App.Log("AUTOEXPORT/ERR", ex, null);
                TxtHint.Text += "\n自动导出失败：" + ex.Message;
            }
        }

        /// <summary>首次自动导出时写一份导入说明（只写一次，不覆盖用户改过的）。</summary>
        void EnsureExportReadme(string dir)
        {
            string path = System.IO.Path.Combine(dir, "导入说明.txt");
            if (System.IO.File.Exists(path))
            {
                // 旧版说明有三处已过时的口径：
                //   1) "点云 = 与本目录同名的 .ply（本工具不改点云）"——手眼生效时那是相机系，差 1.5 米；
                //   2) 只提一份 *.plan.json —— 现在拆成 _traj / _safe 两份；
                //   3) 说违规点会被删掉 —— 现在一点不删。
                // 说明只写一次、不会自己更新，留着比没有更糟。只替换旧版本，用户改过的不动。
                try
                {
                    string old = System.IO.File.ReadAllText(path);
                    bool stale = !old.Contains("_traj.plan.json") || !old.Contains("FlaggedPoints");
                    if (!stale) return;
                    App.Log("README", null, "检测到旧版导入说明，已更新（两份 plan.json + 不拦策略）");
                }
                catch { return; }   // 读不出来就别动它
            }
            System.IO.File.WriteAllText(path,
                "本目录由 SprayTrajViewer 自动生成，文件名取自点云名。\r\n\r\n" +
                "【可导入文件】\r\n" +
                "  *_base.ply         基座系点云（手眼矩阵变换 + 挪到可达位）\r\n" +
                "                     仅在启用手眼标定时生成，是下面几份轨迹的配套工件\r\n" +
                "  *_traj.plan.json   算法原始轨迹 → DensoArmKine 可直接导入\r\n" +
                "  *_safe.plan.json   修正后轨迹   → DensoArmKine 可直接导入\r\n" +
                "                     没跑安全评估时这一份不存在\r\n" +
                "  *_traj.csv         算法原始轨迹（X,Y,Z,I,J,K,speed,spray）\r\n" +
                "  *_safe.csv         修正后轨迹（含 tilt/clearance/关节角/违规位）\r\n" +
                "                     CSV 需要 traj_to_denso_plan.py 转换后才能导入\r\n\r\n" +
                "【安全层违规不删点（report-only）】\r\n" +
                "  两份 plan.json 与对应 CSV **点数完全一致、几何原样**。\r\n" +
                "  GunOpen 严格等于 CSV 的 spray 列，不因违规关枪。\r\n" +
                "  IsValid 恒为 true —— 与 traj_to_denso_plan.py 生成的口径一致。\r\n" +
                "  违规信息只作统计与诊断，落在：\r\n" +
                "    · Meta.FlaggedPoints（被判违规的点数，这些点仍在 Points 里）；\r\n" +
                "    · Meta.BlockPolicy = \"report-only\"；\r\n" +
                "    · *_safe.csv 的 violation 列（逐点位掩码）。\r\n" +
                "  这样做是因为违规**大概率是判定代码的问题**而不是真撞件：\r\n" +
                "  先把完整几何交给模拟平台跑出来看，再回头定位是谁算错了。\r\n\r\n" +
                "【坐标系】\r\n" +
                "  启用手眼标定 → 点云与两份轨迹全在机器人基座系，可直接比对。\r\n" +
                "  未启用手眼   → 两份轨迹在点云采集的相机系，此时没有 _base.ply，\r\n" +
                "                 要对照就用原始那份点云（同一个系）。\r\n" +
                "  别把原始点云 .ply 和基座系 CSV 放在一起看：两者差约 1.5 米平移\r\n" +
                "  加 90° 级旋转，看起来像\"轨迹飞到工件外面\"，其实是坐标系没对上。\r\n" +
                "  每份 plan.json 的 Meta 都有 Frame（base/cloud）与 CloudFile（配套点云名）。\r\n\r\n" +
                "【DensoArmKine 加载口径（默认值正好对，别改）】\r\n" +
                "  轨迹页 → 轨迹口径 = TCP\r\n" +
                "  坐标系   = 基座系\r\n" +
                "  姿态翻转 = 不翻\r\n" +
                "  点云     = 本目录下的 *_base.ply（有就优先用它，别用原始点云）\r\n\r\n" +
                "【两个常见误判】\r\n" +
                "  1) 报「点云不在轨迹位置、最近 ~100mm」= 正常，那就是喷距。\r\n" +
                "  2) --traj-cloud 报一排 NEAR = 工具长口径不同（本仓默认 150mm，\r\n" +
                "     DensoArmKine settings.json 里 VS-087 是 453.08mm）。\r\n",
                new System.Text.UTF8Encoding(true));
        }

        /// <summary>
        /// Meta 里的倾角口径：实际生效的系数（δ=90° 处给多少度）与封顶值。
        /// 同一件用 k=30 / k=90 导出的 json，点数与 GunOpen 逐点相同，**只有枪的朝向差几十度**，
        /// 事后想判断"这份该不该重导"只能靠这两个数。还没评估过时返回空串（绝不能写 NaN 把 JSON 弄坏）。
        /// </summary>
        string TiltMetaText()
        {
            var ci = System.Globalization.CultureInfo.InvariantCulture;
            double k = double.IsNaN(_effTiltK) ? _tiltCoefDeg : _effTiltK;
            double mx = double.IsNaN(_effTiltMaxDeg) ? ParseDouble(TxtTiltMax.Text, double.NaN) : _effTiltMaxDeg;
            var s = new System.Text.StringBuilder();
            if (double.IsFinite(k)) s.Append("\"TiltK\": ").Append(k.ToString("F0", ci)).Append(", ");
            if (double.IsFinite(mx)) s.Append("\"TiltMax\": ").Append(mx.ToString("F0", ci)).Append(", ");
            return s.ToString();
        }

        /// <summary>
        /// 写 DensoArmKine 的 spray_plan JSON。
        ///
        /// 口径（已在 DensoArmKine/Trajectory/SprayPlan.cs 与 Python 版
        /// traj_to_denso_plan.py 两侧核对过，别凭感觉改）：
        ///   · Cartesian = TCP 位姿原样（本仓 CSV 的 X,Y,Z 就是 TCP，**不要再退一个喷距**）；
        ///   · 工具 +Z = (I,J,K)；工具 +X 取"上一点 X 投影到新切平面"，保证相邻姿态最小旋转
        ///     （不做这件事，枪轴接近竖直时参考向量会在 +X/+Z 之间跳，逆解跳到另一支 → J6 自旋）；
        ///   · 欧拉 = Rz·Ry·Rx 分解，单位**度**，Fig=0。
        /// </summary>
        void WriteDensoPlanJson(string path, string name, bool useSafe)
        {
            var ci = System.Globalization.CultureInfo.InvariantCulture;
            // 落盘前搬到基座系：DensoArmKine 的 Cartesian 是基座系位姿。
            // 手眼未启用时 Export* 返回原对象（不拷贝），行为与旧版逐字节一致。
            var trajBase = ExportTraj();
            var safeBase = ExportSafetyPts();
            int n = useSafe ? safeBase.Length : trajBase.Count;
            // 违规标记：useSafe 时按安全层结果打标，但**不删点**。
            // （2026-11-10 用户要求"不安全先不拦"：违规大概率是判定代码的问题，
            //  需要把几何原样交给 DensoArmKine 看清楚，再回头定位是谁算错了。）
            bool[]? ok = useSafe ? BuildExecutableMask(safeBase) : null;
            int nDrop = 0;      // 判了违规的**喷涂**点数（几何仍保留，只计个数）
            if (ok != null) for (int i = 0; i < n; ++i) if (safeBase[i].spray != 0 && !ok[i]) ++nDrop;
            int nSpray = 0;     // 会开枪的喷涂点数（只看 spray 位，**不看**违规）
            for (int i = 0; i < n; i++)
                if ((useSafe ? safeBase[i].spray : trajBase[i].spray) != 0) nSpray++;

            var sb = new System.Text.StringBuilder(n * 190 + 400);
            sb.Append("{\n  \"Name\": \"").Append(JsonEsc(name)).Append("\",\n");
            // Meta：这份 json 到底是**哪一份**轨迹。原来只能靠"旁边有没有 _safe.csv"去猜，
            // 垃圾桶就踩过——重新生成后安全评估没跑成，json 静默留在算法轨迹上，
            // 模拟平台里看到的还是没修正的枪向，还以为是修正算错了。
            // DensoArmKine 的 SprayPlan.Load 用 JsonDocument.TryGetProperty 逐字段取，
            // 顶层多余字段直接忽略，加 Meta 不会影响它加载。
            sb.Append("  \"Meta\": {")
              .Append("\"Source\": \"").Append(useSafe ? "safe" : "raw").Append("\", ")
              .Append("\"SourceText\": \"").Append(useSafe ? "修正后轨迹（安全层输出）" : "算法原始轨迹（未修正）").Append("\", ")
              .Append("\"TiltMode\": \"").Append(ChkTiltPerPoint.IsChecked == true ? "per-point" : "face").Append("\", ")
              .Append("\"Pose\": \"").Append(ChkFastPose.IsChecked == true ? "onepass" : "search").Append("\", ")
              // 倾角口径的三个数必须落盘：同一件用 k=30 和 k=90 导出的 json，点数/GunOpen 完全一致，
              // 只有枪的朝向差几十度 —— 不写这三个数，事后就没法判断该不该重导。
              // 还没评估过（生成即导出的 raw）时两者都是 NaN，此时**不写**，不能吐 NaN 把 JSON 弄坏。
              .Append(TiltMetaText())
              // Points = json 里实际有的几何点数（**含**被判违规的）。
              // 以前这里是 n-nDrop，于是"全拦"时写成 0、Points 数组也空，
              // 打开文件看不出轨迹在哪 —— 排查时最需要的就是这条几何。
              .Append("\"Points\": ").Append(n).Append(", ")
              // SprayPoints = GunOpen=true 的点数，与 CSV 的 spray 列严格一致
              .Append("\"SprayPoints\": ").Append(nSpray).Append(", ")
              // FlaggedPoints = 安全层判违规的点数。**这些点仍然在 Points 里**，
              // 只是统计用。名字刻意不叫 DropPoints：那是旧语义（点被删了），
              // 现在一点都没删，用 Drop 会让人以为 1495 里少了 508 个，反而更糊涂。
              .Append("\"FlaggedPoints\": ").Append(nDrop).Append(", ")
              // 明确写出策略，避免下游猜"到底拦没拦"
              .Append("\"BlockPolicy\": \"report-only\", ")
              .Append("\"Verdict\": \"").Append(JsonEsc(useSafe ? TxtSafetyVerdict.Text : "未评估")).Append("\", ")
              // 坐标系与配套点云：CSV 里的数字是基座系还是相机系，事后从数字本身**看不出来**
              // （两者都是一千多毫米量级的正数）。不写这两个键，拿到一份 json 只能靠
              // "旁边有没有 _base.ply"去猜，猜错就是整套结论反过来。
              // CloudFile 指向与本 json 同目录的基座系点云；未启用标定时为空串。
              .Append("\"Frame\": \"").Append(ExportPose() == null ? "cloud" : "base").Append("\", ")
              // CloudFile 要指向**点云**名（不含 _traj / _safe 后缀）：
              // 这个函数现在被调两次（算法/修正各一份），name 分别带这两个后缀，
              // 直接拼 name+"_base.ply" 会得到不存在的文件名。
              .Append("\"CloudFile\": \"")
              .Append(JsonEsc(ExportPose() == null ? "" : CloudBaseName() + "_base.ply")).Append("\", ")
              .Append("\"GeneratedAt\": \"").Append(DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss")).Append("\"},\n");
            sb.Append("  \"DefaultSpeedPercent\": 30,\n  \"Points\": [\n");

            double px = 0, py = 0, pz = 0;
            bool hasPrevX = false;
            var rows = new List<string>(n);
            int nGunOn = 0;              // 实际开枪点数（供日志统计，写进 Meta 的是 nSpray）
            for (int i = 0; i < n; i++)
            {
                double x, y, z, ax, ay, az, speed;
                int spray;
                bool gunOn = true;
                if (useSafe)
                {
                    var q = safeBase[i];
                    // 2026-11-10 用户要求「不安全先不拦」：违规点**原样写进去**，
                    // GunOpen 也按修正层原本的开关枪意图给，**不因违规关枪**。
                    // 理由：违规大概率是判定代码的问题而不是真撞件，先把完整几何
                    // 交给 DensoArmKine 跑出来看，再回头定位是谁算错了。
                    // 违规信息不丢，落在 *_safe.csv 的 violation 列和 Meta.FlaggedPoints 计数里。
                    x = q.x; y = q.y; z = q.z; ax = q.i; ay = q.j; az = q.k;
                    speed = q.speed; spray = q.spray;
                    gunOn = spray != 0;
                }
                else
                {
                    var q = trajBase[i];
                    x = q.x; y = q.y; z = q.z; ax = q.i; ay = q.j; az = q.k;
                    speed = q.speed; spray = q.spray;
                    gunOn = spray != 0;
                }
                if (gunOn) ++nGunOn;

                // 枪轴归一化（CSV 的 IJK 已经是单位向量，这里只是防脏数据）
                double al = Math.Sqrt(ax * ax + ay * ay + az * az);
                if (al < 1e-12) { ax = 0; ay = 0; az = -1; al = 1; }
                ax /= al; ay /= al; az /= al;

                // 工具 +X：上一点 X 投影到新切平面；首点用"与 +Z 叉乘"起头
                double xx, xy, xz;
                if (!hasPrevX)
                {
                    double tx = 0, ty = 0, tz = 0;
                    if (Math.Abs(az) < 0.9) { tz = 1; } else { tx = 1; }
                    xx = ty * az - tz * ay; xy = tz * ax - tx * az; xz = tx * ay - ty * ax;
                    double l = Math.Sqrt(xx * xx + xy * xy + xz * xz);
                    if (l < 1e-12) { xx = 1; xy = 0; xz = 0; } else { xx /= l; xy /= l; xz /= l; }
                }
                else
                {
                    double d = px * ax + py * ay + pz * az;
                    xx = px - d * ax; xy = py - d * ay; xz = pz - d * az;
                    double l = Math.Sqrt(xx * xx + xy * xy + xz * xz);
                    // 退化阈值 1e-3（≈0.06°）：上一点 X 与新枪轴几乎平行时（盒子这类件在
                    // 面边界处枪轴一次转 90° 就会这样），投影向量的方向完全由浮点噪声决定，
                    // 不同实现/不同编译器会给出相反的滚转 → J6 在相邻两点间翻 180°。
                    // 这种点必须走确定性参考，不能拿噪声当方向。
                    if (l < 1e-3)
                    {
                        double tx = 0, ty = 0, tz = 0;
                        if (Math.Abs(az) < 0.9) { tz = 1; } else { tx = 1; }
                        xx = ty * az - tz * ay; xy = tz * ax - tx * az; xz = tx * ay - ty * ax;
                        l = Math.Sqrt(xx * xx + xy * xy + xz * xz);
                    }
                    if (l < 1e-12) { xx = 1; xy = 0; xz = 0; } else { xx /= l; xy /= l; xz /= l; }
                }
                // Y = Z × X，再正交化 X = Y × Z
                double yx = ay * xz - az * xy, yy = az * xx - ax * xz, yz = ax * xy - ay * xx;
                {
                    double l = Math.Sqrt(yx * yx + yy * yy + yz * yz);
                    if (l < 1e-12) { yx = 0; yy = 1; yz = 0; } else { yx /= l; yy /= l; yz /= l; }
                }
                {
                    double nx = yy * az - yz * ay, ny = yz * ax - yx * az, nz = yx * ay - yy * ax;
                    double l = Math.Sqrt(nx * nx + ny * ny + nz * nz);
                    if (l > 1e-12) { xx = nx / l; xy = ny / l; xz = nz / l; }
                }
                px = xx; py = xy; pz = xz; hasPrevX = true;

                // R = [X Y Z]（列），按 Rz·Ry·Rx 分解。
                // 注意 r20 = **R[2,0] = 工具 X 轴的 z 分量**（不是枪轴 Z 的 x 分量）：
                // 列向量矩阵里 R[2,0] 是第 0 列的 z 分量。写成 ax 会让整条轨迹的姿态偏几十度，
                // 而且位置完全正确 —— 典型的"看不出错、模拟器里枪歪着喷"。
                double r20 = Math.Max(-1.0, Math.Min(1.0, xz));
                double ry = -Math.Asin(r20);
                double cy = Math.Cos(ry);
                double rx, rz;
                if (Math.Abs(cy) < 1e-9)
                {
                    rx = Math.Atan2(-ay, yy);   // R[1,2]=Z.y=ay, R[1,1]=Y.y=yy
                    rz = 0;
                }
                else
                {
                    rz = Math.Atan2(xy, xx);    // R[1,0]=X.y, R[0,0]=X.x
                    rx = Math.Atan2(yz, az);    // R[2,1]=Y.z, R[2,2]=Z.z
                }
                rx = Norm180(rx * 180.0 / Math.PI);
                ry = Norm180(ry * 180.0 / Math.PI);
                rz = Norm180(rz * 180.0 / Math.PI);

                // 逐行收集再 join，不能用 "i == n-1" 决定逗号。
                //
                // 字段格式与 traj_to_denso_plan.py 完全一致（Joints/Cartesian/SpeedPercent/
                // GunOpen/Label/IsValid），不多不少。
                // 用户要求「不安全先不拦」，所以这里**不引入**任何新字段（Risk 之类）：
                // 多了字段下游不认识，反而让人怀疑文件是不是被改坏了。
                // IsValid 恒 true —— 和生成器口径一致，违规信息只留在 *_safe.csv 里诊断。
                rows.Add("    { \"Joints\": null, \"Cartesian\": ["
                  + F6(x, ci) + ", " + F6(y, ci) + ", " + F6(z, ci) + ", "
                  + F6(rx, ci) + ", " + F6(ry, ci) + ", " + F6(rz, ci) + ", 0], "
                  + "\"SpeedPercent\": 30, "
                  + "\"GunOpen\": " + (gunOn ? "true" : "false") + ", "
                  + "\"Label\": \"v" + (int)Math.Round(speed, MidpointRounding.AwayFromZero) + "\", "
                  + "\"IsValid\": true }");
            }

            sb.Append(string.Join(",\n", rows)).Append("\n  ]\n}\n");
            System.IO.File.WriteAllText(path, sb.ToString(), new System.Text.UTF8Encoding(false));
        }

        static double Norm180(double deg)
        {
            double d = (deg + 180.0) % 360.0;
            if (d < 0) d += 360.0;
            return d - 180.0;
        }

        static string F6(double v, System.Globalization.CultureInfo ci)
        {
            if (double.IsNaN(v) || double.IsInfinity(v)) v = 0;
            return Math.Round(v, 6).ToString("0.######", ci);
        }

        static string JsonEsc(string s) =>
            (s ?? "").Replace("\\", "\\\\").Replace("\"", "\\\"");

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
        // 修正后轨迹分两条画：喷涂青、空移紫。LinesVisual3D 只有一个 Color，
        // 合成一条的话空移段就会跟喷涂段同色，"哪一段是空走"在修正层就看不出来了。
        private LinesVisual3D? _fixedSprayVisual;
        private LinesVisual3D? _fixedMoveVisual;
        // 被安全层过滤掉的喷涂线：不上绿（不执行），单独画成暗灰，
        // 让人一眼看出"哪几条线不喷"，而不是在绿线里凭空少一截还找不到原因。
        private LinesVisual3D? _fixedDroppedVisual;
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

        void CmbUpDir_Changed(object sender, System.Windows.Controls.SelectionChangedEventArgs e)
        {
            if (PanelUpCustom == null) return;
            PanelUpCustom.Visibility = CmbUpDir.SelectedIndex == 4
                ? Visibility.Visible : Visibility.Collapsed;
        }

        /// <summary>
        /// 界面「重力向上」→ 点云系里的单位向量。返回 null 表示用库默认（+Z）。
        ///
        /// 这件事本身不难，难在**不能想当然**：倾角口径里"朝上"是相对重力定义的，
        /// 而点云是相机系，+Z 是光轴。所以这里做成显式选项而不是自动从手眼矩阵推 ——
        /// 实测煤气灶背板：手眼 R 给出的 up 会让主面偏离竖直 16.5°，而按「立着挂」
        /// 工艺事实取 ±Y 只有 3.7°（立式背板 860 高 × 400 宽）。矩阵是手工录的，
        /// 工艺摆姿比它可信。
        /// </summary>
        double[]? BuildUpVector()
        {
            switch (CmbUpDir.SelectedIndex)
            {
                case 1: return new double[] { 0, -1, 0 };
                case 2: return new double[] { 0,  1, 0 };
                case 3: return new double[] { 0,  0, -1 };
                case 4:
                    double ux = ParseDouble(TxtUpX.Text, 0);
                    double uy = ParseDouble(TxtUpY.Text, 0);
                    double uz = ParseDouble(TxtUpZ.Text, 1);
                    // 全零 = 没填 ⇒ 交回库默认，别把"朝某个方向"当成有效输入
                    if (Math.Abs(ux) < 1e-9 && Math.Abs(uy) < 1e-9 && Math.Abs(uz) < 1e-9) return null;
                    return new double[] { ux, uy, uz };
                default: return null;   // +Z：库默认，保持旧行为
            }
        }

        /// <summary>
        /// 口径一致性自检：点云已在基座系（恒等）时，重力向上必须是 +Z。
        ///
        /// 这两项是**一对**：手眼变换会把点云搬到基座系，同时把"下肢"的语义一起换掉。
        /// 只改了基座位姿、忘了把「重力向上」从 −Y 改回 +Z，是最容易发生也最难看出来的错 ——
        /// 点云位置全对、可达率也正常，只有倾角按一个歪掉 90° 的"向上"算，
        /// 表现为"明明立着挂的板，算法当成平放"。
        /// 不直接改用户的选项（也许他真的想这么算），而是明确告警。
        /// </summary>
        bool CheckUpDirConsistency()
        {
            if (CmbBaseMode.SelectedIndex != 2) return true;      // 仅恒等口径需要
            if (CmbUpDir.SelectedIndex == 0) return true;         // +Z，一致
            string msg =
                "「基座位姿」选了「点云已在基座系（恒等）」，但「重力向上」不是 +Z。\n\n" +
                "点云已经在机器人基座系里，基座 +Z 就是重力向上，倾角口径必须用 +Z。\n" +
                "现在这样算，点云位置对、可达率也正常，但倾角是相对一个歪掉的方向算的 ——\n" +
                "现象很像「明明立着挂的板，算法当成平放」。\n\n" +
                "要继续吗？（选「否」会先停下，建议把「重力向上」改回 +Z）";
            var r = MessageBox.Show(msg, "口径不一致", MessageBoxButton.YesNo, MessageBoxImage.Warning);
            return r == MessageBoxResult.Yes;
        }
        /// <summary>
        /// 点云 → 基座 的变换矩阵 T_base_in_cloud（4x4 行主序），喂给安全层与导出。
        ///
        /// 优先级：
        ///   1. 手眼标定生效（_handEyeBasePose 非空）→ 直接用标定矩阵 + 可达位平移。
        ///      这是唯一由**标定**给出的真值，界面上的 x/y/z 输入框无法表达它的任意旋转。
        ///   2. 「点云已在基座系（恒等）」→ 单位阵。注意这与返回 null **含义完全不同**：
        ///      null 会让核心库按质心再猜一次原点（baseAuto=1），单位阵才是"别动，已经对齐了"。
        ///      同一份点云该走哪条搞错了代价很大：实测 P0 从 368 变 2611、可达率 99.4%→97.2%，
        ///      看起来像"变换之后算法崩了"，其实是口径用错了。
        ///   3. 「手动输入」→ 按原点 + 朝向偏转拼。约定基座 +Z 竖直向上、+X 指向工件质心。
        ///   4. 「自动假设」→ null，交给核心库按工件质心猜（报告里会标注 baseAuto=1）。
        /// </summary>
        double[]? BuildBasePose()
        {
            // 1) 手眼标定优先：标定值永远比"看着填的数"可信
            if (_handEyeBasePose != null) return _handEyeBasePose;

            // 2) 点云已在基座系：恒等变换，明确告诉核心库"不要猜"
            if (CmbBaseMode.SelectedIndex == 2)
            {
                return new double[]
                {
                    1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1,
                };
            }
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
            // 自动评估不弹窗打断（刚生成完就要出结论），只写日志；
            // 手动点按钮时才拦下来问 —— 那时人是准备好看结论的。
            if (!auto && !CheckUpDirConsistency()) return;

            var p = Native.robot_safety_default_params();
            // 喷距必须与生成时一致 —— 倾斜是绕面点的球面旋转，喷距错了球心就错，
            // 后面「可达 / 碰撞」的结论全部不成立，所以这里强制取界面上的生成值。
            p.standoff = ParseDouble(TxtStandoff.Text, p.standoff);
            p.safety = ParseDouble(TxtSafetyGap.Text, p.safety);
            // 绝对安全距离/回退硬下限：只有低于它才回退，间隙略小但高于它时保留平滑。
            p.hardClearance = Math.Max(0, ParseDouble(TxtHardClearance.Text, p.hardClearance));
            p.toolLen = ParseDouble(TxtToolLen.Text, p.toolLen);
            // 倾角系数：δ=90°（竖直侧壁）时给多少度，其余按 δ/90 线性缩。库默认 60（2026-10-26 由 30 改）。
            // 命令行 --tiltk 可覆盖，用于"系数扫一遍"的 A/B，不必改代码重编。
            if (!double.IsNaN(_tiltCoefDeg)) p.tiltNonUpDeg = _tiltCoefDeg;
            p.tiltMaxDeg = ParseDouble(TxtTiltMax.Text, p.tiltMaxDeg);
            _effTiltK = p.tiltNonUpDeg;
            _effTiltMaxDeg = p.tiltMaxDeg;
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
            // 要不要按"这台机器人够不够得着"判死。取消勾选 = 只查几何（喷距/倾角/枪体），
            // 适用于点云本来就不是给当前机器人型号/站位用的场合。
            p.requireReach = ChkRequireReach.IsChecked == true ? 1 : 0;
            // 姿态求解策略：1 = 单次逆解（每点 1~2 次 IK）。判据一条不少（可达/关节/连杆/自碰撞/
            // 喷角/间隙全查），省掉的只是"为了挑更容易解的倾角而反复试解"。
            // 开口盒子 19082 点实测：完整搜索 30.8s → 单次逆解 1.59s（19×），可达率同为 99.3%；
            // 代价是深腔件的 P0 会报得更多（不再逐点找可避让的位姿），平面件逐值一致。
            // fastPose 是**位掩码**（见 robot_safety.h）：bit0=单次逆解，bit1=倾角逐点连续。
            p.fastPose = (ChkFastPose.IsChecked == true ? 1 : 0)
                       | (ChkTiltPerPoint.IsChecked == true ? 2 : 0);
            // 平放面的倾角方向（基座系）：0=不指定，+1=往基座 +X，−1=往基座 −X。
            // 安全层里 "往上/往下" 是按基座 +Z 定义的重力，平放面没有"往上扳"可言，
            // 只能选往 +X 还是 −X 侧倾。基座位姿未知时核心库会忽略它。
            p.tiltSideX = CmbTiltSide.SelectedIndex switch { 1 => 1.0, 2 => -1.0, _ => 0.0 };

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
                double[]? upVec = BuildUpVector();
                // 安全层在**基座系**里算：基座位姿已知时，竖直向上直接取基座 +Z（核心库内部做），
                // 界面这个「重力向上」只在基座靠"自动假设"时才参与 —— 那时没有基座系可用。
                App.Log("SAFETY", null, basePose != null
                    ? "竖直向上 = 基座 +Z（由基座位姿给出，忽略「重力向上」下拉）"
                    : "竖直向上 = 「重力向上」下拉（基座为自动假设，没有基座系可用）");
                st = await Task.Run(() =>
                    Native.robot_safety_check_up(pts, pts.Length, cloudXyz, _cloud.Count,
                                                 basePose, ref p, upVec));
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
                err = "当前 SprayTraj.dll 不含安全层导出（robot_safety_check_up）。请重新编译 C++ 核心库并复制到程序目录。";
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
            // 可达率只在"要求可达"时才是真测出来的；关掉运动学后绝不能打 100%，
            // 那只是记账口径（R.reachable 被填成 n），会让人以为测过了。
            if (p.requireReach != 0)
                sb.AppendLine($"可达 {sum.reachable}/{sum.count}" +
                              $"（{(sum.count > 0 ? 100.0 * sum.reachable / sum.count : 0):F1}%）");
            else
                sb.AppendLine($"口径：仅几何（未查可达 / 关节 / 连杆 / 自碰撞）· {sum.count} 点");
            sb.AppendLine($"最小表面距离 {sum.minClearanceMm:F2} mm" +
                          (sum.minClearanceMm < 0 ? "  ← 穿透" : ""));
            sb.AppendLine($"　枪体 {ToolDesc()}　连杆 {LinkDesc()}");
            sb.AppendLine($"喷角 均值 {sum.tiltMeanDeg:F1}° · p95 {sum.tiltP95Deg:F1}° · 最大 {sum.tiltMaxDeg:F1}°");
            // 姿态口径：单次逆解是"快速预览"，深腔/窄槽件可能多报 P0（不再逐点找可避让位姿）。
            // 这一行必须写在结论卡上，否则"这份 REJECT 该不该当真"没法判断。
            if ((p.fastPose & 1) != 0)
                sb.AppendLine("姿态：单次逆解（快速预览 · 正式结论请取消勾选，跑完整搜索）");
            else
                sb.AppendLine("姿态：完整搜索（面级择优 + 逐点救援）");
            // 倾角口径：逐点连续 = 倾角由该点法向唯一决定（tilt = 系数×δ/90，δ=法向与 +Z 夹角），
            // 没有面级常量。面级常量会把连续曲面切成台阶（垃圾桶实测 0°/25° 两块）。
            // 系数一并打出来：它同时决定"侧壁给多少度"，是工艺与可达之间的旋钮，看不到就没法判断结论。
            sb.AppendLine((p.fastPose & 2) != 0
                ? $"倾角：逐点连续 {p.tiltNonUpDeg:F0}°×δ/90（δ=法向与正位置夹角，上限 {p.tiltMaxDeg:F0}°）"
                : $"倾角：面级常量（侧壁基准 {p.tiltNonUpDeg:F0}°×δ/90，面内不跳）");
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

            // 可喷性：安全层判完之后，真正下发的是什么。
            // 这一行是"结论 → 执行"的关键：P0 点不再进 plan.json，剩余碎段按长度清掉。
            var execOk = BuildExecutableMask(got);
            var (nKeep, nDropPts, nKeepRuns) = MaskSummary(got, execOk);
            sb.AppendLine(nDropPts == 0
                ? $"可喷 {nKeep} 点（无过滤，全部下发）"
                : $"可喷 {nKeep} 点 · 过滤 {nDropPts} 点（可喷仅剩 {nKeepRuns} 段）← 枪体/连杆/自碰撞/不可达");
            if (p.requireReach != 0)
            {
                // 三种口径必须分得清，否则"这份结论能不能当真"没法判断：
                // 自动假设 = 原点按质心猜的，可达率只能当预判；
                // 恒等 = 点云已经手眼变换到基座系，是三者里唯一"位置可信"的；
                // 手动 = 现场标定值，可信但要知道填的是什么。
                if (sum.baseAuto != 0)
                    sb.AppendLine($"基座：自动假设 ({sum.basePose[3]:F0}, {sum.basePose[7]:F0}, {sum.basePose[11]:F0}) —— 未标定，可达率仅作预判");
                else if (CmbBaseMode.SelectedIndex == 2)
                    sb.AppendLine("基座：点云已在基座系（恒等）—— 手眼标定后口径");
                else
                    sb.AppendLine($"基座：手动 ({sum.basePose[3]:F0}, {sum.basePose[7]:F0}, {sum.basePose[11]:F0})");
            }
            else
                sb.AppendLine("基座：不用（仅几何口径与机器人型号/站位无关）");
            // 算法库构建时间：改完 C++ 却忘了同步 SprayTraj.dll 时，界面会安静地跑旧算法，
            // 现象就是"明明改了还是老样子"。把这一行摆出来，一眼能对上是哪一版。
            sb.AppendLine($"算法库：{LibStamp()}");

            ShowSafetyCard(verdict, vColor, $"耗时 {sw.ElapsedMilliseconds} ms", sb.ToString().TrimEnd());

            _safetyPts = got;
            // 是否画修正层交给「显示修正后轨迹」复选框，不再用 tiltedPoints>0 判。
            // 原来只有"有点被倾斜"才画，可空移点**只抬枪不改倾角**，于是"只抬了枪"的件
            // 修正层整层不出现，用户会以为修正没生效。
            // ⚠ 不能写 IsChecked.GetValueOrDefault(true)：WPF 的 ToggleButton.IsChecked
            // 默认值是 **false 不是 null**，普通未勾选复选框会老实返回 false —— 曾因此
            // 修正层默认整层不画。XAML 里已给 ChkShowFixed 设 IsChecked="True"，
            // 这里只需 != false（null 三态也按显示算）。
            RenderSafetyOverlay(got, ChkShowFixed.IsChecked != false);

            // 修正后的轨迹自动落盘（点云同名 + 模拟平台 JSON）。
            // 放在渲染之后：导出失败不影响已经算出来的结论展示。
            AutoExport(true);

            string safetyHint = p.requireReach != 0
                ? $"安全评估：{verdict} · 可达 {sum.reachable}/{sum.count} · 最小间隙 {sum.minClearanceMm:F1}mm"
                : $"安全评估（仅几何）：{verdict} · 枪体最小间隙 {sum.minClearanceMm:F1}mm";
            if (nDropPts > 0) safetyHint += $" · 可喷 {nKeep} 点（过滤 {nDropPts} 点）";
            // 导出结果不能被评估结论顶掉：失败/回退信息是用户唯一能看见的线索
            // （2026-10-28 就是因为被顶掉，用户以为"轨迹没生成"）。
            string exportLine = TxtHint.Text.Split('\n').FirstOrDefault(l => l.Contains("自动导出")) ?? "";
            TxtHint.Text = exportLine.Length > 0 ? safetyHint + "\n" + exportLine : safetyHint;
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
            if (_fixedSprayVisual != null) { Viewport.Children.Remove(_fixedSprayVisual); _fixedSprayVisual = null; }
            if (_fixedMoveVisual != null) { Viewport.Children.Remove(_fixedMoveVisual); _fixedMoveVisual = null; }
            if (_fixedDroppedVisual != null) { Viewport.Children.Remove(_fixedDroppedVisual); _fixedDroppedVisual = null; }
            if (arr.Length == 0) return;

            const uint MaskLocal = RobotSafetyRule.Clearance | RobotSafetyRule.JointLimit;

            var red = new Point3DCollection();
            var orange = new Point3DCollection();
            var yellowAll = new List<Point3D>();

            foreach (var q in arr)
            {
                if (q.violation == 0) continue;
                // 显示坐标系：风险点必须和点云/轨迹在同一个系里，否则会聚成"工件外面的一团"
                var dq = DispP(q.x, q.y, q.z);
                var pt = new Point3D(dq.x, dq.y, dq.z);
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

            // 修正后轨迹：喷涂（绿）/ 空移（紫）分组画，空移段必须能和喷涂段区分开。
            // 直接用安全结果 arr 自己的点画 —— 它就是完整的修正轨迹（含步长插补点）。
            // 旧版要求 _lastTraj.Count == arr.Length 才画：那是以「倾斜只挪位置、不改点数」
            // 为前提的下标一一对应假设。2026-10-27 加了步长硬上限插点后，安全结果会比
            // 原轨迹多出插补点（垃圾桶 734→752），点数永远对不上，整层修正轨迹被静默
            // 跳过 —— 现象就是"明明评估过了，修正轨迹却总没有"。插补点的 spray 在
            // 核心库里就是按两端开关枪给的，这里直接消费即可。
            // 是否连线看修正后相邻点距（>25mm 不连）：斜坡+插点之后正常步长 ≤18mm，
            // 超过 25mm 只剩"插点被闸门放弃"的真跳跃 —— 不连是对的，那正是要暴露的段。
            var sprayPts = new Point3DCollection();
            var movePts = new Point3DCollection();
            var droppedPts = new Point3DCollection();
            // 可执行掩码：被过滤的喷涂线画成暗灰，不画进绿色喷涂线。
            var drawOk = BuildExecutableMask(arr);
            double dStepDbg = 0.0;
            if (showFixed && arr.Length > 1)
            {
                for (int i = 1; i < arr.Length; i++)
                {
                    var o0 = arr[i - 1];
                    var o1 = arr[i];
                    double dx = o1.x - o0.x, dy = o1.y - o0.y, dz = o1.z - o0.z;
                    double d = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                    if (d > dStepDbg) dStepDbg = d;
                    if (d > 25.0) continue;
                    // 两端都喷才算喷涂段；只要有一端是空移（枪关）就整段按空移画，
                    // 这样"下枪/抬枪"那一小截也不会被涂成喷涂色。
                    // 但被安全层拦下的喷涂段要单独归到 dropped —— 它不会被执行，
                    // 染成绿色就是在展示一条机器人不会走的线。
                    var dst = (o0.spray == 1 && o1.spray == 1)
                        ? ((drawOk[i - 1] && drawOk[i]) ? sprayPts : droppedPts)
                        : movePts;
                    // 显示坐标系：连线的端点也要搬，否则修正层会整层飘在点云外面
                    var d0 = DispP(o0.x, o0.y, o0.z);
                    var d1 = DispP(o1.x, o1.y, o1.z);
                    dst.Add(new Point3D(d0.x, d0.y, d0.z));
                    dst.Add(new Point3D(d1.x, d1.y, d1.z));
                }
                if (sprayPts.Count > 0)
                {
                    _fixedSprayVisual = new LinesVisual3D
                    {
                        Points = sprayPts,
                        // 绿：修正后喷涂。早先用的是青色 #22D3EE —— 和原轨迹的空移蓝 #40B0F0
                        // 是一对邻居色，远看分不出来（这正是"空移段看着还是同一个颜色"的来源之一）。
                        // 绿在这套配色里没被占用（红=原喷涂、蓝=原空移、紫=修正空移、黄/橙=风险），
                        // 四个色相互不挨着，一眼就能分开。
                        Color = Color.FromArgb(215, 74, 222, 128),
                        Thickness = 1.2
                    };
                    Viewport.Children.Add(_fixedSprayVisual);
                }
                if (movePts.Count > 0)
                {
                    _fixedMoveVisual = new LinesVisual3D
                    {
                        Points = movePts,
                        Color = Color.FromArgb(220, 192, 132, 252),  // 紫：修正后空移
                        Thickness = 1.6
                    };
                    Viewport.Children.Add(_fixedMoveVisual);
                }
                if (droppedPts.Count > 0)
                {
                    _fixedDroppedVisual = new LinesVisual3D
                    {
                        Points = droppedPts,
                        Color = Color.FromArgb(190, 120, 120, 120),  // 暗灰：被安全层过滤，不下发
                        Thickness = 1.0
                    };
                    Viewport.Children.Add(_fixedDroppedVisual);
                }
            }

            ApplyRiskVisibility();
            ApplyFixedVisibility();
            App.Log("SAFETY/DRAW", null,
                $"风险点 红={red.Count} 橙={orange.Count} 黄={nYellowDraw}/{yellowAll.Count}（抽稀）；" +
                $"修正轨迹线段 喷涂={(_fixedSprayVisual?.Points.Count ?? 0) / 2}" +
                $" 空移={(_fixedMoveVisual?.Points.Count ?? 0) / 2}" +
                $" 过滤={droppedPts.Count / 2}" +
                $" · showFixed={showFixed} arrLen={arr.Length} 收集喷={sprayPts.Count / 2} 空={movePts.Count / 2} max步={dStepDbg:F1}");
        }

        /// <summary>当前加载的 SprayTraj.dll 的构建时间（旧库没有这个导出 → 明确告诉用户）。</summary>
        static string LibStamp()
        {
            try
            {
                IntPtr p = Native.robot_safety_build_stamp();
                return p == IntPtr.Zero ? "未知" : (Marshal.PtrToStringAnsi(p) ?? "未知");
            }
            catch (EntryPointNotFoundException) { return "旧版库（无构建戳，请同步 SprayTraj.dll）"; }
            catch { return "未知"; }
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
            if (Viewport == null) return;
            bool show = ChkShowFixed == null || ChkShowFixed.IsChecked.GetValueOrDefault(true);
            // 修正层的空移段受两道开关管：修正轨迹总开关 + 「显示空移段」。
            // 关掉空移后修正层也只留喷涂段，跟原轨迹层口径一致，避免"关了还在"。
            bool showMove = show && (ChkShowMove == null || ChkShowMove.IsChecked.GetValueOrDefault(true));
            Sync(_fixedSprayVisual, show);
            Sync(_fixedMoveVisual, showMove);
            // 被过滤的喷涂线跟喷涂段同开关：它替代的是喷涂段，关掉修正层就该一起消失
            Sync(_fixedDroppedVisual, show);

            void Sync(Visual3D? v, bool want)
            {
                if (v == null) return;
                bool inTree = Viewport.Children.Contains(v);
                if (want && !inTree) Viewport.Children.Add(v);
                else if (!want && inTree) Viewport.Children.Remove(v);
            }
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

        /// <summary>
        /// 写修正后轨迹 CSV（无对话框版本，供命令行验收复用）。
        /// 手眼生效时这里会搬到基座系；关节角/间隙/倾角是**相对无关**的标量，原样保留。
        /// </summary>
        // ══════════════════════════════════════════════════════════════════
        //  可喷性判定（安全层 → 执行口径的过滤）
        //
        //  安全层只会「判」：P0 点标 violation。但标了不等于拦住了 —— 早先 plan.json
        //  里 IsValid 恒 true、GunOpen 只看 spray，于是枪体撞件、IK 解不出的点照样
        //  下发给机器人（背板_X+800 工具长 450 时有 128 点枪体撞件、2272 点不可达，
        //  全部原样进 json）。这里把「判」变成「拦」。
        //
        //  为什么按「连续喷涂段」整段拦，而不是逐点剔：
        //  删掉一行中间的几个点，前后两点仍是 GunOpen=true，机器人会关不掉枪地
        //  直线穿过去 —— 从「撞一点」变成「撞一整条」。整段不喷是唯一不需要重新
        //  规划抬枪空移就能保证不撞的做法；漏一条 4mm 的线，比撞坏工件或枪好。
        // ══════════════════════════════════════════════════════════════════

        /// <summary>P0 = 不可执行：IK 解不出、连杆/枪体撞件、自碰撞。</summary>
        static bool IsBlocked(uint violation) => (violation & MaskP0) != 0;

        /// <summary>
        /// 可执行喷涂点的最少连续点数。约 4mm/点，5 点 ≈ 20mm：
        /// 短于这个长度的碎段，机器人光是开关枪就过去了，实际喷不上什么，
        /// 却要额外一次抬枪/关枪动作，且极易在起停处补枪过量。
        /// </summary>
        const int MinRunPoints = 5;

        /// <summary>
        /// 逐点可执行掩码。
        /// 1) P0 点（IK 解不出 / 枪体撞件 / 连杆撞件 / 自碰撞）直接标 false；
        /// 2) 连续可喷点不足 <see cref="MinRunPoints"/> 的碎段一并标 false。
        ///
        /// 为什么不再"整条喷涂线一起否"：背板_X+800 在 450mm 枪长下 87% 的喷涂点
        /// 不可达，整线否的结果是 6 条线全灭、可喷 0 点——把"过滤"变成了"什么都不做"。
        /// 逐点过滤能如实反映"哪段真的喷得到"，剩下的碎段再按长度清掉，
        /// 既不放过撞件点，也不把能喷的地方一起否掉。
        /// </summary>
        static bool[] BuildExecutableMask(RobotSafetyPoint[] pts)
        {
            var ok = new bool[pts.Length];
            for (int i = 0; i < pts.Length; ++i)
                ok[i] = pts[i].spray != 0 && !IsBlocked(pts[i].violation);

            // 清掉过短的连续可喷碎段
            int i0 = 0;
            while (i0 < pts.Length)
            {
                if (!ok[i0]) { ++i0; continue; }
                int i1 = i0;
                while (i1 < pts.Length && ok[i1]) ++i1;
                if (i1 - i0 < MinRunPoints)
                    for (int k = i0; k < i1; ++k) ok[k] = false;
                i0 = i1;
            }
            return ok;
        }

        /// <summary>过滤统计：可喷点 / 被拦点 / 保留的独立喷涂段数。</summary>
        (int keep, int drop, int dropSegs) MaskSummary(RobotSafetyPoint[] pts, bool[] ok)
        {
            int keep = 0, drop = 0, runs = 0;
            for (int i = 0; i < pts.Length; ++i)
            {
                if (pts[i].spray == 0) continue;
                if (ok[i]) ++keep; else ++drop;
            }
            for (int i = 0; i < pts.Length; ++i)
                if (ok[i] && (i == 0 || !ok[i - 1])) ++runs;
            return (keep, drop, runs);
        }

        void WriteFixedCsv(string path)
        {
            var pts = ExportSafetyPts();
            var ok = BuildExecutableMask(pts);
            var sb = new System.Text.StringBuilder();
            // 末列 sprayable：1=这条喷涂段能喷（也就能走），0=被安全层拦下。
            // 保留全部点是为了诊断——想知道"哪一片漏喷、为什么漏"，必须看得见被拦的点；
            // 真正下发给机器人的是 plan.json，那里会按 IsValid/GunOpen 把它们关掉。
            sb.Append("X,Y,Z,I,J,K,speed,spray,tilt_deg,clearance_mm,self_gap_mm,margin_deg,")
              .Append("standoff_err_mm,lifted_mm,cond,reverted,J1_deg,J2_deg,J3_deg,J4_deg,J5_deg,J6_deg,violation,sprayable\n");
            for (int i = 0; i < pts.Length; i++)
            {
                var q = pts[i];
                sb.Append($"{q.x:F3},{q.y:F3},{q.z:F3},{q.i:F4},{q.j:F4},{q.k:F4},{q.speed:F1},{q.spray},");
                sb.Append($"{q.tiltDeg:F2},{q.clearanceMm:F2},{q.selfGapMm:F2},{q.marginDeg:F2},");
                sb.Append($"{q.standoffErrMm:E6},{q.liftedMm:F1},{q.cond:F1},{q.reverted}");
                var th = q.theta ?? new double[6];
                for (int k = 0; k < 6; k++)
                    sb.Append($",{(k < th.Length ? th[k] : 0.0):F3}");
                sb.Append($",{RuleText(q.violation)},{ok[i] switch { true => 1, false => 0 }}\n");
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
            if ((v & RobotSafetyRule.StepJump) != 0) parts.Add("STEP_JUMP");
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

            // 点云与面片都要按**显示坐标系**累加：轨迹视觉对象的 Path 已经是显示系的了，
            // 这里若拿相机系的原始点云/面顶点，包围盒会被撑成两个相距 1.5 米的团，
            // 取景后看起来就是"工件和轨迹各在一头"。
            foreach (var p in _cloud)
            {
                var d = DispP(p.x, p.y, p.z);
                Acc(d.x, d.y, d.z);
            }
            foreach (var plRaw in _planes)
            {
                var pl = DispPlane(plRaw);
                int nv = Math.Min(pl.polyVertCount, pl.polyVerts?.Length / 3 ?? 0);
                for (int k = 0; k < nv; k++)
                    Acc(pl.polyVerts![k * 3], pl.polyVerts[k * 3 + 1], pl.polyVerts[k * 3 + 2]);
            }
            foreach (var t in _trajVisuals)
                foreach (var p in t.Path) Acc(p.X, p.Y, p.Z);
            foreach (var t in _moveVisuals)
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
                // 面0也要按显示坐标系取，否则手眼启用时相机会对着"相机系里的那块面"，
                // 而画面上的几何全在基座系 —— 视角直接跑到工件外面去。
                var p0 = DispPlane(_planes[0]);
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

        void ChkShowMove_Changed(object sender, RoutedEventArgs e)
        {
            ApplyTrajVisibility();
            ApplyFixedVisibility();   // 修正层也分喷/空，开关要同步生效
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
            // XAML 里 CheckBox 的 IsChecked="True" 在**解析阶段**就会触发 Checked 事件，
            // 那时 x:Name 字段还没连上、Viewport 也还没创建（它在本文件靠后的位置）。
            // 少了这道判断，启动直接 NullReferenceException（TargetInvocationException 包着，看不出真因）。
            if (Viewport == null || ChkShowTraj == null || ChkShowMove == null) return;
            bool show = ChkShowTraj.IsChecked.GetValueOrDefault(true);
            // 空移段可以单独关掉：只想看"真正喷涂的那几行"时，蓝线会挡视线
            bool showMove = show && ChkShowMove.IsChecked.GetValueOrDefault(true);
            foreach (var t in _trajVisuals)
            {
                bool inTree = Viewport.Children.Contains(t);
                if (show && !inTree) Viewport.Children.Add(t);
                else if (!show && inTree) Viewport.Children.Remove(t);
            }
            foreach (var t in _moveVisuals)
            {
                bool inTree = Viewport.Children.Contains(t);
                if (showMove && !inTree) Viewport.Children.Add(t);
                else if (!showMove && inTree) Viewport.Children.Remove(t);
            }
        }
    }
}
