using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Windows;
using System.Windows.Media;
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
        private readonly List<ModelVisual3D> _planeVisuals = new();
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
                    if (System.IO.File.Exists(TxtPath.Text))
                        LoadCloud(TxtPath.Text);
                }
                catch { }
            };
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
            }
            catch (Exception ex)
            {
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
                $"  连续 S 形喷涂轨迹（全程 spray=1）";

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

            // 实体管渲染：网格实体，任何缩放级别都不会丢段/虚线
            // 去除连续重复点（行与弧接点处），避免零长度段产生 NaN 网格
            // 多平面分割时，不同面之间的转移是空程跳跃——按距离阈值断开成多条管，
            // 避免把跨面转移画成斜穿的连线。
            double breakDist = 15.0;   // 相邻点距约 4mm，超过 15mm 视为跨段跳跃
            var tubes = new List<TubeVisual3D>();
            Point3DCollection path = new();
            Point3D last = default;
            void FlushTube()
            {
                if (path.Count >= 2)
                {
                    tubes.Add(new TubeVisual3D
                    {
                        Path = path,
                        Diameter = 2.5,
                        Fill = new SolidColorBrush(Color.FromRgb(255, 48, 32)),
                        ThetaDiv = 6
                    });
                }
                path = new Point3DCollection();
            }
            foreach (var p in traj)
            {
                // 空移点（spray=0）：断开当前管并跳过该点——空移段不显示、不画线
                if (p.spray == 0)
                {
                    FlushTube();
                    last = new Point3D(p.x, p.y, p.z);
                    continue;
                }
                var cur = new Point3D(p.x, p.y, p.z);
                if (path.Count > 0 && (cur - last).Length > breakDist)
                    FlushTube();          // 跨段：断开，另起一条管
                if (path.Count == 0 || (cur - last).Length > 1e-6)
                {
                    path.Add(cur);
                    last = cur;
                }
            }
            FlushTube();

            foreach (var t in tubes)
            {
                _trajVisuals.Add(t);
                Viewport.Children.Add(t);
            }
            Viewport.ZoomExtents(0);   // 立即取景到全部内容
            ApplyTrajVisibility();
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
            n.Normalize();
            var c = new Point3D(p.cx, p.cy, p.cz);

            // 直接用 C++ 算好的局部基 + 内点包围盒（不做 PCA、不做距离近似）
            Vector3D u = new Vector3D(p.ux, p.uy, p.uz);
            Vector3D v = new Vector3D(p.vx, p.vy, p.vz);
            if (u.Length < 1e-6 || v.Length < 1e-6)
            {
                // 兜底：旧 DLL 或未填充时用 cross(z,n)
                u = Vector3D.CrossProduct(new Vector3D(0, 0, 1), n);
                if (u.Length < 1e-6) u = new Vector3D(1, 0, 0);
                u.Normalize();
                v = Vector3D.CrossProduct(n, u);
                v.Normalize();
            }
            else { u.Normalize(); v.Normalize(); }

            double minU = p.minU, maxU = p.maxU, minV = p.minV, maxV = p.maxV;
            if (maxU - minU < 1e-6 || maxV - minV < 1e-6) return;

            var mesh = new MeshGeometry3D();
            int seg = 20;
            for (int i = 0; i <= seg; i++)
                for (int j = 0; j <= seg; j++)
                {
                    double fu = minU + (maxU - minU) * i / seg;
                    double fv = minV + (maxV - minV) * j / seg;
                    var pt = c + u * fu + v * fv;
                    mesh.Positions.Add(pt);
                    mesh.TextureCoordinates.Add(new System.Windows.Point((double)i / seg, (double)j / seg));
                }
            for (int i = 0; i < seg; i++)
                for (int j = 0; j < seg; j++)
                {
                    int a = i * (seg + 1) + j;
                    int b = a + 1;
                    int d = (i + 1) * (seg + 1) + j;
                    int e = d + 1;
                    mesh.TriangleIndices.Add(a); mesh.TriangleIndices.Add(d); mesh.TriangleIndices.Add(b);
                    mesh.TriangleIndices.Add(b); mesh.TriangleIndices.Add(d); mesh.TriangleIndices.Add(e);
                }

            var mat = new DiffuseMaterial(new SolidColorBrush(Color.FromArgb(90, 64, 156, 255)));
            var model = new GeometryModel3D(mesh, mat) { BackMaterial = mat };
            var visual = new ModelVisual3D { Content = model };
            _planeVisuals.Add(visual);
            Viewport.Children.Add(visual);
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
