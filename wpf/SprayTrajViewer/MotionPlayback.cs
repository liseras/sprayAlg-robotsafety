using System;
using System.Diagnostics;
using System.Windows;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Media3D;
using System.Windows.Threading;
using HelixToolkit.Wpf;

namespace SprayTrajViewer
{
    public partial class MainWindow
    {
        private RobotMotionPoint[] _motionProfile = Array.Empty<RobotMotionPoint>();
        private RobotMotionLimits _evaluatedMotionLimits;
        private DispatcherTimer? _motionTimer;
        private readonly Stopwatch _motionClock = new();
        private double _motionStartTime;
        private PointsVisual3D? _motionTcp;
        private LinesVisual3D? _motionGunAxis;

        RobotMotionLimits ReadMotionLimits() => new()
        {
            StraightSpeed = ParseDouble(TxtMotionStraight.Text, double.NaN),
            TurnSpeed = ParseDouble(TxtMotionTurn.Text, double.NaN),
            TcpAcceleration = 500,
            TcpJerk = 5000,
            GunAxisSpeed = 40,
            JointSpeed = 90,
            JointAcceleration = 300,
            JointJerk = 2000,
        };

        void MotionTimingOption_Changed(object sender, RoutedEventArgs e)
        {
            if (!IsLoaded || _safetyPts.Length == 0) return;
            var limits = ReadMotionLimits();
            if (_motionProfile.Length != 0 && ChkMotionTiming.IsChecked == true &&
                limits.StraightSpeed == _evaluatedMotionLimits.StraightSpeed &&
                limits.TurnSpeed == _evaluatedMotionLimits.TurnSpeed) return;
            SetMotionProfile(Array.Empty<RobotMotionPoint>());
            CardSafety.Visibility = Visibility.Collapsed;
            TxtSafetyInline.Text = "运动时间参数已修改，请重新评估";
        }

        void SetMotionProfile(RobotMotionPoint[] profile)
        {
            StopMotionPlayback();
            _motionProfile = profile;
            if (_motionTcp != null) Viewport.Children.Remove(_motionTcp);
            if (_motionGunAxis != null) Viewport.Children.Remove(_motionGunAxis);
            _motionTcp = null;
            _motionGunAxis = null;
            MotionPlaybackPanel.Visibility = profile.Length > 1 ? Visibility.Visible : Visibility.Collapsed;
            MotionTimeSlider.Value = 0;
            MotionTimeSlider.Maximum = profile.Length > 1 ? profile[^1].Time : 0;
            if (profile.Length > 1) TxtMotionCursor.Text = $"0.00 / {profile[^1].Time:F2} s";
        }

        void StopMotionPlayback()
        {
            _motionTimer?.Stop();
            _motionClock.Stop();
            BtnMotionPlay.Content = "播放轨迹";
        }

        void BtnMotionPlay_Click(object sender, RoutedEventArgs e)
        {
            if (_motionProfile.Length < 2) return;
            if (_motionClock.IsRunning)
            {
                StopMotionPlayback();
                return;
            }
            if (MotionTimeSlider.Value >= MotionTimeSlider.Maximum) MotionTimeSlider.Value = 0;
            _motionStartTime = MotionTimeSlider.Value;
            RenderMotionFrame(_motionStartTime);
            if (_motionTimer == null)
            {
                _motionTimer = new DispatcherTimer(DispatcherPriority.Render)
                {
                    Interval = TimeSpan.FromMilliseconds(1000d / 60)
                };
                _motionTimer.Tick += (_, _) =>
                {
                    MotionTimeSlider.Value = Math.Min(MotionTimeSlider.Maximum,
                        _motionStartTime + _motionClock.Elapsed.TotalSeconds);
                    if (MotionTimeSlider.Value >= MotionTimeSlider.Maximum) StopMotionPlayback();
                };
            }
            BtnMotionPlay.Content = "暂停";
            _motionClock.Restart();
            _motionTimer.Start();
        }

        void MotionTimeSlider_PreviewMouseLeftButtonDown(object sender, MouseButtonEventArgs e)
            => StopMotionPlayback();

        void MotionTimeSlider_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            if (_motionProfile.Length == 0 || !IsLoaded) return;
            RenderMotionFrame(e.NewValue);
        }

        void RenderMotionFrame(double time)
        {
            if (_motionProfile.Length < 2 || _motionProfile.Length != _safetyPts.Length) return;
            int left = 0, right = _motionProfile.Length - 1;
            while (left + 1 < right)
            {
                int middle = left + (right - left) / 2;
                if (_motionProfile[middle].Time <= time) left = middle;
                else right = middle;
            }
            var a = _safetyPts[left];
            var b = _safetyPts[right];
            var t0 = _motionProfile[left];
            var t1 = _motionProfile[right];
            double duration = t1.Time - t0.Time;
            var start = new Point3D(a.x, a.y, a.z);
            var delta = new Vector3D(b.x - a.x, b.y - a.y, b.z - a.z);
            double length = delta.Length;
            double u = Math.Clamp((time - t0.Time) / duration, 0, 1);
            var (progress, speed) = MotionProgress(u, duration, length, t0.Speed, t1.Speed);
            var p = start + delta * progress;
            var gun = InterpolateGunAxis(new Vector3D(a.i, a.j, a.k), new Vector3D(b.i, b.j, b.k), progress);
            var pos = DispP(p.X, p.Y, p.Z);
            var dir = DispV(gun.X, gun.Y, gun.Z);
            var tcp = new Point3D(pos.x, pos.y, pos.z);
            bool spray = a.spray != 0 && b.spray != 0;
            bool flagged = IsBlocked(a.violation | b.violation) ||
                (t1.Flags & RobotMotionPoint.JointBranchJump) != 0;
            Color color = flagged ? Colors.Orange : spray ? Colors.LawnGreen : Colors.DeepSkyBlue;
            if (_motionTcp == null || _motionGunAxis == null)
            {
                _motionTcp = new PointsVisual3D { Size = 12 };
                _motionGunAxis = new LinesVisual3D { Thickness = 3 };
                Viewport.Children.Add(_motionTcp);
                Viewport.Children.Add(_motionGunAxis);
            }
            _motionTcp.Color = color;
            _motionTcp.Points = new Point3DCollection { tcp };
            _motionGunAxis.Color = color;
            _motionGunAxis.Points = new Point3DCollection { tcp,
                tcp + 35 * new Vector3D(dir.x, dir.y, dir.z) };
            TxtMotionCursor.Text = $"{time:F2} / {_motionProfile[^1].Time:F2} s · {speed:F1} mm/s · " +
                (flagged ? "风险段" : spray ? "喷涂" : "空移");
        }

        // Same quintic as the native timing bound. The stationary case uses
        // normalized travel to animate a change in orientation at a fixed TCP.
        static (double Progress, double Speed) MotionProgress(double u, double duration, double length,
            double startSpeed, double endSpeed)
        {
            double a = length > 1e-6 ? startSpeed * duration / (5 * length) : 0;
            double b = length > 1e-6 ? endSpeed * duration / (5 * length) : 0;
            double v = 1 - u, u2 = u * u, v2 = v * v;
            double progress = 5 * a * u * v2 * v2 + 20 * a * u2 * v2 * v +
                10 * (1 - 2 * b) * u2 * u * v2 + 5 * (1 - b) * u2 * u2 * v + u2 * u2 * u;
            double speed = 5 * length / duration * (a * v2 * (v2 + 4 * u * v) +
                6 * (1 - 2 * (a + b)) * u2 * v2 + b * u2 * (4 * u * v + u2));
            return (Math.Clamp(progress, 0, 1), Math.Max(0, speed));
        }

        static Vector3D InterpolateGunAxis(Vector3D a, Vector3D b, double fraction)
        {
            a.Normalize();
            b.Normalize();
            double cosine = Math.Clamp(Vector3D.DotProduct(a, b), -1, 1);
            if (cosine > 0.999999)
            {
                var interpolated = a * (1 - fraction) + b * fraction;
                interpolated.Normalize();
                return interpolated;
            }
            var tangent = b - a * cosine;
            if (tangent.LengthSquared < 1e-12)
                tangent = Vector3D.CrossProduct(a, Math.Abs(a.Z) < 0.9 ?
                    new Vector3D(0, 0, 1) : new Vector3D(1, 0, 0));
            tangent.Normalize();
            double angle = fraction * Math.Acos(cosine);
            return a * Math.Cos(angle) + tangent * Math.Sin(angle);
        }
    }
}
