using System.Globalization;
using System.Text;
using System.Text.Json;
using DensoArmKine.Kinematics;
using DensoArmKine.Models;

namespace DensoIkProbe;

/// <summary>
/// 离线逆解探针。把主城序（DensoArmKine）的 DLS 逆解接到喷涂轨迹上，逐点求解，
/// 回答唯一一个问题：<b>这条轨迹，VS-087 能不能连续走通</b>。
///
/// ─────────── 与主城序的关系：只读 ───────────
/// 逆解/正解/连杆库全部来自主城序源码的<b>链接</b>引用（见 csproj），本工程不含任何算法实现，
/// 也不修改主城序任何一行。主城序改了，本探针下次编译自动跟随。
///
/// ─────────── 硬约束：绝不动机器人 ───────────
/// 只链接纯数学与纯数据文件；Hardware/（b-CAP / 串口）、Motion/（运动下发）、
/// IkCompareCli（Run 需要 ip，会连真机）一律不链接。编译产物物理上不具备通信能力。
///
/// ─────────── 为什么这里要「种子链」 ───────────
/// DLS 迭代从 seed 出发收敛到<b>离 seed 最近</b>的解支，这与控制器 P2J 的"最近优先"一致。
/// 沿轨迹把上一点解作下一点 seed ⇒ 整条轨迹落在同一解支上 ⇒ 天然避免解支翻转（翻腕）。
/// 若某点收敛到离 seed 很远的解，那本身就是「要翻腕」的硬信号，必须单独计数。
/// </summary>
public static class Program
{
    // ================================================================ 入口

    public static int Main(string[] args)
    {
        Console.OutputEncoding = Encoding.UTF8;

        // ---------- 步骤 0：解析命令行 ----------
        var opt = ProbeOptions.Parse(args);
        if (opt.Help) { ProbeOptions.PrintUsage(); return 0; }

        if (!File.Exists(opt.TrajPath))
        {
            Console.Error.WriteLine($"轨迹文件不存在：{opt.TrajPath}");
            return 2;
        }

        // ---------- 步骤 1：加载轨迹点 ----------
        List<TrajPoint> traj = TrajLoader.Load(opt.TrajPath);
        if (traj.Count == 0)
        {
            Console.Error.WriteLine("轨迹为空或格式无法识别。");
            return 2;
        }
        Console.WriteLine($"轨迹：{Path.GetFileName(opt.TrajPath)}  共 {traj.Count} 点"
                          + (traj[0].HasEuler ? "（含欧拉角）" : "（仅枪轴，欧拉角由枪轴反算）"));

        // ---------- 步骤 2：构建机型与工具 ----------
        RobotModel model = DensoLibrary.Get(opt.ModelId);
        ToolConfig tool = BuildTool(opt, model);
        Console.WriteLine($"机型：{model.Id}（reach {model.SpecReachMm}mm）  "
                          + $"工具：Offset=({tool.OffsetX:F3}, {tool.OffsetY:F3})  LengthZ={tool.LengthZ:F3}");

        // ---------- 步骤 3：工具参数自校验（可选） ----------
        if (opt.CalibPath != null) CalibrateTool(opt, model, tool);

        // ---------- 步骤 4：逐点逆解（种子链） ----------
        SolveReport rep = SolveChain(model, tool, traj, opt);

        // ---------- 步骤 5：体检 ----------
        Diagnose(rep, model, opt);

        // ---------- 步骤 6：输出 ----------
        rep.Print();

        if (opt.OutPath != null)
        {
            rep.Save(opt.OutPath, opt, model, tool, traj);
            Console.WriteLine();
            Console.WriteLine($"报告已写出：{opt.OutPath}");
        }
        return rep.Ok ? 0 : 1;
    }

    // ================================================================ 步骤 2：工具

    /// <summary>
    /// 组装工具参数。命令行显式给了就用给定的，否则用探针内置的默认值。
    ///
    /// 默认值 (0.038, 2.214, 450.756) 不是拍脑袋来的：它由 <c>_safe.csv</c> 里的关节角
    /// 做 FK 反推 <c>target − flange</c> 在法兰系下的位移，对 1467 个点取其恒定分量得到
    /// （dx 抖动 0.07mm、dy 0.016mm，只剩关节角三位小数的舍入误差）。
    /// 与 <see cref="DensoLibrary.Vs087"/> 的默认工具 (-2.21, 6.77, 453.08) 不同 ——
    /// 后者是主城序按现场标定 T4 记的值，而轨迹这道工序实际用的是 450 长枪。
    /// 用错工具会让逆解在网络里"位置对了但腕部差 450mm"，所以这必须是显式可覆盖的参数。
    /// </summary>
    private static ToolConfig BuildTool(ProbeOptions opt, RobotModel model)
    {
        if (opt.ToolExplicit)
            return new ToolConfig
            {
                OffsetX = opt.ToolOffX,
                OffsetY = opt.ToolOffY,
                OffsetZ = 0,
                LengthZ = opt.ToolLen,
                TwistDeg = 0,
            };

        return new ToolConfig { OffsetX = 0.038, OffsetY = 2.214, OffsetZ = 0, LengthZ = 450.756, TwistDeg = 0 };
    }

    // ================================================================ 步骤 3：校准

    /// <summary>
    /// 用一份「已知关节角」的轨迹（通常是安全层输出的 _safe.csv）反查工具参数：
    /// 对每个点算 <c>d = R_flangeᵀ · (p_target − p_flange)</c>，若工具参数正确，d 应恒等于
    /// (OffsetX, OffsetY, LengthZ)。取均值即得，同时报出抖动看是否可信。
    ///
    /// 为什么值得做这一步：工具长度错 400mm 时，逆解会"收敛得很好"却解出完全错误的腕部姿态
    /// ——残差为零、姿态漂亮、真机撞枪。这种错必须在校准阶段拦住，不能靠肉眼。
    /// </summary>
    private static void CalibrateTool(ProbeOptions opt, RobotModel model, ToolConfig tool)
    {
        if (!File.Exists(opt.CalibPath)) { Console.WriteLine($"（校准文件不存在，跳过）{opt.CalibPath}"); return; }
        List<TrajPoint> refs = TrajLoader.Load(opt.CalibPath);
        if (refs.Count == 0 || !refs[0].HasJoints) { Console.WriteLine("（校准文件无关节角，跳过）"); return; }

        var dxs = new List<double>(); var dys = new List<double>(); var dzs = new List<double>();
        foreach (TrajPoint p in refs)
        {
            FkResult fk = FkSolver.Solve(model, p.Joints, tool);
            Mat4 r = fk.Flange.RotationPart;
            Vec3 d = p.Pos - fk.Flange.Translation;
            // 法兰系下表达：Rᵀ · d
            dxs.Add(r.M00 * d.X + r.M10 * d.Y + r.M20 * d.Z);
            dys.Add(r.M01 * d.X + r.M11 * d.Y + r.M21 * d.Z);
            dzs.Add(r.M02 * d.X + r.M12 * d.Y + r.M22 * d.Z);
        }

        Console.WriteLine();
        Console.WriteLine("── 工具参数自校验（用已知关节角反推 target−flange 在法兰系下的位移）──");
        Report("dx", dxs);
        Report("dy", dys);
        Report("dz", dzs);

        static void Report(string nm, List<double> v)
        {
            double mean = v.Average();
            double sd = Math.Sqrt(v.Select(x => (x - mean) * (x - mean)).Average());
            Console.WriteLine($"  {nm}: mean={mean,9:F4}  std={sd,7:F4}  "
                              + $"min={v.Min(),9:F4}  max={v.Max(),9:F4}  {(sd < 0.5 ? "✓ 恒定" : "⚠ 不恒定")}");
        }
    }

    // ================================================================ 步骤 4：逆解

    /// <summary>
    /// 沿轨迹逐点求解，上一点解作下一点 seed（种子链）。
    ///
    /// 首点没有 seed，用「零位 + 一组几何上合理的候选」试，取第一个限位内收敛解。
    /// 之后每点只从上一解出发解一次 —— 这不是为了省时间，而是<b>刻意</b>：
    /// 多起点会掩盖「这点必须翻腕才能到」的事实，而那正是我们要统计的东西。
    /// 只有当链条断了（上一点没解出来）才回退到多起点，尽量续上。
    /// </summary>
    private static SolveReport SolveChain(RobotModel model, ToolConfig tool, List<TrajPoint> traj, ProbeOptions opt)
    {
        var rep = new SolveReport { Count = traj.Count };
        var sw = System.Diagnostics.Stopwatch.StartNew();

        double[] seed = null;
        for (int i = 0; i < traj.Count; i++)
        {
            TrajPoint p = traj[i];
            // 逆解需要完整 6 维位姿。plan.json 自带欧拉角（首选路径）；
            // csv-only 输入只有 I/J/K（枪轴），缺绕轴自转角 ⇒ 无法唯一确定姿态，必须明确报错而不是瞎猜。
            if (!p.HasEuler) { rep.AddFail(i, p, "该点只有枪轴无欧拉角，姿态不唯一（请用 *.plan.json）"); continue; }
            Vec3 euler = p.Euler;

            IkResult r;
            int tried = 1;
            if (seed == null)
            {
                // 首点（或链条断开后）：多起点找第一个限位内收敛解
                r = IkSeeds.SolveBest(model, tool, p.Pos, euler, IkTargetFrame.Tcp,
                                      ZeroSeed(model), out int used, out tried, preferInLimit: true);
                if (r is { Converged: true }) seed = (double[])r.JointsDeg.Clone();
            }
            else
            {
                r = IkSeeds.SolveOne(model, tool, p.Pos, euler, IkTargetFrame.Tcp, (double[])seed.Clone());
                if (r is { Converged: true }) seed = (double[])r.JointsDeg.Clone();
                else seed = null;   // 断链：下一点重新找起点
            }

            if (r == null) { rep.AddFail(i, p, "逆解返回空"); continue; }
            if (!r.Converged) { rep.AddFail(i, p, r.Error ?? "未收敛"); continue; }

            rep.AddOk(i, p, r, tried);
        }

        sw.Stop();
        rep.ElapsedMs = sw.Elapsed.TotalMilliseconds;
        rep.MaxIter = rep.Oks.Count == 0 ? 0 : rep.Oks.Max(x => x.Result.Iterations);
        return rep;
    }

    /// <summary>找出「枪轴」在模型里的表达，用于从 I/J/K 反算欧拉角。</summary>
    private static double[] ZeroSeed(RobotModel m) => new double[m.Joints.Count];

    // ================================================================ 步骤 5：体检

    /// <summary>
    /// 五项硬指标 —— 只问「机器人敢不敢走」，不问余量好不好：
    ///   ① 逆解失败（物理上摆不出来）
    ///   ② 关节超软限位
    ///   ③ 腕部奇异 |θ5| &lt; 6°（此处雅可比退化，真机会抖/失控）
    ///   ④ 相邻步长 max|Δθ| 过大（轨迹不连续，真机走不到）
    ///   ⑤ 解支翻转（单关节相邻 Δ &gt; 120°，= 翻腕）
    /// </summary>
    private static void Diagnose(SolveReport rep, RobotModel model, ProbeOptions opt)
    {
        const double kJ5Sing = 6.0;
        const double kStepWarn = 5.0;
        const double kStepBad = 30.0;
        const double kFlip = 120.0;

        for (int k = 0; k < rep.Oks.Count; k++)
        {
            Solved s = rep.Oks[k];
            s.LimitMask = s.Result.OverLimit;
            s.AnyLimit = s.Result.OverLimit.Any(b => b);
            s.Theta5 = s.Result.JointsDeg[4];
            s.J5Sing = Math.Abs(s.Theta5) < kJ5Sing;
        }

        for (int k = 1; k < rep.Oks.Count; k++)
        {
            double[] a = rep.Oks[k - 1].Result.JointsDeg;
            double[] b = rep.Oks[k].Result.JointsDeg;
            int arg = 0; double mx = 0;
            for (int j = 0; j < 6; j++)
            {
                double d = Math.Abs(b[j] - a[j]);
                if (d > mx) { mx = d; arg = j; }
            }
            rep.Oks[k].PrevStepDeg = mx;
            rep.Oks[k].PrevStepAxis = arg;
            rep.Oks[k].Flip = mx > kFlip;
        }

        rep.StatLimit = rep.Oks.Count(x => x.AnyLimit);
        rep.StatJ5Sing = rep.Oks.Count(x => x.J5Sing);
        rep.StatStepWarn = rep.Oks.Count(x => x.PrevStepDeg > kStepWarn);
        rep.StatStepBad = rep.Oks.Count(x => x.PrevStepDeg > kStepBad);
        rep.StatFlip = rep.Oks.Count(x => x.Flip);
        rep.MaxStep = rep.Oks.Count == 0 ? 0 : rep.Oks.Max(x => x.PrevStepDeg);

        // 逐轴超限分布（归因用：集中一轴 = 系统性需求；散在多轴 = 分支选择问题）
        rep.LimitByAxis = new int[6];
        foreach (Solved s in rep.Oks)
            for (int j = 0; j < 6; j++)
                if (s.LimitMask[j]) rep.LimitByAxis[j]++;
    }

    // ================================================================ 报告

    public sealed class SolveReport
    {
        public int Count;
        public double ElapsedMs;
        public int MaxIter;
        public List<Solved> Oks = new();
        public List<Failed> Fails = new();

        public int StatLimit, StatJ5Sing, StatStepWarn, StatStepBad, StatFlip;
        public double MaxStep;
        public int[] LimitByAxis = new int[6];

        public void AddOk(int i, TrajPoint p, IkResult r, int tried) =>
            Oks.Add(new Solved { Index = i, Point = p, Result = r, Tried = tried });

        public void AddFail(int i, TrajPoint p, string why) =>
            Fails.Add(new Failed { Index = i, Point = p, Reason = why });

        /// <summary>能否走通 = 无失败、无超限、无奇异、无翻转、步长不超过危险阈值。</summary>
        public bool Ok => Fails.Count == 0 && StatLimit == 0 && StatJ5Sing == 0
                          && StatFlip == 0 && StatStepBad == 0;

        public void Print()
        {
            Console.WriteLine();
            Console.WriteLine("── 结果 ──");
            Console.WriteLine($"  点总数          {Count}");
            Console.WriteLine($"  逆解收敛        {Oks.Count}  ({100.0 * Oks.Count / Math.Max(1, Count):F1}%)");
            Console.WriteLine($"  逆解失败        {Fails.Count}");
            Console.WriteLine();
            Console.WriteLine("── 能否走通（五项硬指标）──");
            Line("① 逆解失败（摆不出来）", Fails.Count);
            Line("② 关节超软限位", StatLimit);
            Line("③ 腕部奇异 |θ5|<6°", StatJ5Sing);
            Line("④ 相邻步长 >5°", StatStepWarn);
            Line("④' 相邻步长 >30°（危险）", StatStepBad);
            Line("⑤ 解支翻转 Δ>120°（翻腕）", StatFlip);
            Console.WriteLine();
            Console.WriteLine($"  相邻步长 max    {MaxStep:F2}°");
            Console.WriteLine($"  单点迭代 max    {MaxIter} 步");
            Console.WriteLine($"  平均耗时        {ElapsedMs / Math.Max(1, Count):F2} ms/点（总 {ElapsedMs:F0} ms）");
            Console.WriteLine();
            Console.WriteLine(Ok ? "  ✅ 结论：整条轨迹可连续执行" : "  ❌ 结论：存在阻断项，见上表");
            Console.WriteLine();

            if (StatLimit > 0)
            {
                Console.WriteLine("  超限分布（归因用）：");
                for (int j = 0; j < 6; j++)
                    if (LimitByAxis[j] > 0)
                    {
                        double ex = 0;
                        foreach (Solved s in Oks)
                            if (s.LimitMask[j]) ex = Math.Max(ex, Exceed(s.Result.JointsDeg[j], j));
                        Console.WriteLine($"    J{j + 1}  {LimitByAxis[j],5} 点   最大超出 {ex:F2}°");
                    }
                Console.WriteLine();
            }

            if (Fails.Count > 0)
            {
                Console.WriteLine("  失败前 10 个：");
                foreach (Failed f in Fails.Take(10))
                    Console.WriteLine($"    #{f.Index,5}  {f.Point.Pos}  {f.Reason}");
                Console.WriteLine();
            }

            if (StatFlip > 0)
            {
                Console.WriteLine("  翻转前 10 处：");
                foreach (Solved s in Oks.Where(x => x.Flip).Take(10))
                    Console.WriteLine($"    #{s.Index,5}  J{s.PrevStepAxis + 1} 跳 {s.PrevStepDeg:F1}°");
                Console.WriteLine();
            }
        }

        private static void Line(string name, int n) =>
            Console.WriteLine($"  {(n == 0 ? "✓" : "✗")} {name,-28} {n,5}");

        private static double Exceed(double deg, int axis)
        {
            (double lo, double hi)[] lim =
            {
                (-170,170),(-100,135),(-133,150),(-270,270),(-120,120),(-360,360)
            };
            double c = (lim[axis].lo + lim[axis].hi) / 2;
            double cand = deg + 360 * Math.Round((c - deg) / 360);
            return Math.Max(0, Math.Max(lim[axis].lo - cand, cand - lim[axis].hi));
        }

        public void Save(string path, ProbeOptions opt, RobotModel model, ToolConfig tool, List<TrajPoint> traj)
        {
            var sb = new StringBuilder();
            sb.AppendLine("# DensoArmKine 逆解 → 喷涂轨迹离线体检");
            sb.AppendLine();
            sb.AppendLine($"> 生成时间 {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
            sb.AppendLine($"> **本报告全程离线**：未连接控制器、未下发任何运动指令。");
            sb.AppendLine();
            sb.AppendLine("## 输入");
            sb.AppendLine();
            sb.AppendLine($"- 轨迹 `{Path.GetFileName(opt.TrajPath)}`，{Count} 点");
            sb.AppendLine($"- 机型 `{model.Id}`，reach {model.SpecReachMm} mm");
            sb.AppendLine($"- 工具 Offset=({tool.OffsetX:F3}, {tool.OffsetY:F3}) LengthZ={tool.LengthZ:F3}");
            sb.AppendLine($"- 逆解：`InverseKinematics.Solve`（DLS，种子链，目标口径 TCP）");
            sb.AppendLine();
            sb.AppendLine("## 结论");
            sb.AppendLine();
            sb.AppendLine(Ok ? "**✅ 整条轨迹可连续执行**" : "**❌ 存在阻断项**");
            sb.AppendLine();
            sb.AppendLine("| 指标 | 数量 |");
            sb.AppendLine("|---|---|");
            sb.AppendLine($"| 点总数 | {Count} |");
            sb.AppendLine($"| 逆解收敛 | {Oks.Count} |");
            sb.AppendLine($"| 逆解失败 | {Fails.Count} |");
            sb.AppendLine($"| 关节超软限位 | {StatLimit} |");
            sb.AppendLine($"| 腕部奇异 \\|θ5\\|<6° | {StatJ5Sing} |");
            sb.AppendLine($"| 相邻步长 >5° | {StatStepWarn} |");
            sb.AppendLine($"| 相邻步长 >30° | {StatStepBad} |");
            sb.AppendLine($"| 解支翻转 >120° | {StatFlip} |");
            sb.AppendLine($"| 相邻步长 max | {MaxStep:F2}° |");
            sb.AppendLine();

            if (StatLimit > 0)
            {
                sb.AppendLine("## 超限分布");
                sb.AppendLine();
                sb.AppendLine("| 轴 | 点数 | 最大超出 |");
                sb.AppendLine("|---|---|---|");
                for (int j = 0; j < 6; j++)
                    if (LimitByAxis[j] > 0)
                    {
                        double ex = 0;
                        foreach (Solved s in Oks)
                            if (s.LimitMask[j]) ex = Math.Max(ex, Exceed(s.Result.JointsDeg[j], j));
                        sb.AppendLine($"| J{j + 1} | {LimitByAxis[j]} | {ex:F2}° |");
                    }
                sb.AppendLine();
            }

            sb.AppendLine("## 逐点明细");
            sb.AppendLine();
            sb.AppendLine("| # | X | Y | Z | J1 | J2 | J3 | J4 | J5 | J6 | 步长 | 标记 |");
            sb.AppendLine("|---|---|---|---|---|---|---|---|---|---|---|---|");
            foreach (Solved s in Oks)
            {
                double[] t = s.Result.JointsDeg;
                var marks = new List<string>();
                if (s.AnyLimit) marks.Add("超限");
                if (s.J5Sing) marks.Add("奇异");
                if (s.Flip) marks.Add("翻腕");
                if (s.PrevStepDeg > 30) marks.Add("跳变");
                sb.AppendLine($"| {s.Index} | {s.Point.Pos.X:F2} | {s.Point.Pos.Y:F2} | {s.Point.Pos.Z:F2} "
                              + $"| {t[0]:F2} | {t[1]:F2} | {t[2]:F2} | {t[3]:F2} | {t[4]:F2} | {t[5]:F2} "
                              + $"| {s.PrevStepDeg:F2} | {string.Join(" ", marks)} |");
            }

            File.WriteAllText(path, sb.ToString(), new UTF8Encoding(true));
        }
    }

    public sealed class Solved
    {
        public int Index;
        public TrajPoint Point;
        public IkResult Result;
        public int Tried;
        public bool[] LimitMask;
        public bool AnyLimit;
        public double Theta5;
        public bool J5Sing;
        public double PrevStepDeg;
        public int PrevStepAxis;
        public bool Flip;
    }

    public sealed class Failed
    {
        public int Index;
        public TrajPoint Point;
        public string Reason;
    }
}

/// <summary>一个轨迹点。位置 + 姿态（欧拉角或枪轴二选一）+ 可选的已知关节角。</summary>
public sealed class TrajPoint
{
    public int Index;
    public Vec3 Pos;
    public Vec3 Euler;          // DENSO 约定 R = Rz·Ry·Rx，度
    public Vec3 GunAxis;        // 枪轴方向（单位）
    public bool HasEuler;
    public double[] Joints;     // 已知关节角（来自 _safe.csv），可为 null
    public bool HasJoints => Joints != null && Joints.Length >= 6;
}

/// <summary>命令行选项。</summary>
public sealed class ProbeOptions
{
    public string TrajPath;
    public string CalibPath;
    public string OutPath;
    public string ModelId = "VS-087";
    public bool ToolExplicit;
    public double ToolOffX, ToolOffY, ToolLen;
    public bool Help;

    public static ProbeOptions Parse(string[] a)
    {
        var o = new ProbeOptions();
        for (int i = 0; i < a.Length; i++)
        {
            switch (a[i])
            {
                case "-h": case "--help": o.Help = true; break;
                case "-o": case "--out": o.OutPath = Next(a, ref i); break;
                case "--calib": o.CalibPath = Next(a, ref i); break;
                case "--model": o.ModelId = Next(a, ref i); break;
                case "--tool-len":
                    o.ToolLen = double.Parse(Next(a, ref i), CultureInfo.InvariantCulture);
                    o.ToolExplicit = true; break;
                case "--tool-off":
                    {
                        string[] s = Next(a, ref i).Split(',');
                        o.ToolOffX = double.Parse(s[0], CultureInfo.InvariantCulture);
                        o.ToolOffY = double.Parse(s[1], CultureInfo.InvariantCulture);
                        o.ToolExplicit = true; break;
                    }
                default:
                    if (!a[i].StartsWith("-")) o.TrajPath = a[i];
                    break;
            }
        }
        return o;
    }

    private static string Next(string[] a, ref int i) => i + 1 < a.Length ? a[++i] : null;

    public static void PrintUsage()
    {
        Console.WriteLine("""
        DensoIkProbe —— 用主城序（DensoArmKine）的 DLS 逆解对喷涂轨迹做离线体检。
        ⚠ 全程离线：不连控制器、不下发运动。

        用法：
          DensoIkProbe <轨迹文件> [选项]

        轨迹文件：*.plan.json（推荐，含 DENSO 欧拉角）或 *_traj.csv / *_safe.csv

        选项：
          -o, --out <路径>       把详细报告写成 Markdown
              --calib <路径>     用含关节角的 csv 反查工具参数（自校验）
              --model <机型id>   默认 VS-087
              --tool-len <mm>    工具长度（默认 450.756，由本算例反推）
              --tool-off <x,y>   工具横向偏置（默认 0.038,2.214）
          -h, --help             显示本帮助

        返回码：0 = 整条可走通；1 = 存在阻断项；2 = 输入错误
        """);
    }
}

/// <summary>
/// 轨迹加载。支持两种输入：
///   ① <c>*.plan.json</c>：Points[].Cartesian = [x,y,z,rx,ry,rz,extra]，直接有欧拉角。
///   ② <c>*.csv</c>：列见 <see cref="LoadCsv"/>；若含 J1_deg..J6_deg 则同时带关节角。
/// </summary>
public static class TrajLoader
{
    public static List<TrajPoint> Load(string path)
    {
        string ext = Path.GetExtension(path).ToLowerInvariant();
        if (ext == ".json") return LoadPlanJson(path);
        if (ext == ".csv") return LoadCsv(path);
        // .ply 是点云（碰撞检测用），不是轨迹。这里明确拒绝，避免"喂错文件却算出貌似有理的数"。
        throw new NotSupportedException(
            $"{ext} 不是轨迹格式。轨迹在 *.plan.json / *_traj.csv / *_safe.csv 里；.ply 是点云。");
    }

    private static List<TrajPoint> LoadPlanJson(string path)
    {
        var list = new List<TrajPoint>();
        using JsonDocument doc = JsonDocument.Parse(File.ReadAllText(path));
        if (!doc.RootElement.TryGetProperty("Points", out JsonElement pts)) return list;

        int idx = 0;
        foreach (JsonElement p in pts.EnumerateArray())
        {
            if (!p.TryGetProperty("Cartesian", out JsonElement c)) continue;
            double[] v = c.EnumerateArray().Select(x => x.GetDouble()).ToArray();
            if (v.Length < 6) continue;

            var tp = new TrajPoint
            {
                Index = idx++,
                Pos = new Vec3(v[0], v[1], v[2]),
                Euler = new Vec3(v[3], v[4], v[5]),
                HasEuler = true,
            };
            if (p.TryGetProperty("Joints", out JsonElement j) && j.ValueKind == JsonValueKind.Array)
            {
                double[] jj = j.EnumerateArray().Select(x => x.GetDouble()).ToArray();
                if (jj.Length >= 6) tp.Joints = jj;
            }
            list.Add(tp);
        }
        return list;
    }

    private static List<TrajPoint> LoadCsv(string path)
    {
        string[] lines = File.ReadAllLines(path);
        if (lines.Length < 2) return new List<TrajPoint>();

        string[] head = lines[0].Replace("\uFEFF", "").Split(',');
        var col = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
        for (int i = 0; i < head.Length; i++) col[head[i].Trim()] = i;

        int cX = Col(col, "X"), cY = Col(col, "Y"), cZ = Col(col, "Z");
        int cI = Col(col, "I"), cJ = Col(col, "J"), cK = Col(col, "K");
        if (cX < 0 || cY < 0 || cZ < 0) throw new FormatException("csv 缺少 X/Y/Z 列");

        var jc = new int[6];
        bool hasJ = true;
        for (int k = 0; k < 6; k++)
        {
            jc[k] = Col(col, $"J{k + 1}_deg");
            if (jc[k] < 0) { hasJ = false; break; }
        }

        var list = new List<TrajPoint>();
        for (int ln = 1; ln < lines.Length; ln++)
        {
            if (string.IsNullOrWhiteSpace(lines[ln])) continue;
            string[] f = lines[ln].Split(',');
            if (f.Length < head.Length) continue;

            var tp = new TrajPoint
            {
                Index = list.Count,
                Pos = new Vec3(D(f[cX]), D(f[cY]), D(f[cZ])),
            };
            if (cI >= 0 && cJ >= 0 && cK >= 0)
                tp.GunAxis = new Vec3(D(f[cI]), D(f[cJ]), D(f[cK])).Normalized();
            if (hasJ) tp.Joints = Enumerable.Range(0, 6).Select(k => D(f[jc[k]])).ToArray();
            list.Add(tp);
        }
        return list;

        static int Col(Dictionary<string, int> c, string n) => c.TryGetValue(n, out int i) ? i : -1;
        static double D(string s) => double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);
    }
}
