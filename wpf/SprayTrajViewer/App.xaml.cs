using System;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Threading;

namespace SprayTrajViewer;

/// <summary>
/// 应用入口。除了常规启动，负责两件事：
///   1) 记录全局未处理异常到日志文件（WinExe 没有控制台，崩溃时否则一片空白，无法定位）；
///   2) 把启动参数原样落到日志，配合命令行自运行模式可复现界面验证。
/// </summary>
public partial class App : Application
{
    /// <summary>诊断日志路径（%TEMP%\SprayTrajViewer.log）。同时写 stdout，便于重定向捕获。</summary>
    public static readonly string LogPath =
        Path.Combine(Path.GetTempPath(), "SprayTrajViewer.log");

    protected override void OnStartup(StartupEventArgs e)
    {
        // 每次启动清空日志，避免历史内容干扰本次诊断
        try { File.WriteAllText(LogPath, ""); } catch { }

        AppDomain.CurrentDomain.UnhandledException += (_, a) =>
            Log("FATAL/AppDomain", a.ExceptionObject as Exception, "IsTerminating=" + a.IsTerminating);
        DispatcherUnhandledException += (_, a) =>
        {
            Log("FATAL/Dispatcher", a.Exception);
            // 记录后仍交给系统处理，不掩盖真实崩溃
            a.Handled = false;
        };
        TaskScheduler.UnobservedTaskException += (_, a) =>
            Log("FATAL/Task", a.Exception.GetBaseException());

        Log("BOOT", null, "args: " + string.Join(" ", e.Args));
        Log("BOOT", null, "base=" + AppContext.BaseDirectory);
        Log("BOOT", null, "runtime=" + Environment.Version);

        base.OnStartup(e);
    }

    protected override void OnExit(ExitEventArgs e)
    {
        Log("EXIT", null, "code=" + e.ApplicationExitCode);
        base.OnExit(e);
    }

    /// <summary>写一行诊断日志：同时输出到 stdout（可被重定向）与日志文件。</summary>
    public static void Log(string tag, Exception? ex = null, string? msg = null)
    {
        var sb = new StringBuilder();
        sb.Append(DateTime.Now.ToString("HH:mm:ss.fff")).Append(" [").Append(tag).Append(']');
        if (!string.IsNullOrEmpty(msg)) sb.Append(' ').Append(msg);
        if (ex != null)
        {
            sb.Append(' ').Append(ex.GetType().FullName).Append(": ").Append(ex.Message);
            if (!string.IsNullOrEmpty(ex.StackTrace)) sb.Append(" @ ").Append(ex.StackTrace);
        }
        string line = sb.ToString();

        try { Console.WriteLine(line); Console.Out.Flush(); } catch { }
        try { File.AppendAllText(LogPath, line + Environment.NewLine); } catch { }
    }
}
