using System;
using System.Collections.Generic;
using System.IO;

namespace SprayTrajViewer;

/// <summary>
/// 界面状态的轻量持久化（key=value 纯文本，存 %LOCALAPPDATA%\SprayTrajViewer\ui.txt）。
///
/// 存在的理由：WPF 在 .NET Core/5+ 上不再沿用系统给上一个进程记的浏览目录，
/// OpenFileDialog 不给 InitialDirectory 就固定落回"文档"，每次开都要重新点一遍路径。
/// 目前只存 lastCloudDir（上次打开的 PLY 所在目录），后续窗口尺寸/参数也可挂这里，
/// 所以做成通用的 key-value 而不是单个文本文件。
/// </summary>
static class UiState
{
    static readonly string StateDir = Path.Combine(
        Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "SprayTrajViewer");

    static string StateFile => Path.Combine(StateDir, "ui.txt");

    static Dictionary<string, string> Load()
    {
        var d = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        try
        {
            if (!File.Exists(StateFile)) return d;
            foreach (var line in File.ReadAllLines(StateFile))
            {
                int eq = line.IndexOf('=');
                if (eq <= 0) continue;
                d[line.Substring(0, eq).Trim()] = line.Substring(eq + 1).Trim();
            }
        }
        catch { /* 状态文件坏掉不该影响主流程，按空处理 */ }
        return d;
    }

    public static string? Get(string key)
        => Load().TryGetValue(key, out var v) && v.Length > 0 ? v : null;

    /// <summary>读-改-写，保留其它键（不会因为加一个键把已有记忆清掉）。</summary>
    public static void Set(string key, string value)
    {
        try
        {
            var d = Load();
            d[key] = value;
            Directory.CreateDirectory(StateDir);
            var lines = new List<string>(d.Count);
            foreach (var kv in d) lines.Add(kv.Key + "=" + kv.Value);
            File.WriteAllLines(StateFile, lines);
        }
        catch { }
    }
}
