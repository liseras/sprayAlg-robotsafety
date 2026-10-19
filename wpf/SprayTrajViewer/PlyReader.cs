using System;
using System.Buffers.Text;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace SprayTrajViewer
{
    // 简单 PLY 点云读取（ascii + binary_little_endian，含 x/y/z）
    //
    // 性能说明（大点云加载慢的根因）：
    //   旧实现用 BinaryReader.Read() 一个字节一个字节地读行，60 万行 × 50 字节
    //   就是 3000 万次虚调用；每行还要 Split 出 string[] 再 double.Parse。
    //   现场大平面件(60 万点)光读文件就要好几秒，比算法本身还慢。
    //
    // 现在：整文件一次读进 byte[]，ASCII 用 Utf8Parser.TryParse 直接在字节上
    // 解 double（零分配、且是正确舍入的，与 double.Parse 结果一致）。
    public static class PlyReader
    {
        public static List<(double x, double y, double z)> Read(string path)
        {
            byte[] buf = File.ReadAllBytes(path);
            int vertexCount = 0;
            bool ascii = true;
            var props = new List<(string type, string name)>();   // 顶点属性（按声明顺序）
            bool inVertex = false;
            int dataStart = buf.Length;

            // ---- 头部：行数很少，用字符串处理即可 ----
            int p = 0;
            while (p < buf.Length)
            {
                int eol = p;
                while (eol < buf.Length && buf[eol] != (byte)'\n') eol++;
                if (eol > p && buf[eol - 1] == (byte)'\r') eol--;
                string line = Encoding.ASCII.GetString(buf, p, eol - p);
                int next = eol;
                // 跳过 '\r'
                while (next < buf.Length && buf[next] != (byte)'\n') next++;
                next++;

                if (line.StartsWith("format")) ascii = line.Contains("ascii");
                else if (line.StartsWith("element "))
                {
                    inVertex = line.StartsWith("element vertex");
                    if (inVertex)
                    {
                        var t = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                        if (t.Length >= 3) vertexCount = int.Parse(t[t.Length - 1]);
                    }
                }
                else if (line.StartsWith("property") && inVertex)
                {
                    var t = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    if (t.Length >= 3) props.Add((t[1], t[2]));
                }
                else if (line.StartsWith("end_header"))
                {
                    dataStart = next;
                    break;
                }
                p = next;
            }

            var pts = new List<(double, double, double)>(Math.Max(0, vertexCount));

            if (ascii)
            {
                var span = new ReadOnlySpan<byte>(buf);
                int i = dataStart;
                int n = buf.Length;
                for (int k = 0; k < vertexCount && i < n; k++)
                {
                    if (!Utf8Parser.TryParse(span.Slice(i), out double x, out int u1)) break;
                    i += u1; SkipBlanks(span, ref i);
                    if (!Utf8Parser.TryParse(span.Slice(i), out double y, out int u2)) break;
                    i += u2; SkipBlanks(span, ref i);
                    if (!Utf8Parser.TryParse(span.Slice(i), out double z, out int u3)) break;
                    i += u3;
                    pts.Add((x, y, z));
                    // 跳过该行剩余列（强度等）
                    while (i < n && span[i] != (byte)'\n') i++;
                    i++;
                }
            }
            else
            {
                if (props.Count == 0) return pts;
                // 记录 x/y/z 在每点字节块里的偏移与类型。
                // 注意：所有属性都要计入步长——旧实现把未识别的属性跳过不读，
                // 一旦文件里多一列（如 intensity），后面所有点都会错位。
                int stride = 0;
                var off = new Dictionary<string, (int off, string type)>();
                foreach (var (type, name) in props)
                {
                    if (name == "x" || name == "y" || name == "z") off[name] = (stride, type);
                    stride += SizeOf(type);
                }
                for (int k = 0; k < vertexCount; k++)
                {
                    int baseOff = dataStart + k * stride;
                    if (baseOff + stride > buf.Length) break;
                    double x = ReadField(buf, baseOff, off, "x");
                    double y = ReadField(buf, baseOff, off, "y");
                    double z = ReadField(buf, baseOff, off, "z");
                    pts.Add((x, y, z));
                }
            }
            return pts;
        }

        static void SkipBlanks(ReadOnlySpan<byte> s, ref int i)
        {
            while (i < s.Length && (s[i] == (byte)' ' || s[i] == (byte)'\t')) i++;
        }

        static int SizeOf(string type) => type switch
        {
            "char" or "uchar" or "int8" or "uint8" => 1,
            "short" or "ushort" or "int16" or "uint16" => 2,
            "int" or "uint" or "int32" or "uint32" or "float" or "float32" => 4,
            "double" or "float64" => 8,
            _ => 1,
        };

        static double ReadField(byte[] b, int baseOff, Dictionary<string, (int off, string type)> off, string name)
        {
            if (!off.TryGetValue(name, out var f)) return 0;
            int o = baseOff + f.off;
            return f.type switch
            {
                "float" or "float32" => BitConverter.ToSingle(b, o),
                "double" or "float64" => BitConverter.ToDouble(b, o),
                "int" or "int32" => BitConverter.ToInt32(b, o),
                "uint" or "uint32" => BitConverter.ToUInt32(b, o),
                "short" or "int16" => BitConverter.ToInt16(b, o),
                "ushort" or "uint16" => BitConverter.ToUInt16(b, o),
                "char" or "int8" => (sbyte)b[o],
                "uchar" or "uint8" => b[o],
                _ => 0,
            };
        }
    }
}
