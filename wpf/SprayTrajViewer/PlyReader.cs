using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace SprayTrajViewer
{
    // 简单 PLY 点云读取（ascii + binary_little_endian，含 x/y/z）
    public static class PlyReader
    {
        public static List<(double x, double y, double z)> Read(string path)
        {
            var pts = new List<(double, double, double)>();
            using var fs = File.OpenRead(path);
            using var br = new BinaryReader(fs);
            string line;
            int vertexCount = 0;
            bool ascii = true;
            var props = new List<(string type, string name)>();

            while ((line = ReadLine(br)) != null)
            {
                if (line.StartsWith("format")) ascii = line.Contains("ascii");
                else if (line.StartsWith("element vertex"))
                    vertexCount = int.Parse(line.Split(' ').Last());
                else if (line.StartsWith("property"))
                {
                    var t = line.Split(' ');
                    props.Add((t[1], t[2]));
                }
                else if (line.StartsWith("end_header")) break;
            }

            if (ascii)
            {
                for (int i = 0; i < vertexCount; i++)
                {
                    var l = ReadLine(br);
                    if (l == null) break;
                    var v = l.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    if (v.Length >= 3)
                        pts.Add((double.Parse(v[0]), double.Parse(v[1]), double.Parse(v[2])));
                }
            }
            else
            {
                // binary_little_endian，按 property 顺序读，取 x/y/z
                var fmt = new List<(string type, string name)>();
                foreach (var (type, name) in props)
                {
                    if (name == "x" || name == "y" || name == "z" || name == "nx" || name == "ny" || name == "nz"
                        || name == "red" || name == "green" || name == "blue")
                        fmt.Add((type, name));
                }
                for (int i = 0; i < vertexCount; i++)
                {
                    double x = 0, y = 0, z = 0;
                    foreach (var (type, name) in fmt)
                    {
                        if (type == "float")
                        {
                            float val = br.ReadSingle();
                            if (name == "x") x = val;
                            else if (name == "y") y = val;
                            else if (name == "z") z = val;
                        }
                        else if (type == "uchar") br.ReadByte();
                    }
                    pts.Add((x, y, z));
                }
            }
            return pts;
        }

        static string? ReadLine(BinaryReader br)
        {
            var bytes = new List<byte>();
            int b;
            while ((b = br.Read()) != -1)
            {
                if (b == '\n') break;
                if (b != '\r') bytes.Add((byte)b);
            }
            if (bytes.Count == 0 && b == -1) return null;
            return System.Text.Encoding.ASCII.GetString(bytes.ToArray());
        }
    }
}
