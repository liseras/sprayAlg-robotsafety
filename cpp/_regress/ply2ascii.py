# -*- coding: utf-8 -*-
"""PLY -> ASCII PLY 转换（driver 只吃 ASCII）。

支持 ascii / binary_little_endian / binary_big_endian，只保留 x,y,z
（其余属性丢掉，driver 也不用）。
用法: python ply2ascii.py 输入.ply 输出.ply
"""
import sys, struct

T = {
    "char": ("b", 1), "uchar": ("B", 1), "int8": ("b", 1), "uint8": ("B", 1),
    "short": ("h", 2), "ushort": ("H", 2), "int16": ("h", 2), "uint16": ("H", 2),
    "int": ("i", 4), "uint": ("I", 4), "int32": ("i", 4), "uint32": ("I", 4),
    "float": ("f", 4), "float32": ("f", 4), "double": ("d", 8), "float64": ("d", 8),
}

def main():
    src, dst = sys.argv[1], sys.argv[2]
    f = open(src, "rb")
    fmt, nvert = "ascii", 0
    props = []
    while True:
        line = f.readline().decode("ascii", "ignore").strip()
        if not line:
            continue
        t = line.split()
        if t[0] == "format":
            fmt = t[1]
        elif t[0] == "element":
            if t[1] == "vertex":
                nvert = int(t[2])
                props = []
            else:
                break
        elif t[0] == "property" and nvert:
            props.append((t[2], t[1]))
        elif t[0] == "end_header":
            break
    # 跳过非 vertex 的 element
    if fmt == "ascii":
        pts = []
        for _ in range(nvert):
            v = f.readline().split()
            pts.append((float(v[0]), float(v[1]), float(v[2])))
    else:
        en = "<" if fmt == "binary_little_endian" else ">"
        names = [p[0] for p in props]
        idx = {}
        codes = []
        for nm, ty in props:
            c, sz = T[ty]
            codes.append(c)
            idx[nm] = len(codes) - 1
        st = struct.Struct(en + "".join(codes))
        ix, iy, iz = idx.get("x", 0), idx.get("y", 1), idx.get("z", 2)
        buf = f.read(st.size * nvert)
        pts = []
        for k in range(nvert):
            v = st.unpack_from(buf, k * st.size)
            pts.append((v[ix], v[iy], v[iz]))
    f.close()
    with open(dst, "w", newline="\n") as o:
        o.write("ply\nformat ascii 1.0\nelement vertex %d\n" % len(pts))
        o.write("property float x\nproperty float y\nproperty float z\nend_header\n")
        for p in pts:
            o.write("%.4f %.4f %.4f\n" % p)
    print("wrote", dst, len(pts), "pts (from", fmt + ")")

main()
