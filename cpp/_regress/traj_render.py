# -*- coding: utf-8 -*-
"""把点云 + 轨迹渲染成 PNG（纯标准库，无 matplotlib），用来肉眼核对轨迹形态。

用法:
  python traj_render.py box.ply _reg.csv out.png [--view 1] [--w 760] [--h 560]

投影与查看器一致用**等轴测**风格：方位角 az、仰角 el，正交投影后按包围盒缩放。
点云画成暗灰小点（稀疏采样），轨迹按 spray 分色（1=红 喷，0=黄 空移）。
"""
import sys, zlib, struct, math

def read_ply(path):
    f = open(path, "rb")
    n = 0
    fmt = "ascii"
    while True:
        line = f.readline().decode("ascii", "ignore")
        if line.startswith("format"):
            fmt = line.split()[1]
        if line.startswith("element vertex"):
            n = int(line.split()[-1])
        if line.strip() == "end_header":
            break
    if fmt != "ascii":
        raise SystemExit("仅支持 ascii PLY")
    pts = []
    for _ in range(n):
        v = f.readline().split()
        pts.append((float(v[0]), float(v[1]), float(v[2])))
    f.close()
    return pts

def read_csv(path):
    pts = []
    for i, line in enumerate(open(path)):
        if i == 0:
            continue
        v = line.strip().split(",")
        if len(v) < 8:
            continue
        pts.append((float(v[0]), float(v[1]), float(v[2]), int(float(v[7]))))
    return pts

def write_png(path, w, h, pix):
    raw = b"".join(b"\x00" + bytes(pix[y * w * 3:(y + 1) * w * 3]) for y in range(h))
    def chunk(tag, data):
        c = struct.pack(">I", len(data)) + tag + data
        return c + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)
    png = b"\x89PNG\r\n\x1a\n"
    png += chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
    png += chunk(b"IDAT", zlib.compress(raw, 6))
    png += chunk(b"IEND", b"")
    open(path, "wb").write(png)

def main():
    ply, csv, out = sys.argv[1], sys.argv[2], sys.argv[3]
    W = int(sys.argv[sys.argv.index("--w") + 1]) if "--w" in sys.argv else 760
    H = int(sys.argv[sys.argv.index("--h") + 1]) if "--h" in sys.argv else 560
    az = math.radians(float(sys.argv[sys.argv.index("--az") + 1])) if "--az" in sys.argv else math.radians(-28.0)
    el = math.radians(float(sys.argv[sys.argv.index("--el") + 1])) if "--el" in sys.argv else math.radians(24.0)

    cloud = read_ply(ply)
    traj = read_csv(csv)
    allp = [p[:3] for p in cloud] + [p[:3] for p in traj]
    mn = [min(p[i] for p in allp) for i in range(3)]
    mx = [max(p[i] for p in allp) for i in range(3)]
    ctr = [(mn[i] + mx[i]) * 0.5 for i in range(3)]
    rad = max(max(mx[i] - mn[i] for i in range(3)), 1e-6)

    ca, sa = math.cos(az), math.sin(az)
    ce, se = math.cos(el), math.sin(el)
    def proj(p):
        x, y, z = p[0] - ctr[0], p[1] - ctr[1], p[2] - ctr[2]
        # 先绕 Z 转 az（分开 X/Y 两个方向），再按仰角 el 把 Y 压扁、Z 抬正
        xr = ca * x - sa * y
        yr = sa * x + ca * y
        return xr, yr * se - z * ce
    sc = (min(W, H) * 0.82) / rad
    def to_px(p):
        sx, sy = proj(p)
        return int(W * 0.5 + sx * sc), int(H * 0.5 - sy * sc)

    pix = bytearray(W * H * 3)
    for i in range(0, W * H * 3, 3):
        pix[i], pix[i + 1], pix[i + 2] = 16, 16, 20
    def dot(x, y, r, g, b, rad=0):
        for dy in range(-rad, rad + 1):
            for dx in range(-rad, rad + 1):
                xx, yy = x + dx, y + dy
                if 0 <= xx < W and 0 <= yy < H:
                    k = (yy * W + xx) * 3
                    pix[k], pix[k + 1], pix[k + 2] = r, g, b

    step = max(1, len(cloud) // 40000)
    for p in cloud[::step]:
        x, y = to_px(p)
        dot(x, y, 70, 70, 78)
    for p in traj:
        x, y = to_px(p)
        if p[3] == 1:
            dot(x, y, 235, 60, 55, 1)
        else:
            dot(x, y, 210, 190, 60, 0)
    # 轨迹折线
    prev = None
    for p in traj:
        cur = to_px(p)
        if prev is not None and (p[3] == 1):
            x0, y0 = prev
            n = max(abs(cur[0] - x0), abs(cur[1] - y0))
            for k in range(1, n + 1):
                t = k / n
                dot(int(x0 + (cur[0] - x0) * t), int(y0 + (cur[1] - y0) * t), 255, 70, 64)
        prev = cur
    write_png(out, W, H, pix)
    print("wrote", out, "traj pts", len(traj))

main()
