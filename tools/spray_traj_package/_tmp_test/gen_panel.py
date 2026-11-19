#!/usr/bin/env python3
# 生成回归用"扁平类工件"点云：平板 + 两个大孔 + 右侧斜面（喷粉工艺的典型件）
# 用法: python gen_panel.py out.ply [pass_spacing_hint]
import math, random, sys

def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "panel.ply"
    W, H = 620.0, 360.0          # 板尺寸 mm
    step = 3.0                   # 扫描点距
    holes = [(-150.0, 60.0, 55.0), (60.0, -50.0, 45.0)]
    slope_u0, slope_deg = 215.0, 35.0
    t = math.tan(math.radians(slope_deg))
    random.seed(7)
    pts = []
    nu = int(W / step) + 1
    nv = int(H / step) + 1
    for j in range(nv):
        for i in range(nu):
            u = -W / 2 + i * step
            v = -H / 2 + j * step
            ok = True
            for (cu, cv, r) in holes:
                if math.hypot(u - cu, v - cv) < r:
                    ok = False
                    break
            if not ok:
                continue
            z = 0.0
            if u > slope_u0:                      # 右侧斜面（抬起）
                z = (u - slope_u0) * t
            # 少量起伏 + 噪声（模拟真实扫描)
            z += 0.8 * math.sin(u / 85.0) * math.cos(v / 70.0)
            z += random.gauss(0.0, 0.25)
            pts.append((u + random.gauss(0, 0.15),
                        v + random.gauss(0, 0.15),
                        z))
    with open(out, "w", encoding="ascii") as f:
        f.write("ply\nformat ascii 1.0\n")
        f.write("comment generated panel with 2 holes + slope\n")
        f.write("element vertex %d\n" % len(pts))
        f.write("property float x\nproperty float y\nproperty float z\n")
        f.write("end_header\n")
        for p in pts:
            f.write("%.3f %.3f %.3f\n" % p)
    print("wrote", out, len(pts), "pts")

main()
