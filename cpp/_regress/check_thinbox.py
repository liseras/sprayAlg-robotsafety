# -*- coding: utf-8 -*-
"""薄壁开口盒轨迹体检：覆盖范围 + 是否穿进材料。

用法: python check_thinbox.py <traj.csv> [...]
盒体: 400x300x300，壁厚 4，顶面开口。
判据:
  · 材料 = 外盒 [0,400]x[0,300]x[0,300] 减去 空腔 [4,396]x[4,296]x[4,300]
  · 轨迹点（含相邻点连线上的采样点）若落在"材料膨胀 1mm"之内 → 计入穿材料
"""
import sys
import csv

W, D, H, T = 400.0, 300.0, 300.0, 4.0
TOL = 1.0   # 允许贴面 1mm


def in_material(x, y, z):
    # 外盒外 → 不算
    if not (-TOL <= x <= W + TOL and -TOL <= y <= D + TOL and -TOL <= z <= H + TOL):
        return False
    # 空腔（含向上开口）→ 不算
    if (T - TOL <= x <= W - T + TOL and T - TOL <= y <= D - T + TOL and z >= T - TOL):
        return False
    return True


def main(path):
    pts = []
    with open(path) as f:
        for row in csv.DictReader(f):
            pts.append((float(row["X"]), float(row["Y"]), float(row["Z"]),
                        int(row["spray"])))
    bad = []
    for a, b in zip(pts, pts[1:]):
        d = ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5
        n = max(1, int(d / 2.0))
        for i in range(n + 1):
            t = i / n
            p = (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t)
            if in_material(*p):
                bad.append((p, a[3], d))
                break
    print("%-28s 点=%-6d 段穿材料=%-4d" % (path.split("/")[-1], len(pts), len(bad)))
    for p, sp, d in bad[:5]:
        print("   穿入 (%.1f,%.1f,%.1f) spray=%d 段长%.1f" % (p[0], p[1], p[2], sp, d))


if __name__ == "__main__":
    for p in sys.argv[1:]:
        main(p)
