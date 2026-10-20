# -*- coding: utf-8 -*-
"""按面统计：点数 / 行数 / 行距 / 覆盖率 / 枪位贴面性
   用 plane*.ply 的 4 顶点反推面局部基，只留"距该面 standoff 且落在外接矩形内"的轨迹点。
   行数用 b 轴 1mm 直方图峰值计数（蛇形转弯会填满 b 的中间值，不能按间隔聚类）。
"""
import numpy as np, glob, os, sys

pre = sys.argv[1]
standoff = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
csvp = os.path.splitext(pre)[0] + ".csv"


def rdply(p):
    L = open(p, errors="replace").read().split("\n")
    nv = 0
    for i, l in enumerate(L):
        if l.startswith("element vertex"):
            nv = int(l.split()[2])
        elif l.strip() == "end_header":
            h = i + 1
            break
    return np.array([[float(x) for x in L[h + k].split()[:3]] for k in range(nv)])


T = np.loadtxt(csvp, delimiter=",", skiprows=1)
P = T[:, :3]
meshes = sorted(glob.glob(pre + ".plane*.ply"),
                key=lambda s: int(s.rsplit("plane", 1)[1].split(".")[0]))
print("%-22s %6s %5s %8s %7s %8s" % ("face", "pts", "rows", "spacing", "覆盖", "枪位sd"))
for mp in meshes:
    V = rdply(mp)
    c0 = V.mean(0)
    _, S, vt = np.linalg.svd(V - c0, full_matrices=False)
    n = vt[2]
    u = vt[0] - (vt[0] @ n) * n; u /= np.linalg.norm(u)
    v = np.cross(n, u)
    U = (V - c0) @ u; W = (V - c0) @ v
    lu, lv = U.max() - U.min(), W.max() - W.min()
    c = c0 + u * (U.max() + U.min()) / 2 + v * (W.max() + W.min()) / 2
    Q = P - c
    hh = Q @ n
    uu = Q @ u; vv = Q @ v
    m = (np.abs(np.abs(hh) - standoff) < 1.5) & (np.abs(uu) <= lu / 2 + 2) & (np.abs(vv) <= lv / 2 + 2)
    if m.sum() < 20:
        print("%-22s %6d  (无轨迹点)" % (os.path.basename(mp), m.sum()))
        continue
    a, b = uu[m], vv[m]

    def peaks_of(x):
        """1mm 直方图 + 局部极大（间距 >=5mm，高度 >=0.3 峰值）→ 行中心"""
        edges = np.arange(np.floor(x.min()) - 1, np.ceil(x.max()) + 2, 1.0)
        h, e = np.histogram(x, bins=edges)
        if h.max() <= 0:
            return np.array([])
        cut = 0.30 * h.max()
        c = [(e[k] + e[k + 1]) / 2 for k in range(1, len(h) - 1)
             if h[k] >= h[k - 1] and h[k] >= h[k + 1] and h[k] >= cut]
        if not c:
            return np.array([])
        keep = [c[0]]
        for t in c[1:]:
            if t - keep[-1] < 5.0:
                keep[-1] = t if abs(t) > abs(keep[-1]) else keep[-1]
            else:
                keep.append(t)
        return np.array(keep)

    pa, pb = peaks_of(a), peaks_of(b)
    pk = pb if len(pb) >= len(pa) else pa
    d = np.diff(pk)
    print("%-22s %6d %5d %8.2f %7d %8.3f" % (
        os.path.basename(mp), m.sum(), len(pk), np.median(d) if len(d) else 0,
        m.sum(), np.std(np.abs(hh[m]))))
