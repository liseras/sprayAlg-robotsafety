# -*- coding: utf-8 -*-
"""重建曲面网格平滑度体检（第十三轮新增）

用法:
    python mesh_smoothness.py <planeN.ply>

输出相邻三角面二面角、相邻顶点法向夹角的分布。判定口径：
真实曲面的二面角 ≈ 格距 / 曲率半径。曲型板拱面 R≈400mm、格距 18.2mm → 期望 ~2.6°。
p95 明显大于这个值 = 网格有棱（曲率不连续），不是"曲率大"。
"""
import sys
import numpy as np
from collections import defaultdict


def load_ply(path):
    V, F = [], []
    nv = 0
    cnt = 0
    mode = None
    with open(path) as f:
        for ln in f:
            a = ln.split()
            if not a:
                continue
            if a[0] == "element" and a[1] == "vertex":
                nv = int(a[2])
            if a[0] == "end_header":
                mode = "v"
                continue
            if mode == "v":
                if cnt < nv:
                    V.append([float(x) for x in a[:3]])
                    cnt += 1
                else:
                    F.append([int(a[1]), int(a[2]), int(a[3])])
    return np.array(V), np.array(F)


def report(name, arr, unit="deg"):
    print("  %-22s n=%4d  p50=%6.2f  p90=%6.2f  p95=%6.2f  p99=%6.2f  max=%6.2f %s"
          % (name, len(arr), np.percentile(arr, 50), np.percentile(arr, 90),
             np.percentile(arr, 95), np.percentile(arr, 99), arr.max(), unit))


def main(path):
    V, F = load_ply(path)
    print("mesh: verts=%d tris=%d  (%s)" % (len(V), len(F), path))
    e1 = V[F[:, 1]] - V[F[:, 0]]
    e2 = V[F[:, 2]] - V[F[:, 0]]
    Nf = np.cross(e1, e2)
    area = 0.5 * np.linalg.norm(Nf, axis=1)
    Nf /= np.maximum(np.linalg.norm(Nf, axis=1), 1e-12)[:, None]

    edges = defaultdict(list)
    for fi, (a, b, c) in enumerate(F):
        for u, v in ((a, b), (b, c), (c, a)):
            edges[(min(u, v), max(u, v))].append(fi)
    dif = []
    for _, fs in edges.items():
        if len(fs) == 2:
            dif.append(np.degrees(np.arccos(np.clip(Nf[fs[0]].dot(Nf[fs[1]]), -1, 1))))
    dif = np.array(dif)

    Vn = np.zeros_like(V)
    for fi, (a, b, c) in enumerate(F):
        for i in (a, b, c):
            Vn[i] += Nf[fi] * area[fi]
    Vn /= np.maximum(np.linalg.norm(Vn, axis=1), 1e-12)[:, None]
    adj = defaultdict(list)
    for (u, v) in edges:
        adj[u].append(v)
        adj[v].append(u)
    vd = []
    for i in range(len(V)):
        for j in adj[i]:
            vd.append(np.degrees(np.arccos(np.clip(Vn[i].dot(Vn[j]), -1, 1))))
    vd = np.array(vd)

    print("\ndihedral angle between adjacent faces (curvature consistency):")
    report("face dihedral", dif)
    print("\nangle between adjacent vertex normals:")
    report("vertex normal", vd)
    print("\nrule of thumb: p95 should be near grid_spacing / curvature_radius;")
    print("for the curved panel (R~400mm, cell 18.2mm) expect ~2.6 deg.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1])
