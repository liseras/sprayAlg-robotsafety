#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""3D 视图：原始点云 + 重建面片 + 轨迹，按 WPF 观感近似。"""
import sys, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "SimSun"]
plt.rcParams["axes.unicode_minus"] = False


def rd_ply(path):
    with open(path, "r", errors="replace") as f:
        L = f.read().split("\n")
    nv = nf = 0
    for i, l in enumerate(L):
        if l.startswith("element vertex"):
            nv = int(l.split()[2])
        elif l.startswith("element face"):
            nf = int(l.split()[2])
        elif l.strip() == "end_header":
            h = i + 1
            break
    V = np.array([[float(x) for x in L[h + k].split()[:3]] for k in range(nv)])
    F = [[int(x) for x in L[h + nv + k].split()[1:]] for k in range(nf)] if nf else []
    return V, F


def read_csv(p):
    out = []
    with open(p, "r", errors="replace") as f:
        f.readline()
        for ln in f:
            s = ln.strip().split(",")
            if len(s) >= 8:
                out.append([float(v) for v in s[:6]] + [float(s[6]), int(s[7])])
    return np.array(out)


def main():
    cloud, prefix, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    elev = float(sys.argv[4]) if len(sys.argv) > 4 else 28
    azim = float(sys.argv[5]) if len(sys.argv) > 5 else -60

    fig = plt.figure(figsize=(13, 10))
    ax = fig.add_subplot(111, projection="3d")

    if cloud and os.path.exists(cloud):
        P, _ = rd_ply(cloud)
        s = max(1, len(P) // 12000)
        ax.scatter(P[::s, 0], P[::s, 1], P[::s, 2], s=0.7, c="#c8d4e2", alpha=0.55,
                   label="原始点云")

    files = sorted(glob.glob(prefix + ".plane*.ply"))
    if not files:
        files = sorted(glob.glob(os.path.splitext(prefix)[0] + ".csv.plane*.ply"))
    for f in files:
        V, F = rd_ply(f)
        if not F:
            continue
        tris = [[V[i] for i in t] for t in F]
        pc = Poly3DCollection(tris, alpha=0.42, facecolor="#3f7fd0",
                             edgecolor="#1c4f8c", linewidths=0.8)
        ax.add_collection3d(pc)

    csvp = prefix if prefix.lower().endswith(".csv") else (os.path.splitext(prefix)[0] + ".csv")
    t = read_csv(csvp)
    on = t[:, 7] > 0.5
    ax.plot(t[~on, 0], t[~on, 1], t[~on, 2], "-", lw=0.8, color="#9e9e9e", label="空走 spray=0")
    ax.plot(t[on, 0], t[on, 1], t[on, 2], "-", lw=1.0, color="#d62728", label="喷涂 spray=1")

    allp = np.vstack([t[:, :3]])
    if cloud and os.path.exists(cloud):
        P, _ = rd_ply(cloud)
        allp = np.vstack([allp, P])
    c = allp.mean(0)
    r = 0.6 * max(allp.max(0) - allp.min(0))
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title(os.path.basename(prefix))
    fig.tight_layout()
    fig.savefig(outp, dpi=115)
    print("saved", outp)


if __name__ == "__main__":
    main()
