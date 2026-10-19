#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""面片 + 该面上全部轨迹点（按点到面距离筛选，不依赖 CSV 分段）—— 最可信的规整度检查图。"""
import sys, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False


def rd_ply(p):
    L = open(p, errors="replace").read().split("\n")
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
    with open(p, errors="replace") as f:
        f.readline()
        for ln in f:
            s = ln.strip().split(",")
            if len(s) >= 8:
                out.append([float(v) for v in s[:6]] + [float(s[6]), int(s[7])])
    return np.array(out)


def basis_minrect(V):
    c = V.mean(0)
    _, _, vt = np.linalg.svd(V - c, full_matrices=False)
    n = vt[2]
    b0 = vt[0]
    b1 = vt[1]
    W = np.stack([(V - c) @ b0, (V - c) @ b1], 1)
    best = (1e18, 0.0)
    for deg in range(0, 180):
        th = np.radians(deg)
        R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
        Qr = W @ R.T
        e = Qr.max(0) - Qr.min(0)
        if e[0] * e[1] < best[0]:
            best = (e[0] * e[1], th)
    th = best[1]
    u = np.cos(th) * b0 + np.sin(th) * b1
    v = -np.sin(th) * b0 + np.cos(th) * b1
    return c, u, v, n


def main():
    cloud, prefix, outf = sys.argv[1], sys.argv[2], sys.argv[3]
    meshes = sorted(glob.glob(prefix + ".plane*.ply"),
                    key=lambda p: int(p.rsplit("plane", 1)[1].split(".")[0]))
    csvp = prefix if prefix.lower().endswith(".csv") else (os.path.splitext(prefix)[0] + ".csv")
    T = read_csv(csvp)
    P = T[:, :3]
    on = T[:, 7] > 0

    n = len(meshes)
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))
    fig = plt.figure(figsize=(8.0 * ncol, 5.4 * nrow))
    for k, mp in enumerate(meshes):
        V, F = rd_ply(mp)
        c, u, v, nr = basis_minrect(V)
        Q = np.stack([(V - c) @ u, (V - c) @ v], 1)
        ax = fig.add_subplot(nrow, ncol, k + 1)
        for tri in F:
            ax.fill(Q[[t for t in tri], 0], Q[[t for t in tri], 1],
                    color="#3f7fd0", alpha=0.13, zorder=1, ec="#1c4f8c", lw=0.5)
        # 该面上的轨迹点：轨迹是喷头位置（沿法向抬了 standoff），所以只按面内
        # 坐标是否落在面片范围内筛选，不用点面距离。
        Qp = np.stack([(P[on] - c) @ u, (P[on] - c) @ v], 1)
        lo = Q.min(0) - 4.0
        hi = Q.max(0) + 4.0
        m = np.all((Qp >= lo) & (Qp <= hi), axis=1)
        Qp = Qp[m]
        ax.plot(Qp[:, 0], Qp[:, 1], ".", ms=1.0, color="#d62728", zorder=3)
        # 行数/行距（按短边方向统计）
        ext = Q.max(0) - Q.min(0)
        if len(Qp) > 20:
            # 行方向 = 用户配置的扫描方向；用点云直方图找真实的“行进轴”
            best = None
            for deg in range(0, 180, 1):
                th = np.radians(deg)
                a = np.array([np.cos(th), np.sin(th)])
                t = Qp @ a
                b = -np.sin(th) * Qp[:, 0] + np.cos(th) * Qp[:, 1]
                # 行的“横向”投影应该有明显的周期性 -> 用直方图峰值数量衡量
                hh, edges = np.histogram(b, bins=max(8, int((b.max() - b.min()) / 2)))
                peaks = np.sum((hh[1:-1] > hh[:-2]) & (hh[1:-1] >= hh[2:]) & (hh[1:-1] > hh.max() * 0.3))
                if best is None or peaks < best[1]:
                    best = (b, peaks)
            b = best[0]
            o = np.sort(b)
            grp = []
            cur = [o[0]]
            for x in o[1:]:
                if x - cur[-1] <= 2.5:
                    cur.append(x)
                else:
                    grp.append(np.mean(cur))
                    cur = [x]
            grp.append(np.mean(cur))
            dy = np.diff(grp)
            txt = "行数 %d  行距中位 %.2f mm（%.2f~%.2f）" % (
                len(grp), np.median(dy) if len(dy) else 0,
                dy.min() if len(dy) else 0, dy.max() if len(dy) else 0)
        else:
            txt = "轨迹点 %d" % len(Qp)
        ax.set_title("%s  顶点=%d 三角=%d  面尺寸 %.0f×%.0f\n%s"
                     % (os.path.basename(mp), len(V), len(F), ext[0], ext[1], txt), fontsize=9)
        ax.set_aspect("equal")
        ax.grid(True, ls=":", lw=0.4)
    fig.suptitle(os.path.basename(prefix) + " —— 重建面（蓝）与落在其上的喷涂轨迹（红）", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outf, dpi=115)
    print("saved", outf)


if __name__ == "__main__":
    main()
