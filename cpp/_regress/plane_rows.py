#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""每个重建面单独铺平显示：面片轮廓 + 该面的轨迹行（最能看出规整度）。"""
import sys, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter
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
    prefix, outf = sys.argv[1], sys.argv[2]
    meshes = sorted(glob.glob(prefix + ".plane*.ply"),
                    key=lambda p: int(p.rsplit("plane", 1)[1].split(".")[0]))
    if not meshes:
        meshes = sorted(glob.glob(os.path.splitext(prefix)[0] + ".csv.plane*.ply"),
                        key=lambda p: int(p.rsplit("plane", 1)[1].split(".")[0]))
    csvp = prefix if prefix.lower().endswith(".csv") else (os.path.splitext(prefix)[0] + ".csv")
    T = read_csv(csvp)
    # 按面间哨兵拆分：哨兵 = spray==0 且坐标与前一点完全相同
    cuts = []
    for i in range(1, len(T)):
        if T[i, 7] == 0 and np.allclose(T[i, :3], T[i - 1, :3]):
            cuts.append(i)
    blocks, s = [], 0
    for c in cuts:
        blocks.append(T[s:c]); s = c
    blocks.append(T[s:])

    n = len(meshes)
    ncol = min(2, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.4 * ncol, 5.6 * nrow), squeeze=False)

    for k, mp in enumerate(meshes):
        V, F = rd_ply(mp)
        c = V.mean(0)
        # 面内基：最小面积矩形的两条边（不能用 V[1]-V[0]，那只是轮廓上的一段）
        nn = np.zeros(3)
        for tri in F:
            nn += np.cross(V[tri[1]] - V[tri[0]], V[tri[2]] - V[tri[0]])
        nn /= (np.linalg.norm(nn) + 1e-12)
        ref = np.array([1.0, 0.0, 0.0])
        if abs(ref @ nn) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])
        b0 = ref - nn * (ref @ nn); b0 /= np.linalg.norm(b0)
        b1 = np.cross(nn, b0)
        W = np.stack([(V - c) @ b0, (V - c) @ b1], 1)
        best = (1e18, 0.0)
        for deg in range(0, 180):
            th = np.radians(deg)
            R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
            Qr = W @ R.T
            ext = Qr.max(0) - Qr.min(0)
            if ext[0] * ext[1] < best[0]:
                best = (ext[0] * ext[1], th)
        th = best[1]
        u = np.cos(th) * b0 + np.sin(th) * b1
        v = -np.sin(th) * b0 + np.cos(th) * b1
        Q = np.stack([(V - c) @ u, (V - c) @ v], 1)
        ax = axes[k // ncol][k % ncol]

        ec = Counter()
        for tri in F:
            for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                ec[(min(a, b), max(a, b))] += 1
        for (a, b), cnt in ec.items():
            if cnt == 1:
                ax.plot(Q[[a, b], 0], Q[[a, b], 1], color="#1c4f8c", lw=2.2, zorder=3)
        ax.fill(Q[:, 0], Q[:, 1], color="#3f7fd0", alpha=0.10, zorder=1)

        if k < len(blocks):
            B = blocks[k][blocks[k][:, 7] > 0.5]
            if len(B) > 1:
                nrm = np.cross(u, v); nrm /= np.linalg.norm(nrm)
                if nrm @ (B[:, :3].mean(0) - c) < 0:
                    nrm = -nrm
                Bp = B[:, :3] - nrm * 100.0      # 退掉 standoff
                BQ = np.stack([(Bp - c) @ u, (Bp - c) @ v], 1)
                ax.plot(BQ[:, 0], BQ[:, 1], "-", color="#d62728", lw=1.1, zorder=4)
                d = np.diff(BQ, axis=0)
                seg = np.hypot(d[:, 0], d[:, 1])
                ang = np.degrees(np.arctan2(d[:, 1], d[:, 0]))
                ang = ((ang + 90) % 180) - 90
                main = np.median(ang[seg > 2]) if (seg > 2).any() else 0
                # 行数：按垂直于行的投影分簇
                th = np.radians(main)
                Vp = -np.sin(th) * BQ[:, 0] + np.cos(th) * BQ[:, 1]
                o = np.argsort(Vp); ys = []; cur = [o[0]]
                for i in o[1:]:
                    if Vp[i] - Vp[cur[-1]] <= 2.0:
                        cur.append(i)
                    else:
                        ys.append(np.median(Vp[cur])); cur = [i]
                ys.append(np.median(Vp[cur]))
                dy = np.diff(ys)
                dy = dy[dy > 3]        # 去掉转弯点造成的假行
                info = "行数≈%d  行距 %.2f~%.2f mm  行方向偏差 %.1f°" % (
                    len(dy) + 1, dy.min() if len(dy) else 0, dy.max() if len(dy) else 0,
                    abs(main))
                ax.set_title("%s\n顶点=%d 三角=%d\n%s" % (os.path.basename(mp), len(V), len(F), info),
                             fontsize=9)
        ax.set_aspect("equal")
        ax.grid(True, ls=":", lw=0.4)
        ax.set_xlabel("u (mm)"); ax.set_ylabel("v (mm)")

    for k in range(len(meshes), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")

    fig.suptitle(os.path.basename(prefix) + "   —— 每个重建面铺平（蓝框=重建轮廓，红线=轨迹行）",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(outf, dpi=115)
    print("saved", outf)


if __name__ == "__main__":
    main()
