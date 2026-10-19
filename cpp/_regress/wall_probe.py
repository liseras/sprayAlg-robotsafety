#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""剖面 + 陡壁点分布，判断侧壁是环形柱面还是长条平面。"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "SimSun"]
plt.rcParams["axes.unicode_minus"] = False
from numpy.linalg import svd


def read_ply(path):
    with open(path, "r", errors="replace") as f:
        lines = f.read().split("\n")
    nv = 0
    for i, ln in enumerate(lines):
        if ln.startswith("element vertex"):
            nv = int(ln.split()[2])
        elif ln.strip() == "end_header":
            hdr = i + 1
            break
    return np.array([[float(x) for x in lines[hdr + k].split()[:3]] for k in range(nv)])


def main():
    P = read_ply(sys.argv[1])
    c = P.mean(0)
    _, _, vt = svd(P - c, full_matrices=False)
    nrm, u, v = vt[2], vt[0], vt[1]
    if nrm[2] < 0:
        nrm = -nrm
    H = (P - c) @ nrm
    U, V = (P - c) @ u, (P - c) @ v

    from scipy.spatial import cKDTree
    _, idx = cKDTree(P).query(P, k=25)
    N = np.array([svd(P[idx[i]] - P[idx[i]].mean(0), full_matrices=False)[2][2] for i in range(len(P))])
    ang = np.degrees(np.arccos(np.clip(np.abs(N @ nrm), 0, 1)))
    steep = ang > 40

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0][0]
    ax.scatter(U[~steep], V[~steep], s=0.5, c="#b0c4de", label="平缓")
    ax.scatter(U[steep], V[steep], s=1.2, c="#d62728", label="陡壁(>40°)")
    ax.set_aspect("equal"); ax.legend(fontsize=8, markerscale=6)
    ax.set_title("俯视：陡壁点分布（环形=柱面，直条=平面侧壁）")
    ax.grid(True, ls=":", lw=0.4)

    # 剖面：沿 u 的若干窄带
    ax = axes[0][1]
    for v0, col in zip(np.percentile(V, [20, 40, 60, 80]), ["#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd"]):
        b = np.abs(V - v0) < 1.5
        o = np.argsort(U[b])
        ax.plot(U[b][o], H[b][o], ".", ms=1.2, color=col, label="v≈%.0f" % v0)
    ax.set_xlabel("u (mm)"); ax.set_ylabel("高度 (mm)")
    ax.set_title("沿 u 的剖面束（竖向段=侧壁）")
    ax.legend(fontsize=7, markerscale=6); ax.grid(True, ls=":", lw=0.4)

    # 剖面：沿 v 的若干窄带
    ax = axes[1][0]
    for u0, col in zip(np.percentile(U, [20, 40, 60, 80]), ["#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd"]):
        b = np.abs(U - u0) < 1.5
        o = np.argsort(V[b])
        ax.plot(V[b][o], H[b][o], ".", ms=1.2, color=col, label="u≈%.0f" % u0)
    ax.set_xlabel("v (mm)"); ax.set_ylabel("高度 (mm)")
    ax.set_title("沿 v 的剖面束（竖向段=侧壁）")
    ax.legend(fontsize=7, markerscale=6); ax.grid(True, ls=":", lw=0.4)

    # 陡壁点做平面拟合：RANSAC 看能拟出几个面
    ax = axes[1][1]
    Q = P[steep]
    Qn = N[steep]
    # 按法向聚类（简单 k-means）
    from scipy.cluster.vq import kmeans2
    K = 4
    cent, lab = kmeans2(Qn, K, minit="++", seed=3)
    print("陡壁点 %d 个，按法向聚成 %d 组：" % (len(Q), K))
    for k in range(K):
        m = lab == k
        if m.sum() < 5:
            continue
        # 组内平面拟合残差
        cc = Q[m].mean(0)
        _, s, vt2 = svd(Q[m] - cc, full_matrices=False)
        rms = s[2] / np.sqrt(m.sum())
        a = np.degrees(np.arccos(np.clip(abs(cent[k] @ nrm), 0, 1)))
        print("  组%d: %4d 点 法向与地面夹角=%5.1f°  平面拟合残差RMS=%6.2f mm" % (k, m.sum(), a, rms))
        ax.scatter((Q[m] - cc) @ vt2[0], (Q[m] - cc) @ vt2[1], s=1.0,
                   label="组%d %d点 %.0f°" % (k, m.sum(), a))
    ax.set_aspect("equal"); ax.legend(fontsize=7, markerscale=6)
    ax.set_title("陡壁点按法向分组（每组铺平后）")
    ax.grid(True, ls=":", lw=0.4)

    fig.suptitle(sys.argv[1], fontsize=11)
    fig.tight_layout()
    out = sys.argv[1].rsplit(".", 1)[0] + "_wall.png"
    fig.savefig(out, dpi=110)
    print("saved", out)


if __name__ == "__main__":
    main()
