# -*- coding: utf-8 -*-
"""把审计结果画成图：漏喷位置 + 喷距偏差位置

左图：点云按"到轨迹折线的距离"着色（灰=覆盖到，红=空档 > 0.65×行距）
右图：喷枪点按"沿枪轴实际喷距"着色（蓝=偏近，灰=正常，红=偏远）

用法：python audit_plot.py <cloud.ply> <traj.csv> <out.png> [--row 12] [--so 100]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from audit_traj import load_cloud, load_traj, polyline_dist, spray_dist, Grid


def main():
    ply, csvp, out = sys.argv[1], sys.argv[2], sys.argv[3]
    L, so = 12.0, 100.0
    if "--row" in sys.argv:
        L = float(sys.argv[sys.argv.index("--row") + 1])
    if "--so" in sys.argv:
        so = float(sys.argv[sys.argv.index("--so") + 1])

    Q = load_cloud(ply)
    P, A, sp = load_traj(csvp)
    C = (P + so * A)[sp == 1]

    rs = np.random.default_rng(3)
    Qs = Q[rs.choice(len(Q), min(40000, len(Q)), replace=False)]
    d = polyline_dist(Qs, C, L, max(4 * 0.65 * L, so))
    gap = d > 0.65 * L

    fig = plt.figure(figsize=(13, 5.6))
    ax = fig.add_subplot(121, projection="3d")
    ax.computed_zorder = False
    ax.scatter(Qs[~gap, 0], Qs[~gap, 1], Qs[~gap, 2], s=1, c="#c8ccd4", label="覆盖到", zorder=1)
    ax.scatter(Qs[gap, 0], Qs[gap, 1], Qs[gap, 2], s=4, c="#d81e06", label="空档 > 0.65L", zorder=3)
    ax.plot(C[:, 0], C[:, 1], C[:, 2], lw=0.6, c="#1f4e9c", alpha=0.9, zorder=4)
    ax.set_title("漏喷位置（红 = 到最近轨迹线 > %.1fmm）" % (0.65 * L), fontsize=10)
    ax.legend(loc="upper right", fontsize=7, markerscale=4)
    ax.set_box_aspect([1, 1, 0.7])
    ax.view_init(elev=28, azim=-58)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_ticklabels([])

    # 喷距
    grid = Grid(Q, cell=8.0)
    idx = np.where(sp == 1)[0]
    sel = idx[np.linspace(0, len(idx) - 1, min(1500, len(idx))).astype(int)]
    dv, pv = [], []
    for k in sel:
        p, a = P[k], A[k]
        hi = p + a * 3 * so
        c = grid.cand(np.minimum(p, hi) - 6, np.maximum(p, hi) + 6)
        if c is None:
            continue
        q = Q[c] - p
        t = q @ a
        perp = np.linalg.norm(q - t[:, None] * a, axis=1)
        ok = (t > 1) & (t < 3 * so) & (perp <= 6)
        if ok.any():
            dv.append(t[ok].min())
            pv.append(p)
    dv = np.array(dv)
    pv = np.array(pv)

    ax2 = fig.add_subplot(122, projection="3d")
    ax2.computed_zorder = False
    sc = ax2.scatter(pv[:, 0], pv[:, 1], pv[:, 2], s=9, c=dv, cmap="coolwarm",
                     vmin=so * 0.6, vmax=so * 1.4, zorder=3)
    fig.colorbar(sc, ax=ax2, shrink=0.6, label="实际喷距 mm（设定 %.0f）" % so)
    ax2.set_title("喷距偏差（红=偏远 蓝=偏近）", fontsize=10)
    ax2.set_box_aspect([1, 1, 0.7])
    ax2.view_init(elev=28, azim=-58)
    for a in (ax2.xaxis, ax2.yaxis, ax2.zaxis):
        a.set_ticklabels([])

    plt.tight_layout()
    plt.savefig(out, dpi=110, facecolor="white")
    print("saved", out, "| 空档占比 %.1f%%" % (100.0 * gap.mean()))


if __name__ == "__main__":
    sys.exit(main())
