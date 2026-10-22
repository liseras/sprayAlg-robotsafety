# -*- coding: utf-8 -*-
"""修复前/后 枪位安全距离对比。

左图：每个轨迹点到工件点云的最近距离；右图：每段折线中点到工件的最近距离。
修复前（outputs/_trace_curved.csv）行间直连穿过曲型板，折线最小净空 0.09mm；
修复后（cpp/_regress/regress_out/curved.csv）同一用例最小净空 19.65mm。
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "SimSun"]
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "cpp", "_regress"))
from traj_clear_scan import read_ply, read_traj, Grid, pt_seg_dist  # noqa: E402

CLEAR = 20.0
CAP = 80.0


def clearance(ply, csv):
    Q = read_ply(ply)
    G, sp = read_traj(csv)
    grid = Grid(Q)
    dp = np.full(len(G), CAP)
    for k in range(len(G)):
        c = grid.cand(G[k] - CAP, G[k] + CAP)
        if c is not None:
            dp[k] = min(CAP, float(np.linalg.norm(Q[c] - G[k], axis=1).min()))
    ds = np.full(max(0, len(G) - 1), CAP)
    for k in range(len(G) - 1):
        A, B = G[k], G[k + 1]
        c = grid.cand(np.minimum(A, B) - CAP, np.maximum(A, B) + CAP)
        if c is not None:
            ds[k] = min(CAP, float(pt_seg_dist(Q[c], A, B).min()))
    return dp, ds, sp


def main():
    ply = os.path.join(HERE, "..", "cpp", "_regress", "ucurved.ply")
    old_csv = os.path.join(HERE, "_trace_curved.csv")
    new_csv = os.path.join(HERE, "..", "cpp", "_regress", "regress_out", "curved.csv")
    old_p, old_s, _ = clearance(ply, old_csv)
    new_p, new_s, _ = clearance(ply, new_csv)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 7.6), sharex=False)
    for ax, old, new, ylab, ttl in (
        (ax1, old_p, new_p, "点到工件最近距离 (mm)", "枪位点级净空"),
        (ax2, old_s, new_s, "折线到工件最近距离 (mm)", "相邻轨迹点折线中段净空（穿模指标）"),
    ):
        ax.plot(np.arange(len(old)), old, color="#d64545", lw=1.0,
                label="修复前  min=%.2fmm" % old.min())
        ax.plot(np.arange(len(new)), new, color="#2f6fd0", lw=1.0,
                label="修复后  min=%.2fmm" % new.min())
        ax.axhline(CLEAR, color="#444", ls="--", lw=0.9, label="安全半径 20mm")
        ax.set_ylim(0, 60)
        ax.set_ylabel(ylab)
        ax.set_title(ttl, fontsize=11)
        ax.grid(alpha=0.25)
        ax.legend(loc="lower right", fontsize=9)

    # 标出修复前最危险的区段（195mm 法向翻转 + 648mm 跨行直连）
    ax2.annotate("行内法向翻转 195mm\n折线穿板 0.09mm",
                 xy=(1, old_s[1]), xytext=(60, 18),
                 arrowprops=dict(arrowstyle="->", color="#d64545"),
                 color="#d64545", fontsize=9)
    ax1.annotate("648mm 跨行直连\n枪位贴板 1.71mm",
                 xy=(139, old_p[139]), xytext=(300, 10),
                 arrowprops=dict(arrowstyle="->", color="#d64545"),
                 color="#d64545", fontsize=9)

    fig.suptitle("曲型板默认行距12：点云重建修好后暴露的行间穿模，已改为安全空移", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(HERE, "safety_compare.png")
    fig.savefig(out, dpi=150)
    print("saved", out)


if __name__ == "__main__":
    main()
