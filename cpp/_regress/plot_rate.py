# -*- coding: utf-8 -*-
"""图：各样本的枪姿角速度峰值 vs 工业机器人手腕常见限速"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

CASES = ["box", "tray", "zhe", "cao", "complex", "curved", "bin"]
OUT = "syn/plot_rate.png"
DIR = "regress_out"


def main():
    p95s, maxs = [], []
    for name in CASES:
        with open("%s/%s.csv" % (DIR, name)) as f:
            rd = csv.reader(f)
            h = next(rd)
            idx = {k: i for i, k in enumerate(h)}
            rows = [r for r in rd]
        P = np.array([[float(r[idx[c]]) for c in "XYZ"] for r in rows])
        A = np.array([[float(r[idx[c]]) for c in "IJK"] for r in rows])
        sp = np.array([int(float(r[idx["spray"]])) for r in rows])
        spd = np.array([float(r[idx["speed"]]) for r in rows])

        step = np.linalg.norm(np.diff(P, axis=0), axis=1)
        ang = np.degrees(np.arccos(np.clip((A[1:] * A[:-1]).sum(1), -1, 1)))
        m = (sp[:-1] == 1) & (sp[1:] == 1) & (step < 8.0) & (step > 1e-6)
        rad = ang[m] / step[m] * np.minimum(spd[:-1][m], spd[1:][m])
        p95s.append(np.percentile(rad, 95))
        maxs.append(rad.max())

    x = np.arange(len(CASES))
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    ax.bar(x - 0.19, p95s, width=0.36, color="#4c78a8", label="95% 分位角速度")
    ax.bar(x + 0.19, maxs, width=0.36, color="#d81e06", label="峰值角速度")
    ax.axhline(180, ls="--", lw=1.2, c="#666")
    ax.text(-0.45, 195, "手腕常见限速 180°/s", ha="left", fontsize=8, color="#666")
    ax.text(-0.45, 375, "上限 360°/s", ha="left", fontsize=8, color="#999")
    ax.axhline(360, ls=":", lw=1.2, c="#999")
    for i, v in enumerate(maxs):
        if v > 500:
            ax.text(i + 0.19, v * 1.02, "%.0f" % v, ha="center", fontsize=8, color="#d81e06")
    ax.set_xticks(x)
    ax.set_xticklabels(CASES)
    ax.set_ylabel("枪姿角速度 °/s")
    ax.set_title("同行相邻点之间的枪轴转动速率（按速度列折算）", fontsize=11)
    ax.legend(fontsize=8)
    ax.set_ylim(0, max(maxs) * 1.12)
    plt.tight_layout()
    plt.savefig(OUT, dpi=120, facecolor="white")
    print("saved", OUT)
    for n, a, b in zip(CASES, p95s, maxs):
        print("%-9s p95=%7.0f  max=%7.0f °/s" % (n, a, b))


if __name__ == "__main__":
    main()
