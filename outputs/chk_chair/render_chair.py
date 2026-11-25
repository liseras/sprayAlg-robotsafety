#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""椅子点云的轨迹出图（临时脚手架）：前视 X-Y 叠落点 + 侧视 Z-Y 看喷距。

落点 = TCP + 枪轴×standoff，与 chk_two/render.py 同一口径。
用法: python render_chair.py <cloud.ply> <traj.csv> <out.png> [standoff]
"""
import sys, os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPRAY = "#dc2626"
MOVE = "#2563eb"
XFACE = "#f97316"
PTS = "#94a3b8"


def read_ply_xyz(path):
    xs = []
    with open(path, "r", errors="ignore") as f:
        hdr = True
        for ln in f:
            if hdr:
                if ln.strip() == "end_header":
                    hdr = False
                continue
            p = ln.split()
            if len(p) >= 3:
                xs.append((float(p[0]), float(p[1]), float(p[2])))
    return np.array(xs)


def read_csv(path):
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as f:
        r = csv.reader(f)
        head = next(r, None) or []
        rows = [x for x in r if x and len(x) >= 8]
    low = [h.strip().lower() for h in head]
    xyz = np.array([[float(t[0]), float(t[1]), float(t[2])] for t in rows])
    gun = np.array([[float(t[3]), float(t[4]), float(t[5])] for t in rows])
    gun /= np.maximum(np.linalg.norm(gun, axis=1), 1e-12)[:, None]
    spray = np.array([0 if t[7].strip() in ("0", "0.0", "") else 1 for t in rows], int)
    kind = np.zeros(len(rows), int)
    if "kind" in low:
        ki = low.index("kind")
        kind = np.array([int(float(t[ki])) if ki < len(t) and t[ki].strip() else 0 for t in rows], int)
    return xyz, gun, spray, kind


def main():
    cloud_p, csv_p, out_p = sys.argv[1], sys.argv[2], sys.argv[3]
    so = float(sys.argv[4]) if len(sys.argv) > 4 else 150.0
    cloud = read_ply_xyz(cloud_p)
    tcp, gun, spray, kind = read_csv(csv_p)
    hit = tcp + gun * so

    fig = plt.figure(figsize=(13.4, 6.2), dpi=120)
    ax = fig.add_subplot(1, 2, 1)
    ax.scatter(cloud[:, 0], cloud[:, 1], s=2.2, c=PTS, linewidths=0, label="cloud")
    m = spray == 0
    ax.scatter(hit[m, 0], hit[m, 1], s=5, c=MOVE, linewidths=0, label="move")
    m = (spray != 0) & (kind == 0)
    ax.scatter(hit[m, 0], hit[m, 1], s=5, c=SPRAY, linewidths=0, label="spray")
    m = (spray != 0) & (kind != 0)
    ax.scatter(hit[m, 0], hit[m, 1], s=7, c=XFACE, linewidths=0, label="cross-face")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(f"front view X-Y, landing pts (standoff {so:.0f}mm), {len(hit)} pts", fontsize=10)
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.legend(markerscale=3, fontsize=8, loc="upper right")
    ax.grid(True, linewidth=0.3, alpha=0.4)

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.scatter(cloud[:, 2], cloud[:, 1], s=2.2, c=PTS, linewidths=0)
    ax2.scatter(hit[:, 2], hit[:, 1], s=4, c=SPRAY, linewidths=0)
    ax2.set_aspect("equal")
    ax2.set_title("side view Z-Y (nearest = camera side)", fontsize=10)
    ax2.set_xlabel("Z (mm)")
    ax2.set_ylabel("Y (mm)")
    ax2.grid(True, linewidth=0.3, alpha=0.4)

    fig.tight_layout()
    fig.savefig(out_p, facecolor="white")
    print("saved", out_p)


if __name__ == "__main__":
    main()
