# -*- coding: utf-8 -*-
"""重建面片对比图：修复前/后 × 行距20/60。

上排是修复前：行距20 正常 3 面片，行距60 被拆成 13 个碎矩形；
下排是修复后：两种行距都稳定成同样的 3 个面片（逐字节一致）。
"""
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
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
    F = []
    for k in range(nf):
        s = L[h + nv + k].split()
        F.append([int(x) for x in s[1:1 + int(s[0])]])
    return V, F


def load_cloud(path, maxn=12000):
    V, _ = rd_ply(path)
    step = max(1, len(V) // maxn)
    return V[::step]


def draw(ax, cloud, prefix, title, color_seed=0):
    ax.scatter(cloud[:, 0], cloud[:, 1], cloud[:, 2], s=0.5, c="#d8dee8",
               alpha=0.35, depthshade=False)
    files = sorted(glob.glob(prefix + ".ply.plane*.ply"))
    palette = (plt.cm.tab20(np.linspace(0, 1, 20)) if hasattr(plt.cm, "tab20")
               else plt.cm.get_cmap("tab20")(np.linspace(0, 1, 20)))
    for k, f in enumerate(files):
        V, F = rd_ply(f)
        if not F:
            continue
        tris = [[V[i] for i in t] for t in F]
        pc = Poly3DCollection(tris, alpha=0.85,
                              facecolor=palette[(k + color_seed) % 20],
                              edgecolor="#243447", linewidths=0.5)
        ax.add_collection3d(pc)
    allp = cloud
    c = allp.mean(0)
    r = 0.62 * max(allp.max(0) - allp.min(0))
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    ax.view_init(elev=30, azim=-58)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_zlabel("Z (mm)")
    ax.tick_params(labelsize=7)


def main():
    cloud_path = r"cpp/_regress/ucurved.ply"
    cloud = load_cloud(cloud_path)
    panels = [
        (r"outputs/rowtest20/out", "修复前 · 行距20\n3 个面片（正常）", 0),
        (r"outputs/rowtest60/out", "修复前 · 行距60\n13 个碎面片（重建被行距带偏）", 3),
        (r"outputs/fix12d_row20/out", "修复后 · 行距20\n3 个面片，网格哈希 7E59A025…", 0),
        (r"outputs/fix12d_row60/out", "修复后 · 行距60\n3 个面片，与行距20逐字节一致", 0),
    ]
    fig = plt.figure(figsize=(13.5, 11.5))
    for i, (prefix, title, seed) in enumerate(panels, 1):
        ax = fig.add_subplot(2, 2, i, projection="3d")
        draw(ax, cloud, prefix, title, color_seed=seed)
    fig.suptitle("曲型板 ucurved.ply：扫描行距只影响布行，不应影响点云重建", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = r"outputs/recon_compare_row20_row60.png"
    fig.savefig(out, dpi=150)
    print("saved", out)


if __name__ == "__main__":
    main()
