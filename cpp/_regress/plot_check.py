#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""渲染重建面片 + 轨迹，用于判定"平面是否规则、轨迹是否规整"。"""
import sys, os, math, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "SimSun"]
plt.rcParams["axes.unicode_minus"] = False


def read_ply_ascii(path):
    with open(path, "r", errors="replace") as f:
        lines = f.read().split("\n")
    nv = nf = 0
    for i, ln in enumerate(lines):
        if ln.startswith("element vertex"):
            nv = int(ln.split()[2])
        elif ln.startswith("element face"):
            nf = int(ln.split()[2])
        elif ln.strip() == "end_header":
            hdr = i + 1
            break
    V = np.array([[float(x) for x in lines[hdr + k].split()[:3]] for k in range(nv)])
    F = []
    for k in range(nf):
        p = lines[hdr + nv + k].split()
        F.append([int(x) for x in p[1:1 + int(p[0])]])
    return V, np.array(F, dtype=int)


def read_csv(path):
    rows = []
    with open(path, "r", errors="replace") as f:
        f.readline()
        for ln in f:
            p = ln.strip().split(",")
            if len(p) >= 8:
                rows.append([float(v) for v in p[:6]] + [float(p[6]), int(p[7])])
    return np.array(rows)


def plane_basis(V):
    c = V.mean(0)
    _, _, vt = np.linalg.svd(V - c, full_matrices=False)
    u, v = vt[0], vt[1]
    return c, u, v


def to_uv(P, c, u, v):
    D = P - c
    return np.stack([D @ u, D @ v], axis=1)


def main():
    prefix = sys.argv[1]           # 例如 _regress/tray_f30.ply
    outpng = sys.argv[2]
    title = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(prefix)

    meshes = sorted(glob.glob(prefix + ".plane*.ply"),
                    key=lambda p: int(p.rsplit("plane", 1)[1].split(".")[0]))
    csvp = prefix[:-4] + ".csv"
    traj = read_csv(csvp) if os.path.exists(csvp) else np.zeros((0, 8))

    n = max(1, len(meshes) + (1 if len(traj) else 0))
    ncol = min(2, n)
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2 * ncol, 6.0 * nrow), squeeze=False)

    for idx, mp in enumerate(meshes):
        V, F = read_ply_ascii(mp)
        c, u, v = plane_basis(V)
        UV = to_uv(V, c, u, v)
        ax = axes[idx // ncol][idx % ncol]

        # 三角面片填充
        ax.triplot(UV[:, 0], UV[:, 1], F, color="#9ec5ee", lw=0.8, zorder=1)
        # 外轮廓加粗（只画被 1 个三角形使用的边）
        from collections import Counter
        ec = Counter()
        for tri in F:
            for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                ec[(min(a, b), max(a, b))] += 1
        for (a, b), cnt in ec.items():
            if cnt == 1:
                ax.plot(UV[[a, b], 0], UV[[a, b], 1], color="#1c4f8c", lw=2.0, zorder=3)

        # 轨迹（只显示属于本面的点：用 UV 范围粗筛 + 同法向）
        if len(traj):
            npl = np.array([0.0, 0.0, 0.0])
        info = "%s\n顶点=%d 三角=%d" % (os.path.basename(mp), len(V), len(F))
        ax.set_title(info, fontsize=9)
        ax.set_aspect("equal")
        ax.grid(True, ls=":", lw=0.4)
        ax.set_xlabel("u (mm)")
        ax.set_ylabel("v (mm)")

    # 单独一张：所有面 + 轨迹（俯视，使用全局最佳拟合方向）
    if len(traj):
        ax2 = axes.flat[len(meshes)] if len(meshes) < nrow * ncol else None
        if ax2 is not None:
            T = traj[:, :3]
            c, u, v = plane_basis(T)
            TUV = to_uv(T, c, u, v)
            on = traj[:, 7] > 0.5
            ax2.plot(TUV[~on, 0], TUV[~on, 1], ".", ms=1.0, color="#bdbdbd", label="spray=0")
            ax2.plot(TUV[on, 0], TUV[on, 1], "-", lw=0.7, color="#d62728", label="spray=1")
            ax2.set_aspect("equal")
            ax2.grid(True, ls=":", lw=0.4)
            ax2.legend(fontsize=7, markerscale=4)
            ax2.set_title("轨迹俯视（点=%d，喷=%d）" % (len(traj), int(on.sum())), fontsize=9)
            for k in range(len(meshes)):
                V, F = read_ply_ascii(meshes[k])
                UVm = to_uv(V, c, u, v)
                from collections import Counter
                ec = Counter()
                for tri in F:
                    for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                        ec[(min(a, b), max(a, b))] += 1
                for (a, b), cnt in ec.items():
                    if cnt == 1:
                        ax2.plot(UVm[[a, b], 0], UVm[[a, b], 1], color="#1c4f8c", lw=1.6)

    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpng, dpi=110)
    print("saved", outpng)


if __name__ == "__main__":
    main()
