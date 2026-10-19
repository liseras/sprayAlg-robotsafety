#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""分析点云几何：bbox、Z 直方图、剖面，判断有没有"接近垂直的侧壁"。"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "SimSun"]
plt.rcParams["axes.unicode_minus"] = False


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
    P = np.array([[float(x) for x in lines[hdr + k].split()[:3]] for k in range(nv)])
    return P


def main():
    path = sys.argv[1]
    P = read_ply(path)
    print("点数=%d" % len(P))
    lo, hi = P.min(0), P.max(0)
    print("bbox min=(%.1f, %.1f, %.1f) max=(%.1f, %.1f, %.1f)" % (*lo, *hi))
    print("尺寸=(%.1f, %.1f, %.1f) mm" % tuple(hi - lo))

    # 拟合主平面
    c = P.mean(0)
    _, _, vt = np.linalg.svd(P - c, full_matrices=False)
    nrm, u, v = vt[2], vt[0], vt[1]
    if nrm[2] < 0:
        nrm = -nrm
    H = (P - c) @ nrm
    print("主法向=(%.3f, %.3f, %.3f) 高度范围=[%.1f, %.1f]" % (*nrm, H.min(), H.max()))

    # 用法向差判断"陡壁"：局部 PCA 法向
    from numpy.linalg import svd
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(P)
        _, idx = tree.query(P, k=25)
    except Exception:
        rng = np.random.default_rng(0)
        idx = np.empty((len(P), 25), dtype=int)
        for i in range(0, len(P), 500):
            d = ((P[None, :, :] - P[i:i + 500, None, :]) ** 2).sum(2)
            idx[i:i + 500] = np.argsort(d, axis=1)[:, :25]
    N = np.zeros_like(P)
    for i in range(len(P)):
        Q = P[idx[i]] - P[idx[i]].mean(0)
        N[i] = svd(Q, full_matrices=False)[2][2]
    cosang = np.abs(N @ nrm)
    ang = np.degrees(np.arccos(np.clip(cosang, 0, 1)))
    print("局部法向与主法向夹角分位数: "
          + ", ".join("P%d=%.1f°" % (q, np.percentile(ang, q)) for q in (50, 75, 90, 95, 99)))
    steep = ang > 45
    print("陡壁点(>45°)占比=%.2f%% (%d 点)" % (100.0 * steep.mean(), steep.sum()))

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))

    # 1) Z 直方图
    ax = axes[0]
    ax.hist(H, bins=120, color="#4a90d9")
    ax.set_xlabel("沿主法向高度 (mm)")
    ax.set_ylabel("点数")
    ax.set_title("高度直方图（多峰=有台阶/侧壁）")

    # 2) 俯视散点（按夹角着色）
    ax = axes[1]
    U, V = (P - c) @ u, (P - c) @ v
    sc = ax.scatter(U, V, c=ang, s=0.6, cmap="turbo", vmin=0, vmax=90)
    plt.colorbar(sc, ax=ax, label="局部法向与主法向夹角 (°)")
    ax.set_aspect("equal")
    ax.set_title("俯视：夹角>45° 即陡壁（红/黄）")

    # 3) 剖面：沿某一条线取点，看横向轮廓
    ax = axes[2]
    # 沿 V 方向每 2mm 一个切片，画 U-Z 剖面
    band = np.abs(U - np.median(U)) < 3.0
    Q = P[band]
    Qu, Qv = (Q - c) @ u, (Q - c) @ v
    ax.scatter(Qu, H[band], s=0.8, c="#d62728")
    ax.set_xlabel("u (mm)")
    ax.set_ylabel("高度 (mm)")
    ax.set_title("中间纵剖面（可见竖直壁则 U 方向有近垂直段）")
    ax.grid(True, ls=":", lw=0.4)

    fig.suptitle(path, fontsize=11)
    fig.tight_layout()
    out = path.rsplit(".", 1)[0] + "_geom.png"
    fig.savefig(out, dpi=110)
    print("saved", out)


if __name__ == "__main__":
    main()
