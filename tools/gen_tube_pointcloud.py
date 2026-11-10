#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""圆柱形管子点云生成器（单位：mm，输出 ASCII PLY / CSV）。

典型用法：
    # 外表面圆柱（直径 50、长 500，轴线沿 X，点间距 1mm）
    python tools/gen_tube_pointcloud.py

    # 壁厚 3mm 的开口管：内表面 + 两端环形端面
    python tools/gen_tube_pointcloud.py --wall 3 --caps

    # 加点云噪声、轴线改成 Z、整体平移到机器人坐标系附近
    python tools/gen_tube_pointcloud.py --noise 0.05 --axis Z --translate 534 -451 1457
"""

from __future__ import annotations

import argparse
import math
import os
import numpy as np

# 局部坐标里 z 为管轴方向；按目标轴做循环置换（右旋正交变换）
AXIS_ORDER = {"X": (2, 0, 1), "Y": (1, 2, 0), "Z": (0, 1, 2)}
# 预览视角：让管轴在画面里大致水平，便于看全长
VIEW_AZIM = {"X": -70, "Y": 20, "Z": -60}
# 垂直于管轴的两个坐标分量，用于横截面视图
CROSS_PLANE = {"X": (1, 2), "Y": (0, 2), "Z": (0, 1)}
AXIS_NAME = ("X", "Y", "Z")
# 工程根目录（tools 的上一级），默认输出仍落在 outputs/ 下
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------- 步骤 1：解析命令行参数 ----------
def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="生成圆柱形管子点云（mm）")
    parser.add_argument("--length", type=float, default=500.0, help="管长 mm，默认 500")
    parser.add_argument("--diameter", type=float, default=50.0, help="外径 mm，默认 50")
    parser.add_argument("--wall", type=float, default=0.0,
                        help="壁厚 mm；0 表示只生成外表面（实心管外观）")
    parser.add_argument("--res", type=float, default=1.0, help="点间距 mm，默认 1.0")
    parser.add_argument("--axis", choices=["X", "Y", "Z"], default="X", help="管轴方向，默认 X")
    parser.add_argument("--noise", type=float, default=0.0, help="各向同性高斯噪声标准差 mm")
    parser.add_argument("--caps", action="store_true", help="生成两端端面（实心盘或环形环）")
    parser.add_argument("--translate", type=float, nargs=3, default=(0.0, 0.0, 0.0),
                        metavar=("TX", "TY", "TZ"), help="整体平移 mm")
    parser.add_argument("--seed", type=int, default=20261109, help="随机种子")
    parser.add_argument("--out", default=None, help="输出 .ply 路径（默认 outputs/tube_D50_L500.ply）")
    parser.add_argument("--csv", action="store_true", help="同时输出同名 .csv")
    parser.add_argument("--preview", action="store_true", default=True, help="输出预览 PNG")
    parser.add_argument("--no-preview", dest="preview", action="store_false")
    return parser.parse_args(argv)


# ---------- 步骤 2：生成圆柱侧面（外表面 / 内表面共用） ----------
def build_lateral(radius: float, length: float, res: float, rng, noise: float) -> np.ndarray:
    """沿圆周与轴向均匀采样，得到近似等密度的圆柱面点云。"""
    # 圆周点数取偶数，保证 theta=0 / pi 处落在 ±R，包围盒左右对称
    n_theta = max(12, int(round(2.0 * math.pi * radius / res)))
    n_theta += n_theta % 2
    n_axial = max(2, int(round(length / res)) + 1)
    theta = np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=False)
    axial = np.linspace(-length / 2.0, length / 2.0, n_axial)
    theta_grid, axial_grid = np.meshgrid(theta, axial, indexing="ij")
    x = radius * np.cos(theta_grid)
    y = radius * np.sin(theta_grid)
    pts = np.stack([x, y, axial_grid], axis=-1).reshape(-1, 3)
    return apply_noise(pts, noise, rng)


# ---------- 步骤 3：生成端面（实心圆盘或圆环） ----------
def build_cap(outer_r: float, inner_r: float, length: float, res: float,
              z_sign: float, rng, noise: float) -> np.ndarray:
    """从内径到外径一圈圈采样，每圈点数随半径变化以保持密度一致。"""
    radii = np.arange(inner_r, outer_r + 1e-9, res)
    if radii.size == 0 or radii[-1] < outer_r - 1e-6:
        radii = np.append(radii, outer_r)

    chunks = []
    for radius in radii:
        if radius <= 1e-9:                      # 实心盘圆心单独一个点
            count, phase = 1, 0.0
        else:
            count = max(6, int(round(2.0 * math.pi * radius / res)))
            phase = rng.uniform(0.0, 2.0 * math.pi / count)   # 每圈错开，避免放射状条纹
        theta = np.linspace(0.0, 2.0 * math.pi, count, endpoint=False) + phase
        ring = np.stack([radius * np.cos(theta),
                         radius * np.sin(theta),
                         np.full(count, z_sign * length / 2.0)], axis=-1)
        chunks.append(ring)
    return apply_noise(np.concatenate(chunks, axis=0), noise, rng)


# ---------- 步骤 4：加噪声 ----------
def apply_noise(pts: np.ndarray, noise: float, rng) -> np.ndarray:
    if noise <= 0.0:
        return pts
    return pts + rng.normal(0.0, noise, pts.shape)


# ---------- 步骤 5：姿态与位置变换 ----------
def transform(pts: np.ndarray, axis: str, translate) -> np.ndarray:
    pts = pts[:, AXIS_ORDER[axis]]              # 局部 z 轴 -> 目标轴
    return pts + np.asarray(translate, dtype=float)


# ---------- 步骤 6：写出 PLY / CSV ----------
def write_ply(path: str, pts: np.ndarray, rgb=(200, 205, 212)) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    colors = np.tile(np.asarray(rgb, dtype=np.float64), (pts.shape[0], 1))
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("ply\nformat ascii 1.0\n")
        fh.write("comment synthetic cylinder tube point cloud (unit: mm)\n")
        fh.write(f"element vertex {pts.shape[0]}\n")
        fh.write("property float x\nproperty float y\nproperty float z\n")
        fh.write("property uchar red\nproperty uchar green\nproperty uchar blue\n")
        fh.write("end_header\n")
        np.savetxt(fh, np.hstack([pts, colors]),
                   fmt="%.3f %.3f %.3f %d %d %d")


def write_csv(path: str, pts: np.ndarray, rgb=(200, 205, 212)) -> None:
    colors = np.tile(np.asarray(rgb, dtype=np.float64), (pts.shape[0], 1))
    np.savetxt(path, np.hstack([pts, colors]),
               fmt=["%.3f", "%.3f", "%.3f", "%d", "%d", "%d"],
               delimiter=",", header="x,y,z,r,g,b", comments="")


# ---------- 步骤 7：预览图（3D 全貌 + 横截面） ----------
def save_preview(path: str, pts: np.ndarray, axis: str, title: str,
                 diameter: float, inner_diameter: float = 0.0) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    view = subsample(pts, 40000)
    plane = CROSS_PLANE[axis]                    # 垂直于管轴的两个坐标

    fig = plt.figure(figsize=(12.6, 4.9), dpi=130)
    grid = fig.add_gridspec(1, 2, width_ratios=[2.0, 1.0], wspace=0.18)

    # 左图：3D 全貌
    ax = fig.add_subplot(grid[0, 0], projection="3d")
    ax.scatter(view[:, 0], view[:, 1], view[:, 2], s=0.35, c="#4a7ebb",
               marker=".", depthshade=False, linewidths=0)
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    center, extent = (lo + hi) / 2.0, np.maximum(hi - lo, 1e-6)
    half = extent * 0.54
    ax.set_xlim(center[0] - half[0], center[0] + half[0])
    ax.set_ylim(center[1] - half[1], center[1] + half[1])
    ax.set_zlim(center[2] - half[2], center[2] + half[2])
    # 按真实比例设盒体长宽比，细长管就该看起来是细长的
    ax.set_box_aspect(extent / extent.max())
    # 3D 面板只保留管轴向刻度，Y/Z 交给右侧截面图，避免 3D 轴标签互相压字
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_title(title, fontsize=10)
    ax.view_init(elev=26, azim=VIEW_AZIM[axis])

    # 右图：垂直管轴的截面，叠加理论圆做尺寸对照
    ax2 = fig.add_subplot(grid[0, 1])
    ax2.scatter(view[:, plane[0]] - center[plane[0]], view[:, plane[1]] - center[plane[1]],
                s=2.5, c="#2f6fb5", marker=".", linewidths=0, zorder=3)
    theta = np.linspace(0.0, 2.0 * math.pi, 361)
    for dia, style in ((diameter, "--"), (inner_diameter, ":")):
        if dia <= 0:
            continue
        ax2.plot(dia / 2.0 * np.cos(theta), dia / 2.0 * np.sin(theta),
                 style, color="#d9534f", linewidth=1.0, zorder=2)
    ax2.set_aspect("equal")
    ax2.set_xlabel(f"{AXIS_NAME[plane[0]]} (mm)")
    ax2.set_ylabel(f"{AXIS_NAME[plane[1]]} (mm)")
    ax2.set_title(f"Cross-section (nominal D{diameter:g} mm)", fontsize=10)
    ax2.grid(True, linewidth=0.4, alpha=0.5)

    fig.subplots_adjust(left=0.02, right=0.97, top=0.90, bottom=0.09)
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def subsample(pts: np.ndarray, limit: int) -> np.ndarray:
    if pts.shape[0] <= limit:
        return pts
    idx = np.random.default_rng(0).choice(pts.shape[0], limit, replace=False)
    return pts[np.sort(idx)]


# ---------- 步骤 8：打印统计校验 ----------
def report(pts: np.ndarray, args, out_ply: str) -> None:
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    print(f"点数        : {pts.shape[0]}")
    print(f"包围盒 min  : {lo[0]:.3f}, {lo[1]:.3f}, {lo[2]:.3f}")
    print(f"包围盒 max  : {hi[0]:.3f}, {hi[1]:.3f}, {hi[2]:.3f}")
    print(f"尺寸 XYZ    : {hi[0]-lo[0]:.3f} x {hi[1]-lo[1]:.3f} x {hi[2]-lo[2]:.3f} mm")
    outer_r = args.diameter / 2.0
    print(f"Spec        : 长 {args.length:.1f} mm, 外径 {args.diameter:.1f} mm (R={outer_r:.2f}), "
          f"壁厚 {args.wall:.1f} mm, 点间距 {args.res:.2f} mm")
    print(f"输出        : {out_ply}")


# ---------- 步骤 9：主流程调度 ----------
def main(argv=None) -> int:
    args = parse_args(argv)
    if args.diameter <= 0 or args.length <= 0:
        raise SystemExit("直径和长度必须为正数")
    if args.wall >= args.diameter / 2.0:
        raise SystemExit("壁厚必须小于半径")

    rng = np.random.default_rng(args.seed)
    outer_r = args.diameter / 2.0
    inner_r = outer_r - args.wall if args.wall > 0.0 else 0.0

    parts = [build_lateral(outer_r, args.length, args.res, rng, args.noise)]
    if inner_r > 0.0:                                   # 有壁厚 -> 补内表面
        parts.append(build_lateral(inner_r, args.length, args.res, rng, args.noise))
    if args.caps:                                       # 端面：环形或实心盘
        for z_sign in (-1.0, 1.0):
            parts.append(build_cap(outer_r, inner_r, args.length, args.res,
                                   z_sign, rng, args.noise))

    pts = transform(np.concatenate(parts, axis=0), args.axis, args.translate)

    tag = f"D{args.diameter:g}_L{args.length:g}"
    if args.wall > 0:
        tag += f"_wall{args.wall:g}"
    out_ply = args.out or os.path.join(PROJECT_ROOT, "outputs", f"tube_{tag}.ply")
    write_ply(out_ply, pts)
    if args.csv:
        write_csv(os.path.splitext(out_ply)[0] + ".csv", pts)

    report(pts, args, out_ply)

    if args.preview:
        png = os.path.splitext(out_ply)[0] + "_preview.png"
        title = (f"Cylinder tube  D{args.diameter:g} x L{args.length:g} mm"
                 + (f"  wall {args.wall:g} mm" if args.wall > 0 else "")
                 + f"  (axis {args.axis})")
        save_preview(png, pts, args.axis, title,
                     diameter=args.diameter, inner_diameter=2.0 * inner_r)
        print(f"预览        : {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
