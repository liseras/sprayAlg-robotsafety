#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""访客椅点云生成器（网布靠背 + 一体坐垫 + 塑料扶手 + 悬臂钢管腿，单位 mm）。

建模思路和"像不像"的关键：
  1. 先按实物尺寸参数化建模（坐垫挤出体、弯板靠背、扫掠钢管）；
  2. 再用**虚拟深度相机**正投影 + z-buffer 取点 —— 只保留最先被打到的表面。
     所以产出和真实扫描一致：背面没有点、扶手挡住的地方没有点、坐垫前沿后面是空的。
这与 cloud-dump 里 p0xx-cam.ply 的口径一致：相机系、可喷侧朝 -Z、无颜色只有 intensity。

典型用法：
    python tools/gen_chair_pointcloud.py                     # 默认 → outputs/chair_visitor.ply
    python tools/gen_chair_pointcloud.py --pitch 0 --yaw 0    # 正视图（坐垫顶面会很稀，真实口径）
    python tools/gen_chair_pointcloud.py --pitch 15 --yaw 10 --pitch-step 2.2
"""

from __future__ import annotations

import argparse
import math
import os

import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 相机系放置：椅脚中心落在哪。默认对齐现有工位（相机 +Z = 基座 +X，相机 y ≈ 基座 z）
DEFAULT_ORIGIN = (900.0, -711.0, 1100.0)   # 椅脚中心在相机系的 (x, y, z)


# ---------------------------------------------------------------- 平面轮廓工具
def rounded_rect_boundary(a: float, b: float, r: float, step: float):
    """圆角矩形（半宽 a / 半深 b / 圆角 r）的外轮廓点与外法向。返回 (N,2),(N,2)。"""
    r = max(0.0, min(r, min(a, b) * 0.999))
    P, N = [], []
    corners = ((a - r, b - r, 0.0), (-(a - r), b - r, 0.5 * math.pi),
               (-(a - r), -(b - r), math.pi), (a - r, -(b - r), 1.5 * math.pi))
    n_arc = max(4, int(0.5 * math.pi * r / step))
    for i, (cx, cz, th0) in enumerate(corners):
        th = th0 + np.linspace(0.0, 0.5 * math.pi, n_arc + 1)
        P.append(np.stack([cx + r * np.cos(th), cz + r * np.sin(th)], axis=1))
        N.append(np.stack([np.cos(th), np.sin(th)], axis=1))
        # 连到下一条边的起点：直边外法向就是该边的轴向
        edge_n = [(0.0, 1.0), (-1.0, 0.0), (0.0, -1.0), (1.0, 0.0)][i]
        p0 = P[-1][-1]
        nxt = corners[(i + 1) % 4]
        p1 = np.array([nxt[0] + r * math.cos(nxt[2]), nxt[1] + r * math.sin(nxt[2])])
        seg = float(np.linalg.norm(p1 - p0))
        n_ed = max(2, int(seg / step) + 1)
        t = np.linspace(0.0, 1.0, n_ed)[1:, None]
        P.append(p0[None, :] * (1.0 - t) + p1[None, :] * t)
        N.append(np.tile(np.array(edge_n, dtype=float), (n_ed - 1, 1)))
    return np.concatenate(P, axis=0), np.concatenate(N, axis=0)


def rounded_rect_mask(X, Z, a, b, r):
    """点是否落在圆角矩形内（a,b 为半宽半深，r 为圆角）。"""
    r = max(1e-6, min(r, min(a, b) * 0.999))
    ix = np.abs(X) - (a - r)
    iz = np.abs(Z) - (b - r)
    outside = np.hypot(np.maximum(ix, 0.0), np.maximum(iz, 0.0))
    inside_core = np.minimum(np.maximum(ix, iz), 0.0)
    return (outside + inside_core) <= r


# ---------------------------------------------------------------- 基本体：挤出体
def extrude_rounded_rect(cx, cz, w, d, y_bot, y_top, r_corner, r_fillet, step,
                         dish=0.0, slope=0.0, z_ref=0.0):
    """圆角矩形挤出体（六面 + 上下圆角过渡），用于坐垫、扶手垫、坐垫下的盆。

    dish   : 顶面中央下凹量（mm，正数=中间低）
    slope  : 顶面沿 +Z 的坡（正数=后面低），单位 mm/mm
    """
    a, b = 0.5 * w, 0.5 * d
    r_corner = min(r_corner, min(a, b) * 0.98)
    r_fillet = min(r_fillet, 0.45 * (y_top - y_bot), min(a, b) - r_corner - 1.0)
    r_fillet = max(r_fillet, 0.5)

    def top_y(x, z):
        y = y_top + slope * (z - z_ref)
        if dish and a > 1e-6:
            u = np.clip(x / max(a - r_fillet, 1e-6), -1.0, 1.0)
            v = np.clip((z - z_ref) / max(b - r_fillet, 1e-6), -1.0, 1.0)
            y = y - dish * (1.0 - u * u) * (1.0 - v * v)
        return y

    parts = []
    # 顶面 / 底面：内缩 r_fillet 的圆角矩形网格
    na = int(2 * (a - r_fillet) / step) + 1
    nb = int(2 * (b - r_fillet) / step) + 1
    gx = np.linspace(-(a - r_fillet), (a - r_fillet), max(na, 3))
    gz = np.linspace(-(b - r_fillet), (b - r_fillet), max(nb, 3))
    GX, GZ = np.meshgrid(gx, gz, indexing="ij")
    keep = rounded_rect_mask(GX, GZ, a - r_fillet, b - r_fillet, max(r_corner - r_fillet, 0.0))
    for sgn in (+1.0, -1.0):
        x = GX[keep] + cx
        z = GZ[keep] + cz
        y = top_y(GX[keep], GZ[keep]) if sgn > 0 else np.full(x.shape, y_bot)
        parts.append(np.stack([x, y, z], axis=1))

    # 侧壁 + 上下圆角
    B, Nrm = rounded_rect_boundary(a, b, r_corner, step)
    n_arc = max(3, int(0.5 * math.pi * r_fillet / step))
    alpha = np.linspace(0.0, 0.5 * math.pi, n_arc + 1)
    for sgn in (+1.0, -1.0):
        for al in alpha:
            inset = r_fillet * (1.0 - math.cos(al))
            y = (y_top - r_fillet + r_fillet * math.sin(al)) if sgn > 0 else \
                (y_bot + r_fillet - r_fillet * math.sin(al))
            x = B[:, 0] - Nrm[:, 0] * inset + cx
            z = B[:, 1] - Nrm[:, 1] * inset + cz
            yy = np.full(x.shape, y)
            if sgn > 0:                                   # 顶面圆角要跟着 dish/slope 走
                wgt = math.sin(al)
                yy = yy + wgt * (top_y(B[:, 0], B[:, 1]) - y_top)
            parts.append(np.stack([x, yy, z], axis=1))
    return np.concatenate(parts, axis=0)


# ---------------------------------------------------------------- 基本体：弯板
def curved_panel(spec, step):
    """靠背弯板：横向内凹 + 纵向 S 形（腰托 / 上沿后卷），带厚度与侧边。

    spec 里给 y0/H/width(v)/zc(v)/sag(v)/thick/圆角。
    """
    y0, H = spec["y0"], spec["H"]
    hw_f, zc_f, sag_f = spec["hw"], spec["zc"], spec["sag"]
    thick, rc = spec["thick"], spec["rc"]
    y_top = y0 + H

    def surf(u, y):
        v = np.clip((y - y0) / H, 0.0, 1.0)
        x = u * hw_f(v)
        z = zc_f(v) + sag_f(v) * (1.0 - u * u)
        return x, z

    def mask(u, y):
        hw = hw_f(np.clip((y - y0) / H, 0.0, 1.0))
        x = u * hw
        ok = np.abs(u) <= 1.0
        # 上沿两角倒圆
        corner = y > (y_top - rc)
        if np.any(corner):
            dx = np.abs(x) - (hw - rc)
            dy = y - (y_top - rc)
            ok = ok & (~corner | ((np.maximum(dx, 0.0) ** 2 + np.maximum(dy, 0.0) ** 2) <= rc * rc))
        return ok

    nu = int(2.0 / step * 0.5 * float(hw_f(0.5))) + 3
    nv = int(H / step) + 3
    U, Y = np.meshgrid(np.linspace(-1.0, 1.0, max(nu, 8)),
                       np.linspace(y0, y_top, max(nv, 8)), indexing="ij")
    ok = mask(U, Y)
    u = U[ok]
    y = Y[ok]
    x, z = surf(u, y)

    # 由高度场算外法向，厚度沿法向偏置，边缘才不会出现"双层壳"
    dzdx = np.gradient(z.reshape(-1), 1.0)  # placeholder，真实梯度用解析式
    v = np.clip((y - y0) / H, 0.0, 1.0)
    dzdu = -2.0 * sag_f(v) * u / np.maximum(hw_f(v), 1e-6)
    dzdy = sag_f(v) - 2.0 * sag_f(v) * u * u  # 占位，下面用数值微分覆盖
    _ = (dzdx, dzdu, dzdy)

    # 数值法向（对 (u,y) 的参数面做差分，尺度一致后才安全）
    du = 1.0 / max(nu - 1, 1)
    dy = H / (nv - 1)
    uu = np.clip(u + du, -1.0, 1.0)
    x2, z2 = surf(uu, y)
    y3 = np.clip(y + dy, y0, y_top)
    x3, z3 = surf(u, y3)
    tx = np.stack([x2 - x, np.zeros_like(x), z2 - z], axis=1)
    ty = np.stack([x3 - x, y3 - y, z3 - z], axis=1)
    n = np.cross(tx, ty)
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-9)
    n *= np.sign(n[:, 2:3] + 1e-12)                     # 统一朝向 -Z（可喷侧）

    front = np.stack([x, y, z], axis=1)
    back = front - n * thick
    parts = [front, back]

    # 轮廓侧边：底边 / 两侧 / 上沿两角 / 顶边
    ring = _panel_outline(spec, step)
    ru, ry = ring
    rx, rz = surf(ru, ry)
    rv = np.clip((ry - y0) / H, 0.0, 1.0)
    rdup = -2.0 * sag_f(rv) * ru / np.maximum(hw_f(rv), 1e-6)
    rn = np.stack([-rdup, np.zeros_like(ru), np.ones_like(ru)], axis=1)
    rn /= np.maximum(np.linalg.norm(rn, axis=1, keepdims=True), 1e-9)
    rn *= np.sign(rn[:, 2:3] + 1e-12)
    rf = np.stack([rx, ry, rz], axis=1)
    for t in (0.25, 0.5, 0.75, 1.0):
        parts.append(rf - rn * (thick * t))
    return np.concatenate(parts, axis=0)


def _panel_outline(spec, step):
    """靠背轮廓在 (u, y) 参数面里的闭环，用于生成侧边。"""
    y0, H, rc = spec["y0"], spec["H"], spec["rc"]
    hw_f = spec["hw"]
    y_top = y0 + H
    # 侧边沿 y 走到圆角起点
    ys = np.arange(y0, y_top - rc, step)
    us_r = np.ones_like(ys)
    us_l = -np.ones_like(ys)
    # 上沿圆角（在 x-y 平面里画弧，再折回 u = x/hw）
    th = np.linspace(0.0, 0.5 * math.pi, max(4, int(0.5 * math.pi * rc / step)))
    hw_top = hw_f(1.0)
    ua, ya = [], []
    for sgn in (+1.0, -1.0):
        cx = sgn * (hw_top - rc)
        xs = cx + sgn * rc * np.cos(th)[::-1] if sgn > 0 else cx - rc * np.cos(th)
        yy = (y_top - rc) + rc * np.sin(th)
        if sgn > 0:
            xs = cx + rc * np.cos(th)
        ua.append(xs / max(hw_top, 1e-6))
        ya.append(yy)
    # 顶边（两圆角之间）
    n_top = max(3, int((2.0 * (hw_top - rc)) / step) + 1)
    ut = np.linspace(ua[0][-1], ua[1][-1], n_top)
    yt = np.full_like(ut, y_top)
    # 底边
    n_bot = max(3, int((2.0 * hw_f(0.0)) / step) + 1)
    ub = np.linspace(-1.0, 1.0, n_bot)
    yb = np.full_like(ub, y0)

    u = np.concatenate([ub, us_r, ua[0], ut, ua[1], us_l[::-1]])
    y = np.concatenate([yb, ys, ya[0], yt, ya[1], ys[::-1]])
    return u, y


# ---------------------------------------------------------------- 基本体：钢管
def smooth_path(pts, sigma_mm, passes=3, ds=1.5):
    """折线按 ds 重采样 + 多次滑动平均，角点自然变成圆弧（对 Ø25 管视觉上等同真弯管）。

    ⚠ 窗口宽度必须远小于路径长度：早先按"点数"取窗口，短路径被整段抹平成一点，
      腿直接消失。这里窗口按 mm 换算成样本数，且限制在路径长度的 1/6 以内。
    """
    pts = np.asarray(pts, dtype=float)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    total = float(seg.sum())
    n = max(int(total / ds) + 1, 8)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    t = np.linspace(0.0, total, n)
    out = np.stack([np.interp(t, s, pts[:, i]) for i in range(pts.shape[1])], axis=1)
    k = int(np.clip(round(sigma_mm / ds), 1, max(1, n // 6)))
    for _ in range(passes):
        pad = np.concatenate([np.repeat(out[:1], k, axis=0), out, np.repeat(out[-1:], k, axis=0)])
        ker = np.ones(2 * k + 1) / (2 * k + 1)
        out = np.stack([np.convolve(pad[:, i], ker, mode="same")[k:k + n]
                        for i in range(out.shape[1])], axis=1)
    return out


def sweep_tube(path, radius, step, closed=False):
    """沿折线扫掠圆截面。参考轴取 (1,0,0)，适用于 YZ 平面内的路径。"""
    path = np.asarray(path, dtype=float)
    if closed:
        path = np.vstack([path, path[:1]])
    d = np.gradient(path, axis=0)
    d /= np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-9)
    ref = np.array([1.0, 0.0, 0.0])
    n1 = np.cross(d, ref)
    bad = np.linalg.norm(n1, axis=1) < 1e-6
    n1[bad] = np.cross(d[bad], np.array([0.0, 1.0, 0.0]))
    n1 /= np.maximum(np.linalg.norm(n1, axis=1, keepdims=True), 1e-9)
    n2 = np.cross(d, n1)
    n2 /= np.maximum(np.linalg.norm(n2, axis=1, keepdims=True), 1e-9)
    n_th = max(10, int(2.0 * math.pi * radius / step))
    th = np.linspace(0.0, 2.0 * math.pi, n_th, endpoint=False)
    C, S = np.cos(th), np.sin(th)
    out = path[:, None, :] + radius * (n1[:, None, :] * C[None, :, None] + n2[:, None, :] * S[None, :, None])
    rings = [out.reshape(-1, 3)]
    # 端盖：点云里开口管看起来像个洞，补个圆盘
    for i in (0, -1):
        r_in = np.linspace(0.0, radius * 0.95, max(2, int(radius / step) + 1))
        cap = [path[i][None, :] + n1[i][None, :] * rr * C[:, None] + n2[i][None, :] * rr * S[:, None]
               for rr in r_in]
        rings.append(np.concatenate(cap, axis=0))
    return np.concatenate(rings, axis=0)


# ---------------------------------------------------------------- 椅子装配
SEAT_TOP = 452.0
SEAT_FRONT = -235.0
SEAT_BACK = 222.0


def build_chair(step: float):
    """返回 [(名称, 点(N,3))]；坐标系：椅子正放，地面 y=0，正面朝 -Z，x 左右对称。"""
    parts = []
    cz_seat = 0.5 * (SEAT_FRONT + SEAT_BACK)

    # 1) 坐垫：480×457×60，前低后高（后倾 12mm），顶面中间微凹
    parts.append(("seat", extrude_rounded_rect(
        0.0, cz_seat, 480.0, SEAT_BACK - SEAT_FRONT, 392.0, SEAT_TOP,
        r_corner=52.0, r_fillet=20.0, step=step,
        dish=3.0, slope=-10.0 / (SEAT_BACK - SEAT_FRONT), z_ref=cz_seat)))

    # 2) 坐垫下的薄盆（黑色塑料底座），比坐垫小一圈
    parts.append(("seat_pan", extrude_rounded_rect(
        0.0, cz_seat, 430.0, 420.0, 370.0, 394.0,
        r_corner=45.0, r_fillet=9.0, step=step)))

    # 3) 靠背弯板：底 440 宽，中部略宽，上沿倒圆；腰部前凸 20，整体后仰
    def hw(v):
        return 0.5 * (404.0 + 20.0 * np.sin(math.pi * np.clip(v, 0.0, 1.0)))

    def zc(v):
        v = np.clip(v, 0.0, 1.0)
        return (222.0 + 74.0 * v
                - 21.0 * np.exp(-((v - 0.32) / 0.20) ** 2)
                + 14.0 * np.exp(-((v - 0.96) / 0.10) ** 2))

    def sag(v):
        v = np.clip(v, 0.0, 1.0)
        return 19.0 + 13.0 * v

    parts.append(("back", curved_panel(
        {"y0": 448.0, "H": 582.0, "hw": hw, "zc": zc, "sag": sag,
         "thick": 13.0, "rc": 80.0}, step)))

    # 4) 扶手：黑色塑料垫 + 内衬钢管环路
    for sx in (-1.0, 1.0):
        parts.append((f"arm_pad{sx:+.0f}", extrude_rounded_rect(
            sx * 220.0, 25.0, 60.0, 250.0, 628.0, 660.0,
            r_corner=28.0, r_fillet=13.0, step=step, dish=2.0,
            slope=-0.02, z_ref=25.0)))
        loop = smooth_path([
            (sx * 208.0, 396.0, 148.0),      # 从坐垫侧后方立起来
            (sx * 208.0, 608.0, 152.0),
            (sx * 208.0, 622.0, 60.0),       # 往前（贴住扶手下表面）
            (sx * 208.0, 622.0, -64.0),
            (sx * 206.0, 556.0, -128.0),     # 前段下弯
            (sx * 204.0, 398.0, -164.0),     # 回到坐垫前侧
        ], sigma_mm=26.0)
        parts.append((f"arm_tube{sx:+.0f}", sweep_tube(loop, 11.0, step)))

    # 5) 悬臂腿：一根管从坐垫前侧下到地面、向后走、再从后半段立起来
    for sx in (-1.0, 1.0):
        def X(y, top=198.0, bot=224.0):
            t = np.clip((392.0 - y) / 392.0, 0.0, 1.0)
            return sx * (top + (bot - top) * t)

        yz = [(398.0, -178.0), (250.0, -190.0), (78.0, -196.0),
              (16.0, -140.0), (16.0, 120.0), (78.0, 168.0),
              (250.0, 176.0), (398.0, 170.0)]
        path = smooth_path(np.stack([X(np.array([p[0] for p in yz])),
                                     np.array([p[0] for p in yz]),
                                     np.array([p[1] for p in yz])], axis=1),
                           sigma_mm=24.0)
        parts.append((f"leg{sx:+.0f}", sweep_tube(path, 12.5, step)))

    # 6) 横撑：坐垫前/后各一根，地面最前面一根
    for name, y, z, x0 in (("bar_front", 392.0, -160.0, 212.0),
                           ("bar_rear", 390.0, 178.0, 212.0),
                           ("bar_floor", 16.0, -150.0, 222.0)):
        path = np.array([[sx * x0, y, z] for sx in np.linspace(-1.0, 1.0, 2)])
        path = smooth_path(path, sigma_mm=6.0)
        parts.append((name, sweep_tube(path, 10.0, step)))
    return parts


# ---------------------------------------------------------------- 虚拟深度相机
def rotate_to_camera(pts, pitch_deg, yaw_deg):
    """把正放的椅子转成"相机看到的姿态"：pitch>0 让相机略微俯视（能看到坐垫顶面）。"""
    a = math.radians(-pitch_deg)          # 绕 X：椅背朝相机倾
    b = math.radians(yaw_deg)             # 绕 Y：左右转一点，像手持扫描
    Rx = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    Ry = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
    return pts @ (Ry @ Rx).T


def to_depth_image(pts, pitch_deg, yaw_deg, origin, pitch_px):
    """正投影 + z-buffer：每个像素只留最先被打到的那个点（和真实扫描一个口径）。"""
    pts_cam = rotate_to_camera(pts, pitch_deg, yaw_deg)
    lo = pts_cam.min(axis=0)
    ij = np.floor((pts_cam[:, :2] - lo[:2]) / pitch_px).astype(np.int64)
    nx, ny = ij[:, 0].max() + 1, ij[:, 1].max() + 1
    key = ij[:, 0] * ny + ij[:, 1]
    order = np.lexsort((pts_cam[:, 2], key))     # 同像素按 z 升序 → 取最近的
    key_s = key[order]
    first = np.ones(key_s.shape, dtype=bool)
    first[1:] = key_s[1:] != key_s[:-1]
    idx = order[first]
    return pts_cam[idx], lo, nx, ny


# ---------------------------------------------------------------- 输出
def write_ply(path, pts, intensity=1.0):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("ply\nformat ascii 1.0\n")
        fh.write("comment synthetic visitor chair cloud (unit: mm)\n")
        fh.write(f"element vertex {pts.shape[0]}\n")
        fh.write("property float x\nproperty float y\nproperty float z\n")
        fh.write("property float intensity\n")
        fh.write("end_header\n")
        for p in pts:
            fh.write(f"{p[0]:.3f} {p[1]:.3f} {p[2]:.3f} {intensity:g}\n")


def save_preview(path, pts, parts_info=""):
    """预览：正视 / 3/4 立体 / 侧视，都是"相机看到的那一面"。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lo, hi = pts.min(axis=0), pts.max(axis=0)
    ctr, ext = (lo + hi) / 2.0, np.maximum(hi - lo, 1e-6)
    depth = np.clip((pts[:, 2] - lo[2]) / max(ext[2], 1e-6), 0.0, 1.0)

    fig = plt.figure(figsize=(15.4, 5.4), dpi=140)

    ax1 = fig.add_subplot(1, 3, 1)
    ax1.scatter(pts[:, 0], pts[:, 1], s=2.4, c=depth, cmap="Greys", linewidths=0)
    ax1.set_aspect("equal")
    ax1.set_title("front view  (darker = nearer to camera)", fontsize=10)
    ax1.set_xlabel("X (mm)")
    ax1.set_ylabel("Y (mm)")
    ax1.grid(True, linewidth=0.3, alpha=0.4)

    # 3/4 视角：把点云绕 Y 转 34°、绕 X 抬 14° 后正投影，等价"从斜前方拍一张"
    a, b = math.radians(14.0), math.radians(34.0)
    Rx = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    Ry = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
    q = pts @ (Ry @ Rx).T
    q = q - (q.min(axis=0) + q.max(axis=0)) / 2.0
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.scatter(q[:, 0], q[:, 1], s=2.4, c=depth, cmap="Greys", linewidths=0)
    ax2.set_aspect("equal")
    ax2.set_title("3/4 view", fontsize=10)
    ax2.set_xlabel("(rotated)")
    ax2.grid(True, linewidth=0.3, alpha=0.4)
    ax2.set_xticks([])
    ax2.set_yticks([])

    ax3 = fig.add_subplot(1, 3, 3)
    ax3.scatter(pts[:, 2], pts[:, 1], s=2.4, c=depth, cmap="Greys", linewidths=0)
    ax3.set_aspect("equal")
    ax3.set_title("side view  (Z = depth)", fontsize=10)
    ax3.set_xlabel("Z (mm)")
    ax3.set_ylabel("Y (mm)")
    ax3.grid(True, linewidth=0.3, alpha=0.4)

    fig.suptitle(f"Visitor chair point cloud — {len(pts)} pts  {parts_info}", fontsize=11)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.90, bottom=0.08, wspace=0.26)
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description="生成访客椅点云（mm）")
    ap.add_argument("--pitch-step", type=float, default=2.6, help="深度图像素间距 mm，默认 2.6（≈真实扫描）")
    ap.add_argument("--pitch", type=float, default=12.0, help="虚拟相机俯视角度，默认 15°")
    ap.add_argument("--yaw", type=float, default=10.0, help="虚拟相机左右偏角，默认 10°")
    ap.add_argument("--noise", type=float, default=0.35, help="各向同性噪声 σ mm")
    ap.add_argument("--dropout", type=float, default=0.03, help="随机丢点比例")
    ap.add_argument("--origin", type=float, nargs=3, default=DEFAULT_ORIGIN,
                    metavar=("X", "Y", "Z"), help="椅脚中心在相机系的位置")
    ap.add_argument("--seed", type=int, default=20261124)
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-preview", dest="preview", action="store_false")
    args = ap.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    step = min(1.4, args.pitch_step * 0.5)
    parts = build_chair(step)
    all_pts = np.concatenate([p for _, p in parts], axis=0)
    print(f"建模点数  : {all_pts.shape[0]}（{len(parts)} 个部件，采样步长 {step:.2f}mm）")

    img, lo, nx, ny = to_depth_image(all_pts, args.pitch, args.yaw, args.origin, args.pitch_step)
    print(f"深度图    : {nx} × {ny} px @ {args.pitch_step}mm → {img.shape[0]} 点（单视角，含遮挡）")

    # 深度图 → 相机系坐标：像素中心 + 原定位 + 噪声
    ij = np.floor((img[:, :2] - lo[:2]) / args.pitch_step)
    cam = np.stack([lo[0] + (ij[:, 0] + 0.5) * args.pitch_step,
                    lo[1] + (ij[:, 1] + 0.5) * args.pitch_step,
                    img[:, 2]], axis=1)
    cam += rng.normal(0.0, max(args.noise, 0.0), cam.shape)
    keep = rng.random(cam.shape[0]) >= args.dropout
    cam = cam[keep]

    # 平移到相机系：把"椅脚中心"搬到 --origin（旋转后 (0,0,0) 仍是原点）
    cam = cam + np.asarray(args.origin, dtype=float) - np.array([0.0, 0.0, 0.0])

    out = args.out or os.path.join(PROJECT_ROOT, "outputs", "chair_visitor.ply")
    write_ply(out, cam)

    lo2, hi2 = cam.min(axis=0), cam.max(axis=0)
    print(f"点数      : {cam.shape[0]}")
    print(f"包围盒 min: {lo2[0]:.1f}, {lo2[1]:.1f}, {lo2[2]:.1f}")
    print(f"包围盒 max: {hi2[0]:.1f}, {hi2[1]:.1f}, {hi2[2]:.1f}")
    print(f"尺寸 XYZ  : {hi2[0]-lo2[0]:.1f} x {hi2[1]-lo2[1]:.1f} x {hi2[2]-lo2[2]:.1f} mm")
    print(f"输出      : {out}")

    if args.preview:
        png = os.path.splitext(out)[0] + "_preview.png"
        save_preview(png, cam, f"(pitch {args.pitch:g}°, yaw {args.yaw:g}°, {args.pitch_step:g}mm)")
        print(f"预览      : {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
