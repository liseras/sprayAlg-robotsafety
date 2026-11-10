# -*- coding: utf-8 -*-
"""合成测试件：每个件都有"标准答案"，用来暴露算法在边界条件下的行为

生成 ascii PLY（单层点云，模拟激光只打到可及的一面）。

  rect     平板 300x200            期望：行距=设定值、四边覆盖、行数=宽/行距+1
  hole     平板 300x200 中间挖 100x100 通孔
                -fill 300（默认）期望：跨孔连续（用户的"洞当实心"语义）
                -fill 50          期望：绕开孔，孔内无轨迹
  narrow   窄条 300x6             期望：面宽 < 行距，应出 1 行（窄面单道）
  lbracket L 形件 两片 200x150 成 90°  期望：两面各出轨迹、行方向各自最优
  big      大板 1000x800           期望：行距均匀、无累积漂移
  ring     环形窄面 外径200 内径170  期望：同心圆（环宽30 = 2.5 行距）
"""
import sys
import numpy as np


def write_ply(path, P):
    with open(path, "w") as f:
        f.write("ply\nformat ascii 1.0\n")
        f.write("element vertex %d\n" % len(P))
        f.write("property float x\nproperty float y\nproperty float z\nend_header\n")
        np.savetxt(f, P, fmt="%.3f")


def grid_rect(x0, x1, y0, y1, step=1.0, z=0.0):
    xs = np.arange(x0, x1 + 1e-9, step)
    ys = np.arange(y0, y1 + 1e-9, step)
    X, Y = np.meshgrid(xs, ys)
    return np.c_[X.ravel(), Y.ravel(), np.full(X.size, z)]


def make_rect():
    return grid_rect(-150, 150, -100, 100, 1.0)


def make_hole():
    P = grid_rect(-150, 150, -100, 100, 1.0)
    keep = ~((np.abs(P[:, 0]) < 50) & (np.abs(P[:, 1]) < 50))
    return P[keep]


def make_narrow():
    return grid_rect(-150, 150, -3, 3, 0.5)


def make_big():
    return grid_rect(-500, 500, -400, 400, 2.0)


def make_ring():
    P = []
    for r in np.arange(85.0, 100.0, 0.5):
        n = max(16, int(2 * np.pi * r / 0.8))
        t = np.linspace(0, 2 * np.pi, n, endpoint=False)
        P.append(np.c_[r * np.cos(t), r * np.sin(t), np.zeros(n)])
    return np.vstack(P)


def make_lbracket():
    a = grid_rect(-100, 100, -75, 75, 1.0)          # 底面 z=0
    b = grid_rect(-100, 100, -75, 75, 1.0)          # 立板
    b[:, 2] = b[:, 1] + 75.0                        # 沿 y 折起 90°
    b[:, 1] = -75.0
    return np.vstack([a, b])


CASES = {
    "rect": make_rect,
    "hole": make_hole,
    "narrow": make_narrow,
    "big": make_big,
    "ring": make_ring,
    "lbracket": make_lbracket,
}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    for name, fn in CASES.items():
        P = fn()
        write_ply("%s/syn_%s.ply" % (out, name), P)
        print("%-10s %8d 点  bbox=%s" % (name, len(P),
              np.round(P.max(0) - P.min(0), 1)))


if __name__ == "__main__":
    main()
