#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""轨迹规整度量化：行方向一致性、行距均匀性、行覆盖范围。"""
import sys, os, glob
import numpy as np
from collections import Counter


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
    F = [[int(x) for x in lines[hdr + nv + k].split()[1:]] for k in range(nf)]
    return V, F


def read_csv(path):
    out = []
    with open(path, "r", errors="replace") as f:
        f.readline()
        for ln in f:
            p = ln.strip().split(",")
            if len(p) >= 8:
                out.append([float(v) for v in p[:6]] + [float(p[6]), int(p[7])])
    return np.array(out)


def main():
    prefix = sys.argv[1]
    mesh = glob.glob(prefix + ".plane0.ply")
    traj = read_csv(prefix[:-4] + ".csv")
    T = traj[traj[:, 7] > 0.5][:, :3]

    # 面基准
    if mesh:
        V, F = read_ply_ascii(mesh[0])
        c = V.mean(0)
        _, _, vt = np.linalg.svd(V - c, full_matrices=False)
        u, v = vt[0], vt[1]
        print("面顶点=%d 三角=%d" % (len(V), len(F)))
    else:
        c = T.mean(0)
        _, _, vt = np.linalg.svd(T - c, full_matrices=False)
        u, v = vt[0], vt[1]

    D = T - c
    U, Vv = D @ u, D @ v

    # 行方向：用相邻点位移的主方向
    d = np.diff(np.stack([U, Vv], 1), axis=0)
    seg = np.hypot(d[:, 0], d[:, 1])
    longseg = d[seg > 2.0]
    ang = np.degrees(np.arctan2(longseg[:, 1], longseg[:, 0]))
    # 归一到 [-90,90)
    ang = ((ang + 90) % 180) - 90
    ref = np.median(ang)
    dev = np.abs(((ang - ref + 90) % 180) - 90)
    print("行方向中位数=%.2f°  偏离>5°的点段占比=%.1f%%" %
          (ref, 100.0 * (dev > 5).sum() / max(1, len(dev))))

    # 旋转到行方向
    th = np.radians(ref)
    R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
    P = np.stack([U, Vv], 1) @ R.T
    # 按 v 分簇成行
    order = np.argsort(P[:, 1])
    rows, cur = [], [order[0]]
    for i in order[1:]:
        if P[i, 1] - P[cur[-1], 1] <= 2.0:
            cur.append(i)
        else:
            rows.append(cur)
            cur = [i]
    rows.append(cur)
    ys, spans, npts = [], [], []
    for r in rows:
        ys.append(np.median(P[r, 1]))
        spans.append((P[r, 0].min(), P[r, 0].max(), len(r)))
        npts.append(len(r))
    ys = np.array(ys)
    dy = np.diff(ys)
    print("行数=%d 行距: 均值=%.2f 最小=%.2f 最大=%.2f 标准差=%.2f" %
          (len(rows), dy.mean(), dy.min(), dy.max(), dy.std()))
    print("行内点数: 均值=%.1f 最小=%d 最大=%d" % (np.mean(npts), min(npts), max(npts)))
    u0 = min(s[0] for s in spans); u1 = max(s[1] for s in spans)
    print("行覆盖 u=[%.1f, %.1f] (宽 %.1f)" % (u0, u1, u1 - u0))
    if mesh:
        Dm = V - c
        Um, Vm = Dm @ u, Dm @ v
        Q = np.stack([Um, Vm], 1) @ R.T
        print("面覆盖 u=[%.1f, %.1f] v=[%.1f, %.1f]" %
              (Q[:, 0].min(), Q[:, 0].max(), Q[:, 1].min(), Q[:, 1].max()))
        print("行覆盖 v=[%.1f, %.1f]" % (ys.min(), ys.max()))
        # 左右端缺口
        print("左右端缺口: 左=%.1f 右=%.1f" %
              (min(s[0] for s in spans) - Q[:, 0].min(),
               Q[:, 0].max() - max(s[1] for s in spans)))
    # 打印每行
    for k, (y, s) in enumerate(zip(ys, spans)):
        print("  row%02d v=%7.2f  u=[%7.2f, %7.2f] len=%6.1f n=%d" %
              (k, y, s[0], s[1], s[1] - s[0], s[2]))


if __name__ == "__main__":
    main()
