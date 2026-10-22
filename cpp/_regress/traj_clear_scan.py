# -*- coding: utf-8 -*-
"""轨迹枪位安全距离体检

用法:
    python traj_clear_scan.py <工件.ply> <轨迹.csv> [--clear 20] [--top 25]

三项检查：
  1. 点级：每个轨迹点（枪位）到工件点云的最近距离
  2. 折线级：每两个相邻轨迹点之间的连线到工件的最近距离（**解析化的点到线段距离**）
     —— 两端都安全不代表中间安全，面间过渡/跨段跳跃都是"两点一直线"，
        中间穿过壁时，只看端点完全看不出来，这正是三维里"穿模"的样子。
  3. 连续侵入段：起止序号、长度、最小距离、坐标范围

实现：工件点云建 20mm 网格哈希，查询时只取 AABB 膨胀 clear 范围内的候选点，
再算解析距离。二十七万点的点云也是秒级。
"""
import sys
import csv
import numpy as np

CELL = 20.0


def read_ply(path):
    with open(path, "r", errors="ignore") as f:
        n = 0
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("element vertex"):
                n = int(line.split()[-1])
            if line.strip() == "end_header":
                break
        data = np.loadtxt(f, max_rows=n, dtype=np.float64)
    return np.ascontiguousarray(data[:, :3])


def read_traj(path):
    P = []
    with open(path) as f:
        rd = csv.reader(f)
        head = next(rd)
        idx = {k: i for i, k in enumerate(head)}
        for rec in rd:
            P.append((float(rec[idx["X"]]), float(rec[idx["Y"]]), float(rec[idx["Z"]]),
                      int(float(rec[idx["spray"]]))))
    P = np.array(P)
    return np.ascontiguousarray(P[:, :3]), P[:, 3].astype(int)


class Grid:
    def __init__(self, Q, cell=CELL):
        self.Q = Q
        self.cell = cell
        self.org = Q.min(axis=0)
        ijk = np.floor((Q - self.org) / cell).astype(np.int64)
        self.dim = ijk.max(axis=0) + 1
        key = (ijk[:, 0] * self.dim[1] + ijk[:, 1]) * self.dim[2] + ijk[:, 2]
        order = np.argsort(key, kind="stable")
        self.key = key[order]
        self.idx = order

    def cand(self, lo, hi):
        """返回 AABB[lo,hi] 内所有点的下标"""
        c = self.cell
        i0 = np.floor((lo - self.org) / c).astype(np.int64)
        i1 = np.floor((hi - self.org) / c).astype(np.int64)
        i0 = np.maximum(i0, 0)
        i1 = np.minimum(i1, self.dim - 1)
        if np.any(i1 < i0):
            return None
        out = []
        dy, dz = self.dim[1], self.dim[2]
        for ix in range(i0[0], i1[0] + 1):
            for iy in range(i0[1], i1[1] + 1):
                k0 = (ix * dy + iy) * dz + i0[2]
                k1 = (ix * dy + iy) * dz + i1[2]
                a = np.searchsorted(self.key, k0, "left")
                b = np.searchsorted(self.key, k1, "right")
                if b > a:
                    out.append(self.idx[a:b])
        if not out:
            return None
        return np.concatenate(out)


def pt_seg_dist(P, A, B):
    """P: (n,3)   A,B: (3,)  → 点到线段的最近距离"""
    AB = B - A
    L2 = float(AB @ AB)
    if L2 < 1e-12:
        return np.linalg.norm(P - A, axis=1)
    t = ((P - A) @ AB) / L2
    t = np.clip(t, 0.0, 1.0)
    C = A + t[:, None] * AB
    return np.linalg.norm(P - C, axis=1)


def main():
    clear, topn = 20.0, 25
    pos = []
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--clear":
            clear = float(args[i + 1]); i += 2
        elif args[i] == "--top":
            topn = int(args[i + 1]); i += 2
        else:
            pos.append(args[i]); i += 1
    if len(pos) < 2:
        print(__doc__); sys.exit(1)

    Q = read_ply(pos[0])
    G, sp = read_traj(pos[1])
    grid = Grid(Q)
    print("part=%d 点   traj=%d 点   spray=1 %d / spray=0 %d   安全半径 %.1f mm"
          % (len(Q), len(G), (sp == 1).sum(), (sp == 0).sum(), clear))

    # ---- 1) 点级 ----
    dp = np.empty(len(G))
    for k in range(len(G)):
        c = grid.cand(G[k] - clear, G[k] + clear)
        if c is None:
            dp[k] = 1e9
        else:
            dp[k] = np.linalg.norm(Q[c] - G[k], axis=1).min()
    for tag, m in (("全部", np.ones(len(G), bool)), ("喷涂段", sp == 1), ("空移段", sp == 0)):
        if m.sum() == 0:
            continue
        dd = dp[m]
        print("  %s: mean=%.1f p5=%.1f min=%.1f  <%.0fmm %d 点 (%.2f%%)"
              % (tag, dd.mean(), np.percentile(dd, 5), dd.min(), clear,
                 (dd < clear).sum(), 100.0 * (dd < clear).sum() / len(dd)))

    order = np.argsort(dp)[:topn]
    print("\n-- 最贴件的 %d 个点 --" % topn)
    print("  idx      X        Y        Z     spray   dist")
    for k in sorted(order):
        print("  %5d %8.1f %8.1f %8.1f    %d    %6.2f" % (k, G[k, 0], G[k, 1], G[k, 2], sp[k], dp[k]))

    # ---- 2) 折线级 ----
    seg_len = np.linalg.norm(G[1:] - G[:-1], axis=1)
    dsg = np.full(len(G) - 1, 1e9)
    for k in range(len(G) - 1):
        A, B = G[k], G[k + 1]
        lo = np.minimum(A, B) - clear
        hi = np.maximum(A, B) + clear
        c = grid.cand(lo, hi)
        if c is not None:
            dsg[k] = pt_seg_dist(Q[c], A, B).min()
    bad = dsg < clear
    print("\n-- 折线中段侵入（相邻两点连线 < %.0fmm）--" % clear)
    idxs = np.where(bad)[0]
    if len(idxs):
        print("  共 %d 段，总长 %.1f mm  最差 %.2f mm"
              % (len(idxs), seg_len[idxs].sum(), dsg[idxs].min()))
        for k in idxs[np.argsort(dsg[idxs])][:15]:
            print("  seg#%5d  len=%7.1f  min=%.2f  spray=%d->%d  %s -> %s"
                  % (k, seg_len[k], dsg[k], sp[k], sp[k + 1],
                     np.round(G[k], 1).tolist(), np.round(G[k + 1], 1).tolist()))
    else:
        print("  （无）")

    # ---- 3) 连续侵入段（按点级）----
    badp = dp < clear
    print("\n-- 连续侵入段（点级，>=3 点）--")
    runs = []
    k, n = 0, len(badp)
    while k < n:
        if not badp[k]:
            k += 1; continue
        j = k
        while j + 1 < n and badp[j + 1]:
            j += 1
        if j - k + 1 >= 3:
            runs.append((k, j))
        k = j + 1
    runs.sort(key=lambda r: dp[r[0]:r[1] + 1].min())
    for k, j in runs[:20]:
        seg = G[k:j + 1]
        print("  #%5d..%-5d  n=%-4d min=%.2f  spray=%s  X[%.0f,%.0f] Y[%.0f,%.0f] Z[%.0f,%.0f]"
              % (k, j, j - k + 1, dp[k:j + 1].min(),
                 "1" if sp[k] == 1 else "0",
                 seg[:, 0].min(), seg[:, 0].max(), seg[:, 1].min(), seg[:, 1].max(),
                 seg[:, 2].min(), seg[:, 2].max()))
    if not runs:
        print("  （无）")


if __name__ == "__main__":
    main()
