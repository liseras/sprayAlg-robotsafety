# -*- coding: utf-8 -*-
"""覆盖完整性检查：点云每个点到最近"喷涂"轨迹点的距离分布。

用途——验证"按重建轮廓修正行"这类改动是否真的漏喷。只看轨迹点数会误判：
点数减少可能只是夹掉了面外的空走，也可能是真漏了边沿。
判据：点云各点到最近喷点的距离，取 p50/p90/p99 与"距离 <= 行距/2 的比例"。

用法：
    python tools/traj_coverage.py <cloud.ply> <traj.csv> [--row 60]
"""
import sys

import numpy as np


def read_ply(path):
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.strip() == 'end_header':
                break
        return np.loadtxt(fh)[:, :3]


def read_traj(path):
    a = np.loadtxt(path, delimiter=',', skiprows=1)
    return a[:, :3], a[:, 7].astype(int)


def main():
    cloud, csv = sys.argv[1], sys.argv[2]
    row = 60.0
    if '--row' in sys.argv:
        row = float(sys.argv[sys.argv.index('--row') + 1])
    Q = read_ply(cloud)
    G, sp = read_traj(csv)
    S = G[sp == 1]
    if len(S) == 0:
        print('无喷涂点')
        return
    # 网格加速（1/3 行距的格子），避免 O(N*M)
    cell = max(5.0, row / 3.0)
    origin = np.minimum(Q.min(0), S.min(0)) - cell
    idxS = np.floor((S - origin) / cell).astype(int)
    buckets = {}
    for k, key in enumerate(map(tuple, idxS)):
        buckets.setdefault(key, []).append(k)
    d = np.full(len(Q), 1e9)
    idxQ = np.floor((Q - origin) / cell).astype(int)
    for i, key in enumerate(map(tuple, idxQ)):
        best = 1e9
        r = 1
        while True:
            found = False
            for dx in range(-r, r + 1):
                for dy in range(-r, r + 1):
                    for dz in range(-r, r + 1):
                        if max(abs(dx), abs(dy), abs(dz)) != r:
                            continue
                        lst = buckets.get((key[0] + dx, key[1] + dy, key[2] + dz))
                        if not lst:
                            continue
                        found = True
                        cand = S[lst]
                        best = min(best, float(np.linalg.norm(cand - Q[i], axis=1).min()))
            # 已经找到的最近距离小于"已搜索半径"就可以停（保证不再有更近的）
            if best < (r - 1) * cell or r > 12:
                break
            r += 1
        d[i] = best
    print('点云 %d 点，喷涂点 %d，行距 %.0f' % (len(Q), len(S), row))
    print('到最近喷点距离：p50=%.1f p90=%.1f p99=%.1f max=%.1f mm'
          % tuple(np.percentile(d, [50, 90, 99, 100])))
    for lim in (row / 2.0, row, 1.5 * row):
        print('  距离 <= %.0fmm（%.2f×行距）占比 %.1f%%' % (lim, lim / row, 100.0 * (d <= lim).mean()))


if __name__ == '__main__':
    main()
