# -*- coding: utf-8 -*-
"""曲面枪位干净度体检（第十三轮新增）

用法:
    SPRAYTRAJ_DUMP=1 spray_driver.exe <曲面件>.ply out.csv out.ply -flat 40 -fill 300
    python gun_cleanliness.py surf_rows.csv [--erow x,y,z] [--eidx u,v,w] [--standoff 100]

surf_rows.csv 由 SPRAYTRAJ_DUMP=1 导出（row,k,x,y,z,s,nx,ny,nz），
x,y,z 是**贴面点**，法向已做完行内平滑与跨行一致化；枪位 = 点 + standoff·n。

三个指标（都只与姿态有关，与行位置无关）：
  gun_jitter   枪位偏离"行内 25 点滑动平均基准"的距离 → 弯弯绕绕的量化
  crossing     相邻行在同一沿程位置时，枪位横向坐标是否反序 → "相互缠绕"的量化
  normal_step  相邻行点法向夹角 → 姿态连续性

曲面主轴 erow/eidx 从 `SPRAYTRAJ_DUMP=1` 的 `[曲面dbg] ... eRow=(..) eIdx=(..)` 行拷贝；
默认值是曲型板 u08 的实测值。
"""
import sys
import csv
import numpy as np


def parse_vec(s):
    v = np.array([float(x) for x in s.split(",")])
    return v / np.linalg.norm(v)


def main(csv_path, eRow, eIdx, standoff):
    rows = {}
    with open(csv_path) as f:
        rd = csv.reader(f)
        next(rd)
        for rec in rd:
            rows.setdefault(int(rec[0]), []).append(
                (np.array([float(rec[2]), float(rec[3]), float(rec[4])]),
                 np.array([float(rec[6]), float(rec[7]), float(rec[8])])))
    keys = sorted(rows)
    T, S, N, RID, P = [], [], [], [], []
    for r in keys:
        for p, n in rows[r]:
            P.append(p)
            N.append(n)
            T.append(p.dot(eRow))
            S.append(p.dot(eIdx))
            RID.append(r)
    P, N, T, S, RID = map(np.array, (P, N, T, S, RID))
    sg = (P + standoff * N).dot(eIdx)

    # 1) gun jitter
    off = []
    for r in keys:
        m = RID == r
        k = 25
        ker = np.ones(k) / k
        Ns = np.stack([np.convolve(N[m][:, i], ker, mode="same") for i in range(3)], 1)
        ln = np.linalg.norm(Ns, axis=1)
        ln[ln < 1e-9] = 1
        Ns /= ln[:, None]
        off.append(np.linalg.norm((P[m] + standoff * N[m]) - (P[m] + standoff * Ns),
                                  axis=1)[k // 2:-(k // 2)])
    off = np.concatenate(off)

    # 2) crossing between adjacent rows (row index small should have smaller s)
    cross = tot = 0
    for r0, r1 in zip(keys[:-1], keys[1:]):
        m0, m1 = RID == r0, RID == r1
        t0, t1 = T[m0], T[m1]
        lo, hi = max(t0.min(), t1.min()), min(t0.max(), t1.max())
        if hi - lo < 40:
            continue
        for x in np.arange(lo, hi, 10.0):
            ka, kb = np.abs(t0 - x) < 6, np.abs(t1 - x) < 6
            if ka.sum() < 2 or kb.sum() < 2:
                continue
            tot += 1
            if np.median(sg[m0][ka]) > np.median(sg[m1][kb]):
                cross += 1

    # 3) normal step along rows
    angs = []
    for r in keys:
        n = N[RID == r]
        angs.append(np.degrees(np.arccos(np.clip((n[1:] * n[:-1]).sum(1), -1, 1))))
    angs = np.concatenate(angs)

    gap = np.median(S[RID == keys[1]]) - np.median(S[RID == keys[0]])   # 贴面行距，不是枪位差
    print("rows=%d points=%d  (%s)" % (len(keys), len(P), csv_path))
    print("  gun_jitter   mean=%.1f  p95=%.1f  max=%.1f mm   (compare with row gap %.2f mm)"
          % (off.mean(), np.percentile(off, 95), off.max(), gap))
    print("  crossing     %d/%d = %.0f%%  (adjacent rows where gun s is out of order)"
          % (cross, tot, 100.0 * cross / max(1, tot)))
    print("  normal_step  mean=%.2f  p95=%.2f  max=%.2f deg"
          % (angs.mean(), np.percentile(angs, 95), angs.max()))
    print("\nreference (13th round, curved panel u08):")
    print("  before: jitter p95=17.6 max=39.0, crossing=14%, normal_step p95=5.98")
    print("  after : jitter p95=4.6  max=10.7, crossing=0%,  normal_step p95=2.14")


if __name__ == "__main__":
    eRow = np.array([-0.338, -0.941, 0.026])
    eIdx = np.array([0.755, -0.255, 0.604])
    standoff = 100.0
    args = [a for a in sys.argv[1:]]
    pos = []
    i = 0
    while i < len(args):
        if args[i] == "--erow":
            eRow = parse_vec(args[i + 1]); i += 2
        elif args[i] == "--eidx":
            eIdx = parse_vec(args[i + 1]); i += 2
        elif args[i] == "--standoff":
            standoff = float(args[i + 1]); i += 2
        else:
            pos.append(args[i]); i += 1
    if not pos:
        print(__doc__)
        sys.exit(1)
    eRow = eRow / np.linalg.norm(eRow)
    eIdx = eIdx / np.linalg.norm(eIdx)
    main(pos[0], eRow, eIdx, standoff)
