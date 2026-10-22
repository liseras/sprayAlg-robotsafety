# -*- coding: utf-8 -*-
"""复现 SOR 去噪，检查壁顶/棱边点是否被误删。

用法: python sor_check.py <ply> [k] [z]
"""
import sys
import numpy as np

P = np.loadtxt(sys.argv[1], skiprows=9, max_rows=34580, dtype=np.float32)
k = int(sys.argv[2]) if len(sys.argv) > 2 else 20
z = float(sys.argv[3]) if len(sys.argv) > 3 else 2.5

n = len(P)
mea = np.empty(n, dtype=np.float32)
CH = 512
for s in range(0, n, CH):
    e = min(n, s + CH)
    d2 = ((P[s:e, None, :] - P[None, :, :]) ** 2).sum(-1)
    # 排除自身：把自身距离设为 +inf，取前 k 个
    np.fill_diagonal(d2[s:e], np.inf) if e - s == d2.shape[0] else None
    idx = np.argpartition(d2, k, axis=1)[:, :k]
    mea[s:e] = np.sqrt(np.take_along_axis(d2, idx, 1).mean(1))

mu, sd = mea.mean(), mea.std()
thr = mu + z * sd
print("n=%d  mean=%.4f std=%.4f thr=%.4f  removed=%d" % (n, mu, sd, thr, (mea > thr).sum()))

# 检查关键点是否存活
targets = [("壁顶(400,200,300)", [400, 200, 300]), ("壁顶(0,100,300)", [0, 100, 300]),
           ("壁中(400,200,150)", [400, 200, 150]), ("壁顶(400,208,300)", [400, 208, 300])]
for name, t in targets:
    d = np.linalg.norm(P - np.array(t, np.float32), axis=1)
    j = int(np.argmin(d))
    print("%-20s 最近点 %s  d=%.2f  meanKNN=%.4f  %s"
          % (name, np.round(P[j], 1).tolist(), d[j], mea[j],
             "存活" if mea[j] <= thr else "被删"))

# 被删点的空间分布
rm = P[mea > thr]
if len(rm):
    print("被删点 bbox:", np.round(rm.min(0), 1), np.round(rm.max(0), 1))
    zs = rm[:, 2]
    print("被删点 z 分布: >290 占 %.0f%%, <10 占 %.0f%%"
          % (100 * (zs > 290).mean(), 100 * (zs < 10).mean()))
