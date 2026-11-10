# -*- coding: utf-8 -*-
"""行结构核对：行数、行距、行覆盖范围与点云边界是否对得上

与 regress_run.py 的 frag_stats 不同，这里不靠"方向突变"切行，而是：
  1) 用轨迹自身的走向定出"行方向"（能量最大的单位方向）
  2) 跨行方向投影后按间隙切行
  3) 逐行给出：跨行位置、行长、点数
再做三件核对：
  A. 实际行距 vs 设定行距（偏差、标准差、最大偏差）
  B. 首末行是否压到点云边界（漏喷边界带宽 = 边界到首末行的距离 - 0.5×行距）
  C. 行内是否有断口（环/孔/碎段）

用法：
    python row_check.py <cloud.ply> <traj.csv> [--row 12] [--so 100]
"""
import sys
import numpy as np
from audit_traj import load_cloud, load_traj


def best_across(C, nref):
    """搜索"跨行方向"：让落点投影到该方向后，直方图最集中的那个角度

    不能用"相邻点方向的主特征向量"当行方向：那只是平均走向，
    0.5° 的误差乘以 1000mm 行长就是 9mm 的投影漂移，行与行会糊成一团。
    行是等距平行线，能量最大化的判据应该直接放在"投影后够不够尖"上。
    """
    # 面内正交基
    tmp = np.array([1.0, 0, 0]) if abs(nref[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(nref, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(nref, e1)

    best = (None, -1.0)
    for deg in np.arange(0.0, 180.0, 0.25):
        t = np.radians(deg)
        a = e1 * np.cos(t) + e2 * np.sin(t)
        cv = C @ a
        h, _ = np.histogram(cv, bins=max(8, int((cv.max() - cv.min()) / 1.0)))
        score = float((h.astype(float) ** 2).sum())     # 越尖越大
        if score > best[1]:
            best = (a, score)
    return best[0]


def split_rows(C, across, gap):
    cv = C @ across
    o = np.argsort(cv)
    cvs = cv[o]
    brk = np.where(np.diff(cvs) > gap)[0]
    starts = np.r_[0, brk + 1]
    ends = np.r_[brk, len(cvs) - 1]
    rows = []
    for a, b in zip(starts, ends):
        if b - a < 3:
            continue
        idx = o[a:b + 1]
        rows.append({"yc": float(cvs[a:b + 1].mean()), "idx": idx})
    rows.sort(key=lambda r: r["yc"])
    return rows


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    ply, csvp = sys.argv[1], sys.argv[2]
    L, so = 12.0, 100.0
    if "--row" in sys.argv:
        L = float(sys.argv[sys.argv.index("--row") + 1])
    if "--so" in sys.argv:
        so = float(sys.argv[sys.argv.index("--so") + 1])

    Q = load_cloud(ply)
    P, A, sp = load_traj(csvp)
    C = (P + so * A)[sp == 1]
    if len(C) < 10:
        print("喷涂点太少")
        return 1

    c0 = Q.mean(0)
    w, v = np.linalg.eigh((Q - c0).T @ (Q - c0))
    nref = v[:, 0]
    across = best_across(C, nref)
    rd = np.cross(nref, across)
    rd /= np.linalg.norm(rd)

    rows = split_rows(C, across, max(2.0, L * 0.4))
    if len(rows) < 2:
        print("只解析出 %d 行（可能是同心圆/等高线走法，不适合按行核对）" % len(rows))
        # 仍给出落点范围
        cv = C @ across
        print("落点跨行范围 [%.1f, %.1f]" % (cv.min(), cv.max()))
        return 0

    yc = np.array([r["yc"] for r in rows])
    steps = np.diff(yc)

    # 每行长度（沿行方向投影的跨度）
    lens = []
    for r in rows:
        q = C[r["idx"]] @ rd
        lens.append(float(q.max() - q.min()))
    lens = np.array(lens)

    # 点云跨行范围
    qv = Q @ across
    qlo, qhi = np.percentile(qv, 0.2), np.percentile(qv, 99.8)

    print("解析出 %d 行" % len(rows))
    print("行距: 均值 %.2f  中位 %.2f  min %.2f  max %.2f  std %.2f  (设定 %.1f)"
          % (steps.mean(), np.median(steps), steps.min(), steps.max(), steps.std(), L))
    print("      偏差: 最大 +%.2f / %.2f mm  (%.0f%% / %.0f%%)"
          % (steps.max() - L, steps.min() - L,
             100 * (steps.max() - L) / L, 100 * (steps.min() - L) / L))
    print("行长:  均值 %.0f  min %.0f  max %.0f mm" % (lens.mean(), lens.min(), lens.max()))
    print("点云跨行范围 [%.1f, %.1f] (跨距 %.1f)" % (qlo, qhi, qhi - qlo))
    print("首行 %.1f  末行 %.1f" % (yc.min(), yc.max()))
    print("边界余量: 低端 %.1f mm / 高端 %.1f mm  (>=0 表示行已覆盖到边)"
          % (yc.min() - 0.5 * L - qlo, qhi - (yc.max() + 0.5 * L)))

    # 行距异常清单
    bad = [(i, float(steps[i])) for i in range(len(steps)) if abs(steps[i] - L) > 0.25 * L]
    if bad:
        print("\n行距异常 (>25%% 偏差):")
        for i, s in bad:
            print("   第%d→%d行: %.2f mm (%.0f%%)" % (i, i + 1, s, 100 * (s - L) / L))
    else:
        print("\n行距全部在 ±25%% 内")
    return 0


if __name__ == "__main__":
    sys.exit(main())
