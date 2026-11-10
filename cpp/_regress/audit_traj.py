# -*- coding: utf-8 -*-
"""轨迹质量独立审计

不复用算法自身的覆盖率/残差口径，只吃两样东西：
  1) 原始点云 PLY
  2) 导出轨迹 CSV（X,Y,Z,I,J,K,speed,spray）

反算四项与工艺直接相关的指标：

  覆盖率  把每个喷枪点沿枪轴下降 standoff，得到"落点中心"（理论上落回工件表面）。
          表面点到最近落点中心的距离 d：
              d <= 0.5×行距 → 无搭接也必被覆盖
              d >  0.5×行距 → 该处存在 >= 一个行距的空档
          报告：整体覆盖率(0.5L)、局部覆盖率(仅统计枪口曾经过的区域)、
                距离 p50/p95/max（max/1.0 = 最大漏喷半宽，单位 mm）

  空走率  喷涂段沿程每 5mm 取一点，若该处 Rf 内没有任何工件表面点，
          说明这一段是在空气里喷（穿洞、越界、空行程）。
          报告：空走长度占喷涂总长的比例 + 最长的连续空走段

  点距    相邻喷枪点的间距分布。设定 4mm 却在某处 30mm，说明该段是
          "两点之间靠控制器插补"，喷膜厚会在这里变薄。

  姿态    相邻喷枪轴线的夹角（抖动），以及枪轴与最近表面点局部法向的夹角
          （偏斜角：偏了喷膜就偏一边、厚度不匀）。

用法：
    python audit_traj.py <cloud.ply> <traj.csv> [--row 12] [--so 100] [--rf 0.5]
"""
import sys
import os
import csv
import numpy as np

CELL = 20.0


# ---------- 步骤 1：读点云 ----------
def load_cloud(path):
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


# ---------- 步骤 2：读轨迹（含枪轴） ----------
def load_traj(path):
    with open(path) as f:
        rd = csv.reader(f)
        head = next(rd)
        idx = {k: i for i, k in enumerate(head)}
        rows = []
        for rec in rd:
            rows.append([float(rec[idx["X"]]), float(rec[idx["Y"]]), float(rec[idx["Z"]]),
                         float(rec[idx["I"]]), float(rec[idx["J"]]), float(rec[idx["K"]]),
                         int(float(rec[idx["spray"]]))])
    A = np.array(rows, dtype=np.float64)
    return A[:, :3], A[:, 3:6], A[:, 6].astype(int)


# ---------- 步骤 3：空间网格（统一给点云和落点用） ----------
class Grid:
    def __init__(self, Q, cell=CELL):
        self.Q = Q
        self.cell = float(cell)
        self.org = Q.min(axis=0)
        ijk = np.floor((Q - self.org) / self.cell).astype(np.int64)
        self.dim = ijk.max(axis=0) + 1
        key = (ijk[:, 0] * self.dim[1] + ijk[:, 1]) * self.dim[2] + ijk[:, 2]
        order = np.argsort(key, kind="stable")
        self.key = key[order]
        self.idx = order

    def cand(self, lo, hi):
        c = self.cell
        i0 = np.maximum(np.floor((lo - self.org) / c).astype(np.int64), 0)
        i1 = np.minimum(np.floor((hi - self.org) / c).astype(np.int64), self.dim - 1)
        if np.any(i1 < i0):
            return None
        out = []
        dy, dz = self.dim[1], self.dim[2]
        # 网格跨格时逐格归并，避免一次 huge alloc
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


def nearest_dist(grid, P, rmax, chunk=4000):
    """每个 P 点到 grid.Q 的最近距离（超出 rmax 记 rmax）"""
    Q = grid.Q
    out = np.full(len(P), rmax)
    for s in range(0, len(P), chunk):
        blk = P[s:s + chunk]
        for k in range(len(blk)):
            p = blk[k]
            c = grid.cand(p - rmax, p + rmax)
            if c is None:
                continue
            d = np.linalg.norm(Q[c] - p, axis=1).min()
            out[s + k] = min(d, rmax)
    return out


def nearest_center_dist(grid, P, rmax, chunk=4000):
    """每个 P 点到最近落点中心的距离"""
    return nearest_dist(grid, P, rmax, chunk)


# ---------- 步骤 4：覆盖率（到折线段的距离） ----------
def polyline_dist(Q, C, L, rmax, chunk=3000):
    """每个 Q 点到"落点折线"的最近距离

    不能拿表面点到**离散落点**的距离当覆盖判据：行距 12、点距 4 时，
    即使行铺得完全正确，两个落点对角线上的表面点到最近落点的距离也是
    sqrt(6²+2²)=6.3mm，会被误判成"漏喷"。喷涂覆盖的是**一条条线**，
    所以必须算到折线段的垂直距离。
    """
    A = C[:-1]
    B = C[1:]
    v = B - A
    ok = np.linalg.norm(v, axis=1) < 3.0 * max(L, 1.0)   # 行内相邻才成段
    A, B, v = A[ok], B[ok], v[ok]
    mid = 0.5 * (A + B)
    grid = Grid(mid, cell=max(20.0, 2.0 * L))

    out = np.full(len(Q), rmax)
    for s in range(0, len(Q), chunk):
        blk = Q[s:s + chunk]
        for k in range(len(blk)):
            p = blk[k]
            c = grid.cand(p - rmax, p + rmax)
            if c is None:
                continue
            a, b = A[c], B[c]
            ab = b - a
            L2 = np.maximum((ab * ab).sum(axis=1), 1e-12)
            t = np.clip(((p - a) * ab).sum(axis=1) / L2, 0.0, 1.0)
            pr = a + t[:, None] * ab
            d = np.linalg.norm(pr - p, axis=1).min()
            out[s + k] = min(d, rmax)
    return out


def coverage(Q, C, L, so, rf=0.5, rf2=0.65):
    """Q 表面点，C 落点中心（已沿枪轴前进 standoff）"""
    Rf = rf * L
    rmax = max(Rf * 4.0, so)
    d = polyline_dist(Q, C, L, rmax)
    near = d <= so * 1.5
    res = {
        "cover_all": float((d <= Rf).mean()) * 100.0,
        "cover_065": float((d <= rf2 * L).mean()) * 100.0,
        "cover_near": float((d[near] <= Rf).mean()) * 100.0 if near.any() else 0.0,
        "near_frac": float(near.mean()) * 100.0,
        "d_p50": float(np.percentile(d, 50)),
        "d_p95": float(np.percentile(d, 95)),
        "d_max": float(d.max()),
        "gap_gt_0.9L": float((d > 0.9 * L).mean()) * 100.0,
    }
    return res, d


# ---------- 步骤 5：空走率 ----------
def airrun(P, sp, Q, L, step=5.0, rf=0.5, C=None):
    """喷涂段沿程采样，看落点脚下 Rf 内有没有工件点

    P 是枪位（离面 standoff），直接测 P 到点云当然永远是 standoff 那么远。
    必须用落点中心 C（= P + standoff·枪轴）来测"这段是不是在空气里喷"。
    """
    if C is None:
        C = P
    grid = Grid(Q, cell=max(20.0, 2.0 * L))
    Rf = rf * L
    air = 0.0
    total = 0.0
    longest = 0.0
    cur = 0.0
    for k in range(len(P) - 1):
        if sp[k] != 1 or sp[k + 1] != 1:
            continue
        A, B = C[k], C[k + 1]
        seg = float(np.linalg.norm(B - A))
        if seg < 1e-9:
            continue
        total += seg
        n = max(1, int(seg / step))
        for t in range(n):
            p = A + (B - A) * ((t + 0.5) / n)
            c = grid.cand(p - Rf, p + Rf)
            # 注意：cand 只按网格单元取候选，单元内还有 Rf 之外的远点，
            # 必须再算一次真实距离；否则"孔洞里只要有同格的边缘点"就会被判成非空。
            hit = c is not None and float(np.linalg.norm(Q[c] - p, axis=1).min()) <= Rf
            if not hit:
                air += seg / n
                cur += seg / n
                longest = max(longest, cur)
            else:
                cur = 0.0
    return {"air_pct": (100.0 * air / total) if total > 0 else 0.0,
            "air_len": air, "spray_len": total, "air_longest": longest}


# ---------- 步骤 5.5：实际离面距离 ----------
# ---------- 步骤 5.5：实际喷距（沿枪轴射线打到表面） ----------
def spray_dist(P, A, sp, Q, so, nsample=800, rad=6.0):
    """沿枪轴射线找"枪口正下方最近的表面点"，得到真正的喷涂距离

    不能拿"枪位到最近点云点的距离"当喷距：开口盒内壁布行时，最近的往往是
    旁边的侧壁（31mm），而枪轴正下方那块面其实还是 100mm。
    这两件事必须分开——前者是碰撞/安全，后者才是膜厚。
    """
    grid = Grid(Q, cell=8.0)
    idx = np.where(sp == 1)[0]
    if len(idx) == 0:
        return {}
    sel = idx[np.linspace(0, len(idx) - 1, min(nsample, len(idx))).astype(int)]
    out = []
    tmax = 3.0 * so
    for k in sel:
        p, a = P[k], A[k]
        hi = p + a * tmax
        lo = np.minimum(p, hi) - rad
        hh = np.maximum(p, hi) + rad
        c = grid.cand(lo, hh)
        if c is None:
            continue
        q = Q[c] - p
        t = q @ a
        perp = np.linalg.norm(q - t[:, None] * a, axis=1)
        ok = (t > 1.0) & (t < tmax) & (perp <= rad)
        if not ok.any():
            continue
        out.append(t[ok].min())
    if not out:
        return {}
    d = np.array(out)
    return {
        "sd_p05": float(np.percentile(d, 5)),
        "sd_p50": float(np.percentile(d, 50)),
        "sd_p95": float(np.percentile(d, 95)),
        "sd_min": float(d.min()),
        "sd_max": float(d.max()),
        "sd_out20": float((np.abs(d - so) > 0.2 * so).mean()) * 100.0,
        "sd_out50": float((np.abs(d - so) > 0.5 * so).mean()) * 100.0,
        "sd_hit": float(len(d)) / max(1, len(sel)) * 100.0,
    }


def standoff_stats(P, sp, Q, so, near=12.0):
    """枪位到工件表面的实际最近距离分布

    工艺假定"枪位 = 面点 + standoff·法向"，即恒为 standoff。
    实测距离一旦偏离，膜厚就跟着变（喷距偏差 ±20% 通常已超出工艺窗口）。
    """
    grid = Grid(Q, cell=8.0)
    G = P[sp == 1]
    rmax = 3.0 * so
    d = nearest_dist(grid, G, rmax)
    dev = d - so
    return {
        "so_p05": float(np.percentile(d, 5)),
        "so_p50": float(np.percentile(d, 50)),
        "so_p95": float(np.percentile(d, 95)),
        "so_min": float(d.min()),
        "so_max": float(d.max()),
        "so_out20": float((np.abs(dev) > 0.2 * so).mean()) * 100.0,
        "so_out50": float((np.abs(dev) > 0.5 * so).mean()) * 100.0,
    }


# ---------- 步骤 6：点距 ----------
def steps(P, sp):
    d = np.linalg.norm(np.diff(P, axis=0), axis=1)
    d = d[d > 1e-9]
    if len(d) == 0:
        return {}
    return {"step_p50": float(np.percentile(d, 50)),
            "step_p99": float(np.percentile(d, 99)),
            "step_max": float(d.max()),
            "step_zero": int((d < 0.05).sum())}


# ---------- 步骤 7：姿态 ----------
def pca_normals(Qsub, Q, tree, k=16, rad=15.0):
    """局部 PCA 法向 + 平面度

    平面度 = λmin / (λ0+λ1+λ2)。薄壁/棱边/两片点混在一球邻域里时，
    平面度会很大（形状根本不是一张面），这类点的法向没有参考价值，
    必须剔除，否则"枪轴偏斜 63°"全是估计噪声，不是算法问题。
    """
    N = np.tile(np.array([0.0, 0.0, 1.0]), (len(Qsub), 1))
    flat = np.ones(len(Qsub))
    for i, p in enumerate(Qsub):
        c = tree.cand(p - rad, p + rad)
        if c is None or len(c) < 6:
            continue
        nb = Q[c]
        d = np.linalg.norm(nb - p, axis=1)
        nb = nb[np.argsort(d)[:k]]
        X = nb - nb.mean(axis=0)
        w, v = np.linalg.eigh(X.T @ X)
        w = np.clip(w, 1e-12, None)
        N[i] = v[:, 0]
        flat[i] = w[0] / w.sum()
    return N, flat


def pose_stats(P, A, sp, Q, tree, L, so=100.0):
    """相邻轴线夹角 + 轴线与表面法向夹角"""
    idx = np.where(sp == 1)[0]
    if len(idx) < 3:
        return {}
    Aa = A[idx]
    # 相邻轴线夹角（只在同一行内看，跨行会跳到另一条行，没有可比性）
    d = np.linalg.norm(np.diff(Aa, axis=0), axis=1)
    cosang = np.clip((Aa[1:] * Aa[:-1]).sum(axis=1), -1, 1)
    ang = np.degrees(np.arccos(cosang))
    # 跨行不算：位置跳 > 3×点距 的相邻对子剔除
    pos = P[idx]
    jump = np.linalg.norm(np.diff(pos, axis=0), axis=1)
    same = jump < 3.0 * max(L, 1.0)
    # 枪轴 vs 局部法向
    Ps = pos + so * Aa          # 落点（P = S - so·A）
    sel = np.arange(0, len(Ps), max(1, len(Ps) // 3000))
    Ns, flat = pca_normals(Ps[sel], Q, tree)
    good = flat < 0.05                       # 只统计"确实在一张平面上"的落点
    dotn = np.clip((Aa[sel] * Ns).sum(axis=1), -1, 1)
    skew = np.degrees(np.arccos(np.abs(dotn)))
    sk = skew[good] if good.any() else skew
    return {"axis_jump_p95": float(np.percentile(ang[same], 95)) if same.any() else 0.0,
            "axis_jump_max": float(ang[same].max()) if same.any() else 0.0,
            "skew_p50": float(np.percentile(sk, 50)),
            "skew_p95": float(np.percentile(sk, 95)),
            "skew_max": float(sk.max()),
            "skew_eval_frac": float(good.mean()) * 100.0}


# ---------- 步骤 7.5：可喷点子集 ----------
def sprayable_subset(Q, cap=12000, win=130.0):
    """按算法自己的角度窗口筛"本来就该喷"的点

    覆盖率如果拿"全部点云"当分母，朝下的底面（angZ>130°，算法明确不喷）
    会白占一大块，看上去像漏喷，其实是按规则排除的。
    这里用局部 PCA 法向估计每个点朝向，只保留 angZ<=win 且局部足够平的点。
    """
    if len(Q) > cap:
        rs = np.random.default_rng(11)
        Qs = Q[rs.choice(len(Q), cap, replace=False)]
    else:
        Qs = Q
    tree = Grid(Q, cell=8.0)
    N, flat = pca_normals(Qs, Q, tree)
    angZ = np.degrees(np.arccos(np.clip(np.abs(N[:, 2]), -1, 1)))
    # |nz| 取绝对值：法向无符号，朝下的面 angZ 会算成 180-，这里直接按"与Z轴夹角"判
    az = np.where(N[:, 2] >= 0, angZ, 180.0 - angZ)
    keep = (az <= win) & (flat < 0.08)
    return Qs[keep], float(keep.mean()) * 100.0


# ---------- 步骤 8：主流程 ----------
def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    ply, csvp = sys.argv[1], sys.argv[2]
    row = 12.0
    so = 100.0
    rf = 0.5
    if "--row" in sys.argv:
        row = float(sys.argv[sys.argv.index("--row") + 1])
    if "--so" in sys.argv:
        so = float(sys.argv[sys.argv.index("--so") + 1])
    if "--rf" in sys.argv:
        rf = float(sys.argv[sys.argv.index("--rf") + 1])

    Q = load_cloud(ply)
    # 覆盖率对点云做确定性抽样：20 万点逐点算最近距离太慢，统计意义不变
    if len(Q) > 30000:
        rs = np.random.default_rng(7)
        Qcov = Q[rs.choice(len(Q), 30000, replace=False)]
    else:
        Qcov = Q
    P, A, sp = load_traj(csvp)
    # 落点中心 = 枪位沿枪轴下降 standoff（i,j,k 指向表面）
    C = P + so * A
    Cspray = C[sp == 1]

    print("点云 %d 点 | 轨迹 %d 点（喷 %d / 空移 %d）| 行距 %.1f 离面 %.1f"
          % (len(Q), len(P), (sp == 1).sum(), (sp == 0).sum(), row, so))

    Qsp, sp_frac = sprayable_subset(Q)
    cov, d = coverage(Qsp, Cspray, row, so, rf)
    print("\n[覆盖率]（到轨迹折线的垂距；分母只算可喷点：局部法向 angZ<=130° 且平面度<0.08，占全点云 %.0f%%）" % sp_frac)
    print("  可喷区覆盖   %6.2f%%  (Rf=%.1f=0.5×行距) / %6.2f%% (Rf=%.1f=0.65×行距 即 30%% 搭接)"
          % (cov["cover_all"], rf * row, cov["cover_065"], 0.65 * row))
    print("  枪口所及区域 %6.2f%%  (占可喷区 %.1f%%)"
          % (cov["cover_near"], cov["near_frac"]))
    print("  表面点到最近落点中心  p50=%.2f  p95=%.2f  max=%.2f mm"
          % (cov["d_p50"], cov["d_p95"], cov["d_max"]))
    print("  空档 > 0.9×行距 的占比 %.2f%%" % cov["gap_gt_0.9L"])

    ar = airrun(P, sp, Q, row, rf=rf, C=C)
    print("\n[空走率]")
    print("  喷涂总长 %.0fmm  空中喷涂 %.0fmm  (%.1f%%)  最长连续 %.0fmm"
          % (ar["spray_len"], ar["air_len"], ar["air_pct"], ar["air_longest"]))

    sd = standoff_stats(P, sp, Q, so)
    sd2 = spray_dist(P, A, sp, Q, so)
    print("\n[实际喷距]（沿枪轴射线打到表面；工艺设定 %.0fmm）" % so)
    if sd2:
        print("  p05=%.1f  p50=%.1f  p95=%.1f  min=%.1f  max=%.1f mm"
              % (sd2["sd_p05"], sd2["sd_p50"], sd2["sd_p95"], sd2["sd_min"], sd2["sd_max"]))
        print("  偏离设定值 >20%% 的点 %.2f%%   >50%% 的点 %.2f%%   (射线命中 %.0f%%)"
              % (sd2["sd_out20"], sd2["sd_out50"], sd2["sd_hit"]))
    print("  [参考] 枪位到最近点云点（含侧壁，安全口径）: p05=%.1f p50=%.1f p95=%.1f min=%.1f mm"
          % (sd["so_p05"], sd["so_p50"], sd["so_p95"], sd["so_min"]))

    st = steps(P, sp)
    print("\n[点距] 设定应为较小值")
    print("  p50=%.2f  p99=%.2f  max=%.2f mm  零长段 %d 个"
          % (st.get("step_p50", 0), st.get("step_p99", 0),
             st.get("step_max", 0), st.get("step_zero", 0)))

    tree = Grid(Q, cell=8.0)
    po = pose_stats(P, A, sp, Q, tree, row, so)
    print("\n[姿态]")
    print("  同行相邻轴线夹角  p95=%.2f°  max=%.2f°"
          % (po.get("axis_jump_p95", 0), po.get("axis_jump_max", 0)))
    print("  枪轴 vs 局部法向偏斜  p50=%.1f°  p95=%.1f°  max=%.1f°  (可评估落点 %.0f%%)"
          % (po.get("skew_p50", 0), po.get("skew_p95", 0), po.get("skew_max", 0),
             po.get("skew_eval_frac", 0)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
