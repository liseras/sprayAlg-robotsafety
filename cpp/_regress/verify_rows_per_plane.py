#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""逐面验收：把轨迹点按"到哪个面正好等于离面距离"归属到各个重建面，
再在该面的局部 2D 坐标里量化 行数 / 行距均匀性 / 是否溢出重建轮廓。

用法：
    python verify_rows_per_plane.py <out_csv_prefix> [standoff] [row]

例：
    python verify_rows_per_plane.py /tmp/q 100 12

输入：
    <prefix>.csv             轨迹点（X,Y,Z,I,J,K,speed,spray）
    <prefix>.ply.planeK.ply  每个重建面的三角网格（顶点+面）

输出：
    每个面：顶点数/三角数、面内边长、行数、行距（均值/最小/最大/标准差）、
    轨迹足迹相对重建轮廓的四个方向余量（负值 = 溢出）。
"""
import sys
import glob
import numpy as np


def read_ply_ascii(path):
    with open(path, "r", errors="replace") as f:
        lines = f.read().split("\n")
    nv = nf = hdr = 0
    for i, ln in enumerate(lines):
        if ln.startswith("element vertex"):
            nv = int(ln.split()[2])
        elif ln.startswith("element face"):
            nf = int(ln.split()[2])
        elif ln.strip() == "end_header":
            hdr = i + 1
            break
    V = np.array([[float(x) for x in lines[hdr + k].split()[:3]] for k in range(nv)])
    F = []
    for k in range(nf):
        p = lines[hdr + nv + k].split()
        idx = [int(x) for x in p[1:1 + int(p[0])]]
        for t in range(1, len(idx) - 1):
            F.append((idx[0], idx[t], idx[t + 1]))
    return V, F


def boundary_edges(V, F):
    """只被一个三角形使用的边 = 轮廓边。用来确定面的真实边长方向。"""
    cnt = {}
    for tri in F:
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            k = (a, b) if a < b else (b, a)
            cnt[k] = cnt.get(k, 0) + 1
    edges = [(a, b) for (a, b), c in cnt.items() if c == 1]
    return edges


def read_csv(path):
    rows = []
    with open(path, "r", errors="replace") as f:
        f.readline()
        for ln in f:
            p = ln.strip().split(",")
            if len(p) >= 8:
                rows.append([float(v) for v in p[:3]] + [int(float(p[7]))])
    a = np.array(rows)
    return a[:, :3], a[:, 3]


def plane_frame(V, F):
    """由面片求 (质心, 法向, u, v)。
    u 取【轮廓边】里最长的一条 —— 不能用任意顶点对的最长距离，那会取到矩形对角线（372mm），
    坐标系会斜 45°，面内尺寸被算大、行分簇也会失效。
    """
    c = V.mean(0)
    edges = boundary_edges(V, F)
    if edges:
        best, blen = None, -1.0
        for a, b in edges:
            L = np.linalg.norm(V[a] - V[b])
            if L > blen:
                blen, best = L, (a, b)
        u = V[best[1]] - V[best[0]]
    else:
        u = V[1] - V[0]
    u = u / (np.linalg.norm(u) + 1e-12)

    _, _, vt = np.linalg.svd(V - c, full_matrices=False)
    n = vt[2]
    if np.dot(n, c) < 0:   # 只统一符号，便于打印
        n = -n
    # 把 u 投影到面内并正交化
    u = u - np.dot(u, n) * n
    u = u / (np.linalg.norm(u) + 1e-12)
    v = np.cross(n, u)
    return c, n, u, v


def main():
    prefix = sys.argv[1]
    standoff = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
    row_expected = float(sys.argv[3]) if len(sys.argv) > 3 else 12.0

    P, spray = read_csv(prefix + ".csv")
    P = P[spray > 0]
    meshes = sorted(glob.glob(prefix + ".ply.plane*.ply"))
    if not meshes:
        print("未找到 %s.ply.plane*.ply，请先跑 spray_driver 并确认导出了面片" % prefix)
        return
    if len(P) == 0:
        print("CSV 中没有 spray=1 的点")
        return

    frames = []
    for m in meshes:
        V, F = read_ply_ascii(m)
        c, n, u, v = plane_frame(V, F)
        U = (V - c) @ u
        Vv = (V - c) @ v
        frames.append({
            "file": m, "V": V, "F": F, "c": c, "n": n, "u": u, "v": v,
            "u0": U.min(), "u1": U.max(), "v0": Vv.min(), "v1": Vv.max(),
        })

    # 归属：轨迹点沿法向到该面的距离应当正好等于 standoff。
    # 注意两点：
    #  ① 面法向符号是随机的（PCA/SVD 决定的），必须用 min(|d-s|,|d+s|) 做符号无关判定，
    #     否则法向朝反的面会一点都分不到；
    #  ② 只比距离不够——多个面的距离可能恰好都接近 standoff，故加入"投影是否落在本面矩形内"
    #     的惩罚项，让"离得近且落在面上"的面胜出。
    npl = len(frames)
    score = np.zeros((len(P), npl))
    for k, f in enumerate(frames):
        d = (P - f["c"]) @ f["n"]
        resid = np.minimum(np.abs(d - standoff), np.abs(d + standoff))
        Uq = (P - f["c"]) @ f["u"]
        Vq = (P - f["c"]) @ f["v"]
        outside = np.maximum.reduce([
            Uq - f["u1"], f["u0"] - Uq, Vq - f["v1"], f["v0"] - Vq,
            np.zeros_like(Uq)])
        score[:, k] = resid + 10.0 * outside
    who = np.argmin(score, axis=1)
    best = score[np.arange(len(P)), who]

    print("轨迹点 %d，归属得分（残差+10×出界距离）: 中位数=%.3f 最大=%.3f mm"
          % (len(P), np.median(best), best.max()))

    for k, f in enumerate(frames):
        sel = who == k
        n = int(sel.sum())
        print("\n=== 面%d  %s ===" % (k, f["file"].split("/")[-1].split("\\")[-1]))
        print("顶点=%d  面内尺寸 u=%.1f x v=%.1f mm  法向=(%.3f,%.3f,%.3f)"
              % (len(f["V"]), f["u1"] - f["u0"], f["v1"] - f["v0"], *f["n"]))
        if n == 0:
            print("该面没有轨迹点")
            continue

        Q = P[sel] - f["c"]
        U = Q @ f["u"]
        Vv = Q @ f["v"]

        # ── 行方向：取自轨迹自身的直线段，而不是假定它等于矩形边 ────────────
        # 这样"行是否平行矩形边"才是一个能被检出的量（夹角 > 1° 就说明没对齐）。
        Qall = P - f["c"]
        Uall, Vall = Qall @ f["u"], Qall @ f["v"]
        idx = np.where(sel)[0]

        def step_angle(a, b):
            du = Uall[b] - Uall[a]
            dv = Vall[b] - Vall[a]
            L = np.hypot(du, dv)
            if not (2.0 < L < 12.0):     # 只认行内正常步长（≈点距 4mm）
                return None
            return np.arctan2(dv, du)

        angs = [a for a in (step_angle(p, q) for p, q in zip(idx[:-1], idx[1:])
                            if q == p + 1 and who[q] == k) if a is not None]
        if not angs:
            print("该面没有可用的行内步长")
            continue
        ang = np.array(angs)
        theta = np.median((ang + np.pi) % np.pi)   # 180° 周期下的中位
        dev_deg = (np.degrees(theta) + 90) % 180 - 90
        print("行方向相对矩形长边偏离 = %.2f°" % dev_deg)

        # 旋转到行方向
        ct, st = np.cos(-theta), np.sin(-theta)
        A = Uall * ct - Vall * st      # 沿行
        B = Uall * st + Vall * ct      # 跨行

        # ── 分行：蛇形轨迹每一行是一段"方向恒定"的连续段，按方向反转切分 ──────
        # 比按跨行坐标聚簇稳：不会因为行内轻微漂移把一个行拆成两段，也不会把相邻行并起来。
        rows = []
        cur = [idx[0]]
        for a, b in zip(idx[:-1], idx[1:]):
            th = step_angle(a, b) if (b == a + 1 and who[b] == k) else None
            if th is None:
                # 折返弧/空移：行结束
                if len(cur) >= 3:
                    rows.append(cur)
                cur = [b] if who[b] == k else []
                continue
            aligned = abs(((th - theta + np.pi / 2) % np.pi) - np.pi / 2) < np.radians(25)
            if not aligned and len(cur) >= 3:
                rows.append(cur)
                cur = []
            cur.append(b)
        if len(cur) >= 3:
            rows.append(cur)

        # 只保留够长的（真正的喷涂行，而非短碎片）
        rows = [r for r in rows if A[r].max() - A[r].min() >= 30.0]
        if len(rows) < 2:
            print("有效行数=%d（不足 2 行，无法评估行距）" % len(rows))
            continue

        ys = np.array([np.median(B[r]) for r in rows])
        dy = np.abs(np.diff(ys))   # 行按 B 递减排列，取绝对值避免出现负行距
        print("行数=%d  行距: 均值=%.2f 最小=%.2f 最大=%.2f 标准差=%.2f (设定 %.1f)"
              % (len(rows), dy.mean(), dy.min(), dy.max(), dy.std(), row_expected))

        lens = np.array([A[r].max() - A[r].min() for r in rows])
        print("行长度: 均值=%.1f 最短=%.1f 最长=%.1f mm" % (lens.mean(), lens.min(), lens.max()))

        # 溢出检查：把矩形四角也转到行坐标系，比较足迹
        cu = np.array([f["u0"], f["u1"], f["u1"], f["u0"]])
        cv = np.array([f["v0"], f["v0"], f["v1"], f["v1"]])
        Aq = cu * ct - cv * st
        Bq = cu * st + cv * ct
        ra = np.concatenate([A[r] for r in rows])
        rb = np.concatenate([B[r] for r in rows])
        print("足迹 A=[%.1f, %.1f] B=[%.1f, %.1f]；轮廓 A=[%.1f, %.1f] B=[%.1f, %.1f]"
              % (ra.min(), ra.max(), rb.min(), rb.max(), Aq.min(), Aq.max(), Bq.min(), Bq.max()))
        print("边缘余量: 左=%.1f 右=%.1f 下=%.1f 上=%.1f mm  (负值=溢出)"
              % (ra.min() - Aq.min(), Aq.max() - ra.max(),
                 rb.min() - Bq.min(), Bq.max() - rb.max()))
        print("行距一致性: 最大/最小 = %.2f" % (dy.max() / max(1e-6, dy.min())))


if __name__ == "__main__":
    main()
