# -*- coding: utf-8 -*-
"""按面拆分的覆盖诊断

faces.csv（driver 导出）给出每个面的质心/法向/内点数。
先把点云按"到各面平面的垂距最近"归属到面，再对每个面统计：
    点云跨度(mm)  该面在自身局部基下的尺寸
    覆盖率        该面点云中，落点中心 Rf 内有落点的比例
    缺口          未覆盖点到最近落点的距离 p95/max
    V 向行距      x 向行距  （按落点聚类得到，看是否等于设定值）

用法：
    python face_cover.py <cloud.ply> <traj.csv> <faces.csv> [--row 12] [--so 100] [--rf 0.5]
"""
import sys
import csv
import numpy as np
from audit_traj import load_cloud, load_traj, Grid


def load_faces(path):
    out = []
    with open(path, newline="") as f:
        rd = csv.DictReader(f)
        for r in rd:
            out.append({
                "idx": int(r["face"]),
                "count": int(r["count"]),
                "n": np.array([float(r["nx"]), float(r["ny"]), float(r["nz"])]),
                "n0": np.array([float(r["n0x"]), float(r["n0y"]), float(r["n0z"])]),
                "angZ": float(r["angZ"]),
                "isInner": int(r["isInner(1内/0外/-1未分)"]),
                "c": np.array([float(r["cx"]), float(r["cy"]), float(r["cz"])]),
                "inWin": int(r["inWin(1=夹角<=130)"]),
            })
    return out


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    ply, csvp, fcsv = sys.argv[1], sys.argv[2], sys.argv[3]
    row, so, rf = 12.0, 100.0, 0.5
    for k in ("--row", "--so", "--rf"):
        if k in sys.argv:
            val = float(sys.argv[sys.argv.index(k) + 1])
            if k == "--row":
                row = val
            elif k == "--so":
                so = val
            else:
                rf = val

    Q = load_cloud(ply)
    P, A, sp = load_traj(csvp)
    C = (P + so * A)[sp == 1]
    faces = load_faces(fcsv)

    ### 归属：到各面平面的垂距最近
    D = np.stack([np.abs((Q - f["c"]) @ (f["n"] / np.linalg.norm(f["n"]))) for f in faces], axis=1)
    lab = D.argmin(axis=1)

    gridC = Grid(C, cell=max(20.0, 2.0 * row))
    Rf = rf * row
    rmax = max(4.0 * Rf, so)

    print("%-4s %8s %10s %10s %8s %7s %7s %6s  %s"
          % ("面", "点数", "跨度mm", "面积mm2", "覆盖%", "缺口p95", "缺口max", "内壁", "行数(V向)"))
    for fi, f in enumerate(faces):
        sel = lab == fi
        n = int(sel.sum())
        if n < 10:
            print("%-4d %8d  (点太少，跳过)" % (fi, n))
            continue
        q = Q[sel]
        c = f["c"]
        nn = f["n"] / np.linalg.norm(f["n"])
        # 局部基
        tmp = np.array([0, 0, 1.0]) if abs(nn[2]) < 0.9 else np.array([1.0, 0, 0])
        u = np.cross(nn, tmp)
        u /= np.linalg.norm(u)
        v = np.cross(nn, u)
        qu = (q - c) @ u
        qv = (q - c) @ v
        area = (qu.max() - qu.min()) * (qv.max() - qv.min())

        d = np.full(n, rmax)
        for k in range(n):
            cc = gridC.cand(q[k] - rmax, q[k] + rmax)
            if cc is not None:
                d[k] = min(np.linalg.norm(C[cc] - q[k], axis=1).min(), rmax)
        cov = float((d <= Rf).mean()) * 100.0
        # 行数：落点在该面的 v 向聚类
        cs = C[np.linalg.norm(C - c, axis=1) < 1e6]
        cv = (cs - c) @ v
        cu = (cs - c) @ u
        # 只取落在该面范围内的落点
        inb = ((cu > qu.min() - row) & (cu < qu.max() + row) &
               (cv > qv.min() - row) & (cv < qv.max() + row))
        cvv = np.sort(cv[inb])
        nrow = 0
        if len(cvv) > 2:
            brk = np.where(np.diff(cvv) > max(2.0, row * 0.4))[0]
            nrow = len(brk) + 1
        print("%-4d %8d %10s %10.0f %8.1f %7.1f %7.1f %6d  %d"
              % (fi, n, "%.0fx%.0f" % (qu.max() - qu.min(), qv.max() - qv.min()),
                 area, cov, np.percentile(d, 95), d.max(), f["isInner"], nrow))
    return 0


if __name__ == "__main__":
    sys.exit(main())
