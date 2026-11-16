# -*- coding: utf-8 -*-
"""量化"轨迹出框"：喷涂点落在所有重建面多边形之外的比例与外溢距离。

用法：
    python tools/traj_outside_mesh.py <traj.csv> <out.ply> [--standoff 100]

traj.csv 列：x,y,z,i,j,k,speed,spray（drv 输出口径，spray 在最后一列）
out.ply   : spray_driver 的轨迹 ply（含 .planeN.ply 面片文件，读同级目录）

输出：
    每个面：面片尺寸 / 落在该面平面带内的喷点数 / 其中真正在多边形内的点数
    全局  ：认定给某面但落在多边形外的喷点数、占比、最大外溢距离（到多边形边界的距离）
"""
import glob
import os
import sys

import numpy as np


def read_ply_verts(path):
    nv = 0
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('element vertex'):
                nv = int(line.split()[-1])
            if line.strip() == 'end_header':
                break
        return np.loadtxt(fh, max_rows=nv)[:, :3]


def read_ply_mesh(path):
    nv = nf = 0
    for line in open(path, 'r', encoding='utf-8', errors='ignore'):
        if line.startswith('element vertex'):
            nv = int(line.split()[-1])
        elif line.startswith('element face'):
            nf = int(line.split()[-1])
        elif line.strip() == 'end_header':
            break
    V = np.loadtxt(path, skiprows=0, max_rows=0)
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        lines = fh.readlines()
    hdr_end = next(i for i, l in enumerate(lines) if l.strip() == 'end_header')
    V = np.array([[float(x) for x in l.split()] for l in lines[hdr_end + 1:hdr_end + 1 + nv]])[:, :3]
    T = []
    if nf:
        for l in lines[hdr_end + 1 + nv:hdr_end + 1 + nv + nf]:
            p = l.split()
            if len(p) >= 4:
                T.append([int(p[1]), int(p[2]), int(p[3])])
    return V, np.array(T, dtype=int) if T else np.zeros((0, 3), dtype=int)


def boundary_loops(V, T):
    """三角网格的边界环（外环 + 孔洞环），2D 坐标在面内 u/v 平面上给出。"""
    if len(T) == 0:
        return []
    cnt = {}
    for a, b, c in T:
        for e in ((a, b), (b, c), (c, a)):
            k = (min(e), max(e))
            cnt[k] = cnt.get(k, 0) + 1
    adj = {}
    for (a, b), n in cnt.items():
        if n != 1:
            continue
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    loops, used = [], set()
    for s in list(adj):
        if s in used:
            continue
        loop, cur, prev = [], s, None
        while cur is not None and cur not in used:
            used.add(cur)
            loop.append(cur)
            nxt = [q for q in adj[cur] if q != prev and q not in used]
            prev, cur = cur, (nxt[0] if nxt else None)
        if len(loop) >= 3:
            loops.append(loop)
    return [V[lp] for lp in loops]


def inside_2d(P, V, T):
    loops = boundary_loops(V, T)
    if not loops:
        return np.zeros(len(P), dtype=bool)
    inside = np.zeros(len(P), dtype=bool)
    for lp in loops:
        n = len(lp)
        for i in range(n):
            ax, ay = lp[i]
            bx, by = lp[(i + 1) % n]
            dy = by - ay
            dy = dy if abs(dy) > 1e-12 else 1e-12
            cond = ((ay > P[:, 1]) != (by > P[:, 1])) & (P[:, 0] < (bx - ax) * (P[:, 1] - ay) / dy + ax)
            inside ^= cond
    return inside


def dist_to_loops_2d(P, loops):
    d = np.full(len(P), 1e18)
    for lp in loops:
        n = len(lp)
        for i in range(n):
            a, b = lp[i], lp[(i + 1) % n]
            ab = b - a
            L2 = float(ab @ ab)
            if L2 < 1e-12:
                t = np.zeros(len(P))
            else:
                t = np.clip(((P - a) @ ab) / L2, 0.0, 1.0)
            q = a + ab * t[:, None]
            d = np.minimum(d, np.linalg.norm(P - q, axis=1))
    return d


def main():
    csv, oply = sys.argv[1], sys.argv[2]
    standoff = 100.0
    if '--standoff' in sys.argv:
        standoff = float(sys.argv[sys.argv.index('--standoff') + 1])

    t = np.loadtxt(csv, delimiter=',', skiprows=1)
    P = t[:, :3]
    sp = t[:, 7].astype(int)
    faces = sorted(glob.glob(oply + '.plane*.ply'))
    if not faces:
        print('没有面片文件：' + oply + '.plane*.ply')
        return

    assigned = np.full(len(P), -1)
    bestd = np.full(len(P), 1e18)
    info = []
    for k, f in enumerate(faces):
        V, T3 = read_ply_mesh(f)
        if len(T3) == 0:
            continue
        n = np.cross(V[T3[0][1]] - V[T3[0][0]], V[T3[0][2]] - V[T3[0][0]])
        n = n / max(np.linalg.norm(n), 1e-12)
        c = V.mean(0)
        u = np.array([1.0, 0.0, 0.0])
        if abs(float(n @ u)) > 0.9:
            u = np.array([0.0, 1.0, 0.0])
        u = u - (n @ u) * n
        u = u / np.linalg.norm(u)
        v = np.cross(n, u)
        d = (P - c) @ n
        sel = np.abs(d - standoff) < 8.0
        idx = np.flatnonzero(sel)
        Q = np.stack([(P[idx] - c) @ u, (P[idx] - c) @ v], axis=1)
        loops = boundary_loops(V @ np.stack([u, v], 1), T3) if False else None
        # 面片顶点本身已经世界坐标，需要投到 (u,v)
        V2 = np.stack([(V - c) @ u, (V - c) @ v], axis=1)
        lp = boundary_loops(V2, T3)
        ins = inside_2d(Q, V2, T3)
        size = np.round(V.max(0) - V.min(0), 1)
        info.append((os.path.basename(f), size, len(idx), int(ins.sum())))
        dd = np.abs(d - standoff)
        better = sel & (dd < bestd)
        assigned[better] = k
        bestd[better] = dd[better]

    print('轨迹点 %d（喷涂 %d）  面片 %d 个' % (len(P), int((sp == 1).sum()), len(faces)))
    print('%-28s %-22s %8s %8s' % ('面片', '尺寸(mm)', '平面带点', '多边形内'))
    for nm, size, tot, ins in info:
        print('%-28s %-22s %8d %8d' % (nm, 'x'.join('%.0f' % s for s in size), tot, ins))

    # 全局：喷点归属到最近面，判是否落在多边形内
    out_n = 0
    out_max = 0.0
    for k, f in enumerate(faces):
        V, T3 = read_ply_mesh(f)
        if len(T3) == 0:
            continue
        n = np.cross(V[T3[0][1]] - V[T3[0][0]], V[T3[0][2]] - V[T3[0][0]])
        n = n / max(np.linalg.norm(n), 1e-12)
        c = V.mean(0)
        u = np.array([1.0, 0.0, 0.0])
        if abs(float(n @ u)) > 0.9:
            u = np.array([0.0, 1.0, 0.0])
        u = u - (n @ u) * n
        u = u / np.linalg.norm(u)
        v = np.cross(n, u)
        idx = np.flatnonzero((assigned == k) & (sp == 1))
        if len(idx) == 0:
            continue
        Q = np.stack([(P[idx] - c) @ u, (P[idx] - c) @ v], axis=1)
        V2 = np.stack([(V - c) @ u, (V - c) @ v], axis=1)
        lp = boundary_loops(V2, T3)
        ins = inside_2d(Q, V2, T3)
        bad = ~ins
        out_n += int(bad.sum())
        if bad.any():
            dm = dist_to_loops_2d(Q[bad], lp)
            out_max = max(out_max, float(dm.max()))
            print('  %-24s 多边形外喷点 %4d / %4d  最远外溢 %.1f mm'
                  % (os.path.basename(f), int(bad.sum()), len(idx), float(dm.max())))

    nsp = int((sp == 1).sum())
    print('汇总：喷涂点 %d，落在所属面多边形外 %d 个（%.1f%%），最远外溢 %.1f mm'
          % (nsp, out_n, 100.0 * out_n / max(nsp, 1), out_max))


if __name__ == '__main__':
    main()
