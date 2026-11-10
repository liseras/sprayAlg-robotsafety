# -*- coding: utf-8 -*-
"""整件完整轨迹计划：喷涂段（P 型逐点）+ 空移段（J 型关节点），一条序列可整条模拟。

输出格式（给 DensoArmKine 的「J 型空移模拟」面板读）：
{
  "Name": "...",
  "Summary": {...},
  "Segments": [
    {"Type":"spray", "Points":[{"p":[x,y,z],"th":[j1..j6]}, ...]},
    {"Type":"jmove", "Points":[{"p":[...],"th":[...],"g":[...],"tilt":0}, ...]}
  ]
}
喷涂段：每个规划点一个关节角（连续 seed），播放时点间关节线性插值。
空移段：关节点（与 move_joint 同源：J5 超限即弃、倾斜回退）。
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys

import numpy as np

from . import robot, geom
from .move_plan import CloudGrid, read_ply_xyz
from .move_joint import _dir, _pos, _segments, solve_joint, solve_any_tilt, _tilt
from .move_pose import build_dirs

J5_LIM = 120.0


def solve_spray(p, g, tool, seed):
    """喷涂点：用多初值 IK（限位优先），并返回是否超限。"""
    R = geom.pose_from_z(g)[:3, :3]
    sol = robot.ik(p, R, tool, seed=seed)
    if not sol.ok:
        return None, True
    th = sol.theta.copy()
    if seed is not None:
        th = th + 360.0 * np.round((seed - th) / 360.0)
    return th, bool(sol.over)


def run(csv_in, cloud_in, out_json, step_mm=50.0):
    rows = list(csv.DictReader(open(csv_in, encoding='utf-8-sig')))
    grid = CloudGrid(read_ply_xyz(cloud_in))
    tool = robot.Tool(length=453.08)
    spray = [int(r['spray']) for r in rows]
    segs = _segments(spray)

    out = {'Name': csv_in.rsplit('\\', 1)[-1].rsplit('.', 1)[0] + '_full',
           'Summary': {}, 'Segments': []}
    ik_fail = 0
    over_limit = 0
    seed = None

    # 按原顺序：喷涂段与空移段交替
    cursor = 0
    for (i0, i1) in segs:
        # 前面的喷涂点 [cursor, i0-1]
        if i0 > cursor:
            pts = []
            for k in range(cursor, i0):
                p = _pos(rows[k]); g = _dir(rows[k])
                th, over = solve_spray(p, g, tool, seed)
                if th is None:
                    ik_fail += 1
                    th = seed if seed is not None else np.zeros(6)
                elif over:
                    over_limit += 1
                seed = th
                pts.append({'p': p.tolist(), 'th': th.tolist()})
            out['Segments'].append({'Type': 'spray', 'Points': pts})
        # 空移段
        g0 = _dir(rows[i0 - 1]) if i0 > 0 else _dir(rows[i1 + 1])
        g1 = _dir(rows[i1 + 1]) if i1 + 1 < len(rows) else g0
        P = np.array([_pos(rows[k]) for k in range(i0, i1 + 1)])
        L = float(np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1)))
        pole = np.array([0.0, 0.0, -1.0])
        dirs = build_dirs(g0, g1, pole, len(P), 5.0)
        s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
        n_wp = max(2, int(math.ceil(L / step_mm)) + 1)
        centroid = grid.centroid
        wp = []
        for k in range(n_wp):
            t = L * k / (n_wp - 1)
            idx = min(int(np.searchsorted(s, t)), len(P) - 1)
            pnt, g = P[idx], dirs[idx]
            hd = pnt - centroid; hd[2] = 0.0
            nh = np.linalg.norm(hd)
            hdir = hd / nh if nh > 1e-6 else np.array([1.0, 0.0, 0.0])
            seam = (k == 0 or k == n_wp - 1)
            th = None; used = 0.0
            if seam:
                th = solve_joint(pnt, g, tool, seed)
                if th is None:
                    th, used, g2 = solve_any_tilt(pnt, pole, g, tool, seed, hdir, 15.0)
                    if th is not None:
                        g = g2
            else:
                th, used, g2 = solve_any_tilt(pnt, pole, g, tool, seed, hdir, 45.0)
                if th is not None:
                    g = g2
            if th is None:
                ik_fail += 1
                th = seed if seed is not None else np.zeros(6)
            else:
                seed = th
            wp.append({'p': pnt.tolist(), 'th': th.tolist(), 'g': g.tolist(), 'tilt': used})
        out['Segments'].append({'Type': 'jmove', 'Points': wp})
        cursor = i1 + 1

    if cursor < len(rows):
        pts = []
        for k in range(cursor, len(rows)):
            p = _pos(rows[k]); g = _dir(rows[k])
            th, over = solve_spray(p, g, tool, seed)
            if th is None:
                ik_fail += 1
                th = seed if seed is not None else np.zeros(6)
            elif over:
                over_limit += 1
            seed = th
            pts.append({'p': p.tolist(), 'th': th.tolist()})
        out['Segments'].append({'Type': 'spray', 'Points': pts})

    out['Summary'] = {
        'points': len(rows),
        'segments': len(out['Segments']),
        'ik_fail': ik_fail,
        'over_limit': over_limit,
        'step_mm': step_mm,
    }
    open(out_json, 'w', encoding='utf-8').write(json.dumps(out, ensure_ascii=False))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('csv'); ap.add_argument('cloud'); ap.add_argument('out')
    ap.add_argument('--step-mm', type=float, default=50.0)
    a = ap.parse_args(argv)
    o = run(a.csv, a.cloud, a.out, a.step_mm)
    print(json.dumps(o['Summary'], ensure_ascii=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
