# -*- coding: utf-8 -*-
"""喷涂段可达性修复：对 IK 失败 / 关节超限的喷涂点，绕水平方向倾斜工具轴，
在保持 TCP 位置（喷距）不变的前提下，把关节拉回软限位内。

只改这些点的 IJK（姿态），XYZ 不动；喷涂角度会偏离原法向 tilt 度（记录在报告里）。
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys

import numpy as np

from . import robot, geom

J5_LIM = 120.0


def _dir(r):
    v = np.array([float(r['I']), float(r['J']), float(r['K'])])
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else np.array([0.0, 0.0, -1.0])


def _pos(r):
    return np.array([float(r['X']), float(r['Y']), float(r['Z'])])


def solve(p, g, tool, seed):
    R = geom.pose_from_z(g)[:3, :3]
    # 快路径：单初值局部解（~5ms）。解不出再退回多初值。
    if seed is not None:
        sol = robot._ik_once(p, R, tool, seed, 50, 0.05, 0.05)
        if sol.ok and not sol.over:
            th = sol.theta + 360.0 * np.round((seed - sol.theta) / 360.0)
            return th, False
    sol = robot.ik(p, R, tool, seed=seed)
    if not sol.ok:
        return None, True
    th = sol.theta.copy()
    if seed is not None:
        th = th + 360.0 * np.round((seed - th) / 360.0)
    return th, bool(sol.over)


def _tilt(g, axis, deg):
    r = math.radians(deg)
    v = math.cos(r) * g + math.sin(r) * axis
    return v / np.linalg.norm(v)


def run(csv_in, csv_out, report_out, max_tilt=40.0, seed_from_csv=True):
    rows = list(csv.DictReader(open(csv_in, encoding='utf-8-sig')))
    fields = list(rows[0].keys())
    tool = robot.Tool(length=453.08)
    spray = [int(r['spray']) for r in rows]
    xs = [_pos(r) for r in rows]
    ijk = [_dir(r) for r in rows]

    # 初始 seed：用第一个能解出来的点
    seed = None
    fixed = 0
    failed = 0
    over0 = 0
    details = []

    for k in range(len(rows)):
        if spray[k] != 1:
            continue
        p, g = xs[k], ijk[k]
        th, over = solve(p, g, tool, seed)
        if th is not None and not over:
            seed = th
            continue
        if th is None:
            failed += 1
        else:
            over0 += 1
        # 倾斜搜索：四个水平方向 × 5..max_tilt
        base_h = np.cross(g, np.array([0.0, 0.0, 1.0]))
        if np.linalg.norm(base_h) < 1e-6:
            base_h = np.array([1.0, 0.0, 0.0])
        base_h /= np.linalg.norm(base_h)
        perp = np.cross(g, base_h)
        perp /= np.linalg.norm(perp)
        best = None
        for axis in (perp, -perp, base_h, -base_h):
            for deg in (5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0):
                g2 = _tilt(g, axis, float(deg))
                th2, over2 = solve(p, g2, tool, seed)
                if th2 is not None and not over2:
                    best = (g2, th2, float(deg))
                    break
            if best:
                break
        if best:
            ijk[k] = best[0]
            seed = best[1]
            fixed += 1
            details.append({'idx': k, 'tilt': best[2], 'was_fail': th is None})
        else:
            failed += 1 if th is None else 0
            over0 += 0 if th is None else 1

    with open(csv_out, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(fields)
        for i, r in enumerate(rows):
            row = [r.get(c, '') for c in fields]
            row[fields.index('I')] = '%.6f' % ijk[i][0]
            row[fields.index('J')] = '%.6f' % ijk[i][1]
            row[fields.index('K')] = '%.6f' % ijk[i][2]
            w.writerow(row)

    rep = {'csv': csv_in, 'out': csv_out, 'max_tilt': max_tilt,
           'initial_fail': failed, 'initial_over': over0, 'fixed': fixed,
           'unfixed': (failed + over0) - fixed,
           'tilt_used': [d['tilt'] for d in details]}
    if report_out:
        open(report_out, 'w', encoding='utf-8').write(json.dumps(rep, ensure_ascii=False, indent=1))
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('csv'); ap.add_argument('out')
    ap.add_argument('--json', default=None)
    ap.add_argument('--max-tilt', type=float, default=40.0)
    a = ap.parse_args(argv)
    r = run(a.csv, a.out, a.json, a.max_tilt)
    print(json.dumps({k: v for k, v in r.items() if k != 'tilt_used'}, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
