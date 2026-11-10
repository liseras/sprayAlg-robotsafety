# -*- coding: utf-8 -*-
"""Constrained idle-segment pose planner (constrained back-off from -Z).

Rules (2026-10-29, user-confirmed defaults):
  * Target pose of an idle segment = tool axis along -Z (nozzle down).
  * Three gates per point, all must pass:
      1) gun-body clearance to cloud  >= --clear-need  (default 15 mm)
      2) robot self-collision gap     >= --self-need   (default 10 mm)
      3) J5 within +-120 deg (soft limit) and IK converged
  * Pick the SMALLEST tilt angle (0 = exactly -Z) whose segment passes all gates.
    Never fall back to "horizontal" unless that is literally the only feasible one.
  * Spray points (spray=1) are never touched; seam points keep the neighbour spray pose.
  * Per-step tool-axis rotation <= --max-step deg.

Usage:
  python -m robot_safety.move_pose <traj.csv> <cloud.ply> <out.csv> [--json rep.json]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys

import numpy as np

from . import robot, geom
from .move_plan import CloudGrid, read_ply_xyz, SEARCH_R

GUN_R = 20.0
J5_LIM = 120.0


def _dir(r):
    v = np.array([float(r['I']), float(r['J']), float(r['K'])])
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else np.array([0.0, 0.0, -1.0])


def _slerp(a, b, t):
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    d = float(np.clip(np.dot(a, b), -1.0, 1.0))
    if d > 0.9995:
        v = a * (1 - t) + b * t
        return v / np.linalg.norm(v)
    th = math.acos(d); s = math.sin(th)
    return a * math.sin((1 - t) * th) / s + b * math.sin(t * th) / s


def _smooth(t):
    return t * t * (3.0 - 2.0 * t)


def _segments(spray):
    out, n, i = [], len(spray), 0
    while i < n:
        if spray[i] != 0:
            i += 1
            continue
        j = i
        while j + 1 < n and spray[j + 1] == 0:
            j += 1
        out.append((i, j))
        i = j + 1
    return out


def build_dirs(g0, g1, pole, n, max_step_deg):
    """g0 -> pole -> g1, smoothstep in t, then cap per-step angle."""
    dirs = []
    for k in range(n):
        t = _smooth(k / max(1, n - 1))
        if t < 0.5:
            dirs.append(_slerp(g0, pole, t * 2))
        else:
            dirs.append(_slerp(pole, g1, (t - 0.5) * 2))
    dirs[0] = g0.copy(); dirs[-1] = g1.copy()
    # forward/backward cap so no step exceeds max_step_deg
    cap = math.radians(max_step_deg)
    for _ in range(3):
        for k in range(1, n):
            d = math.acos(float(np.clip(np.dot(dirs[k - 1], dirs[k]), -1, 1)))
            if d > cap:
                dirs[k] = _slerp(dirs[k - 1], dirs[k], cap / d)
        for k in range(n - 2, -1, -1):
            d = math.acos(float(np.clip(np.dot(dirs[k + 1], dirs[k]), -1, 1)))
            if d > cap:
                dirs[k] = _slerp(dirs[k + 1], dirs[k], cap / d)
        dirs[0] = g0.copy(); dirs[-1] = g1.copy()
    return dirs


def _align(a, b):
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); s = float(np.linalg.norm(v)); c = float(np.clip(np.dot(a, b), -1, 1))
    if s < 1e-9:
        if c > 0:
            return np.eye(3)
        axis = np.cross(a, [1.0, 0.0, 0.0])
        if np.linalg.norm(axis) < 1e-6:
            axis = np.cross(a, [0.0, 1.0, 0.0])
        axis = axis / np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def eval_segment(points, dirs, grid, tool, sizing, pairs, seed0, clear_need, self_need,
                 R0=None):
    res = dict(clear_min=math.inf, self_min=math.inf, j5_max=0.0, ik_fail=0,
               fail_clear=0, fail_self=0, fail_j5=0, ok=True, ik_fail_at=[], self_fail_at=[], j5_fail_at=[])
    seed = seed0
    R = R0 if R0 is not None else geom.pose_from_z(dirs[0])[:3, :3]
    for k, (p, g) in enumerate(zip(points, dirs)):
        clear = grid.seg_clear(p, p - g * tool.length, GUN_R)
        res['clear_min'] = min(res['clear_min'], clear)
        if clear < clear_need:
            res['fail_clear'] += 1
        if k > 0:
            R = _align(dirs[k - 1], dirs[k]) @ R
        sol = robot.ik(p, R, tool, seed=seed)
        if not sol.ok:
            res['ik_fail'] += 1
            res['ik_fail_at'].append(k)
            continue
        seed = sol.theta
        j5 = abs(float(sol.theta[4]))
        res['j5_max'] = max(res['j5_max'], j5)
        if j5 > J5_LIM:
            res['fail_j5'] += 1
            res['j5_fail_at'].append(k)
        fk = robot.fk(sol.theta, tool)
        segs = robot.arm_segments(fk, sizing)
        gap = min(geom.capsule_capsule_dist(segs[a].a, segs[a].b, segs[b].a, segs[b].b)
                  - segs[a].r - segs[b].r for a, b in pairs)
        res['self_min'] = min(res['self_min'], gap)
        if gap < self_need:
            res['fail_self'] += 1
            res['self_fail_at'].append(k)
    res['ok'] = (res['ik_fail'] == 0 and res['fail_clear'] == 0
                 and res['fail_self'] == 0 and res['fail_j5'] == 0)
    return res


def plan(csv_in, cloud_in, csv_out, clear_need=15.0, self_need=10.0, max_step=5.0,
         tool_len=453.08, tilts=None):
    pts = read_ply_xyz(cloud_in)
    grid = CloudGrid(pts)
    rows = list(csv.DictReader(open(csv_in, encoding='utf-8-sig')))
    fieldnames = list(rows[0].keys())
    spray = [int(r['spray']) for r in rows]
    xs = np.array([[float(r['X']), float(r['Y']), float(r['Z'])] for r in rows])
    ijk = np.array([_dir(r) for r in rows])
    speed = [r['speed'] for r in rows]

    tool = robot.Tool(length=tool_len)
    sizing = robot.ArmSizing.of()
    pairs = robot.self_collision_pairs(
        robot.arm_segments(robot.fk(np.zeros(6), tool), sizing))

    if tilts is None:
        tilts = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

    report = {'csv': csv_in, 'cloud': cloud_in,
              'params': {'clear_need': clear_need, 'self_need': self_need,
                         'max_step': max_step, 'tool_len': tool_len},
              'segments': []}

    for (i0, i1) in _segments(spray):
        if i1 - i0 < 1:
            continue
        g0 = ijk[i0 - 1] if i0 > 0 else ijk[i1 + 1]
        g1 = ijk[i1 + 1] if i1 + 1 < len(rows) else g0
        P = xs[i0:i1 + 1]
        n = i1 - i0 + 1

        h = g0 + g1
        h[2] = 0.0
        nh = np.linalg.norm(h)
        h = h / nh if nh > 1e-6 else np.array([1.0, 0.0, 0.0])

        cands = []
        seed0 = None
        # R0 must have its z column equal to g0 (CSV IJK), never the FK joint pose
        # (safe.csv's J columns do not reproduce the IJK: dot(z_fk, ijk) ~ 0.68 measured).
        R0 = geom.pose_from_z(g0)[:3, :3]
        for tilt in tilts:
            rad = math.radians(tilt)
            pole = -math.cos(rad) * np.array([0.0, 0.0, 1.0]) + math.sin(rad) * h
            dirs = build_dirs(g0, g1, pole, n, max_step)
            ev = eval_segment(P, dirs, grid, tool, sizing, pairs, seed0, clear_need, self_need, R0)
            ev['tilt'] = tilt
            ev['dirs'] = dirs
            cands.append(ev)

        good = [c for c in cands if c['ok']]
        if good:
            chosen = min(good, key=lambda c: c['tilt'])
            feasible = True
        else:
            # best effort: fewest gate failures, then smallest tilt
            chosen = min(cands, key=lambda c: (c['fail_clear'] + c['fail_self'] + c['fail_j5']
                                               + c['ik_fail'], c['tilt']))
            feasible = False

        for k in range(n):
            ijk[i0 + k] = chosen['dirs'][k]

        report['segments'].append({
            'idx': [int(i0), int(i1)], 'n': int(n), 'feasible': feasible,
            'tilt': chosen['tilt'],
            'clear_min': round(float(chosen['clear_min']), 2),
            'self_min': round(float(chosen['self_min']), 2),
            'j5_max': round(float(chosen['j5_max']), 1),
            'fail_clear': int(chosen['fail_clear']),
            'fail_self': int(chosen['fail_self']),
            'fail_j5': int(chosen['fail_j5']),
            'ik_fail': int(chosen['ik_fail']),
            'candidates': [{k2: (round(float(v), 2) if isinstance(v, (int, float)) else v)
                            for k2, v in c.items() if k2 != 'dirs'} for c in cands],
        })

    with open(csv_out, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(fieldnames)
        for i, r in enumerate(rows):
            row = [r.get(c, '') for c in fieldnames]
            row[fieldnames.index('I')] = '%.6f' % ijk[i][0]
            row[fieldnames.index('J')] = '%.6f' % ijk[i][1]
            row[fieldnames.index('K')] = '%.6f' % ijk[i][2]
            w.writerow(row)
    return report


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('csv'); ap.add_argument('cloud'); ap.add_argument('out')
    ap.add_argument('--clear-need', type=float, default=15.0)
    ap.add_argument('--self-need', type=float, default=10.0)
    ap.add_argument('--max-step', type=float, default=5.0)
    ap.add_argument('--tool-len', type=float, default=453.08)
    ap.add_argument('--json', default=None)
    a = ap.parse_args(argv)
    rep = plan(a.csv, a.cloud, a.out, a.clear_need, a.self_need, a.max_step, a.tool_len)
    txt = json.dumps(rep, ensure_ascii=False, indent=2)
    if a.json:
        open(a.json, 'w', encoding='utf-8').write(txt)
    print(txt)
    return 0


if __name__ == '__main__':
    sys.exit(main())
