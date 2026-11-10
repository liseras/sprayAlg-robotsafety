# -*- coding: utf-8 -*-
"""J 型空移：把空移段折叠成「关节点列」，控制器按关节插补走，中间不生成 P 型点。

设计（2026-10-29 用户拍板：关节运动 + 多路点 + 先仿真碰撞）：
  1. 空移段 → 一连串关节点（起点接缝、竖直路点、终点接缝），相邻点间距约 --step-mm。
  2. 每个关节点求解关节角；**用上一点解做唯一初值**保证解支连续（不跳分支）。
  3. 接缝点（段首/段末）姿态允许向竖直方向偏移最多 --seam-tilt 度，用于救「接缝姿态本身不可达」。
  4. 仿真：相邻关节点之间按关节线性插值（控制器 @J 的行为）采样，逐点算
     枪身净空 / 连杆自碰撞 / J5，输出最小值和是否碰撞。

用法：
  python -m robot_safety.move_joint <traj.csv> <cloud.ply> <out.json> [--step-mm 50]
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
from .move_pose import build_dirs, _slerp

GUN_R = 20.0
J5_LIM = 120.0
MAX_JUMP = 30.0


def _dir(r):
    v = np.array([float(r['I']), float(r['J']), float(r['K'])])
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else np.array([0.0, 0.0, -1.0])


def _pos(r):
    return np.array([float(r['X']), float(r['Y']), float(r['Z'])])


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


def _tilt(vertical, horiz, deg):
    """从竖直方向朝水平方向 horiz 倾斜 deg（deg=0 即竖直）。"""
    r = math.radians(deg)
    v = math.cos(r) * vertical + math.sin(r) * horiz
    return v / np.linalg.norm(v)


def solve_joint(p, g, tool, seed, R_prev=None, j5_lim=J5_LIM):
    """解一个关节点的关节角。seed=None 时用多初值；否则用唯一初值（连续解支）。"""
    R = geom.pose_from_z(g)[:3, :3] if R_prev is None else R_prev
    if seed is None:
        sol = robot.ik(p, R, tool)
    else:
        sol = robot._ik_once(p, R, tool, seed, 60, 0.02, 0.02)
    if not sol.ok:
        return None
    if j5_lim is not None and abs(float(sol.theta[4])) > j5_lim:
        return None
    th = sol.theta.copy()
    if seed is not None:
        th = th + 360.0 * np.round((seed - th) / 360.0)
    return th


def solve_any_tilt(p, g_vert, g_seam, tool, seed, horiz, max_tilt=45.0):
    """竖直解不出/超限时，向「接缝方向」或「离件方向」逐步倾斜，直到解出且 J5 合规。

    返回 (theta, used_deg, g)。全都不行则 (None, None, g_vert)。
    """
    for deg in np.arange(0.0, max_tilt + 0.1, 5.0):
        g = _tilt(g_vert, horiz, float(deg))
        th = solve_joint(p, g, tool, seed)
        if th is not None:
            return th, float(deg), g
    return None, None, g_vert


def simulate(pairs, pts_list, tool, sizing, grid, interp=20):
    worst_clear = math.inf
    worst_self = math.inf
    worst_j5 = 0.0
    max_jump = 0.0
    samples = 0
    valid = [p for p in pts_list if p is not None]
    for j in range(len(valid) - 1):
        a, b = valid[j], valid[j + 1]
        max_jump = max(max_jump, float(np.max(np.abs(b - a))))
        steps = max(interp, int(float(np.max(np.abs(b - a))) / 2.0))
        for t in range(steps + 1):
            th = a + (b - a) * (t / steps)
            fk = robot.fk(th, tool)
            tcp = fk.tcp[:3, 3]
            gd = fk.tcp[:3, 2] / np.linalg.norm(fk.tcp[:3, 2])
            worst_clear = min(worst_clear, grid.seg_clear(tcp, tcp - gd * tool.length, GUN_R))
            segs = robot.arm_segments(fk, sizing)
            worst_self = min(worst_self, min(
                geom.capsule_capsule_dist(segs[x].a, segs[x].b, segs[y].a, segs[y].b)
                - segs[x].r - segs[y].r for x, y in pairs))
            worst_j5 = max(worst_j5, abs(float(th[4])))
            samples += 1
    return dict(clear_min=round(worst_clear, 2), self_min=round(worst_self, 2),
                j5_max=round(worst_j5, 1), max_jump=round(max_jump, 1),
                samples=samples,
                collision=bool(worst_clear < 0 or worst_self < 0))


def run(csv_in, cloud_in, out_json, step_mm=50.0, seam_tilt=15.0, max_jump=MAX_JUMP):
    rows = list(csv.DictReader(open(csv_in, encoding='utf-8-sig')))
    grid = CloudGrid(read_ply_xyz(cloud_in))
    tool = robot.Tool(length=453.08)
    sizing = robot.ArmSizing.of()
    pairs = robot.self_collision_pairs(robot.arm_segments(robot.fk(np.zeros(6), tool), sizing))
    spray = [int(r['spray']) for r in rows]

    report = {'csv': csv_in, 'cloud': cloud_in,
              'params': dict(step_mm=step_mm, seam_tilt=seam_tilt, max_jump=max_jump),
              'segments': []}

    for (i0, i1) in _segments(spray):
        if i1 - i0 < 1:
            continue
        g0 = _dir(rows[i0 - 1]) if i0 > 0 else _dir(rows[i1 + 1])
        g1 = _dir(rows[i1 + 1]) if i1 + 1 < len(rows) else g0
        P = np.array([_pos(rows[k]) for k in range(i0, i1 + 1)])
        L = float(np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1)))
        h = g0 + g1; h[2] = 0.0
        nh = np.linalg.norm(h); h = h / nh if nh > 1e-6 else np.array([1.0, 0.0, 0.0])
        pole = np.array([0.0, 0.0, -1.0])
        dirs = build_dirs(g0, g1, pole, len(P), 5.0)
        s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
        n_wp = max(2, int(math.ceil(L / step_mm)) + 1)

        # ---- 选接缝姿态的允许偏移量：先试原姿态，解不出/跳支就向竖直偏 ----
        def seam_ok(pos, g, seed):
            th = solve_joint(pos, g, tool, seed)
            if th is None:
                return None, 'ik'
            if seed is not None and float(np.max(np.abs(th - seed))) > max_jump:
                return th, 'jump'
            return th, ''

        wp = []
        seed = None
        centroid = grid.centroid
        for k in range(n_wp):
            t = L * k / (n_wp - 1)
            idx = min(int(np.searchsorted(s, t)), len(P) - 1)
            pnt, g = P[idx], dirs[idx]
            hd = pnt - centroid
            hd[2] = 0.0
            nhd = np.linalg.norm(hd)
            hdir = hd / nhd if nhd > 1e-6 else np.array([1.0, 0.0, 0.0])
            seam = (k == 0 or k == n_wp - 1)
            th = None
            used = 0.0
            if seam:
                th = solve_joint(pnt, g, tool, seed)
                if th is None:
                    th, used, g2 = solve_any_tilt(pnt, pole, g, tool, seed, hdir, seam_tilt)
                    if th is not None:
                        g = g2
            else:
                th, used, g2 = solve_any_tilt(pnt, pole, g, tool, seed, hdir, 45.0)
                if th is not None:
                    g = g2
            wp.append(dict(i=idx, p=pnt.tolist(), g=g.tolist(),
                           th=None if th is None else th.tolist(),
                           tilt_deg=used))
            if th is not None:
                seed = th

        ths = [None if w['th'] is None else np.array(w['th']) for w in wp]
        sim = simulate(pairs, ths, tool, sizing, grid)
        report['segments'].append(dict(
            idx=[int(i0), int(i1)], len_mm=round(L, 1), waypoints=len(wp),
            ik_fail=int(sum(1 for t in ths if t is None)),
            **sim,
            points=wp))

    open(out_json, 'w', encoding='utf-8').write(json.dumps(report, ensure_ascii=False, indent=1))
    return report


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('csv'); ap.add_argument('cloud'); ap.add_argument('out')
    ap.add_argument('--step-mm', type=float, default=50.0)
    ap.add_argument('--seam-tilt', type=float, default=15.0)
    ap.add_argument('--max-jump', type=float, default=30.0)
    a = ap.parse_args(argv)
    rep = run(a.csv, a.cloud, a.out, a.step_mm, a.seam_tilt, a.max_jump)
    for s in rep['segments']:
        print('seg %d-%d len %.0f waypoints %d ik_fail %d clear %.1f self %.1f j5 %.1f max_jump %.1f collision %s'
              % (s['idx'][0], s['idx'][1], s['len_mm'], s['waypoints'], s['ik_fail'],
                 s['clear_min'], s['self_min'], s['j5_max'], s['max_jump'], s['collision']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
