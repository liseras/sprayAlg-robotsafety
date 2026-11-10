# -*- coding: utf-8 -*-
"""空移段姿态重规划（move_plan）。

病灶（2026-10-28 实测，见 outputs/diag_move/）：
  生成器给空移段的姿态是「面法向的零阶保持」——一段 900mm 的空移从头到尾一个朝向，
  在换面处**原地**从这一面的法向跳到另一面的法向。实测：
    * 开口盒子 #972→#973：位置一点没动（Δpos=0.000mm），枪轴翻 180.0°；
    * 槽形件 #995→#996：走 4.00mm，枪轴翻 146.3°（等效 7315°/s）；
    * 曲型板 #1068→#1069：走 4.00mm，枪轴翻 101.2°。
  真机上这些点要么报「轨迹不可执行」，要么控制器自己插补出一个穿过工件的翻转弧。

本模块只改**空移点（spray=0）的 IJK**，不动位置、不动喷涂点，规则三条：

  1. 姿态变化只发生在「离工件最空」的地方。用每点到点云的净空 d 做权重
     w=(d-floor)²，把姿态的进度 s 按 w 累计——d 小的点 s 几乎不动，枪在低处不转身。
  2. 两端姿态向相邻喷涂点对齐（首点=上一喷涂点朝向、末点=下一喷涂点朝向），
     所以喷涂段一个数都不变，接缝处也没有跳变。
  3. 中间经过「中性姿态」gn = 从 TCP 指向工件质心的方向（枪体朝外、喷嘴朝内），
     分两段 slerp：g0→gn→g1。180° 的翻转因此变成两次 <90° 的转身，
     且中间时刻枪体是竖起来的，而不是横着扫过薄板。

用法:
    python -m robot_safety.move_plan <traj.csv> <cloud.ply> <out.csv> \
        [--standoff 100] [--tool-len 453.08] [--gun-radius 20] [--floor 30] [--json out.json]
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from . import traj as traj_mod
from . import geom


# 碰撞查询的搜索半径。枪身净空判定的阈值只有 15mm，80mm 的搜索窗足够，
# 超出窗口一律返回窗口值（下界）—— 绝不当成实测值。
SEARCH_R = 80.0


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])


def read_ply_xyz(path: str) -> np.ndarray:
    """ASCII PLY 的顶点（样本点云都是 ascii，与 cpp/_regress 的 read_ply 同口径）。"""
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


class CloudGrid:
    """点云最近距离查询（网格哈希，零 scipy 依赖）。

    安全层的 cloud.CloudCollider 用 scipy 的 cKDTree；本机托管 python 没装 scipy，
    而原型要能在现场随手跑起来，所以这里自带一个够用的等价物：
    20mm 网格 + 逐格候选点，查询半径外的点直接跳过（不会把下界当实测值）。
    """

    def __init__(self, pts: np.ndarray, cell: float = 20.0):
        self.pts = np.ascontiguousarray(pts, dtype=np.float64)
        self.cell = float(cell)
        self.centroid = self.pts.mean(axis=0)
        self.org = self.pts.min(axis=0)
        ijk = np.floor((self.pts - self.org) / self.cell).astype(np.int64)
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
        for ix in range(i0[0], i1[0] + 1):
            for iy in range(i0[1], i1[1] + 1):
                k0 = (ix * dy + iy) * dz + i0[2]
                k1 = (ix * dy + iy) * dz + i1[2]
                a = np.searchsorted(self.key, k0, "left")
                b = np.searchsorted(self.key, k1, "right")
                if b > a:
                    out.append(self.idx[a:b])
        return np.concatenate(out) if out else None

    def nearest(self, p: np.ndarray, reach: float = SEARCH_R) -> float:
        idx = self.cand(p - reach, p + reach)
        if idx is None:
            return reach
        return float(np.linalg.norm(self.pts[idx] - p, axis=1).min())

    def nearest_many(self, P: np.ndarray, reach: float = SEARCH_R) -> np.ndarray:
        return np.array([self.nearest(p, reach) for p in P])

    def seg_clear(self, a: np.ndarray, b: np.ndarray, radius: float,
                  reach: float = SEARCH_R) -> float:
        """线段/胶囊到点云的表面距离（axis_dist - radius）。"""
        idx = self.cand(np.minimum(a, b) - reach, np.maximum(a, b) + reach)
        if idx is None:
            return reach
        Q = self.pts[idx]
        AB = b - a
        L2 = float(AB @ AB)
        if L2 < 1e-12:
            d = np.linalg.norm(Q - a, axis=1)
        else:
            t = np.clip(((Q - a) @ AB) / L2, 0.0, 1.0)
            d = np.linalg.norm(Q - (a + t[:, None] * AB), axis=1)
        return float(d.min()) - radius


def _smoothstep(t: float) -> float:
    t = max(0.0, min(1.0, t))
    return t * t * (3.0 - 2.0 * t)


def _slerp(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    """最短弧插值。调用方保证 a,b 夹角 < 170°（中间必过中性姿态）。"""
    d = float(np.clip(np.dot(a, b), -1.0, 1.0))
    if d > 0.9995:
        return _unit(a * (1.0 - t) + b * t)
    th = np.arccos(d)
    s = np.sin(th)
    return _unit(a * (np.sin((1.0 - t) * th) / s) + b * (np.sin(t * th) / s))


def _segments(spray: np.ndarray) -> list[tuple[int, int]]:
    """spray=0 的连续段 [i0, i1]。"""
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


def _angle_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    c = np.clip(np.einsum("ij,ij->i", a, b), -1.0, 1.0)
    return np.degrees(np.arccos(c))


def _gun_body_min(collider, pos: np.ndarray, gun: np.ndarray, tool_len: float,
                  radius: float, search: float = 200.0) -> np.ndarray:
    """枪身胶囊（TCP → TCP-tool_len*gun ... 注意枪体在 -gun 方向）到点云净空。"""
    out = np.full(len(pos), np.inf)
    for i in range(len(pos)):
        a = pos[i]
        b = pos[i] - gun[i] * tool_len
        h = collider.seg_clear(a, b, radius, search)
        out[i] = h
    return out


def stats(traj: traj_mod.Traj, label: str) -> dict:
    p, g, sp = traj.pos, traj.gun, traj.spray
    d = np.linalg.norm(np.diff(p, axis=0), axis=1)
    ang = _angle_deg(g[1:], g[:-1])
    spd = float(np.median(traj.speed[traj.speed > 0])) if (traj.speed > 0).any() else 200.0
    rate = ang / np.maximum(d, 1e-6) * spd
    mv = (sp[1:] == 0) | (sp[:-1] == 0)
    zero_move = (d < 0.05)
    big = ang > 45.0
    return {
        "label": label,
        "n": int(len(p)),
        "move_pts": int((sp == 0).sum()),
        "jump45": int(big.sum()),
        "jump90": int((ang > 90.0).sum()),
        "jump_in_move": int((big & mv).sum()),
        "zero_move_jump": int((big & zero_move).sum()),
        "max_step_deg": float(ang.max()) if len(ang) else 0.0,
        "max_rate_dps": float(np.max(rate[np.isfinite(rate)])) if len(ang) else 0.0,
        "p95_rate_dps": float(np.percentile(rate[np.isfinite(rate)], 95)) if len(ang) else 0.0,
        "move_max_rate_dps": float(np.max(rate[mv])) if mv.any() else 0.0,
        "move_p95_rate_dps": float(np.percentile(rate[mv], 95)) if mv.any() else 0.0,
    }


def _nearest_point(collider, p: np.ndarray, reach: float = SEARCH_R):
    idx = collider.cand(p - reach, p + reach)
    if idx is None:
        return p + np.array([0.0, 0.0, 1.0]), reach
    Q = collider.pts[idx]
    d = np.linalg.norm(Q - p, axis=1)
    k = int(np.argmin(d))
    return Q[k], float(d[k])


def _sph_deg(v: np.ndarray) -> tuple[float, float]:
    """方向 → (仰角, 方位角) 度。仰角以 +Z 为 +90、水平为 0、-Z 为 −90。"""
    v = _unit(v)
    return float(np.degrees(np.arcsin(np.clip(v[2], -1.0, 1.0)))), \
        float(np.degrees(np.arctan2(v[1], v[0])))


def _from_sph_deg(el: float, az: float) -> np.ndarray:
    e, a = np.radians(el), np.radians(az)
    return np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])


def _wrap_deg(a: float) -> float:
    return (a + 180.0) % 360.0 - 180.0


def _solve_segment(P, g0, g1, dclr, centroid, collider, tool_len, gun_radius,
                   clear_floor, neutral, max_step_deg, move_dir: str = "z-"):
    """给一段空移点解姿态（位置不动）。返回 gun 数组 + 诊断。

    转向路径（2026-10-28 用户口径）：**空移段枪轴要沿 Z 轴，从浅到深**，即 −Z。

    之前两版都错了：
      · 第一版"中性姿态选最不贴件"→ 枪被抬到竖直再落回，画面上仰角一路 −90~0；
      · 第二版"等高插值"→ 侧壁↔侧壁的空移变成**在水平面内转 180°**，枪身横着扫一大圈。
    两版的共同毛病是"空移段自身的朝向没有规定"，只保证两端接得上。

    现在明确规定：空移段的主体朝向 = −Z（喷枪对着工件、枪体朝上，从 Z 浅指向 Z 深）。
      · 前 1/2 段：仰角从 el0 抬到 −90，方位角保持 az0（沿经线升到极点）；
      · 后 1/2 段：从极点落到 el1，方位角取 az1。
    极点处方位角无定义、但方向向量连续，所以枪轴是一根竖直线，不会横着扫。
      · 两端仍严格等于相邻喷涂点的姿态 ⟹ 接缝零跳变。
      · 若两端本来就是竖直（顶面↔顶面），整段一点都不转。
    """
    n = len(P)
    el0, az0 = _sph_deg(g0)
    el1, az1 = _sph_deg(g1)
    sign = -1.0 if move_dir == "z-" else 1.0
    pole = sign * 90.0
    theta = max(abs(el0 - pole) + abs(pole - el1), 1e-3)

    # ---- 转身进度 s：按净空加权，低处几乎不转 ----
    d = dclr
    w = np.maximum(0.0, d - clear_floor) + 1e-6
    s = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]))])
    s = s / s[-1] if s[-1] > 1e-9 else np.linspace(0.0, 1.0, n)

    # ---- 每步转角限幅 ----
    dd = np.maximum(0.0, d - clear_floor)
    gain = 0.25 + 0.75 * (dd / dd.max() if dd.max() > 1e-9 else np.ones_like(dd))
    cap0 = (max_step_deg / theta) * gain
    s_out = s.copy()
    for attempt in range(25):
        cap = cap0 * (1.25 ** attempt)
        x = s.copy()
        x[0] = 0.0
        for _ in range(2):
            for i in range(1, n):
                x[i] = min(x[i], x[i - 1] + cap[i])
            for i in range(n - 2, -1, -1):
                x[i] = max(x[i], x[i + 1] - cap[i + 1])
        x[0] = 0.0
        s_out = x
        if s_out[-1] >= 1.0 - 1e-9:
            break
    s = np.maximum.accumulate(np.minimum(s_out, 1.0))
    s[0], s[-1] = 0.0, 1.0

    G = np.zeros((n, 3))
    if move_dir == "level":
        # 保持仰角（在两端仰角之间线性过渡），只在水平面内转方位角。
        # 两端同向 ⟹ 整段不转；两端反向 ⟹ 在同一个水平面里转半圈（枪身横扫）。
        daz = _wrap_deg(az1 - az0)
        for k in range(n):
            t = _smoothstep(s[k])
            G[k] = _from_sph_deg(el0 + (el1 - el0) * t, az0 + daz * t)
    else:
        for k in range(n):
            t = _smoothstep(s[k])
            if t < 0.5:
                u = t * 2.0
                G[k] = _from_sph_deg(el0 + (pole - el0) * u, az0)
            else:
                u = (t - 0.5) * 2.0
                G[k] = _from_sph_deg(pole + (el1 - pole) * u, az1)

    # 枪身净空兜底：低于 soft 的点保持上一点姿态（宁可少转，不许扫过工件）
    soft = max(6.0, 0.5 * gun_radius)
    n_fix = 0
    for k in range(1, n):
        if collider.seg_clear(P[k], P[k] - G[k] * tool_len, gun_radius, SEARCH_R) < soft:
            G[k] = G[k - 1]
            n_fix += 1
    G[0], G[-1] = g0, g1
    body_min = min(collider.seg_clear(P[k], P[k] - G[k] * tool_len, gun_radius, SEARCH_R)
                   for k in range(n))
    return G, {"theta": float(theta), "s": s, "n_fix": int(n_fix),
               "el0": el0, "el1": el1, "pole": pole,
               "body_min": float(body_min)}


def _raise_segment(P, G, collider, need, tool_len, radius):
    """枪身净空不足的点，沿「离开最近工件点」方向推开。首末点锁定（接缝不动）。"""
    n = len(P)
    disp = np.zeros((n, 3))
    moved = False
    worst = 0.0
    for k in range(1, n - 1):
        cl = collider.seg_clear(P[k], P[k] - G[k] * tool_len, radius, SEARCH_R)
        if cl >= need:
            continue
        np_pt, _ = _nearest_point(collider, P[k])
        u = P[k] - np_pt
        ln = float(np.linalg.norm(u))
        u = u / ln if ln > 1e-6 else np.array([0.0, 0.0, 1.0])
        step = (need - cl + 0.5) * 1.2
        disp[k] = u * step
        worst = max(worst, step)
        moved = True
    if not moved:
        return P, 0.0
    ds = disp.copy()
    for k in range(1, n - 1):
        ds[k] = 0.5 * disp[k] + 0.25 * disp[k - 1] + 0.25 * disp[k + 1]
    ds[0] = 0.0
    ds[-1] = 0.0
    return P + ds, float(np.linalg.norm(ds, axis=1).max())


def plan(traj: traj_mod.Traj, collider, standoff: float = 100.0, tool_len: float = 453.08,
         gun_radius: float = 20.0, clear_floor: float = 30.0,
         neutral: str = "centroid", max_step_deg: float = 5.0, move_dir: str = "z-",
         clear_need: float = 15.0, raise_iters: int = 4) -> tuple[traj_mod.Traj, dict]:
    """空移段重规划：先解姿态，再按枪身净空抬路径，交替迭代。

    只改 spray=0 点的位置与 IJK；喷涂点（spray=1）一个字节都不动。
    首末空移点锁定位置与姿态（= 相邻喷涂点的姿态），所以接缝零跳变。
    """
    n = len(traj)
    pos = traj.pos.copy()
    gun = traj.gun.copy()
    spray = traj.spray
    centroid = collider.centroid
    segs = _segments(spray)

    info = []
    for (i0, i1) in segs:
        if i1 - i0 < 1:
            continue
        g0 = gun[i0 - 1].copy() if i0 > 0 else (gun[i1 + 1].copy() if i1 + 1 < n
                                              else np.array([0.0, 0.0, -1.0]))
        g1 = gun[i1 + 1].copy() if i1 + 1 < n else g0.copy()
        idx = np.arange(i0, i1 + 1)
        P = pos[idx].copy()
        P0 = P.copy()
        L0 = float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())

        G = np.tile(g0, (len(P), 1))
        meta = {}
        raised = 0.0
        infeasible = False
        modes = [move_dir] if move_dir != "auto" else ["level", "z-"]
        for it in range(max(1, raise_iters)):
            dclr = collider.nearest_many(P)
            # 多个候选路径各解一次，取"枪身最不贴件"的那条。
            # box 实测：level 净空 14.0mm、z- 净空 0.4mm（枪体竖起来顶到腔顶），
            # 自动选优会选 level；薄壁开放件反过来会选 z-（枪身不横扫）。
            cands = []
            for md in modes:
                Gx, mx = _solve_segment(P, g0, g1, dclr, centroid, collider, tool_len,
                                        gun_radius, clear_floor, neutral, max_step_deg, md)
                cands.append((mx.get("body_min", -1e9), Gx, mx, md))
            cands.sort(key=lambda z: -z[0])
            _, G, meta, md_used = cands[0]
            meta["mode"] = md_used
            if clear_need <= 0 or it == raise_iters - 1:
                break
            # 几何上就过不去的段（最优姿态下枪身仍穿透 >5mm）不要抬：
            # 开口盒子用 453mm 枪就是这种 —— 腔深 300mm 比枪还短，抬 TCP 只会把
            # 路径拉长 51%（实测 913→1383mm）而净空仍 -19.5mm，纯属毁轨迹。
            # 这类只能改站位/加变位机，如实报"不可行"，让工艺端决策。
            if meta.get("body_min", 0.0) < -5.0:
                infeasible = True
                break
            P2, step = _raise_segment(P, G, collider, clear_need, tool_len, gun_radius)
            raised = max(raised, step)
            if step <= 0.2:
                break
            P = P2
        # 保险：抬完再解一次姿态（位置变了，姿态权重也要跟着变）
        dclr = collider.nearest_many(P)
        cands = []
        for md in modes:
            Gx, mx = _solve_segment(P, g0, g1, dclr, centroid, collider, tool_len,
                                    gun_radius, clear_floor, neutral, max_step_deg, md)
            cands.append((mx.get("body_min", -1e9), Gx, mx, md))
        cands.sort(key=lambda z: -z[0])
        _, G, meta, md_used = cands[0]
        meta["mode"] = md_used
        pos[idx] = P
        gun[idx] = G

        body = np.array([collider.seg_clear(P[k], P[k] - G[k] * tool_len,
                                            gun_radius, SEARCH_R) for k in range(len(P))])
        gstep = _angle_deg(G[1:], G[:-1]) if len(P) > 1 else np.zeros(0)
        shift = np.linalg.norm(P - P0, axis=1)
        L1 = float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())
        info.append({
            "seg": [int(i0), int(i1)],
            "pts": int(len(P)),
            "len_mm": L1,
            "len_before_mm": L0,
            "raise_max_mm": float(raised),
            "shift_max_mm": float(shift.max()) if len(shift) else 0.0,
            "ang_endpoints": float(np.degrees(np.arccos(np.clip(np.dot(g0, g1), -1, 1)))),
            "theta_mid": float(meta.get("theta", 0.0)),
            "gun_max_step": float(gstep.max()) if len(gstep) else 0.0,
            "clear_min": float(dclr.min()),
            "clear_max": float(dclr.max()),
            "body_min": float(body.min()),
            "body_hold": int(meta.get("n_fix", 0)),
            "best_neutral_body": float(meta.get("body_min", 0.0)),
            "el_min": float(min(_sph_deg(G[k])[0] for k in range(len(P)))),
            "el_max": float(max(_sph_deg(G[k])[0] for k in range(len(P)))),
            "pole": float(meta.get("pole", 0.0)),
            "mode": str(meta.get("mode", "")),
            "infeasible": bool(infeasible),
        })

    out = traj_mod.Traj(pos=pos, gun=gun, speed=traj.speed.copy(), spray=spray.copy(),
                        extra={k: v.copy() for k, v in traj.extra.items()},
                        header=list(traj.header))
    return out, {"segments": info, "n_segments": len(info)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="空移段姿态重规划")
    ap.add_argument("csv")
    ap.add_argument("cloud")
    ap.add_argument("out")
    ap.add_argument("--standoff", type=float, default=100.0)
    ap.add_argument("--tool-len", type=float, default=453.08, help="法兰→喷嘴长度；DensoArmKine 实测 453.08，本仓 URDF 默认 150")
    ap.add_argument("--gun-radius", type=float, default=20.0)
    ap.add_argument("--floor", type=float, default=30.0, help="低于此净空的点不参与转身")
    ap.add_argument("--max-step", type=float, default=5.0,
                    help="空移段相邻点允许的最大转角（度/点），默认 5°≈250°/s@200mm/s")
    ap.add_argument("--clear-need", type=float, default=15.0,
                    help="空移段枪身胶囊要求的最小净空 mm（<=0 关闭路径抬升）")
    ap.add_argument("--raise-iters", type=int, default=4, help="姿态/抬升交替迭代次数")
    ap.add_argument("--neutral", default="centroid", choices=["centroid", "down"])
    ap.add_argument("--move-dir", default="z-", choices=["auto", "level", "z-", "z+"],
                    help="空移段枪轴策略：auto=每段在 level / z- 里选枪身净空更好的（默认）；"
                         "level=保持仰角只在水平面内转；z-=绕经 -Z 极点（浅→深）；z+=绕经 +Z")
    ap.add_argument("--json", default=None)
    args = ap.parse_args(argv)

    t0 = traj_mod.load_csv(args.csv)
    coll = CloudGrid(read_ply_xyz(args.cloud))

    before = stats(t0, "before")
    t1, plan_info = plan(t0, coll, args.standoff, args.tool_len, args.gun_radius,
                         args.floor, args.neutral, args.max_step,
                         args.move_dir, args.clear_need, args.raise_iters)
    after = stats(t1, "after")

    # 枪身扫掠净空（改前改后同口径，才能说"没变差"）。
    # 口径分两种，别混：
    #   inner = 纯空移点（前后都是 spray=0）—— 这些点的位置/姿态可以改；
    #   seam  = 接缝点（紧邻喷涂段）—— 姿态与位置都由喷涂段决定，不属于空移段的责任。
    mv = t1.spray == 0
    inner = np.zeros(len(t1), dtype=bool)
    if len(t1) > 2:
        inner[1:-1] = (t1.spray[1:-1] == 0) & (t1.spray[:-2] == 0) & (t1.spray[2:] == 0)
    seam = mv & ~inner
    body0 = np.full(len(t0), np.inf)
    body0[inner] = _gun_body_min(coll, t0.pos[inner], t0.gun[inner], args.tool_len,
                                 args.gun_radius)
    b0 = body0[np.isfinite(body0)]
    body = np.full(len(t1), np.inf)
    body[inner] = _gun_body_min(coll, t1.pos[inner], t1.gun[inner], args.tool_len,
                                args.gun_radius)
    body_ok = body[np.isfinite(body)]
    seam_body = _gun_body_min(coll, t1.pos[seam], t1.gun[seam], args.tool_len,
                              args.gun_radius)
    seam_ok = seam_body[np.isfinite(seam_body)]

    traj_mod.save_csv(args.out, t1)
    sp = t0.spray == 1
    unchanged = bool(np.array_equal(t0.pos[sp], t1.pos[sp]) and
                     np.array_equal(t0.gun[sp], t1.gun[sp]) and
                     np.array_equal(t0.speed[sp], t1.speed[sp]))
    rep = {
        "csv": args.csv, "cloud": args.cloud, "out": args.out,
        "params": {"standoff": args.standoff, "tool_len": args.tool_len,
                   "gun_radius": args.gun_radius, "floor": args.floor,
                   "neutral": args.neutral, "max_step": args.max_step,
                   "clear_need": args.clear_need, "raise_iters": args.raise_iters},
        "spray_segment_unchanged": unchanged,
        "before": before, "after": after, "plan": plan_info,
        "move_body_clear": {
            "min": float(body_ok.min()) if len(body_ok) else None,
            "p1": float(np.percentile(body_ok, 1)) if len(body_ok) else None,
            "p50": float(np.percentile(body_ok, 50)) if len(body_ok) else None,
        },
        "move_body_clear_before": {
            "min": float(b0.min()) if len(b0) else None,
            "p1": float(np.percentile(b0, 1)) if len(b0) else None,
            "p50": float(np.percentile(b0, 50)) if len(b0) else None,
        },
        "seam_body_clear": {
            "n": int(seam.sum()),
            "min": float(seam_ok.min()) if len(seam_ok) else None,
            "p50": float(np.percentile(seam_ok, 50)) if len(seam_ok) else None,
        },
    }
    txt = json.dumps(rep, ensure_ascii=False, indent=2)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            f.write(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
