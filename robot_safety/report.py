"""报告输出：verdict.json / report.md / report.png（三视图）。

约定（沿用 ROBOT_SAFETY 交接说明 5.3）：
  1. `verdict.json` —— 机器可读：`{pass, violations:[{index, rule, value, limit, severity}], stats}`
  2. `trajectory_checked.csv` —— **只加列不改列**，逐点可回查
  3. `report.png` —— 违规点在三视图上的分布（点云灰、安全绿、告警黄、穿透红）

matplotlib 的中文字体在本机不保证存在，图内文字**一律用英文**（避免方框字），
中文结论放 report.md。
"""

from __future__ import annotations

import json
import math
import os

import numpy as np

from .planner import Params, Result
from . import cloud as cloud_mod
from . import traj as traj_mod


def _f(x, nd=2):
    if x is None:
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    return round(v, nd)


def write_verdict(path: str, res: Result, traj: traj_mod.Traj, base_note: str) -> dict:
    vios = []
    for v in res.violations:
        vios.append({
            "index": int(v["index"]),
            "rule": v["rule"],
            "value": _f(v["value"], 3),
            "limit": _f(v["limit"], 3),
            "severity": v["severity"],
            "note": v.get("note", ""),
        })
    p0 = [v for v in vios if v["severity"] == "P0"]
    p1 = [v for v in vios if v["severity"] == "P1"]
    # 三档结论：P0 = 机器人不能执行（碰撞/不可达/自碰撞）；P1 = 能执行但有代价
    # （间隙偏小、喷角超工艺窗口、贴限位）；只有 P2 才算干净通过。
    if p0:
        verdict = "REJECT"
    elif p1:
        verdict = "PASS_WITH_WARNINGS"
    else:
        verdict = "PASS"
    doc = {
        "pass": len(p0) == 0,
        "verdict": verdict,
        "counts": {"P0": len(p0), "P1": len(p1),
                   "P2": sum(1 for v in vios if v["severity"] == "P2")},
        "base_frame_note": base_note,
        "stats": {k: (_f(v, 4) if isinstance(v, float) else v) for k, v in res.stats.items()},
        "cloud": {k: (_f(v, 3) if isinstance(v, float) else v) for k, v in res.cloud_info.items()},
        "params": {k: v for k, v in res.params.__dict__.items()},
        "violations": vios,
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=2)
    return doc


def _cluster_runs(vios: list[dict], max_gap: int = 2) -> list[dict]:
    """把同一规则的**连续点**并成一个区段。

    为什么必须聚合：一条 5000 点的轨迹，喷角或间隙只要在一小片区域超标，
    逐点列表就会刷出几百行，反而把"哪一片有问题"这个真正有用的信息埋掉。
    区段化之后 751 行变成几行，每行给出「点号范围 / 点数 / 最差值」。
    最差值同时保留：判断严重程度要看极值，不能只看平均值。
    """
    by_rule: dict[str, list[dict]] = {}
    for v in vios:
        if v["index"] < 0:          # L0 这类无点号的规则单独处理
            continue
        by_rule.setdefault(v["rule"], []).append(v)
    runs: list[dict] = []
    for rule, vs in by_rule.items():
        vs = sorted(vs, key=lambda d: d["index"])
        cur = [vs[0]]
        for v in vs[1:]:
            if v["index"] - cur[-1]["index"] <= max_gap:
                cur.append(v)
            else:
                runs.append(_run_summary(rule, cur))
                cur = [v]
        runs.append(_run_summary(rule, cur))
    return sorted(runs, key=lambda r: ({"P0": 0, "P1": 1, "P2": 2}[r["severity"]], -r["count"]))


def _run_summary(rule: str, vs: list[dict]) -> dict:
    # 极值方向按规则声明：碰撞类"越小越糟"，奇异度"越大越糟"。
    higher_is_worse = rule in ("NEAR_SINGULAR", "TILT_JUMP")
    worst = None
    for v in vs:
        if v["value"] is None:
            continue
        if worst is None:
            worst = v
            continue
        a, b = float(v["value"]), float(worst["value"])
        if (a > b) if higher_is_worse else (a < b):
            worst = v
    worst = worst or vs[0]
    return {
        "rule": rule,
        "severity": vs[0]["severity"],
        "count": len(vs),
        "from": vs[0]["index"],
        "to": vs[-1]["index"],
        "worst_index": worst["index"],
        "worst_value": worst["value"],
        "limit": worst["limit"],
        "note": worst.get("note", ""),
    }


def write_report_md(path: str, res: Result, traj: traj_mod.Traj, src_traj: str, src_cloud: str,
                    base_note: str, verdict: dict) -> None:
    s = res.stats
    p = res.params
    lines: list[str] = []
    A = lines.append
    A("# 机器人安全评估报告")
    A("")
    A(f"- 轨迹：`{src_traj}`（{len(traj)} 点，其中喷涂 {s['spray_points']}）")
    A(f"- 点云：`{src_cloud}`（原始 {res.cloud_info['raw_points']} 点，判定用 {res.cloud_info['points']} 点，"
      f"采样间距 {res.cloud_info['spacing_mm']:.2f}mm）")
    A(f"- 机型：DENSO VS-087（reach 905mm，工具长 {p.tool_len:.0f}mm）")
    A(f"- 结论：**{verdict['verdict']}**"
      f"（P0 {sum(1 for v in res.violations if v['severity'] == 'P0')} 条 · "
      f"P1 {sum(1 for v in res.violations if v['severity'] == 'P1')} 条 · "
      f"P2 {sum(1 for v in res.violations if v['severity'] == 'P2')} 条）")
    A("")
    A("## 1. 关键数字")
    A("")
    A("| 指标 | 值 | 说明 |")
    A("|---|---|---|")
    A(f"| 可达率 | {s['reach_fraction'] * 100:.1f}% | IK 收敛且残差达标的点数占比 |")
    A(f"| 倾斜点数 | {s['tilted_points']} / {s['points']} | 喷角 > 0.5° 的点 |")
    A(f"| 喷角均值 / p95 / 最大 | {s['tilt_mean_deg']:.1f}° / {s['tilt_p95_deg']:.1f}° / {s['tilt_max_deg']:.1f}° | 喷角 = 倾角 = 枪轴与法向夹角 |")
    A(f"| **喷距保持最大偏差** | {s['standoff_max_err_mm']:.6f} mm | 倾斜后 \\|TCP-face\\| 与 standoff 之差的上限 |")
    A(f"| 最小表面距离 | {s['min_clearance_mm']:.2f} mm | 全部碰撞段（枪体 + 连杆）的最小值 |")
    A(f"| 过渡点抬升 | {s['transit_lifted']} 点 | spray=0 且原本贴太近，沿 +Z 抬开 |")
    A(f"| 平滑后重算点 | {s['recomputed_after_smoothing']} | 倾角被包络抬高、需要重新求解的点 |")
    A(f"| 平滑回退点 | {s.get('pose_reverted', 0)} | 平滑位姿低于绝对安全距离/碰撞/自碰撞/关节超限，才退回本点最优；间隙略小不回退 |")
    A("")
    by = s["violations_by_rule"]
    if by:
        A("### 违规统计")
        A("")
        A("| 规则 | 数量 | 含义 |")
        A("|---|---|---|")
        rules = {
            "COLLISION": ("P0", "连杆/枪体与点云穿透（表面距离 < 0）"),
            "CLEARANCE": ("P1", "表面距离小于安全距离"),
            "SELF_COLLISION": ("P0", "非相邻连杆之间干涉"),
            "JOINT_LIMIT": ("P1", "关节超出软限位或余量不足"),
            "IK_FAIL": ("P0", "逆解不收敛（不可达）"),
            "SPRAY_ANGLE": ("P1", "喷角超工艺上限（膜厚会变薄/不均）"),
            "NEAR_SINGULAR": ("P2", "接近奇异位形（关节速度会放大）"),
            "TILT_JUMP": ("P2", "倾角相对上一点跳变过大（平滑回退的代价，手腕会甩一下）"),
        }
        for k, cnt in sorted(by.items(), key=lambda kv: -kv[1]):
            sev, desc = rules.get(k, ("?", ""))
            A(f"| `{k}` | {cnt} | [{sev}] {desc} |")
        A("")

    # ---- 工艺窗口 vs 机器人需求：把"喷角超标"从散点数量变成一句能决策的话
    ang = np.array([q.spray_angle_to_normal for q in res.points if q.ik_ok], dtype=float)
    if ang.size:
        over = int(np.count_nonzero(ang > p.max_spray_angle_deg + 1e-9))
        A("### 工艺窗口 vs 机器人需求")
        A("")
        A(f"假设的工艺喷角上限 = **{p.max_spray_angle_deg:.0f}°**"
          "（⚠ 这个数不是文档给定的：`ROBOT_SAFETY_交接说明` 把「喷角上限」列为待你拍板的开放问题，"
          "当前取 30° 作为占位假设，改 `--max-spray-angle` 可覆盖）。")
        A("")
        A("| 喷角区间 | 点数 | 占比 |")
        A("|---|---|---|")
        for lo, hi in ((0, 5), (5, 15), (15, 30), (30, 45), (45, 61)):
            c = int(np.count_nonzero((ang >= lo) & (ang < hi)))
            A(f"| {lo}–{hi}° | {c} | {c / ang.size * 100:.1f}% |")
        A("")
        if over:
            A(f"⇒ **{over} 点（{over / ang.size * 100:.1f}%）的喷角超出上述假设上限**。"
              "这些点是「保安全间隙」与「保膜厚均匀」之间的真实冲突，不是算法缺陷：")
            A("")
            A("- 若要保安全：接受这些喷角，需工艺确认膜厚是否可接受（或加大喷距/换短枪/拆片覆盖）。")
            A("- 若要保喷角：把 `--safety` 调小看撞不撞（当前 "
              f"{p.safety:.0f}mm），或缩小 `--tilt-max` 让算法不再靠加大倾角换间隙。")
            A("")
        else:
            A(f"⇒ 全部点的喷角都在假设上限 {p.max_spray_angle_deg:.0f}° 以内，工艺与机器人不冲突。")
            A("")
    A("## 2. 方法与口径")
    A("")
    A("### 2.1 倾斜 = 绕表面点做球面旋转（喷距不变）")
    A("")
    A("```")
    A("d(θ) = Rot(axis, θ) · n            # axis 在 n 的切平面内、指向「上」")
    A("TCP  = face + standoff · d(θ)      # face 是原始面点，不随倾斜移动")
    A("⇒ |TCP - face| ≡ standoff           喷距严格不变，变的只有喷角 θ")
    A("```")
    A("")
    A("每点都实测 `standoff_err`，见第 1 节「喷距保持最大偏差」。"
      f"本轮最大偏差 {s['standoff_max_err_mm']:.6f}mm（浮点级，非近似）。")
    A("")
    A("### 2.2 姿态决策：面级常量 + 逐点最小修正")
    A("")
    A(f"1. **候选池**：基准倾角取 {p.tilt_nonup_deg:.0f}°（法向偏离竖直 >{p.up_threshold_deg:.0f}° 的侧壁/内壁）"
      "或 0°（朝上的面，工艺最优；朝下的悬垂面不往上扳，往上扳会顶穿工件），"
      f"按 {p.tilt_step_deg:.0f}° 步进扩到上限 {p.tilt_max_deg:.0f}°；上翻方向再左右取 0/±30°/±60°"
      "（内壁往上扳可能顶到对面的壁，得左右试）。")
    A("2. **候选比序**：P0 点数（穿透+不可达+自碰撞，带 2% 抽样 + 2 点容忍带）→ 工艺喷角 ≤ "
      f"{p.max_spray_angle_deg:.0f}° → 间隙不足点数 → P0 精确比 → 倾角更小 → 方位趋 0。"
      "所有候选必须在同一采样数下比较，否则会挑出「只在小样本上好看」的姿态。")
    A(f"3. **姿态是面级常量**：法向夹角 ≤2°、沿法向偏移 ≤2mm 的相邻行归成一个面组，"
      "一个面组只取一个 `(倾角, 方位)`，面内逐点不变（同一面喷涂过程中喷角恒定），"
      "只允许在换面时改变。")
    A("4. **回退只受硬约束限制**：倾角与间隙不是单调关系（凹槽里倾角越大间隙反而越小），"
      "整行统一倾斜方位同理。但只要面级姿态仍满足硬约束（无碰撞/自碰撞、间隙 ≥ "
      f"{p.hard_clearance:.1f}mm、关节不超限），就保留，哪怕间隙比本点最优小几毫米；"
      "只有硬约束被破坏时才按最小偏离（先 ±1 步、再 ±2/±3 步）逐点退回"
      f"（见 CSV 的 `tilt_jump_deg` 列）。本轮退回 {s.get('pose_reverted', 0)} 点。")
    A("5. **退回时不许翻方位**：方位是「从哪一侧进枪」，翻一次 TCP 就绕面点转 "
      "2·喷距·sin(Δ方位/2)（60° → 100mm 级），画出来是长刺、落盘是不可插补的突跳。"
      "所以退回只在「面级方位 + 该点最优倾角」里找，且必须真过得了硬约束；"
      "缓存原样退回也只认方位相同的。都不成立就保持面级姿态，违规照实报，"
      "交给站位/基座去解决 —— 不为一个点把整段姿态甩出去。")
    A("")
    A("### 2.3 碰撞判据（与 DensoArmKine 同源）")
    A("")
    A("| 项 | 口径 |")
    A("|---|---|")
    A("| 工件碰撞 | 连杆/枪体 = 胶囊（A-B 轴 + 半径 R）；表面距离 = 点到轴最短距离 − R，**<0 即穿透** |")
    A("| 碰撞体 | 底座(可选)/J1..J6 关节球/6 根连杆/工具握把/细管/喷嘴，半径来源 `ArmCollider.ArmSizing` |")
    A("| 自碰撞 | 非相邻连杆胶囊两两最小距离 − 半径和；排除同一关节与其两根连杆 |")
    A("| 可达 | 阻尼最小二乘 IK，位置残差 ≤ 0.5mm、姿态残差 ≤ 0.5°，否则判不可达 |")
    A("| 限位 | VS-087 软限位（J1 ±170°，J2 −100~135°，J3 −133~150°，J4 ±270°，J5 ±120°，J6 ±360°） |")
    A("| 奇异 | 雅可比条件数 > " + f"{p.cond_max:.0f}" + " 记 P2 告警 |")
    A("")
    A("## 3. 坐标系")
    A("")
    A(base_note)
    A("")
    if res.violations:
        A("## 4. 违规明细（按连续区段聚合）")
        A("")
        A("> 同一规则在相邻点上一旦成立就并成一个区段。逐点列表会把「哪一片有问题」埋掉，")
        A("> 完整逐点清单见 `verdict.json` 的 `violations`。")
        A("")
        A("| 规则 | 级别 | 点数 | 点号范围 | 最差点 | 最差值 | 限值 | 说明 |")
        A("|---|---|---|---|---|---|---|---|")
        for r in _cluster_runs(res.violations):
            wv = "—" if r["worst_value"] is None else f"{float(r['worst_value']):.2f}"
            lim = "—" if r["limit"] is None else f"{float(r['limit']):.2f}"
            span = f"{r['from']}" if r["from"] == r["to"] else f"{r['from']}–{r['to']}"
            A(f"| `{r['rule']}` | {r['severity']} | {r['count']} | {span} | {r['worst_index']} | "
              f"{wv} | {lim} | {r['note']} |")
        A("")
        l0 = [v for v in res.violations if v["index"] < 0]
        if l0:
            A("### 输入层（L0）问题")
            A("")
            for v in l0:
                A(f"- [{v['severity']}] `{v['rule']}`：{v.get('note', '')}")
            A("")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def write_png(path: str, res: Result, collider: cloud_mod.CloudCollider, traj: traj_mod.Traj) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pts = collider.pts[:: max(1, len(collider.pts) // 20000)].astype(float)
    pos = np.array([p.pos for p in res.points], dtype=float)
    clr = np.array([p.clearance for p in res.points], dtype=float)
    spray = np.array([traj.spray[p.idx] for p in res.points])

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), dpi=110)
    names = [("X-Z", 0, 2), ("Y-Z", 1, 2), ("X-Y", 0, 1)]
    # 安全=绿，告警=黄，穿透=红（与 CSV 的 clearance 列一一对应）
    colors = np.where(clr < 0, "#d62728", np.where(clr < res.params.safety, "#e8a33d", "#2e9e5b"))
    for ax, (title, ia, ib) in zip(axes, names):
        ax.scatter(pts[:, ia], pts[:, ib], s=0.6, c="#c8ccd2", linewidths=0, label="cloud")
        ax.scatter(pos[spray == 0, ia], pos[spray == 0, ib], s=3, c="#7f8c9b", linewidths=0, label="transit")
        ax.scatter(pos[spray == 1, ia], pos[spray == 1, ib], s=3.5, c=colors[spray == 1], linewidths=0, label="spray")
        ax.set_title(f"{title}  (green<{res.params.safety:.0f}mm? see colours)")
        ax.set_xlabel("XYZ"[ia] + " (mm)")
        ax.set_ylabel("XYZ"[ib] + " (mm)")
        ax.set_aspect("equal", adjustable="datalim")
        ax.grid(alpha=0.25, lw=0.4)
        ax.legend(markerscale=4, fontsize=7, loc="best")
    fig.suptitle(
        f"robot-safety  pass={'YES' if not any(v['severity'] == 'P0' for v in res.violations) else 'NO'}   "
        f"reach {res.stats['reach_fraction'] * 100:.0f}%   tilt mean {res.stats['tilt_mean_deg']:.1f}deg max {res.stats['tilt_max_deg']:.1f}deg   "
        f"min clearance {res.stats['min_clearance_mm']:.1f}mm   "
        f"colour: green=safe  amber=near  red=penetrate",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path)
    plt.close(fig)


def write_all(out_dir: str, res: Result, traj: traj_mod.Traj, collider: cloud_mod.CloudCollider,
              src_traj: str, src_cloud: str, base_note: str) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    vpath = os.path.join(out_dir, "verdict.json")
    verdict = write_verdict(vpath, res, traj, base_note)
    write_report_md(os.path.join(out_dir, "report.md"), res, traj, src_traj, src_cloud, base_note, verdict)
    try:
        write_png(os.path.join(out_dir, "report.png"), res, collider, traj)
    except Exception as ex:      # 画图失败不该让整个评估失败
        print(f"[warn] 出图失败：{ex}")
    return verdict
