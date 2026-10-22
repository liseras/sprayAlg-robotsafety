"""机器人安全层 CLI。

用法（Windows / Git Bash 均可，路径用引号包住）：

    python -m robot_safety.run_safety ^
        --traj E:/tmp_build/final_tray.csv ^
        --cloud cpp/_regress/tray.ply ^
        --out-dir outputs/robot_safety/tray

    # 给真实基座位姿（点云系 -> 基座系）：x y z rx ry rz（度，ZYX 欧拉）
    python -m robot_safety.run_safety --traj ... --cloud ... --base-pose 0 0 0 0 0 0

    # 只检查、不倾斜（用于给「现状」出一份基线报告）
    python -m robot_safety.run_safety --traj ... --cloud ... --no-tilt

退出码：0 = PASS（无 P0），1 = 有 P0 违规，2 = 输入/运行错误。
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

from . import cloud as cloud_mod
from . import planner, report
from . import robot, traj as traj_mod


def _parse_base(args) -> tuple[np.ndarray | None, str]:
    if args.base_json:
        with open(args.base_json, "r", encoding="utf-8") as f:
            m = np.array(json.load(f), dtype=float)
        if m.shape != (4, 4):
            raise ValueError("--base-json 必须是 4x4 矩阵")
        return m, f"基座位姿来自 `{args.base_json}`（点云系→基座系的逆，直接给出 T_base_in_cloud）"
    if args.base_pose:
        v = [float(x) for x in args.base_pose]
        if len(v) != 6:
            raise ValueError("--base-pose 需要 6 个数：x y z rx ry rz（度）")
        return planner.base_from_euler(*v), f"基座位姿来自命令行 `--base-pose {args.base_pose}`"
    return None, ""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="喷涂轨迹机器人安全层：内壁倾斜角 + 碰撞 + 可达性")
    ap.add_argument("--traj", required=True, help="轨迹 CSV（master spray_driver 输出的格式）")
    ap.add_argument("--cloud", required=True, help="工件点云 PLY（原始扫描点云，越大越准）")
    ap.add_argument("--out-dir", default="robot_safety_out", help="输出目录")
    ap.add_argument("--standoff", type=float, default=100.0, help="喷距 mm（必须与生成轨迹时一致）")
    ap.add_argument("--safety", type=float, default=20.0, help="间隙告警阈值 mm（P1，软约束）")
    ap.add_argument("--hard-clearance", type=float, default=2.0,
                    help="绝对安全距离/回退硬下限 mm（默认 2；0=只要不穿透）")
    ap.add_argument("--search", type=float, default=120.0, help="碰撞判定带 mm（默认 120，与 DensoArmKine 一致）")
    ap.add_argument("--tilt-nonup", type=float, default=20.0, help="非朝上面的基准倾角（度）")
    ap.add_argument("--up-threshold", type=float, default=45.0, help="法向偏离竖直多少度算「非朝上」")
    ap.add_argument("--tilt-max", type=float, default=60.0, help="倾角上限（度）")
    ap.add_argument("--tilt-rate", type=float, default=3.0, help="倾角沿轨迹每点最大变化（度）")
    ap.add_argument("--max-spray-angle", type=float, default=30.0, help="工艺允许的最大喷角（度）")
    ap.add_argument("--tool-len", type=float, default=robot.VS087_TOOL_LEN, help="法兰→喷嘴长度 mm")
    ap.add_argument("--include-base", action="store_true", help="把底座立柱也算碰撞体（点云含地面时用）")
    ap.add_argument("--link-stride", type=int, default=3, help="连杆碰撞每 N 点查一次")
    ap.add_argument("--max-points", type=int, default=25000, help="点云判定用点数上限（体素抽稀）")
    ap.add_argument("--base-pose", nargs=6, metavar=("X", "Y", "Z", "RX", "RY", "RZ"),
                    help="基座在点云系中的位姿（mm / 度 ZYX）")
    ap.add_argument("--base-json", help="4x4 矩阵 JSON（T_base_in_cloud）")
    ap.add_argument("--no-tilt", action="store_true", help="完全不加倾斜角（出「现状」基线）")
    ap.add_argument("--no-fix-transit", action="store_true", help="过渡点碰撞不做抬升修正")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    t0 = time.time()
    try:
        tr = traj_mod.load_csv(args.traj)
    except Exception as ex:
        print(f"❌ 轨迹读取失败：{ex}", file=sys.stderr)
        return 2
    if not os.path.exists(args.cloud):
        print(f"❌ 点云不存在：{args.cloud}", file=sys.stderr)
        return 2
    try:
        co = cloud_mod.load_cloud(args.cloud, max_points=args.max_points)
    except Exception as ex:
        print(f"❌ 点云读取失败：{ex}", file=sys.stderr)
        return 2

    p = planner.Params(
        standoff=args.standoff,
        safety=args.safety,
        hard_clearance=args.hard_clearance,
        search=args.search,
        tilt_nonup_deg=0.0 if args.no_tilt else args.tilt_nonup,
        up_threshold_deg=args.up_threshold,
        tilt_max_deg=args.tilt_max,
        tilt_rate_deg_per_pt=args.tilt_rate,
        max_spray_angle_deg=args.max_spray_angle,
        tool_len=args.tool_len,
        include_base=args.include_base,
        link_stride=max(1, args.link_stride),
        fix_transit=not args.no_fix_transit,
        az_fallback=not args.no_tilt,
    )

    try:
        base, base_note = _parse_base(args)
    except Exception as ex:
        print(f"❌ 基座位姿解析失败：{ex}", file=sys.stderr)
        return 2
    if base is None:
        base = planner.auto_base(co, tr)
        base_note = (
            "⚠ **基座位姿为默认假设**（未提供 `--base-pose` / `--base-json`）：\n"
            f"工件质心被摆到基座系的 `({planner.AUTO_BASE_FWD:.2f}·reach, 0, "
            f"{planner.AUTO_BASE_DOWN:.2f}·reach)` ≈ "
            f"`({planner.AUTO_BASE_FWD * robot.VS087_REACH:.0f}, 0, "
            f"{planner.AUTO_BASE_DOWN * robot.VS087_REACH:.0f})` mm，即基座自工件外侧后退、"
            "并低于工件质心，+X 指向工件。\n"
            f"换算到点云系的基座原点假设值：`[{base[0, 3]:.1f}, {base[1, 3]:.1f}, "
            f"{base[2, 3]:.1f}]` mm。\n"
            "**坐标闭环（P0）尚未完成**：点云系 → 基座系必须来自手眼标定或工位图纸。\n"
            "本报告中的「可达率」与「连杆碰撞」结论依赖这个假设；枪体碰撞与喷距/喷角结论不依赖它。"
        )

    res = planner.plan(tr, co, base, p, verbose=not args.quiet)
    out_tr = planner.apply_result(tr, res)

    os.makedirs(args.out_dir, exist_ok=True)
    csv_path = os.path.join(args.out_dir, "trajectory_checked.csv")
    # 倾角跳变：平滑回退会让个别点退回本点最优倾角，从而偏离限幅。
    # 这一列就是给现场的"姿态会不会突然抖一下"的依据（只在行内相邻点之间量）。
    tilt_seq = [q.spray_angle_to_normal for q in res.points]
    jump = [0.0] * len(tilt_seq)
    for a, b in out_tr.rows():
        for i in range(a + 1, b + 1):
            jump[i] = abs(tilt_seq[i] - tilt_seq[i - 1])
    add = {
        "tilt_deg": [round(q.spray_angle_to_normal, 2) for q in res.points],
        "tilt_jump_deg": [round(v, 2) for v in jump],
        "clearance_mm": [round(q.clearance, 2) if np.isfinite(q.clearance) else -999.0 for q in res.points],
        "self_gap_mm": [round(q.self_collision, 2) if np.isfinite(q.self_collision) else 999.0 for q in res.points],
        "margin_deg": [round(q.margin_min, 2) for q in res.points],
        "lifted_mm": [round(q.lifted, 1) for q in res.points],
        "violation": ["|".join(q.violations) for q in res.points],
        # 关节角（度）也一并落盘：一来给现场对照，二来 verify_against_denso.py 要用它
        # 喂 DensoArmKine 的 --cloud-check 做独立实现对拍。
        **{f"J{k + 1}_deg": [round(float(q.theta[k]), 3) if np.all(np.isfinite(q.theta)) else 0.0
                             for q in res.points] for k in range(6)},
    }
    traj_mod.save_csv(csv_path, out_tr, add)

    verdict = report.write_all(args.out_dir, res, out_tr, co, args.traj, args.cloud, base_note)
    # 基座矩阵单独落盘：对拍脚本与「换个工位重算」都要用真值，不能只有一句提示语
    with open(os.path.join(args.out_dir, "base_pose.json"), "w", encoding="utf-8") as f:
        json.dump({"T_base_in_cloud": base.tolist(), "note": base_note, "tool_len": p.tool_len,
                   "standoff": p.standoff, "safety": p.safety, "search": p.search}, f,
                  ensure_ascii=False, indent=2)
    dt = time.time() - t0

    print()
    print("═" * 62)
    print(f"结论 {verdict['verdict']}   P0={sum(1 for v in res.violations if v['severity'] == 'P0')}  "
          f"P1={sum(1 for v in res.violations if v['severity'] == 'P1')}  "
          f"P2={sum(1 for v in res.violations if v['severity'] == 'P2')}")
    print(f"可达 {res.stats['reach_fraction'] * 100:.1f}% · 喷角 均值 {res.stats['tilt_mean_deg']:.1f}° "
          f"p95 {res.stats['tilt_p95_deg']:.1f}° 最大 {res.stats['tilt_max_deg']:.1f}° · "
          f"喷距偏差 {res.stats['standoff_max_err_mm']:.6f}mm")
    if res.stats.get("pose_reverted"):
        print(f"⚠ 平滑回退 {res.stats['pose_reverted']} 点（仅硬约束被破坏：碰撞/自碰撞/关节超限/间隙低于 "
              f"{p.hard_clearance:.1f}mm）")
    print(f"输出：{args.out_dir}\\trajectory_checked.csv / verdict.json / report.md / report.png")
    print(f"耗时 {dt:.1f}s")
    print("═" * 62)
    return 0 if verdict["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
