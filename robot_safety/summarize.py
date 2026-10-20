"""把一批 `run_safety` 输出目录汇总成一张表（`SUMMARY.md`）。

用途：回归跑完 10 个样本后，一眼看出「哪个样本被判 REJECT、为什么」，
而不是逐个打开 `report.md`。表里每个数字都能在该样本目录的 `verdict.json` 里查到。

用法：
    python -m robot_safety.summarize --root outputs/robot_safety
    python -m robot_safety.summarize --root outputs/robot_safety --md outputs/robot_safety/SUMMARY.md
"""

from __future__ import annotations

import argparse
import json
import os


def _load(d: str) -> dict | None:
    p = os.path.join(d, "verdict.json")
    if not os.path.isfile(p):
        return None
    try:
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def collect(root: str) -> list[tuple[str, dict]]:
    out = []
    for name in sorted(os.listdir(root)):
        d = os.path.join(root, name)
        if not os.path.isdir(d):
            continue
        v = _load(d)
        if v:
            out.append((name, v))
    # 结论从重到轻：REJECT → PASS_WITH_WARNINGS → PASS
    order = {"REJECT": 0, "PASS_WITH_WARNINGS": 1, "PASS": 2}
    out.sort(key=lambda kv: (order.get(kv[1].get("verdict", "PASS"), 3), kv[0]))
    return out


def render(root: str, items: list[tuple[str, dict]]) -> str:
    L: list[str] = []
    A = L.append
    A("# 机器人安全层 —— 全样本汇总")
    A("")
    A("> 由 `python -m robot_safety.summarize` 生成。逐样本详情看各自的 `report.md`。")
    A("")
    A("| 样本 | 结论 | 点数 | 可达 | 最小间隙 | 倾斜点 | 喷角 p95/最大 | 喷距偏差 | 平滑回退 | P0/P1/P2 |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for name, v in items:
        s = v.get("stats", {})
        c = v.get("counts") or {}
        p0 = c.get("P0", sum(1 for x in v.get("violations", []) if x["severity"] == "P0"))
        p1 = c.get("P1", sum(1 for x in v.get("violations", []) if x["severity"] == "P1"))
        p2 = c.get("P2", sum(1 for x in v.get("violations", []) if x["severity"] == "P2"))
        mc = s.get("min_clearance_mm")
        mc_txt = "—" if mc is None else f"{mc:.2f}"
        A(f"| `{name}` | **{v.get('verdict', '?')}** | {s.get('points', 0)} | "
          f"{s.get('reach_fraction', 0) * 100:.1f}% | {mc_txt} | {s.get('tilted_points', 0)} | "
          f"{s.get('tilt_p95_deg', 0):.1f}° / {s.get('tilt_max_deg', 0):.1f}° | "
          f"{s.get('standoff_max_err_mm', 0):.3f}mm | {s.get('pose_reverted', 0)} | {p0}/{p1}/{p2} |")
    A("")
    n_rej = sum(1 for _, v in items if v.get("verdict") == "REJECT")
    n_warn = sum(1 for _, v in items if v.get("verdict") == "PASS_WITH_WARNINGS")
    A(f"合计 {len(items)} 个样本：REJECT {n_rej} · PASS_WITH_WARNINGS {n_warn} · "
      f"PASS {len(items) - n_rej - n_warn}。")
    A("")
    bad = [(n, v) for n, v in items if v.get("verdict") == "REJECT"]
    if bad:
        A("## REJECT 的原因")
        A("")
        for name, v in bad:
            by: dict[str, int] = {}
            for x in v.get("violations", []):
                if x["severity"] == "P0":
                    by[x["rule"]] = by.get(x["rule"], 0) + 1
            A(f"- `{name}`：{by}")
        A("")
    return "\n".join(L)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="汇总 robot_safety 输出")
    ap.add_argument("--root", default="outputs/robot_safety")
    ap.add_argument("--md", default=None, help="写出 markdown 的路径（默认 <root>/SUMMARY.md）")
    args = ap.parse_args(argv)

    items = collect(args.root)
    if not items:
        print(f"没有找到任何 verdict.json：{args.root}")
        return 1
    md = render(args.root, items)
    out = args.md or os.path.join(args.root, "SUMMARY.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write(md)
    print(md)
    print(f"\n已写入 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
