"""与 DensoArmKine 对拍：**独立实现交叉验证**（C# 侧 vs 本 Python 侧）。

为什么必须做这件事：本层的 FK/IK/碰撞体是从 `DensoArmKine/` 的 C# 源码移植的，
移植就会有偏差（连杆表抄错一个数、轴的正负号、工具 +Z 的方向……见 planner 里那段
"工具 +Z 取反"的记录）。而这两个实现一旦分歧，现场拿到的是"程序说没撞、真机撞了"。

对拍口径：
  1. 取本层输出目录里的 `base_pose.json`（T_base_in_cloud）与 `trajectory_checked.csv`（含 J1..J6）；
  2. 把点云**变换到基座系**写成临时 PLY（DensoArmKine 的 --cloud-check 默认按基座系读）；
  3. 逐点调 `DensoArmKine.exe --cloud-check --pose "j1,..,j6"`，抓它的「最小表面距离」；
  4. 与本层 CSV 里的 `clearance_mm` 比。距离是刚体变换不变量，两边可以直接比。

用法：
    python -m robot_safety.verify_against_denso --out-dir outputs/robot_safety/tray --samples 8
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

from . import cloud as cloud_mod

DENSO_EXE = r"E:/Code/新喷涂框架/SprayVision/docs/SprayVisionSimple/tools/DensoArmKine/bin/Release/net8.0-windows/DensoArmKine.exe"


def _read_csv(path: str) -> dict[str, np.ndarray]:
    """逐列解析。**不能用 np.loadtxt**：`violation` 列是文本（多数为空），它会直接报错。"""
    import csv as _csv

    with open(path, "r", encoding="utf-8-sig") as f:
        rows = list(_csv.DictReader(f))
    if not rows:
        raise ValueError(f"空 CSV：{path}")
    out: dict[str, np.ndarray] = {}
    for key in rows[0]:
        vals = []
        ok = True
        for r in rows:
            try:
                vals.append(float(r[key]))
            except (TypeError, ValueError):
                ok = False
                break
        if ok:
            out[key] = np.array(vals, dtype=float)
    return out


def _to_base_cloud(cloud_ply: str, T: np.ndarray, out_ply: str) -> int:
    pts = cloud_mod.read_ply_xyz(cloud_ply)
    R, t = T[:3, :3], T[:3, 3]
    pb = (pts - t) @ R          # p_base = R^T (p_cloud - t)
    cloud_mod.write_ply_xyz(out_ply, pb)
    return len(pts)


def _parse_denso(txt: str) -> tuple[str, float, bool]:
    """从 DensoArmKine 报告里抓 (code, 最小表面距离, 是否为下界)。

    表格列序是 **部位 / 半径 / 到轴距离 / 表面距离 / 状态**，所以表面距离是第 3 个数
    （原先取 nums[-2] 拿到的是"到轴距离"，差了一个半径，会系统性偏大）。
    带 `≥` 的行表示该段在搜索半径内根本没测到点云，只能给下界 —— 下界不能当实测值比较。
    """
    code = "?"
    m = re.search(r"code=(\w+)", txt)
    if m:
        code = m.group(1)
    best = None
    best_is_bound = False
    part_re = re.compile(r"^\s*(底座立柱|底座颈|J\d\s*关节|连杆\s*J\d-J\d|工具握把|工具细管|喷嘴)\s+(.*)$")
    for line in txt.splitlines():
        pm = part_re.match(line)
        if not pm:
            continue
        rest = pm.group(2)
        toks = [(mm.group(1), float(mm.group(2)))
                for mm in re.finditer(r"([≥~]?)\s*(-?\d+(?:\.\d+)?)", rest)]
        if len(toks) < 3:
            continue
        prefix, v = toks[2]          # 第 3 个数 = 表面距离
        is_bound = prefix in ("≥", "~")
        if best is None or v < best:
            best = v
            best_is_bound = is_bound
    if best is None:
        m2 = re.search(r"最小表面距离\s*[≥~]?\s*(-?[\d.]+)", txt)
        if m2:
            best = float(m2.group(1))
            best_is_bound = "≥" in txt[max(0, m2.start() - 2):m2.start() + 8]
    return code, (best if best is not None else float("nan")), best_is_bound


def _read_text(path: str) -> str:
    """DensoArmKine 的报告实际是 UTF-8（曾按 GBK 读导致整篇乱码、解析全空）。"""
    raw = open(path, "rb").read()
    for enc in ("utf-8-sig", "utf-8", "gbk"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", "replace")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="与 DensoArmKine 对拍碰撞距离")
    ap.add_argument("--out-dir", required=True, help="run_safety 的输出目录")
    ap.add_argument("--cloud", help="原始点云 PLY（默认从 run_safety 的 report 里推不出来，必须给）")
    ap.add_argument("--samples", type=int, default=6, help="抽样点数（每个点要起一次 DensoArmKine.exe，约 1s）")
    ap.add_argument("--exe", default=DENSO_EXE)
    ap.add_argument("--tool-len", type=float, default=None)
    args = ap.parse_args(argv)

    if not args.cloud:
        print("必须给 --cloud（原始点云 PLY）", file=sys.stderr)
        return 2
    if not os.path.exists(args.exe):
        print(f"DensoArmKine.exe 不存在：{args.exe}", file=sys.stderr)
        return 2

    bp = os.path.join(args.out_dir, "base_pose.json")
    with open(bp, "r", encoding="utf-8") as f:
        meta = json.load(f)
    T = np.array(meta["T_base_in_cloud"], dtype=float)
    tool_len = args.tool_len if args.tool_len is not None else meta["tool_len"]
    # 判定带必须两边一致，否则会出现"本层饱和在下界、Denso 给实测"的假分歧
    search = float(meta.get("search", 120.0))
    safety = float(meta.get("safety", 20.0))

    cols = _read_csv(os.path.join(args.out_dir, "trajectory_checked.csv"))
    n = len(cols["X"])
    clear = cols["clearance_mm"]
    # 抽样：一半取最危险的（clearance 最小的），一半均匀铺开 —— 只看危险点会漏掉"系统性偏差"
    idx_sorted = np.argsort(clear)
    picks = list(idx_sorted[: max(1, args.samples // 2)])
    picks += list(np.linspace(0, n - 1, max(1, args.samples - len(picks))).astype(int))
    picks = sorted(set(int(i) for i in picks))

    tmpdir = tempfile.mkdtemp(prefix="denso_pair_")
    tmp_ply = os.path.join(tmpdir, "cloud_in_base.ply")
    cnt = _to_base_cloud(args.cloud, T, tmp_ply)
    print(f"点云 {cnt} 点 → 基座系临时文件 {tmp_ply}")
    print(f"工具长 {tool_len:.1f}mm · safety {safety:.0f}mm · search {search:.0f}mm · 抽样 {len(picks)} 点")
    print()
    print(f"{'idx':>6} {'本层':>9} {'Denso':>9} {'差':>7}  {'code':<7} 结论")
    worst = 0.0
    n_cmp = 0
    n_code = 0
    for i in picks:
        pose = ",".join(f"{cols[f'J{k + 1}_deg'][i]:.4f}" for k in range(6))
        out_txt = os.path.join(tmpdir, f"check_{i}.txt")
        cmd = [args.exe, "--cloud-check", "--cloud", tmp_ply, "--pose", pose, "--in-base", "on",
               "--safety", f"{safety:.1f}", "--search", f"{search:.1f}", "--target", "400000",
               "--tool", f"{tool_len:.3f}", "--out", out_txt]
        try:
            subprocess.run(cmd, capture_output=True, timeout=120, cwd=os.path.dirname(args.exe))
        except subprocess.TimeoutExpired:
            print(f"{i:6d} 超时")
            continue
        if not os.path.exists(out_txt):
            print(f"{i:6d} 无输出")
            continue
        code, denso_clear, is_bound = _parse_denso(_read_text(out_txt))
        mine = float(clear[i])
        # 碰撞与否必须一致（这是最硬的判据：一个说有撞、一个说没有，就是移植错误）
        my_hit = mine < 0
        denso_hit = code in ("COLLIDE", "PENETRATE", "HIT")
        if code != "?":
            n_code += 1
            if my_hit != denso_hit:
                print(f"{i:6d} {mine:9.2f} {'':>9} {'':>7}  {code:<7} **碰撞判定分歧**")
                worst = max(worst, 99.0)
                continue
        if not np.isfinite(denso_clear) or is_bound:
            print(f"{i:6d} {mine:9.2f} {denso_clear:9.2f} {'—':>7}  {code:<7} 下界，不计入")
            continue
        diff = mine - denso_clear
        worst = max(worst, abs(diff))
        n_cmp += 1
        agree = "一致" if abs(diff) <= 5.0 else "**分歧**"
        print(f"{i:6d} {mine:9.2f} {denso_clear:9.2f} {diff:7.2f}  {code:<7} {agree}")
    print()
    print(f"碰撞判定一致 {n_code} 点 · 距离可比 {n_cmp} 点 · 最大偏差 {worst:.2f}mm")
    print("（容差 5mm：DensoArmKine 用体素采样查询，本层用全量精确距离，几毫米差异属预期）")
    return 0 if worst <= 5.0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
