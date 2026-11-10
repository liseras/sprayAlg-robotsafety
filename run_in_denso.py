# -*- coding: utf-8 -*-
"""一条命令把工件跑进 DensoArmKine：点云 → 轨迹 CSV → spray_plan JSON → 离线整跑。

链路（三段都是既有程序，这里只做编排，不重算任何几何）：

    spray_driver.exe  cloud.ply → traj.csv          X,Y,Z,I,J,K,speed,spray（X,Y,Z 是 TCP）
    traj_to_denso_plan.py  traj.csv → spray_plan.json   DensoArmKine 的 Points[].Cartesian
    DensoArmKine.exe --traj-run --traj-frame base --frame tcp → 逐点逆解统计

为什么要 --traj-run 而不是开界面：界面跑几千点要几十秒且中途看不出跳分支，
本命令落盘的统计（摆出姿态 / 限位内可执行 / 跳变最大的点号）才是"这条轨迹能不能走"的证据。

例：
    # 已平移过的点云（坐标就是基座系），轨迹直接可用
    python run_in_denso.py --cloud "C:/Users/DYF/Desktop/点云/喷涂点云_DENSO可达/面板.ply"

    # 原始点云（相机系）：轨迹要跟着点云一起平移，交给 --shift-from-cloud 自动算
    python run_in_denso.py --cloud "C:/Users/DYF/Desktop/点云/喷涂点云/面板.ply" --shift-from-cloud
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DRIVER = ROOT / "cpp/SprayTrajCore/build/bin/Release/spray_driver.exe"
DENSO_EXE = Path(r"E:/Code/新喷涂框架/SprayVision/docs/SprayVisionSimple/tools/DensoArmKine/bin/Release/net8.0-windows/DensoArmKine.exe")
OUT_DIR = ROOT / "outputs/denso_link"


def grab(txt: str, key: str) -> int | None:
    m = re.search(re.escape(key) + r"\s+(\d+)", txt)
    return int(m.group(1)) if m else None


def main() -> int:
    ap = argparse.ArgumentParser(description="点云 → DensoArmKine 整跑")
    ap.add_argument("--cloud", required=True)
    ap.add_argument("--name", default=None, help="输出文件前缀（默认取点云文件名）")
    ap.add_argument("--model", default="VS-087")
    ap.add_argument("--frame", default="tcp", choices=["tcp", "flange"])
    ap.add_argument("--shift", nargs=3, type=float, default=None, help="轨迹平移量 dx dy dz")
    ap.add_argument("--shift-from-cloud", action="store_true",
                    help="按 shift_cloud_to_denso.py 的同一口径自动算平移量（原始点云用）")
    ap.add_argument("--stages", default="all", choices=["all", "traj", "convert", "run"],
                    help="只跑某一段（调试用）")
    args = ap.parse_args()

    cloud = Path(args.cloud)
    if not cloud.exists():
        print(f"[错误] 点云不存在：{cloud}")
        return 2
    if not DENSO_EXE.exists():
        print(f"[错误] 找不到 DensoArmKine.exe：{DENSO_EXE}")
        return 2
    name = args.name or cloud.stem
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv, json_p, rep = OUT_DIR / f"{name}.csv", OUT_DIR / f"{name}_plan.json", OUT_DIR / f"trajrun_{name}.txt"

    if args.stages in ("all", "traj"):
        print(f"[1/3] 生成轨迹  {cloud.name}")
        r = subprocess.run([str(DRIVER), str(cloud), str(csv), str(OUT_DIR / f"{name}.traj.ply")],
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        if r.returncode != 0 or not csv.exists():
            print("      失败：\n" + (r.stdout or "")[-800:] + (r.stderr or "")[-400:])
            print("      提示：spray_driver.exe 只认 ascii PLY")
            return 3
        print(f"      → {csv}  ({sum(1 for _ in csv.open()) - 1} 点)")

    if args.stages in ("all", "convert"):
        print("[2/3] 转 DensoArmKine JSON")
        cmd = [sys.executable, str(ROOT / "traj_to_denso_plan.py"), "--csv", str(csv), "--out", str(json_p)]
        if args.shift:
            cmd += ["--shift", *[str(v) for v in args.shift]]
        elif args.shift_from_cloud:
            cmd += ["--cloud", str(cloud)]
        else:
            cmd += ["--no-shift"]
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
        if r.returncode != 0:
            print("      失败：\n" + (r.stdout or "")[-800:] + (r.stderr or "")[-400:])
            return 4
        print("      " + "\n      ".join(l for l in (r.stdout or "").strip().splitlines()[1:]))

    if args.stages in ("all", "run"):
        print("[3/3] DensoArmKine 离线整跑")
        subprocess.run([str(DENSO_EXE), "--traj-run", "--traj", str(json_p), "--model", args.model,
                        "--frame", args.frame, "--traj-frame", "base", "--out", str(rep)],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
        txt = rep.read_text(encoding="utf-8-sig", errors="replace") if rep.exists() else ""
        if not txt:
            print("      没拿到报告")
            return 5
        n = grab(txt, "点总数") or 0
        pose = grab(txt, "摆出姿态") or 0
        inlim = grab(txt, "限位内可执行") or 0
        jump = grab(txt, "疑似换分支")
        print(f"      点数 {n} · 摆出姿态 {pose} ({100.0*pose/max(n,1):.1f}%) · "
              f"限位内可执行 {inlim} ({100.0*inlim/max(n,1):.1f}%) · 疑似换分支 {jump}")
        print(f"      报告 {rep}")
        if inlim < n * 0.9:
            print("      ⚠ 限位内可执行率偏低：先看报告里「超软限位分布」集中在哪根轴，"
                  "集中在单轴 ⇒ 摆放/安装高度问题，散在多轴 ⇒ 分支选择问题")
    return 0


if __name__ == "__main__":
    sys.exit(main())
