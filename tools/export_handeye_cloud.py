# -*- coding: utf-8 -*-
"""把一份相机系点云按手眼标定搬到基座系，输出两版点云 + 校验报告。

为什么单独出这个脚本
--------------------
界面上的「勾手眼」**不会改点云**（ApplyHandEyeIfNeeded 只设 _handEyeBasePose /
_handEyeExportPose 和 up 方向），点云在内存里始终是相机系。
想拿到"经过手眼标定转换后的点云"这个实物文件，只能离线做一次变换。

两版输出
--------
1. <name>_handeye.ply           纯标定：p_base = R·p_cam + t，位置是标定真值
2. <name>_handeye_X450.ply      标定 + 可达位平移（与界面 facex=450 / floorz=-200 一致），
                                 可直接用「点云已在基座系（恒等）」模式喂回界面

⚠ 用途边界
-----------
基座系点云的几何轴是 X=前 / Z=上，而生成层内部约定**可喷侧固定 up=(0,0,-1)**
（那是相机系口径：+Z 是光轴指向景物）。把基座系点云直接喂回生成层，方向判定会全错。
所以这两份文件是给 **安全层 / 可视化 / DENSO 对接 / 人工核对位置** 用的，
不要当生成层的输入。

用法
----
    python tools/export_handeye_cloud.py \
        --src "D:/.../槽型盘悬挂_X+800.ply" \
        --handeye config/handeye_cam2base.txt \
        --outdir outputs \
        --face-x 450 --floor-z -200
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from align_cloud_to_robot import load_handeye, transform_ply  # noqa: E402
from shift_cloud_to_denso import Ply                              # noqa: E402

np.set_printoptions(suppress=True, precision=2)


def fmt(v: np.ndarray) -> str:
    return "[" + ", ".join(f"{x:9.2f}" for x in v) + "]"


def report(tag: str, xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lo, hi = xyz.min(0), xyz.max(0)
    print(f"  {tag}")
    print(f"    bbox min {fmt(lo)}")
    print(f"    bbox max {fmt(hi)}")
    print(f"    尺寸     {fmt(hi - lo)}")
    print(f"    中心     {fmt((lo + hi) / 2.0)}")
    print(f"    离基座   X 最近 {lo[0]:.1f}mm / 最远 {hi[0]:.1f}mm"
          f"  （VS-087 臂展 905）")
    return lo, hi


def main() -> int:
    ap = argparse.ArgumentParser(description="导出经过手眼标定转换的点云")
    ap.add_argument("--src", required=True, help="源 ply（相机系）")
    ap.add_argument("--handeye", default=str(ROOT / "config" / "handeye_cam2base.txt"))
    ap.add_argument("--outdir", default=str(ROOT / "outputs"))
    ap.add_argument("--face-x", type=float, default=450.0)
    ap.add_argument("--floor-z", type=float, default=-200.0)
    ap.add_argument("--no-align", action="store_true", help="只出纯标定版，不做可达位平移")
    args = ap.parse_args()

    src = Path(args.src)
    if not src.exists():
        raise SystemExit(f"源点云不存在：{src}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- 步骤 1：读标定矩阵并校验
    print("【1】手眼标定矩阵  p_base = R · p_cam + t")
    R, t = load_handeye(args.handeye)
    for row in R:
        print("      " + "  ".join(f"{v:8.3f}" for v in row))
    print(f"      t = {fmt(t)}")
    up_cloud = R[2, :]     # 基座 +Z（重力向上）在点云系里的方向 = Rᵀ·ẑ = R 第 3 行
    print(f"  重力向上(点云系) = R 第3行 = {fmt(up_cloud)}   ← 安全层 -up 用这组数")

    # ---- 步骤 2：源点云概况
    print(f"\n【2】源点云 {src.name}（相机系）")
    p = Ply(src)
    xyz = p.read_xyz()
    names = [n for _, n in p.props]
    print(f"    点数 {p.vertex_count}  格式 {p.fmt}  属性 {names}")
    has_n = all(k in names for k in ("nx", "ny", "nz"))
    print(f"    含法向 {'是（会一起旋转）' if has_n else '否（只有 xyz，法向由生成层现算）'}")
    src_lo, src_hi = report("相机系 bbox", xyz)

    # ---- 步骤 3：纯手眼变换 → 基座系真值
    print("\n【3】纯手眼变换（标定真值，不做任何摆位）")
    base_xyz = xyz @ R.T + t
    lo, hi = report("基座系 bbox", base_xyz)
    out1 = outdir / f"{src.stem}_handeye.ply"
    transform_ply(src, out1, R, t)
    print(f"    → {out1}")

    # ---- 步骤 4：可达位平移（与界面 facex/floorz 同口径）
    if not args.no_align:
        print(f"\n【4】可达位平移（face-x={args.face_x:g}, floor-z={args.floor_z:g}）")
        shift = np.array([args.face_x - lo[0],
                          -0.5 * (lo[1] + hi[1]),
                          args.floor_z - lo[2]])
        print(f"    平移 shift = {fmt(shift)}")
        out2 = outdir / f"{src.stem}_handeye_X{args.face_x:g}.ply"
        transform_ply(src, out2, R, t + shift)
        report("基座系 bbox（含可达位）", base_xyz + shift)
        print(f"    → {out2}")

    # ---- 步骤 5：一致性校验（拿输出文件重新读回来验，不信内存）
    print("\n【5】校验（读回文件复核）")
    back = Ply(out1).read_xyz()
    res = np.abs(back - base_xyz).max()
    print(f"    纯标定版 与 手工 R·p+t 的最大偏差 = {res:.6f} mm  {'✓' if res < 1e-3 else '✗'}")
    # 竖直方向：源点云 y 轴（重力向上）应落到基座 Z；X 为薄法向轴
    span_src = src_hi - src_lo
    span_base = hi - lo
    print(f"    源  尺寸(X, Y, Z) = {fmt(span_src)}")
    print(f"    基座尺寸(X, Y, Z) = {fmt(span_base)}   （X=前 / Y=左 / Z=上）")
    print(f"    薄法向轴：源最薄轴 argmin={int(np.argmin(span_src))} → "
          f"基座最薄轴 argmin={int(np.argmin(span_base))}（0=X，说明主面正对机器人）")
    print(f"    竖直跨度：源 Y {span_src[1]:.1f}mm → 基座 Z {span_base[2]:.1f}mm "
          f"（差 {abs(span_src[1] - span_base[2]):.2f}mm）")

    print("\n完成。注意：这两份点云是**基座系**口径，给安全层/可视化/对接用；"
          "生成层内部要求可喷侧 up=(0,0,-1)（相机系），不要直接喂回生成层。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
