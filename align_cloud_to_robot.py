# -*- coding: utf-8 -*-
"""把点云摆成「正面朝向机器人」并移到可达区（旋转 + 平移）。

背景
----
上一版只做了平移，工件保持采集时的朝向 —— 相机水平向前拍，点云里那根
「深度轴」到了机器人坐标里就变成了「竖直」，结果工件平躺、主面朝天，
机器人只能从下面/上面去喷，腕心被顶到臂展边缘。

本脚本把「深度轴」转成机器人的前向 +X，让要喷的那一面正对机器人。

两种源坐标系（自动判断，也可 --mode 强制）
------------------------------------------
camera（现场激光/相机采集，z=前向深度，通常 z 中心 > 800mm）
    实测样例：大平面 481×1000×82（z 只有 82 ⇒ 板正对相机）
    旋转后：主面法向 −z_cam → −X_base，即正对机器人
model（建模坐标，z=高度，通常 z 中心小）
    实测样例：开口盒子 400×300×300（z 上是开口面）
    旋转后：顶面 +z_cam → −X_base，即开口正对机器人

用法
----
    python align_cloud_to_robot.py --src <目录> --dst <目录> --dry-run
    python align_cloud_to_robot.py --src <目录> --dst <目录>
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from shift_cloud_to_denso import Ply  # 复用 PLY 头解析与二进制字段逻辑

# 基座系里的“舒适位置”，沿用 planner.py::auto_base：(0.5·reach, 0, 0.35·reach)
TARGET_DEFAULT = np.array([452.5, 0.0, 316.75])

# p_base = R @ p_src
R_CAMERA = np.array([[0.0, 0.0, 1.0],      # X_base = +z_cam（前向深度 → 前方）
                     [-1.0, 0.0, 0.0],     # Y_base = -x_cam
                     [0.0, -1.0, 0.0]])    # Z_base = -y_cam（图像纵向 → 竖直）

R_MODEL = np.array([[0.0, 0.0, -1.0],      # X_base = -z_cam（顶面朝上 → 朝机器人）
                    [0.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0]])


def pick_mode(xyz: np.ndarray) -> str:
    """z 中心高 ⇒ 相机系（深度）；否则建模系（高度）"""
    zc = 0.5 * (xyz[:, 2].min() + xyz[:, 2].max())
    return "camera" if zc > 800.0 else "model"


def rot_matrix(mode: str, flip_v: bool = False) -> np.ndarray:
    R = R_CAMERA if mode == "camera" else R_MODEL
    if flip_v:
        # 绕 X 轴翻 180°：只影响“上下”，不影响朝向
        Rx = np.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
        return Rx @ R
    return R


def load_handeye(spec: str) -> tuple[np.ndarray, np.ndarray]:
    """读手眼标定矩阵 p_base = R·p_cam + t。

    spec 可以是：
      · 一个文本文件（4 行 4 列，或 12/16 个数连续排列，逗号/空格/分号都行）；
      · 直接内联的数字串，例如 "R 的 9 个数 + t 的 3 个数"。

    相机系点云走到基座系全靠它，所以这里把校验做实：
      · R 必须接近正交（|RRᵀ−I| 小）且 det≈+1，否则说明抄错了或行列颠倒；
      · 第 4 行若存在必须是 (0,0,0,1)。
    只告警不抛异常 —— 现场标定本身就有几度误差（这份文件自称手工滑条对齐），
    硬拦下来反而没法用；但偏差大到 1e-2 以上基本就是抄错，要让人看见。
    """
    p = Path(spec)
    if p.exists():
        raw = p.read_text(encoding="utf-8", errors="replace")
        # 注释行必须剔掉：文件头里常写"860 高 × 399 宽"这类说明，正则会把数字一起抓进来
        raw = "\n".join(ln.split("#", 1)[0] for ln in raw.splitlines())
    else:
        raw = spec
    nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", raw)]
    if len(nums) < 12:
        raise SystemExit(f"手眼矩阵需要至少 12 个数（3x4），只解析到 {len(nums)} 个：{spec}")
    if len(nums) >= 16:
        M = np.asarray(nums[:16], dtype=float).reshape(4, 4)
        if np.abs(M[3] - np.array([0.0, 0.0, 0.0, 1.0])).max() > 1e-6:
            print(f"⚠ 第4行不是 (0,0,0,1)：{M[3]}")
    else:
        M = np.eye(4)
        M[:3, :] = np.asarray(nums[:12], dtype=float).reshape(3, 4)

    R, t = M[:3, :3].copy(), M[:3, 3].copy()
    ortho_err = np.abs(R @ R.T - np.eye(3)).max()
    det = np.linalg.det(R)
    print(f"手眼校验：|RRᵀ−I|max={ortho_err:.4g}  det={det:.6f}")
    if ortho_err > 1e-2:
        print("  ⚠ R 偏离正交超过 1e-2：可能把矩阵抄错/行列颠倒（标定误差通常 <1e-3）")
    if abs(det - 1.0) > 1e-2:
        print(f"  ⚠ det={det:.4f} 不等于 +1：这不是纯旋转（含缩放或镜像）")
    return R, t


def transform_ply(src: Path, dst: Path, R: np.ndarray, delta: np.ndarray) -> int:
    """写一个新 PLY：xyz 与法向一起旋转平移，其余属性原样保留。"""
    p = Ply(src)
    names = [n for _, n in p.props]
    has_n = all(k in names for k in ("nx", "ny", "nz"))
    dst.parent.mkdir(parents=True, exist_ok=True)

    if p.fmt == "ascii":
        text = p.body.decode("ascii", "replace")
        lines = text.splitlines()
        ix, iy, iz = names.index("x"), names.index("y"), names.index("z")
        inn = [names.index(k) for k in ("nx", "ny", "nz")] if has_n else None
        out = []
        for i, ln in enumerate(lines[: p.vertex_count]):
            t = ln.split()
            if len(t) < 3:
                out.append(ln)
                continue
            v = np.array([float(t[ix]), float(t[iy]), float(t[iz])])
            v = R @ v + delta
            t[ix], t[iy], t[iz] = f"{v[0]:.6f}", f"{v[1]:.6f}", f"{v[2]:.6f}"
            if inn is not None and max(inn) < len(t):
                n = np.array([float(t[k]) for k in inn])
                n = R @ n
                for k, val in zip(inn, n):
                    t[k] = f"{val:.6f}"
            out.append(" ".join(t))
        out.extend(lines[p.vertex_count:])
        body = ("\n".join(out) + "\n").encode("ascii", "replace")
        dst.write_bytes(p.header + body)

    else:
        dt = p.vertex_dtype()
        need = p.vertex_count * dt.itemsize
        rec = np.frombuffer(p.body[:need], dtype=dt).copy()
        xyz = np.column_stack([rec["x"], rec["y"], rec["z"]]).astype(np.float64)
        xyz = xyz @ R.T + delta
        for k, v in zip(("x", "y", "z"), xyz.T):
            rec[k] = v.astype(rec[k].dtype)
        if has_n:
            nrm = np.column_stack([rec["nx"], rec["ny"], rec["nz"]]).astype(np.float64)
            nrm = nrm @ R.T
            for k, v in zip(("nx", "ny", "nz"), nrm.T):
                rec[k] = v.astype(rec[k].dtype)
        rest = p.body[need:]
        dst.write_bytes(p.header + rec.tobytes() + rest)
    return p.vertex_count


def main() -> int:
    ap = argparse.ArgumentParser(description="点云旋转成正面朝机器人 + 平移到可达区")
    ap.add_argument("--src", required=True, help="源目录（*.ply）")
    ap.add_argument("--dst", required=True, help="输出目录")
    ap.add_argument("--mode", choices=["auto", "camera", "model"], default="auto")
    ap.add_argument("--handeye", default=None,
                    help="用真实手眼标定矩阵（4x4 文本文件，或 12/16 个数内联）。"
                         "给了它就用 p_base = R·p_cam + t，不再猜 --mode")
    ap.add_argument("--face-x", type=float, default=950.0,
                    help="「面向机器人的那一面」到基座的距离（默认 950mm；"
                         "枪 453 时对应腕心≈400，在舒适区）")
    ap.add_argument("--floor-z", type=float, default=0.0,
                    help="工件最低点的高度（默认 0 = 落地）")
    ap.add_argument("--old-target", nargs=3, type=float, default=None,
                    help="改用旧口径：直接指定包围盒中心")
    ap.add_argument("--flip-vertical", action="store_true", help="上下翻转（若看到工件倒立再开）")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    src, dst = Path(args.src), Path(args.dst)
    files = sorted(src.glob("*.ply"))
    if not files:
        print(f"源目录没有 ply: {src}")
        return 1

    # ---- 手眼矩阵：给了就照它变换，之后只做「移到可达位」的平移
    he_R = he_t = None
    if args.handeye:
        he_R, he_t = load_handeye(args.handeye)
        print(f"手眼矩阵 p_base = R·p_cam + t")
        print("  R =")
        for row in he_R:
            print("      " + "  ".join(f"{v:8.3f}" for v in row))
        print(f"  t = {np.round(he_t, 2)}")
        up = he_R[2, :]          # 基座 +Z(重力向上) 在点云系 = Rᵀ·ẑ = R 的第 3 行
        print(f"  重力向上(点云系) = R 第3行 = {np.round(up, 4)}"
              f"  ⇒ upCloud = ({up[0]:.4f}, {up[1]:.4f}, {up[2]:.4f})")
        print(f"  （安全层 -up 用上面这组数；也可以近似成最近的坐标轴）")
        print()

    print(f"{'文件':<32}{'模式':>7}{'原尺寸':>18}{'旋转后尺寸':>18}{'正面x':>8}{'底部z':>7}")
    for f in files:
        xyz = Ply(f).read_xyz()
        if he_R is not None:
            # 手眼变换后已经在基座系里，只差「挪到可达位」
            mode = "handeye"
            R = he_R
            base_t = he_t.copy()
        else:
            mode = pick_mode(xyz) if args.mode == "auto" else args.mode
            R = rot_matrix(mode, args.flip_vertical)
            base_t = np.zeros(3)

        lo, hi = (xyz @ R.T + base_t).min(0), (xyz @ R.T + base_t).max(0)
        center = (lo + hi) / 2.0

        if args.old_target is not None:
            shift = np.asarray(args.old_target, dtype=float) - center
        else:
            # 正面（min_x 那一侧）对齐到 face-x，y 居中，最低点落到 floor-z
            shift = np.array([args.face_x - lo[0], -center[1], args.floor_z - lo[2]])

        delta = base_t + shift
        final_lo = lo + shift

        s0, s1 = xyz.max(0) - xyz.min(0), hi - lo
        print(f"{f.stem[:30]:<32}{mode:>7}"
              f"{'[%4.0f %4.0f %4.0f]' % tuple(s0):>18}"
              f"{'[%4.0f %4.0f %4.0f]' % tuple(s1):>18}"
              f"{final_lo[0]:>8.1f}{final_lo[2]:>7.1f}")

        if not args.dry_run:
            transform_ply(f, dst / f.name, R, delta)

    if args.dry_run:
        print("\n(dry-run：未写文件)")
    else:
        print(f"\n已写出 {len(files)} 个文件 → {dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
