"""把本仓的轨迹 CSV 转成 DensoArmKine 能直接跑的 spray_plan JSON。

两侧口径（都已在本仓/DensoArmKine 源码里核对过，不是猜的）：

  本仓 CSV：`X,Y,Z,I,J,K,speed,spray`
      · X,Y,Z = **TCP（喷枪位姿点）**，不是面点 —— 实测 `pos + standoff·(I,J,K)` 落在工件表面
        （开口盒子 19082 点：中位 0.83mm、p90 1.77mm、99.9% 在 5mm 内，见 outputs/denso_link/）
      · I,J,K = **工具 +Z = 喷射方向**（= 表面外法向取反），见 planner.py `z_tool = -d`
      · 与 planner.py 的 `face = pos + gun·standoff`、`TCP = face + d·standoff` 完全自洽
      · 轨迹坐标系 = 点云坐标系

  DensoArmKine（Trajectory/SprayPlan.cs）：
      · Points[i].Cartesian = [x,y,z,rx,ry,rz,fig]，**TCP 位姿**
      · 欧拉 R = Rz(rz)·Ry(ry)·Rx(rx)，单位**度**（Kinematics/EulerZyx.cs）
      · --frame tcp 表示按 TCP 求逆解；--traj-frame base 表示坐标已是基座系

  所以换算只有三步：
      1) TCP = CSV 的 (X,Y,Z) 原样（**不要再退一个 standoff** —— 退两次等于把枪整体拉近 100mm，
         症状是 J3/J5 大面积超限，看起来像"工位摆得太近"，实测盒子限位内只有 20.4%）
      2) 工具 +Z = (I,J,K)；工具 +X 取「上一点 X 投影到新切平面」保证相邻姿态最小旋转
         （与 planner.py 同一套做法：不这么做，枪轴接近竖直时参考向量会在 +X/+Z 之间跳，
          逆解会跳到另一个解支，表现出来就是 J6 逐点自旋）
      3) R → (rx,ry,rz)（Rz·Ry·Rx 分支），单位度，Fig 默认 0
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

# 与 robot_safety/run_safety.py 的 standoff 默认值一致
DEFAULT_STANDOFF = 100.0


def unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-12 else np.zeros(3)


def read_traj_csv(path: Path) -> dict:
    # utf-8-sig：界面导出的 CSV 带 BOM（Excel 友好），用 utf-8 读会让表头第一列
    # 变成 "\ufeffX"，然后报"缺少列 X"——踩过一次，别再踩。
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    header = [c.strip() for c in lines[0].split(",")]
    low = [h.lower() for h in header]
    need = {}
    for k in ("x", "y", "z", "i", "j", "k"):
        if k not in low:
            raise ValueError(f"{path.name} 缺少列 {k.upper()}（表头：{header}）")
        need[k] = low.index(k)
    sp = low.index("speed") if "speed" in low else None
    gw = low.index("spray") if "spray" in low else None
    pos, axis, speed, spray = [], [], [], []
    for ln in lines[1:]:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        t = s.split(",")
        if len(t) <= max(need.values()):
            continue
        pos.append([float(t[need["x"]]), float(t[need["y"]]), float(t[need["z"]])])
        axis.append([float(t[need["i"]]), float(t[need["j"]]), float(t[need["k"]])])
        speed.append(int(round(float(t[sp]))) if sp is not None else 0)
        spray.append(int(round(float(t[gw]))) if gw is not None else 1)
    if not pos:
        raise ValueError(f"{path.name} 没有数据行")
    return {"pos": np.asarray(pos), "axis": np.asarray(axis), "speed": speed, "spray": spray}


def euler_zyx_deg(R: np.ndarray) -> tuple[float, float, float]:
    """R = Rz(rz)·Ry(ry)·Rx(rx) 的分解，单位度。与 EulerZyx.ToAnglesDeg 同一分支。"""
    r20 = max(-1.0, min(1.0, R[2, 0]))
    ry = -math.asin(r20)
    cy = math.cos(ry)
    if abs(cy) < 1e-9:                      # 万向锁：令 rz=0 反解 rx
        rx = math.atan2(-R[1, 2], R[1, 1])
        return math.degrees(rx), math.degrees(ry), 0.0
    rz = math.atan2(R[1, 0], R[0, 0])
    rx = math.atan2(R[2, 1], R[2, 2])
    return math.degrees(rx), math.degrees(ry), math.degrees(rz)


def delta_from_cloud(cloud: Path, target=(452.5, 0.0, 316.75)) -> np.ndarray:
    """与 shift_cloud_to_denso.py 同一口径：目标点 - 原云包围盒中心。"""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from shift_cloud_to_denso import Ply  # 局部导入，避免无 cloud 时也要依赖

    xyz = Ply(cloud).read_xyz()
    c = (xyz.min(0) + xyz.max(0)) / 2.0
    return np.asarray(target, dtype=float) - c


def convert(csv_path: Path, standoff: float, shift: np.ndarray, fig: float,
            speed_percent: int | None, name: str, pos_is: str = "tcp") -> dict:
    t = read_traj_csv(csv_path)
    pos = t["pos"] + shift                     # 平移到目标坐标系
    ax = np.array([unit(a) for a in t["axis"]])

    prev_x = None
    points = []
    for i in range(len(pos)):
        z = ax[i]
        if prev_x is None:
            tmp = np.array([0.0, 0.0, 1.0]) if abs(z[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
            x0 = unit(np.cross(tmp, z))
        else:                                   # 上一点 X 投影到新切平面
            x0 = prev_x - z * float(np.dot(prev_x, z))
            # 退化阈值 1e-3（≈0.06°）：上一点 X 与新枪轴几乎平行时，投影方向由浮点噪声决定，
            # 不同实现会给出相反的滚转 → 相邻两点 J6 翻 180°。此时必须走确定性参考。
            x0 = unit(x0) if np.linalg.norm(x0) > 1e-3 else None
            if x0 is None:
                tmp = np.array([0.0, 0.0, 1.0]) if abs(z[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
                x0 = unit(np.cross(tmp, z))
        y0 = unit(np.cross(z, x0))
        x0 = unit(np.cross(y0, z))
        prev_x = x0

        R = np.column_stack([x0, y0, z])
        rx, ry, rz = euler_zyx_deg(R)
        # 归一化到 (-180,180]：只是同一个旋转的另一种写法（矩阵不变），
        # 但能让 JSON 可读、也让 DensoArmKine 的姿态翻转/显示功能不被上万度的数骗到。
        rx, ry, rz = ((a + 180.0) % 360.0 - 180.0 for a in (rx, ry, rz))
        tcp = pos[i] - z * standoff if pos_is == "surface" else pos[i]
        points.append({
            "Joints": None,
            "Cartesian": [round(float(tcp[0]), 6), round(float(tcp[1]), 6), round(float(tcp[2]), 6),
                          round(rx, 6), round(ry, 6), round(rz, 6), float(fig)],
            "SpeedPercent": int(speed_percent) if speed_percent is not None else 30,
            "GunOpen": bool(t["spray"][i]),
            "Label": f"v{t['speed'][i]}",
            "IsValid": True,
        })

    return {
        "Name": name,
        "DefaultSpeedPercent": int(speed_percent) if speed_percent is not None else 30,
        "Points": points,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="轨迹 CSV → DensoArmKine spray_plan JSON")
    ap.add_argument("--csv", required=True, help="轨迹 CSV（X,Y,Z,I,J,K,speed,spray）")
    ap.add_argument("--out", required=True, help="输出 JSON")
    ap.add_argument("--standoff", type=float, default=DEFAULT_STANDOFF, help="喷距 mm，必须与生成轨迹时一致")
    ap.add_argument("--shift", nargs=3, type=float, default=None,
                    help="轨迹平移量 dx dy dz（轨迹生成在原始点云上时用）")
    ap.add_argument("--cloud", default=None,
                    help="原始点云 PLY：自动算出与点云平移相同的 delta（与 shift_cloud_to_denso.py 一致）")
    ap.add_argument("--no-shift", action="store_true", help="轨迹已在目标坐标系（在平移后的点云上生成的）")
    ap.add_argument("--pos-is", choices=["tcp", "surface"], default="tcp",
                    help="CSV 的 X,Y,Z 是 TCP（默认，master 的 driver 口径）还是工件表面点")
    ap.add_argument("--fig", type=float, default=0.0, help="P 型第 7 元素 Fig")
    ap.add_argument("--speed-percent", type=int, default=None, help="统一速度百分比（默认 30）")
    ap.add_argument("--name", default=None)
    args = ap.parse_args()

    if args.no_shift:
        shift = np.zeros(3)
    elif args.shift:
        shift = np.asarray(args.shift, dtype=float)
    elif args.cloud:
        shift = delta_from_cloud(Path(args.cloud))
    else:
        shift = np.zeros(3)
        print("[提示] 未给 --shift/--cloud/--no-shift，按零平移处理", file=sys.stderr)

    name = args.name or Path(args.csv).stem
    plan = convert(Path(args.csv), args.standoff, shift, args.fig, args.speed_percent, name, args.pos_is)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8")

    c = np.array([p["Cartesian"] for p in plan["Points"]])
    eu = c[:, 3:6]
    print(f"写出 {out}")
    print(f"  点数 {len(c)}（喷枪开 {sum(1 for p in plan['Points'] if p['GunOpen'])}）")
    print(f"  平移 dx,dy,dz = {shift[0]:.2f}, {shift[1]:.2f}, {shift[2]:.2f} mm   standoff={args.standoff:.0f}")
    print(f"  TCP 包围盒 X[{c[:,0].min():.1f},{c[:,0].max():.1f}] "
          f"Y[{c[:,1].min():.1f},{c[:,1].max():.1f}] Z[{c[:,2].min():.1f},{c[:,2].max():.1f}]")
    print(f"  欧拉角范围 Rx[{eu[:,0].min():.0f},{eu[:,0].max():.0f}] "
          f"Ry[{eu[:,1].min():.0f},{eu[:,1].max():.0f}] Rz[{eu[:,2].min():.0f},{eu[:,2].max():.0f}] 度")
    return 0


if __name__ == "__main__":
    sys.exit(main())
