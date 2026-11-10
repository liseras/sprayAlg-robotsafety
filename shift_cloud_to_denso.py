"""把点云整体平移到 DENSO VS-087 的可达舒适区（只平移，不改朝向、不改内容）。

为什么是 (452.5, 0, 316.75)：
    与 robot_safety/planner.py::auto_base 同一口径——基座假设为"让工件质心落在基座系
    (0.5·reach, 0, 0.35·reach)"。VS-087 reach=905，故 0.5·905=452.5、0.35·905=316.75。
    这个点不是拍的：J2 轴在基座系 z=395，臂展 445+430，目标点落在 z≈317、水平 452 时
    J2/J3/J5 余量都大。

平移量 = 目标点 - 点云"中间位置"。中间位置默认取**包围盒中心**（几何中点，
不受点的疏密影响）；可选 --anchor centroid 改取面积/点密度加权的质心。

用法：
    python shift_cloud_to_denso.py --src DIR --dst DIR [--anchor bbox|centroid]
只写 dst，不动 src。PLY 保留原有格式（ascii/binary）、法向、颜色、其他自定义属性。
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np

REACH = 905.0
TARGET = np.array([0.5 * REACH, 0.0, 0.35 * REACH])   # (452.5, 0, 316.75)

_NP_TYPES = {
    "char": "i1", "uchar": "u1", "int8": "i1", "uint8": "u1",
    "short": "i2", "ushort": "u2", "int16": "i2", "uint16": "u2",
    "int": "i4", "uint": "u4", "int32": "i4", "uint32": "u4",
    "float": "f4", "float32": "f4", "double": "f8", "float64": "f8",
}


class Ply:
    """最小 PLY 读写：只关心 vertex 元素，其余元素原样搬运。"""

    def __init__(self, path: Path):
        self.path = path
        raw = path.read_bytes()
        end = raw.find(b"end_header")
        if end < 0:
            raise ValueError(f"{path.name}: 找不到 end_header")
        nl = raw.find(b"\n", end)
        self.header = raw[: nl + 1]
        self.body = raw[nl + 1:]
        lines = self.header.decode("ascii", "replace").splitlines()
        self.fmt = ""
        self.vertex_count = 0
        self.props: list[tuple[str, str]] = []   # (type, name)
        self.other_elements: list[tuple[str, int]] = []
        cur = None
        for ln in lines:
            t = ln.split()
            if not t:
                continue
            if t[0] == "format":
                self.fmt = t[1]
            elif t[0] == "element":
                cur = t[1]
                if cur == "vertex":
                    self.vertex_count = int(t[2])
                else:
                    self.other_elements.append((cur, int(t[2])))
            elif t[0] == "property" and cur == "vertex":
                if t[1] == "list":
                    raise ValueError(f"{path.name}: vertex 含 list 属性，不支持")
                self.props.append((t[1], t[2]))

    def vertex_dtype(self) -> np.dtype:
        return np.dtype([(n, _NP_TYPES[t]) for t, n in self.props])

    # ---------------------------------------------------------------- 读
    def read_xyz(self) -> np.ndarray:
        if self.fmt == "ascii":
            arr = np.loadtxt(self.path, skiprows=self._body_line_count(), max_rows=self.vertex_count, usecols=(0, 1, 2), dtype=np.float64, ndmin=2)
            return arr
        dt = self.vertex_dtype()
        need = self.vertex_count * dt.itemsize
        if len(self.body) < need:
            raise ValueError(f"{self.path.name}: 顶点数据长度不足")
        rec = np.frombuffer(self.body[:need], dtype=dt)
        return np.stack([rec["x"], rec["y"], rec["z"]], axis=1).astype(np.float64)

    def _body_line_count(self) -> int:
        return self.header.count(b"\n")

    # ---------------------------------------------------------------- 写
    def write_shifted(self, dst: Path, delta: np.ndarray) -> None:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if self.fmt == "ascii":
            self._write_ascii(dst, delta)
        else:
            self._write_binary(dst, delta)

    def _write_ascii(self, dst: Path, delta: np.ndarray) -> None:
        text = self.body.decode("ascii", "replace")
        lines = text.splitlines()
        out = []
        for i, ln in enumerate(lines[: self.vertex_count]):
            t = ln.split()
            if len(t) < 3:
                out.append(ln)
                continue
            x = float(t[0]) + delta[0]
            y = float(t[1]) + delta[1]
            z = float(t[2]) + delta[2]
            out.append(" ".join([f"{x:.6f}", f"{y:.6f}", f"{z:.6f}"] + t[3:]))
        out.extend(lines[self.vertex_count:])   # 其他元素（若有）原样
        body = ("\n".join(out) + "\n").encode("ascii", "replace")
        dst.write_bytes(self.header + body)

    def _write_binary(self, dst: Path, delta: np.ndarray) -> None:
        dt = self.vertex_dtype()
        need = self.vertex_count * dt.itemsize
        rec = np.frombuffer(self.body[:need], dtype=dt).copy()
        for k, d in zip(("x", "y", "z"), delta):
            rec[k] = (rec[k].astype(np.float64) + d).astype(rec[k].dtype)
        rest = self.body[need:]
        dst.write_bytes(self.header + rec.tobytes() + rest)


def anchor_of(xyz: np.ndarray, mode: str) -> np.ndarray:
    if mode == "centroid":
        return xyz.mean(axis=0)
    lo, hi = xyz.min(axis=0), xyz.max(axis=0)
    return (lo + hi) / 2.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--dst", required=True)
    ap.add_argument("--anchor", choices=["bbox", "centroid"], default="bbox")
    ap.add_argument("--target", nargs=3, type=float, default=TARGET.tolist())
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    src, dst = Path(args.src), Path(args.dst)
    target = np.asarray(args.target, dtype=float)
    files = sorted(src.glob("*.ply"))
    if not files:
        print(f"[ERR] {src} 下没有 ply")
        return 1
    if not args.dry_run:
        dst.mkdir(parents=True, exist_ok=True)

    print(f"目标点 T = ({target[0]:.2f}, {target[1]:.2f}, {target[2]:.2f}) mm   anchor={args.anchor}")
    print(f"{'文件':<44}{'点数':>9}{'包围盒尺寸(mm)':>26}{'平移量 dx,dy,dz':>26}")
    rows = []
    for f in files:
        try:
            ply = Ply(f)
            xyz = ply.read_xyz()
        except Exception as e:      # noqa: BLE001
            print(f"[SKIP] {f.name}: {e}")
            continue
        a = anchor_of(xyz, args.anchor)
        delta = target - a
        lo, hi = xyz.min(axis=0), xyz.max(axis=0)
        size = hi - lo
        print(f"{f.name:<44}{len(xyz):>9}{np.array2string(size, precision=1, suppress_small=True):>26}"
              f"{np.array2string(delta, precision=1, suppress_small=True):>26}")
        if not args.dry_run:
            ply.write_shifted(dst / f.name, delta)
        rows.append((f.name, len(xyz), lo + delta, hi + delta))
    if args.dry_run:
        return 0

    print("\n写出后各文件包围盒（均已移到目标点附近）：")
    for name, n, lo, hi in rows:
        c = (lo + hi) / 2
        far = max(  # 距基座原点最远的点（按包围盒角点保守估计）
            float(np.linalg.norm(np.array([x, y, z])))
            for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])
        )
        print(f"  {name:<44} 中心=({c[0]:7.1f},{c[1]:7.1f},{c[2]:7.1f})  最远角点距基座原点={far:6.1f}mm"
              f"  {'OK' if far <= REACH + 60 else '超臂展'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
