"""轨迹 CSV 读写：`X,Y,Z,I,J,K,speed,spray`（master 的 driver 输出格式，一个字符都不改）。

**语义（这条最容易搞反，搞反了整套碰撞判断全错）**：
    轨迹点 = 枪位（TCP 位置），**不是**工件表面点；
    `(I,J,K)` = **枪轴朝向**，从枪指向工件 → 也就是 `-n`（n 是表面外法向）；
    因此    工件表面点 = TCP + standoff * (I,J,K)        （因为 TCP = 面点 + standoff·n）
    枪身朝外方向（工具系 +Z）= -(I,J,K) = n。

校验办法：panel 的面朝上（点云 z≈1343~1374，枪位 z≈1458），若把 (I,J,K) 当 n 用，
会算出「枪在工件下方 100mm」——与 CSV 里的 Z 直接矛盾。这条口径在 cpp 侧是
`emitPt` 里的 `-nrm`（spray_traj.cpp 约 3277 行），本文件与它一一对应。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import geom


@dataclass
class Traj:
    pos: np.ndarray           # (N,3) TCP 位置 mm
    gun: np.ndarray           # (N,3) 枪轴朝向（指向工件），单位向量
    speed: np.ndarray         # (N,)
    spray: np.ndarray         # (N,) int 0/1
    extra: dict[str, np.ndarray] = field(default_factory=dict)   # 输入里的额外列（原样带出）
    header: list[str] = field(default_factory=lambda: ["X", "Y", "Z", "I", "J", "K", "speed", "spray"])

    def __len__(self) -> int:
        return len(self.pos)

    @property
    def normal(self) -> np.ndarray:
        """表面外法向 n（= -gun）。"""
        return -self.gun

    def face_points(self, standoff: float) -> np.ndarray:
        """反推工件表面点：face = TCP + standoff * gun。"""
        return self.pos + self.gun * standoff

    def rows(self) -> list[tuple[int, int]]:
        """把轨迹切成「行」：spray=1 的连续段。过渡段（spray=0）不属于任何行。"""
        out = []
        n = len(self)
        i = 0
        while i < n:
            if self.spray[i] == 0:
                i += 1
                continue
            j = i
            while j + 1 < n and self.spray[j + 1] == 1:
                j += 1
            out.append((i, j))
            i = j + 1
        return out


def load_csv(path: str) -> Traj:
    with open(path, "r", encoding="utf-8-sig") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    header = [c.strip() for c in lines[0].split(",")]
    lower = [h.lower() for h in header]
    need = ["x", "y", "z", "i", "j", "k"]
    for k in need:
        if k not in lower:
            raise ValueError(f"{path} 缺少列 {k.upper()}（表头：{header}）")
    ix = {k: lower.index(k) for k in need}
    i_speed = lower.index("speed") if "speed" in lower else None
    i_spray = lower.index("spray") if "spray" in lower else None

    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]], dtype=float)
    pos = data[:, [ix["x"], ix["y"], ix["z"]]]
    gun = data[:, [ix["i"], ix["j"], ix["k"]]]
    nrm = np.linalg.norm(gun, axis=1)
    bad = nrm < 1e-9
    if np.any(bad):
        raise ValueError(f"{path} 里有 {int(bad.sum())} 个点的 IJK 是零向量，无法定姿态")
    gun = gun / nrm[:, None]
    speed = data[:, i_speed] if i_speed is not None else np.full(len(data), 200.0)
    spray = data[:, i_spray].astype(int) if i_spray is not None else np.ones(len(data), dtype=int)
    extra = {}
    for j, h in enumerate(header):
        if lower[j] not in need and lower[j] not in ("speed", "spray"):
            extra[h] = data[:, j]
    return Traj(pos=pos, gun=gun, speed=speed, spray=spray, extra=extra, header=header)


def save_csv(path: str, traj: Traj, add_cols: dict[str, np.ndarray] | None = None) -> None:
    """**只加列不改列**：原始 8 列原样按 master 的格式写，新列追加在后面。

    为什么要死守这个顺序：`trajectory_checked.csv` 要能和原始 CSV 逐行对齐比对，
    列序一乱，任何 diff 工具给出的都是「整行都变了」。
    """
    add_cols = add_cols or {}
    cols = ["X", "Y", "Z", "I", "J", "K", "speed", "spray"]
    extra_names = [h for h in traj.header if h not in cols]
    all_names = cols + extra_names + list(add_cols.keys())
    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write(",".join(all_names) + "\n")
        for i in range(len(traj)):
            g = traj.gun[i]
            p = traj.pos[i]
            vals = [
                f"{p[0]:.3f}", f"{p[1]:.3f}", f"{p[2]:.3f}",
                f"{g[0]:.4f}", f"{g[1]:.4f}", f"{g[2]:.4f}",
                f"{traj.speed[i]:.1f}", f"{int(traj.spray[i])}",
            ]
            for h in extra_names:
                vals.append(f"{traj.extra[h][i]:.6g}")
            for v in add_cols.values():
                item = v[i]
                if isinstance(item, str):
                    # 违规标记列（如 "COLLISION|JOINT_LIMIT"）：原样写，不加引号
                    vals.append(item)
                elif isinstance(item, (int, np.integer)):
                    vals.append(str(int(item)))
                else:
                    vals.append(f"{float(item):.2f}")
            f.write(",".join(vals) + "\n")


def set_gun_from_normal(nrm: np.ndarray) -> np.ndarray:
    """由法向得到 CSV 的 (I,J,K)（= -n）。"""
    return -geom.unit(nrm)
