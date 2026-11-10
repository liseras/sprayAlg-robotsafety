"""点云碰撞查询：PLY 读取 + KD-tree 胶囊（连杆）最近距离。

**口径（与 DensoArmKine/Cloud/CloudCheck.cs 完全一致）**：
  零件的每一段碰撞体是「以 A-B 为轴、半径 R 的胶囊」；
  表面距离 = (点云到该轴的最短距离) - R。
  * <0   → 穿透（点云表面已经落进胶囊里）
  * 0~safe → 告警带
  * ≥safe → 安全
这和「点到点云距离」是同一件事，区别只在于先量到**轴**再扣半径 —— 因为连杆是圆柱，
贴着壁面平行走过时，轴到壁面的距离恒等于 radius + 真实间隙，不会像「点到最近点」那样受采样疏密影响。

⚠ 已知局限（必须写进报告，不许当没事）：现场点云是**线扫描**，沿扫描线密、线与线之间疏，
真正的空洞处这段判定会**偏乐观**。报告里给出点云采样间距，用来看乐观程度。
"""

from __future__ import annotations

import struct
from dataclasses import dataclass

import numpy as np
try:                                   # scipy 只用于"抽样估点距"，不是判定链路的一部分
    from scipy.spatial import cKDTree
except ImportError:                    # 现场机/托管 python 常常没装 scipy
    cKDTree = None


def _nn_dists(pts: np.ndarray, query: np.ndarray) -> np.ndarray:
    """每个 query 点到 pts 的最近距离。有 scipy 走 KD 树，没有就分块暴力。"""
    if cKDTree is not None:
        tree = cKDTree(pts)
        d, _ = tree.query(query, k=1)
        return d
    out = np.empty(len(query), dtype=np.float64)
    q = np.asarray(query, dtype=np.float32)
    p = np.asarray(pts, dtype=np.float32)
    for i in range(0, len(q), 64):
        blk = q[i:i + 64]
        d2 = ((blk[:, None, :] - p[None, :, :]) ** 2).sum(axis=2)
        out[i:i + 64] = np.sqrt(d2.min(axis=1))
    return out

from . import geom


# ---------------------------------------------------------------- PLY

def read_ply_xyz(path: str) -> np.ndarray:
    """读 PLY 的 x,y,z（支持 ascii / binary_little_endian / binary_big_endian）。

    刻意不引 plyfile 之类的依赖：本仓库的点云只有 x,y,z,intensity 这种平坦布局，
    手写解析器 40 行，而且**出错时能说清是哪一行**——第三方库抛的异常在 10 万点文件上
    除了堆栈什么信息都没有。
    """
    with open(path, "rb") as f:
        raw = f.read()

    end = raw.find(b"end_header")
    if end < 0:
        raise ValueError(f"不是 PLY 文件（找不到 end_header）：{path}")
    nl = raw.find(b"\n", end)
    header = raw[:nl].decode("ascii", "replace")
    body = raw[nl + 1:]

    fmt = "ascii"
    n_vert = 0
    props: list[tuple[str, str]] = []   # (name, ply_type)
    in_vertex = False
    for line in header.splitlines():
        t = line.strip()
        if t.startswith("format "):
            fmt = t.split()[1]
        elif t.startswith("element "):
            parts = t.split()
            in_vertex = parts[1] == "vertex"
            if in_vertex:
                n_vert = int(parts[2])
        elif t.startswith("property ") and in_vertex:
            parts = t.split()
            if parts[1] != "list":
                props.append((parts[2], parts[1]))

    if n_vert <= 0:
        raise ValueError("PLY 里没有 vertex 元素")
    names = [p[0] for p in props]
    for k in ("x", "y", "z"):
        if k not in names:
            raise ValueError(f"PLY 缺少 {k} 属性")

    if fmt == "ascii":
        rows = body.split(b"\n")
        out = np.empty((n_vert, 3), dtype=float)
        got = 0
        for line in rows:
            s = line.strip()
            if not s:
                continue
            vals = s.split()
            if len(vals) < len(props):
                continue
            out[got, 0] = float(vals[names.index("x")])
            out[got, 1] = float(vals[names.index("y")])
            out[got, 2] = float(vals[names.index("z")])
            got += 1
            if got == n_vert:
                break
        if got < n_vert:
            out = out[:got]
        return out

    # 二进制：按属性表拼 dtype
    typemap = {
        "char": "i1", "int8": "i1", "uchar": "u1", "uint8": "u1",
        "short": "i2", "int16": "i2", "ushort": "u2", "uint16": "u2",
        "int": "i4", "int32": "i4", "uint": "u4", "uint32": "u4",
        "float": "f4", "float32": "f4", "double": "f8", "float64": "f8",
    }
    dt = np.dtype([(n, typemap[t]) for n, t in props])
    bo = "<" if fmt.endswith("little_endian") else ">"
    dt = dt.newbyteorder(bo)
    arr = np.frombuffer(body, dtype=dt, count=n_vert)
    return np.column_stack([arr["x"], arr["y"], arr["z"]]).astype(float)


def write_ply_xyz(path: str, pts: np.ndarray) -> None:
    n = len(pts)
    with open(path, "w", encoding="ascii") as f:
        f.write("ply\nformat ascii 1.0\ncomment robot_safety export\n")
        f.write(f"element vertex {n}\n")
        f.write("property float x\nproperty float y\nproperty float z\n")
        f.write("end_header\n")
        for p in pts:
            f.write(f"{p[0]:.3f} {p[1]:.3f} {p[2]:.3f}\n")


# ---------------------------------------------------------------- 碰撞查询

@dataclass
class SegmentHit:
    surface_dist: float          # <0 穿透；= 轴距 - 半径
    axis_dist: float
    nearest_cloud: np.ndarray
    nearest_axis: np.ndarray
    inside_count: int            # 落在胶囊内部的点云点数（穿透证据，不是估计值）
    saturated: bool              # True = 搜索半径内没有点云，surface_dist 只是下界
    part: str = ""


class CloudCollider:
    """一份点云的碰撞查询器。索引只建一次，可对任意多个姿态反复查询。"""

    def __init__(self, pts: np.ndarray, name: str = ""):
        if len(pts) == 0:
            raise ValueError("点云为空")
        # **float32**：距离判定的量级是 0.1~100mm，float32 在 1000mm 坐标上的分辨率约 6e-5mm，
        # 对判定毫无影响，但内存带宽减半（这是全量计算的主要成本）。
        self.pts = np.ascontiguousarray(pts, dtype=np.float32)
        self.name = name
        self.tree = cKDTree(self.pts.astype(np.float64)) if cKDTree is not None else None
        self.min = self.pts.min(axis=0).astype(np.float64)
        self.max = self.pts.max(axis=0).astype(np.float64)
        self.centroid = self.pts.mean(axis=0).astype(np.float64)
        self.raw_count = len(self.pts)
        if len(self.pts) > 400:
            idx = np.random.default_rng(0).choice(len(self.pts), 400, replace=False)
            if cKDTree is not None:
                d, _ = cKDTree(self.pts.astype(np.float64)).query(
                    self.pts[idx].astype(np.float64), k=2)
                self.spacing = float(np.median(d[:, 1]))
            else:
                out = []
                for j in idx:
                    dd = np.linalg.norm(self.pts - self.pts[j], axis=1)
                    dd[j] = np.inf
                    out.append(dd.min())
                self.spacing = float(np.median(out))
        else:
            self.spacing = 0.0

    # 线段到点云 AABB 的下界距离（剪枝用；过松只会多算不会漏）
    def seg_to_aabb(self, a: np.ndarray, b: np.ndarray) -> float:
        best = np.inf
        for t in np.linspace(0.0, 1.0, 9):
            p = a + (b - a) * t
            d = np.maximum(0.0, np.maximum(self.min - p, p - self.max))
            best = min(best, float(np.linalg.norm(d)))
            if best <= 0:
                return 0.0
        return best

    def segment(self, a: np.ndarray, b: np.ndarray, radius: float, search: float = 30.0,
                part: str = "") -> SegmentHit:
        """线段 A-B、半径 radius 的胶囊与点云的最短表面距离。

        search 是「判定范围」：最近距离超过它就直接报 saturated（下界）——
        绝不把下界当实测值报出去（CloudCheck.cs 同样的纪律）。

        **实现选择：全量向量化，不做 KD-tree 球查询。**
        理由：球查询要按采样点逐个展开，长杆（445mm）要么采样太稀漏掉 5mm 薄壁，
        要么采样太密（40 个球 × 上千点）单次 3ms；而把点云一次性投到线段上算距离，
        14k 点只要 0.2ms、250k 点约 3ms，**没有采样漏检**——碰撞判定里「漏检」是不可接受的错误类型，
        宁可多花一点固定开销换确定性。
        """
        a = np.asarray(a, dtype=np.float32)
        b = np.asarray(b, dtype=np.float32)
        search = float(search)
        radius = float(radius)
        if self.seg_to_aabb(a.astype(np.float64), b.astype(np.float64)) > search + radius:
            return SegmentHit(search, search + radius, a.astype(np.float64), a.astype(np.float64), 0, True, part)

        AB = b - a
        L2 = float(np.dot(AB, AB))
        P = self.pts
        # 用「点到直线距离」的展开式，避免构造 (N,3) 的投影点数组：
        #     d² = |P-a|² - ((P-a)·AB)²/|AB|²
        # 再把 t<0 / t>1 的两端截断补回去。实测这比「先算 proj 再减」快 2~3 倍
        # （proj 的广播乘法在 11 万点上要 5.6ms，是整条链路的最大单项）。
        Q = P - a[None, :]
        q2 = np.einsum("ij,ij->i", Q, Q)
        if L2 < 1e-12:
            d2 = q2
        else:
            dot = Q @ AB
            d2 = q2 - (dot * dot) / L2
            np.maximum(d2, 0.0, out=d2)          # 浮点误差可能给出 -1e-12
            d2 = np.where(dot < 0.0, q2, d2)
            d2 = np.where(dot > L2, q2 - 2.0 * dot + L2, d2)
        k = int(np.argmin(d2))
        axis_dist = float(np.sqrt(max(d2[k], 0.0)))
        if axis_dist - radius > search:
            return SegmentHit(search, search + radius, a.astype(np.float64), a.astype(np.float64), 0, True, part)
        inside = int(np.count_nonzero(d2 < radius * radius))
        tt = float(np.clip(((P[k] - a) @ AB) / L2, 0.0, 1.0)) if L2 > 1e-12 else 0.0
        return SegmentHit(
            surface_dist=axis_dist - radius,
            axis_dist=axis_dist,
            nearest_cloud=P[k].astype(np.float64),
            nearest_axis=(a + AB * tt).astype(np.float64),
            inside_count=inside,
            saturated=False,
            part=part,
        )


def _estimate_spacing(pts: np.ndarray, sample: int = 400) -> float:
    """抽样估计最近邻间距（中位数）。用来看「抽稀后还能不能代表表面」。"""
    if len(pts) < 8:
        return 0.0
    rng = np.random.default_rng(0)
    idx = rng.choice(len(pts), min(sample, len(pts)), replace=False)
    if cKDTree is not None:
        tree = cKDTree(pts)
        d, _ = tree.query(pts[idx], k=2)
        return float(np.median(d[:, 1]))
    out = []
    for j in idx:
        dd = np.linalg.norm(pts - pts[j], axis=1)
        dd[j] = np.inf
        out.append(dd.min())
    return float(np.median(out))


def _decimate(pts: np.ndarray, max_points: int, spacing: float) -> np.ndarray:
    """体素抽稀到大约 max_points 点。

    ⚠ 体素边长**必须按表面点云的密度推**，不能按包围盒对角线推：
    这是表面采样（不是体采样），点数 ≈ 面积/间距²。按对角线除立方根会得到 36mm 的体素，
    11 万点直接掉到 190 点——那已经不是同一件工件了（实测过，最近距离全部饱和）。
    正确做法：目标间距 s' = s·√(N/target)，体素取 s'。
    """
    if len(pts) <= max_points or spacing <= 0:
        return pts
    s_target = spacing * np.sqrt(len(pts) / max_points)
    voxel = max(0.5, s_target)
    idx = np.floor((pts - pts.min(axis=0)) / voxel).astype(np.int64)
    _, keep = np.unique(idx, axis=0, return_index=True)
    return pts[np.sort(keep)]


def load_cloud(path: str, max_points: int = 25_000) -> "CloudCollider":
    pts = read_ply_xyz(path)
    raw_n = len(pts)
    sp = _estimate_spacing(pts)
    pts = _decimate(pts, max_points, sp)
    c = CloudCollider(pts, name=path)
    c.raw_count = raw_n
    return c
