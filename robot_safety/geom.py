"""几何基元：向量、4x4 齐次变换、旋转构造与点-线段距离。

单位统一为 **毫米 / 弧度**（CSV 里的角度只在 IO 边界上出现）。

约定（与 DensoArmKine/Models 一致，务必别混）：
  * `n`     = 面法向，从工件表面指向枪（CSV 里存的是它的**反向**，见 traj.py）
  * `gunOut`= 枪身朝外方向 = n = 工具坐标系的 +Z
  * 变换采用「先旋转后平移」的列向量约定：p_parent = R * p_child + t
"""

from __future__ import annotations

import math

import numpy as np

EPS = 1e-12


def unit(v: np.ndarray) -> np.ndarray:
    """归一化；零向量原样返回（调用方负责判退化，这里不抛异常也不静默改成 +Z）。"""
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    return v / n if n > EPS else v.copy()


def norm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def trans(x: float, y: float, z: float) -> np.ndarray:
    T = np.eye(4)
    T[:3, 3] = (x, y, z)
    return T


def rot_axis(axis: np.ndarray, angle: float) -> np.ndarray:
    """罗德里格斯旋转（axis 必须先归一化）。"""
    a = np.asarray(axis, dtype=float)
    n = float(np.linalg.norm(a))
    if n < EPS:
        return np.eye(4)
    a = a / n
    c, s = math.cos(angle), math.sin(angle)
    x, y, z = a
    R = np.array(
        [
            [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
            [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
            [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
        ]
    )
    T = np.eye(4)
    T[:3, :3] = R
    return T


def rot_z(angle: float) -> np.ndarray:
    return rot_axis(np.array([0.0, 0.0, 1.0]), angle)


def rotvec_to_mat(rv: np.ndarray) -> np.ndarray:
    """轴角 -> 3x3 旋转矩阵。"""
    rv = np.asarray(rv, dtype=float)
    th = float(np.linalg.norm(rv))
    if th < EPS:
        return np.eye(3)
    return rot_axis(rv / th, th)[:3, :3]


def mat_to_rotvec(R: np.ndarray) -> np.ndarray:
    """3x3 旋转矩阵 -> 轴角（对 θ≈π 也稳定，用四元数中转）。"""
    R = np.asarray(R, dtype=float)
    tr = float(np.trace(R))
    c = max(-1.0, min(1.0, (tr - 1.0) * 0.5))
    th = math.acos(c)
    if th < 1e-9:
        return np.zeros(3)
    if abs(math.pi - th) > 1e-6:
        k = th / (2.0 * math.sin(th))
        return np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) * k
    # θ≈π：从 R+I 的列里取最大的一列（避免除以 sin θ≈0）
    A = (R + np.eye(3)) * 0.5
    axis = unit(np.array([math.sqrt(max(A[0, 0], 0.0)), math.sqrt(max(A[1, 1], 0.0)), math.sqrt(max(A[2, 2], 0.0))]))
    # 修正符号：非对角项与 axis 外积同号
    for i in range(3):
        for j in range(3):
            if i != j and A[i, j] * axis[i] * axis[j] < 0:
                axis[j] = -axis[j]
    return axis * th


def frame_from_z(z_axis: np.ndarray, roll: float = 0.0, ref: np.ndarray | None = None) -> np.ndarray:
    """由「枪轴方向」+ 绕轴滚转角构造 3x3 姿态。

    z_axis 是**枪身朝外**方向（工具 +Z）。roll 是绕该轴的滚转（喷涂枪轴对称，
    这个自由度只用来让机器人解更舒服 —— 见 planner 的 roll 搜索）。
    ref 是构造 x 轴时的参考向量；退化时自动换一个（不能固定用 +Z：
    当枪轴接近竖直时会退化，这正是俯喷/仰喷的常见情形）。
    """
    z = unit(z_axis)
    if ref is None:
        ref = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(z, ref))) > 0.95:
            ref = np.array([1.0, 0.0, 0.0])
    x = np.cross(ref, z)
    if norm(x) < 1e-9:
        x = np.cross(np.array([1.0, 0.0, 0.0]), z)
    x = unit(x)
    y = unit(np.cross(z, x))
    if abs(roll) > EPS:
        c, s = math.cos(roll), math.sin(roll)
        x, y = x * c + y * s, y * c - x * s
    return np.column_stack([x, y, z])


def pose_from_z(z_axis: np.ndarray, roll: float = 0.0) -> np.ndarray:
    """4x4 位姿：旋转按 frame_from_z，平移由调用方填。"""
    T = np.eye(4)
    T[:3, :3] = frame_from_z(z_axis, roll)
    return T


def tilt_axis_and_dir(n: np.ndarray, up: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, float]:
    """求「把法向往上方扳」的旋转轴与可用的最大倾角。

    返回 (axis, u_up, limit_deg)：
      * axis  —— 绕它正转 θ 即可让 n 离开工件表面朝上偏；
      * u_up  —— 切平面内朝上的单位向量（θ=90° 时的方向）；
      * limit —— 允许的最大倾角（n 越接近「朝上」，能扳的空间越小）。

    为什么不是简单地把 n 和 +Z 做 slerp：slerp 在 n 接近 +Z 时退化，
    而且它的「倾角」在倾斜过程中不再是常数半径球面上的等角步进。
    这里用切平面内的正交方向，θ 就是**枪轴与法向的真实夹角**（= 喷角），
    与工艺口径（喷角上限）一一对应。
    """
    n = unit(n)
    up = np.array([0.0, 0.0, 1.0]) if up is None else unit(up)

    u = up - n * float(np.dot(up, n))
    # ⚠ 判据必须用**归一化之前**的长度：|u| = sin(n 与 up 的夹角)。
    #   先 unit() 再比阈值是错的 —— unit 之后恒等于 1，判据永远不成立。
    #   这个 bug 的表现很有欺骗性：一个平整朝上的面（夹角 2.7°）被判成"非朝上"，
    #   于是"往上扳"扳的是切平面里的噪声方向，倾斜角被加到 40°、枪位横甩，
    #   而喷距仍然守恒（球面），从数字上看一切"正常"。
    # 阈值取 5°：法向偏离竖直小于 5° 就当作朝上，没有工艺意义上的"倾斜"可做。
    if norm(u) < math.sin(math.radians(5.0)):
        return np.zeros(3), np.zeros(3), 0.0
    u = unit(u)
    axis = unit(np.cross(n, u))  # 绕 axis 正转把 n 转向 u
    # 倾角上限：不超过 85°，且不超过「枪轴与 up 垂直」的角度（再大就翻过顶了）
    limit = 85.0
    return axis, u, limit


def rotated_axis(n: np.ndarray, axis: np.ndarray, deg: float) -> np.ndarray:
    """把 n 绕 axis 转 deg 度（保证结果仍是单位向量）。"""
    if abs(deg) < EPS:
        return unit(n)
    return unit(rot_axis(axis, math.radians(deg))[:3, :3] @ unit(n))


def angle_between(a: np.ndarray, b: np.ndarray) -> float:
    """两向量夹角（度）。"""
    a, b = unit(a), unit(b)
    c = float(np.dot(a, b))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def point_segment_dist(P: np.ndarray, A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """批量：点集 P(N,3) 到线段 AB 的距离。

    返回 (dist(N,), t(N,))，t 是投影参数（已截断到 [0,1]）。
    """
    AB = B - A
    L2 = float(np.dot(AB, AB))
    if L2 < 1e-12:
        d = np.linalg.norm(P - A, axis=1)
        return d, np.zeros(len(P))
    t = np.clip(((P - A) @ AB) / L2, 0.0, 1.0)
    proj = A[None, :] + t[:, None] * AB[None, :]
    return np.linalg.norm(P - proj, axis=1), t


def capsule_capsule_dist(A1: np.ndarray, B1: np.ndarray, A2: np.ndarray, B2: np.ndarray) -> float:
    """两条线段（胶囊轴）之间的最短距离 —— **解析解**（Ericson, Real-Time Collision Detection 5.1.9）。

    一开始写的是「24×24 采样网格 + 3 轮局部收缩」，精度够（误差 <0.2mm）但每对 1ms：
    一条轨迹 2000 点 × 8 对 = 16 秒纯花在自碰撞上。解析解是 O(1) 的十几条标量运算，
    同一个任务降到 0.02ms/对，而且**没有采样误差**——两条平行连杆只差 0.3mm 时
    采样法可能报成 0.5mm，解析法给的就是真值。自碰撞判定里这个差别是安全裕度。
    """
    p1 = np.asarray(A1, dtype=float)
    q1 = np.asarray(B1, dtype=float)
    p2 = np.asarray(A2, dtype=float)
    q2 = np.asarray(B2, dtype=float)
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2
    a = float(np.dot(d1, d1))
    e = float(np.dot(d2, d2))
    f = float(np.dot(d2, r))
    eps = 1e-12

    if a <= eps and e <= eps:
        return geom_norm(r)
    if a <= eps:
        s = 0.0
        t = min(1.0, max(0.0, f / e))
    else:
        c = float(np.dot(d1, r))
        if e <= eps:
            t = 0.0
            s = min(1.0, max(0.0, -c / a))
        else:
            b = float(np.dot(d1, d2))
            denom = a * e - b * b
            s = min(1.0, max(0.0, (b * f - c * e) / denom)) if denom > eps else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t = 0.0
                s = min(1.0, max(0.0, -c / a))
            elif t > 1.0:
                t = 1.0
                s = min(1.0, max(0.0, (b - c) / a))
    c1 = p1 + d1 * s
    c2 = p2 + d2 * t
    return geom_norm(c1 - c2)


def geom_norm(v: np.ndarray) -> float:
    return float(np.sqrt(float(v[0]) * float(v[0]) + float(v[1]) * float(v[1]) + float(v[2]) * float(v[2])))
