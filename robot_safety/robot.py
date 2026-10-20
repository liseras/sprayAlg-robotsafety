"""DENSO VS-087 的连杆模型、正/逆运动学、连杆碰撞体。

**这份文件是 DensoArmKine 的 Python 移植，参数与算法同源，不许另立一套：**
  * 连杆表     <- Models/DensoLibrary.cs::Vs087()      （官方 vs087.urdf，reach 905mm 已用铭牌核实）
  * FK 公式    <- Kinematics/FkSolver.cs               T_i = T_{i-1}·Trans(offset_i)·Rot(axis_i, θ_i)
  * IK         <- 自研阻尼最小二乘（DensoArmKine 的 InverseKinematics.cs 是多初值数值解，语义一致）
  * 碰撞体     <- Cloud/ArmCollider.cs                 ArmSizing.From + Build（半径唯一来源）

为什么要移植而不是直接调 DensoArmKine.exe：
  倾斜角规划要**对每个点试十几个候选姿态**（每个候选都要 IK + 碰撞），
  外部进程一次几秒，一条 2000 点的轨迹要跑几个小时；进程内一次 IK ~2ms，整条 ~1 分钟。
  对拍仍用 DensoArmKine 的 --cloud-check / --traj-cloud 做抽样交叉验证（见 verify_against_denso.py）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import geom


@dataclass(frozen=True)
class Joint:
    name: str
    off: tuple[float, float, float]
    axis: tuple[float, float, float]
    lo: float
    hi: float


# VS-087（DENSO WAVE 官方 vs087.urdf，与 DensoLibrary.Vs087 逐字段一致）
VS087_JOINTS = (
    Joint("J1", (0.0, 0.0, 197.5), (0.0, 0.0, 1.0), -170.0, 170.0),
    Joint("J2", (30.0, 0.0, 197.5), (0.0, 1.0, 0.0), -100.0, 135.0),
    Joint("J3", (0.0, 0.0, 445.0), (0.0, 1.0, 0.0), -133.0, 150.0),
    Joint("J4", (-20.0, 0.0, 287.0), (0.0, 0.0, 1.0), -270.0, 270.0),
    Joint("J5", (0.0, 0.0, 143.0), (0.0, 1.0, 0.0), -120.0, 120.0),
    Joint("J6", (0.0, 0.0, 80.0), (0.0, 0.0, 1.0), -360.0, 360.0),
)

VS087_REACH = 905.0
VS087_TOOL_LEN = 150.0   # DefaultTool.LengthZ
VS087_PAYLOAD = 7.0


@dataclass
class Tool:
    length: float = VS087_TOOL_LEN   # 法兰面 -> 喷嘴出口（TCP）
    twist_deg: float = 0.0
    ox: float = 0.0
    oy: float = 0.0
    oz: float = 0.0


@dataclass
class FkResult:
    frames: np.ndarray      # (7,4,4) frames[0]=基座, frames[i]=J_i 子坐标系
    origins: np.ndarray     # (7,3)  各关节中心
    flange: np.ndarray      # (4,4)
    tcp: np.ndarray         # (4,4)
    axes_world: np.ndarray  # (6,3)  各关节转轴（基座系）


def fk(theta_deg: np.ndarray, tool: Tool) -> FkResult:
    """正运动学。theta_deg 为 6 个关节角（度）。"""
    th = np.asarray(theta_deg, dtype=float)
    n = len(VS087_JOINTS)
    frames = np.zeros((n + 1, 4, 4))
    axes = np.zeros((n, 3))
    frames[0] = np.eye(4)
    for i, j in enumerate(VS087_JOINTS):
        R_pre = frames[i][:3, :3]
        ax = geom.unit(R_pre @ np.array(j.axis))     # 转轴在世界系（Trans 不改旋转）
        axes[i] = ax
        local = geom.trans(*j.off) @ geom.rot_axis(np.array(j.axis), math.radians(th[i]))
        frames[i + 1] = frames[i] @ local
    flange = frames[n]
    tool_local = geom.trans(tool.ox, tool.oy, tool.oz + tool.length) @ geom.rot_z(
        math.radians(tool.twist_deg)
    )
    tcp = flange @ tool_local
    return FkResult(frames=frames, origins=frames[:, :3, 3].copy(), flange=flange, tcp=tcp, axes_world=axes)


def over_limit(theta_deg: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    th = np.asarray(theta_deg, dtype=float)
    return np.array([not (j.lo - tol <= th[i] <= j.hi + tol) for i, j in enumerate(VS087_JOINTS)])


def limit_margin_deg(theta_deg: np.ndarray) -> np.ndarray:
    """每个关节到最近限位的余量（度），负值表示已超限。"""
    th = np.asarray(theta_deg, dtype=float)
    return np.array([min(th[i] - j.lo, j.hi - th[i]) for i, j in enumerate(VS087_JOINTS)])


def jacobian(fk_res: FkResult, tool: Tool) -> np.ndarray:
    """TCP 位姿的 6x6 解析雅可比（线速度在前、角速度在后）。

    TCP 固定在法兰 +Z 上，绕法兰 Z 自转不改变 TCP 位置，所以 J6 的线速度列自然为
    z6 × (p_tcp - o6)；这里不做任何简化，交给公式自己给出。
    """
    p = fk_res.tcp[:3, 3]
    J = np.zeros((6, 6))
    for i in range(6):
        J[:3, i] = np.cross(fk_res.axes_world[i], p - fk_res.origins[i + 1])
        J[3:, i] = fk_res.axes_world[i]
    _ = tool
    return J


# 多初值种子（度）：覆盖肩前/肩后 × 肘上/肘下 × 腕翻/腕正，DensoArmKine 的 IkSeeds.cs 同一思路。
# 这些只是**种子**，不承诺落在哪个解支；最终选解靠 cost（限位余量 + 关节变化量）。
SEEDS = (
    (0, 0, -90, 0, 0, 0),
    (0, 20, -60, 0, 45, 0),
    (0, -20, -60, 0, -45, 0),
    (180, 0, -90, 0, 0, 0),
    (180, 20, -60, 0, 45, 0),
    (180, -20, -60, 0, -45, 0),
    (90, 10, -70, 90, 30, 90),
    (-90, 10, -70, -90, 30, -90),
    (0, 60, -120, 0, 60, 0),
    (180, 60, -120, 0, 60, 0),
    (0, -40, -30, 0, -30, 0),
    (180, -40, -30, 0, -30, 0),
)


@dataclass
class IkSolution:
    theta: np.ndarray
    ok: bool
    pos_err: float
    rot_err_deg: float
    iters: int
    over: bool
    margin: np.ndarray = field(default_factory=lambda: np.zeros(6))
    cond: float = 0.0


def analytic_seeds(p_tcp: np.ndarray, R_tcp: np.ndarray, tool: Tool) -> list[np.ndarray]:
    """用几何关系直接造种子：前三个关节解析求解，腕部三个关节由 ZYZ 分解给初值。

    为什么必须这么干（实测，不是理论顾虑）：纯固定种子 `(0,0,-90,0,0,0)` 在这条轨迹上收敛到
    J2=154°、J3=-102°、J6=-279° 的怪解 —— 位置姿态残差都是 0，但 J2 超限 19°，
    整条小臂从工件底下钻过去（J5-J6 连杆穿透 24.7mm）。DLS 只保证"收敛到附近的解"，
    不保证是**人愿意开的那一支**；种子决定解支，所以种子必须给对。

    推导（VS-087）：
        J2 = (30, 0, 197.5)，J2/J3 都绕 Y ⇒ 臂在 J1 系的 (r,z) 平面内是标准两连杆：
        J3 = J2 + 445·(sinθ2, 0, cosθ2)，腕心 ≈ J3 + 430·(sin(θ2+θ3), 0, cos(θ2+θ3))
        两连杆定理给出 (θ2, θ3) 两支（肘上/肘下）；θ1 由腕心水平方位给出（正向/反向）。
    腕部：R03 = Rz(θ1)·Ry(θ2+θ3)，R36 = R03ᵀ·R_target，再按 ZYZ（Rz·Ry·Rz）分解出 (θ4,θ5,θ6) 两支。

    注意腕心是**近似**：VS-087 的腕不是理想球腕（J4 有 −20mm 横向偏置、J6 相对 J5 前移 80mm），
    所以这里只当种子，末位仍由 DLS 修正到毫米级。
    """
    d = np.asarray(R_tcp, dtype=float)[:, 2]          # 工具 +Z（基座系）
    # 腕心近似：TCP 往回收「工具长 + J5→J6 的 80mm」
    W = np.asarray(p_tcp, dtype=float) - d * (tool.length + 80.0)
    r = math.hypot(W[0], W[1])
    wz = W[2]
    lat, a1, a2, base_h = 30.0, 445.0, 430.0, 197.5
    Rt = np.asarray(R_tcp, dtype=float)
    out: list[np.ndarray] = []
    for s1 in (1.0, -1.0):
        th1 = math.degrees(math.atan2(W[1], W[0])) + (0.0 if s1 > 0 else 180.0)
        dx = r * s1 - lat
        dz = wz - base_h
        D2 = dx * dx + dz * dz
        c3 = (D2 - a1 * a1 - a2 * a2) / (2.0 * a1 * a2)
        if abs(c3) > 1.0:
            continue
        for sgn in (1.0, -1.0):
            th3 = sgn * math.degrees(math.acos(c3))
            th2 = math.degrees(
                math.atan2(dx, dz)
                - math.atan2(a2 * math.sin(math.radians(th3)), a1 + a2 * math.cos(math.radians(th3)))
            )
            R03 = geom.rot_axis(np.array([0.0, 0.0, 1.0]), math.radians(th1))[:3, :3] @ geom.rot_axis(
                np.array([0.0, 1.0, 0.0]), math.radians(th2 + th3)
            )[:3, :3]
            R36 = R03.T @ Rt
            for bsign in (1.0, -1.0):
                sb = math.hypot(R36[0, 2], R36[1, 2]) * bsign
                b = math.atan2(sb, R36[2, 2])
                if abs(sb) < 1e-9:
                    a = 0.0
                    c = math.atan2(-R36[1, 0], R36[0, 0])
                else:
                    a = math.atan2(R36[1, 2], R36[0, 2])
                    c = math.atan2(R36[2, 1], -R36[2, 0])
                out.append(np.array([th1, th2, th3,
                                     math.degrees(a), math.degrees(b), math.degrees(c)]))
    return out


def ik(
    p_target: np.ndarray,
    R_target: np.ndarray,
    tool: Tool,
    seed: np.ndarray | None = None,
    extra_seeds: list[np.ndarray] | None = None,
    iters: int = 90,
    tol_pos: float = 0.02,
    tol_rot_deg: float = 0.02,
    require_no_over: bool = True,
) -> IkSolution:
    """TCP 位姿逆解。返回**最好的一支**（残差最小、且尽量不超限），不做「解不出来就编一个」。

    种子顺序：上一点解（保证连续性）→ 解析种子（保证解支正确）→ 固定种子表（兜底）。
    解不出来（残差降不下去）时 ok=False，调用方必须把它当成**不可达**处理。
    """
    p_target = np.asarray(p_target, dtype=float)
    R_target = np.asarray(R_target, dtype=float)
    seeds: list[np.ndarray] = []
    if seed is not None:
        seeds.append(np.asarray(seed, dtype=float))
    seeds.extend(analytic_seeds(p_target, R_target, tool))
    if extra_seeds:
        seeds.extend(np.asarray(s, dtype=float) for s in extra_seeds)
    seeds.extend(np.array(s, dtype=float) for s in SEEDS)

    best: IkSolution | None = None
    for s in seeds:
        sol = _ik_once(p_target, R_target, tool, s, iters, tol_pos, tol_rot_deg)
        # 快路径：收敛、残差好、**且不超限** —— 直接采纳。
        # require_no_over 是必须的：去掉它之后（曾为了提速去掉过），
        # DLS 会给出 J3=-284°、J6=-279° 这种"数学上正确、机械上要走线"的解，
        # 而且它比正确解收敛得更快，于是 bad case 被稳定地选中。
        if sol.ok and sol.pos_err < 0.05 and sol.rot_err_deg < 0.05:
            if not require_no_over or not sol.over:
                return sol
        if best is None:
            best = sol
            continue

        def rank(x: IkSolution):
            # 收敛 > 不限位 > 残差小 > 单轴超限少 > 余量大
            over_sum = float(np.sum(np.maximum(0.0, -x.margin)))
            return (0 if x.ok else 1, 0 if not x.over else 1, over_sum,
                    x.pos_err + 100.0 * x.rot_err_deg, -float(np.min(x.margin)))

        if rank(sol) < rank(best):
            best = sol
    return best


def _ik_once(
    p_target: np.ndarray,
    R_target: np.ndarray,
    tool: Tool,
    seed: np.ndarray,
    iters: int,
    tol_pos: float,
    tol_rot_deg: float,
) -> IkSolution:
    th = np.array(seed, dtype=float)
    w_rot = 100.0     # 1 rad 姿态误差 ≈ 100mm 位置误差（毫米量纲下把两行拉齐）
    lam = 30.0
    best_err = np.inf
    it = 0
    for it in range(1, iters + 1):
        r = fk(th, tool)
        e_p = p_target - r.tcp[:3, 3]
        e_r = geom.mat_to_rotvec(R_target @ r.tcp[:3, :3].T)
        err = np.concatenate([e_p, w_rot * e_r])
        en = float(np.linalg.norm(err))
        if en < best_err:
            best_err = en
            lam = max(lam * 0.7, 0.5)
        else:
            lam = min(lam * 0.5 + lam, 5e3)   # 变差就加大阻尼，避免发散
        if float(np.linalg.norm(e_p)) < tol_pos and math.degrees(float(np.linalg.norm(e_r))) < tol_rot_deg:
            break
        J = jacobian(r, tool)
        Jw = J.copy()
        Jw[3:, :] *= w_rot
        JJt = Jw @ Jw.T + (lam * lam) * np.eye(6)
        try:
            dth = Jw.T @ np.linalg.solve(JJt, err)
        except np.linalg.LinAlgError:
            break
        # 单步限幅：一步最多 30°，否则线性化失效（在限位附近会剧烈来回）
        mx = float(np.max(np.abs(dth)))
        if mx > math.radians(30.0):
            dth *= math.radians(30.0) / mx
        th = th + np.degrees(dth)

    r = fk(th, tool)
    e_p = float(np.linalg.norm(p_target - r.tcp[:3, 3]))
    e_r = math.degrees(float(np.linalg.norm(geom.mat_to_rotvec(R_target @ r.tcp[:3, :3].T))))
    over = bool(np.any(over_limit(th)))
    margin = limit_margin_deg(th)
    # 条件数必须用**归一化雅可比**：原始 J 的前三行是 mm（量纲 ~500），后三行是 rad（~1），
    # 直接算 SVD 得到的条件数天然上千，与"是否接近奇异"没有关系（实测正常姿态就有 1118）。
    # 位置行除以特征长度 500mm 之后，条件数才是无量纲的姿态灵敏度指标。
    J = jacobian(r, tool)
    Jn = J.copy()
    Jn[:3, :] /= 500.0
    sv = np.linalg.svd(Jn, compute_uv=False)
    cond = float(sv[0] / sv[-1]) if sv[-1] > 1e-9 else float("inf")
    return IkSolution(
        theta=th,
        ok=(e_p <= max(tol_pos * 20, 0.5) and e_r <= max(tol_rot_deg * 20, 0.5)),
        pos_err=e_p,
        rot_err_deg=e_r,
        iters=it,
        over=over,
        margin=margin,
        cond=cond,
    )


# ---------------------------------------------------------------- 碰撞体

@dataclass
class ArmSizing:
    joint_r: float
    joint1_r: float
    link_r: float
    wrist_link_r: float
    base_column_r: float
    base_neck_r: float
    grip_r: float
    tool_pipe_r: float
    nozzle_r: float

    @staticmethod
    def of(reach: float = VS087_REACH) -> "ArmSizing":
        joint_r = min(46.0, max(12.0, reach * 0.038))
        grip_r = joint_r * 0.62
        return ArmSizing(
            joint_r=joint_r,
            joint1_r=joint_r * 1.25,
            link_r=joint_r * 0.62,
            wrist_link_r=joint_r * 0.72,
            base_column_r=joint_r * 2.6,
            base_neck_r=joint_r * 1.15,
            grip_r=grip_r,
            tool_pipe_r=grip_r * 0.45,
            nozzle_r=grip_r * 0.30,
        )


GRIP_FRACTION = 0.72     # 握把占「法兰→TCP」全长比例（与 ArmCollider 一致）
NOZZLE_FRACTION = 0.18   # 喷嘴伸出 TCP 的比例（锥按底半径圆柱保守近似）


@dataclass
class Seg:
    a: np.ndarray
    b: np.ndarray
    r: float
    part: str
    group: str   # 'tool' / 'link' / 'joint' / 'base'

    @property
    def is_tool(self) -> bool:
        return self.group == "tool"


def arm_segments(r: FkResult, sizing: ArmSizing, include_base: bool = False, include_tool: bool = True) -> list[Seg]:
    """由 FK 结果生成碰撞段（胶囊）。与 ArmCollider.Build 逐段对应。"""
    segs: list[Seg] = []
    o = r.origins
    if include_base:
        base_h = max(4.0, o[1][2] * 0.55)
        segs.append(Seg(np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, base_h]), sizing.base_column_r, "底座立柱", "base"))
        segs.append(Seg(np.array([0.0, 0.0, base_h]), o[1].copy(), sizing.base_neck_r, "底座颈", "base"))
    for i in range(1, 7):
        jr = sizing.joint1_r if i == 1 else sizing.joint_r
        segs.append(Seg(o[i].copy(), o[i].copy(), jr, f"J{i} 关节", "joint"))
        if i + 1 < 7:
            rr = sizing.wrist_link_r if i >= 4 else sizing.link_r
            segs.append(Seg(o[i].copy(), o[i + 1].copy(), rr, f"连杆 J{i}-J{i+1}", "link"))
    if include_tool:
        flange = r.flange[:3, 3]
        tcp = r.tcp[:3, 3]
        d = tcp - flange
        L = geom.norm(d)
        if L > 0.5:
            u = d / L
            neck = flange + u * (L * GRIP_FRACTION)
            nozzle_tip = tcp + u * (L * NOZZLE_FRACTION)
            segs.append(Seg(flange.copy(), neck, sizing.grip_r, "工具握把", "tool"))
            segs.append(Seg(neck, tcp.copy(), sizing.tool_pipe_r, "工具细管", "tool"))
            segs.append(Seg(tcp.copy(), nozzle_tip, sizing.nozzle_r, "喷嘴", "tool"))
    return segs


# 自碰撞时必须**排除**的段对：同一关节上的两段、相邻连杆（它们本来就贴在一起）
def _self_collision_pairs(segs: list[Seg]) -> list[tuple[int, int]]:
    """返回需要检查的段对（非相邻、且不是「关节球 + 它的两根连杆」这种天然接触）。"""
    idx = [i for i, s in enumerate(segs) if s.group in ("link", "joint")]
    pairs = []
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            ia, ib = idx[a], idx[b]
            sa, sb = segs[ia], segs[ib]
            # 同属一个关节链节（例如 J2 关节球与 连杆 J2-J3）: 名字里共享 J 号就跳过
            ja = _joint_ids(sa.part)
            jb = _joint_ids(sb.part)
            if ja & jb:
                continue
            pairs.append((ia, ib))
    return pairs


def _joint_ids(part: str) -> set[int]:
    out = set()
    tok = part.replace("连杆", "").replace("关节", "").replace("J", " ")
    for t in tok.replace("-", " ").split():
        if t.isdigit():
            out.add(int(t))
    return out


_PAIRS_CACHE: dict[int, list[tuple[int, int]]] = {}


def self_collision_pairs(segs: list[Seg]) -> list[tuple[int, int]]:
    """段对缓存的键是段数（同一套 include_base/tool 组合段数恒定，拓扑也就恒定）。"""
    key = len(segs)
    if key not in _PAIRS_CACHE:
        _PAIRS_CACHE[key] = _self_collision_pairs(segs)
    return _PAIRS_CACHE[key]
