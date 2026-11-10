"""内壁倾斜角规划 + 逐点安全评估。

本文件解决的三个问题（按用户原话）：

1. **「内壁法向如果不是朝上的，需要加一个倾斜角」**
   非朝上的面（侧壁/内壁，n·Z < cos(45°)）默认给一个倾角 θ0（默认 20°）；
   若这个倾角仍会撞，就按 5° 步进往上加，直到 clearance 达标或到上限 60°。

2. **「法向 90° 朝两侧的，倾角改变了，但机器人末端到那个点之间的距离是不变的」**
   倾斜**不是**「把枪位沿某个方向平移」，而是**绕工件表面点做球面旋转**：
       TCP' = face + standoff · d(θ),   d(θ) = Rot(axis, θ) · n
   于是 |TCP' - face| ≡ standoff，喷距严格不变；变的只有**喷角**（= θ，枪轴与法向夹角）。
   每一点都算 standoff_err = ||TCP'-face| - standoff|，超过 1e-3mm 就打警告，不靠"看起来对"。

3. **「要考虑碰撞：轨迹会不会跟原始点云工件碰撞，关节会不会碰撞」**
   * 工件碰撞：连杆/枪体（胶囊）vs 原始点云（surface = 点到轴距离 - 半径，<0 即穿透）
   * 自碰撞  ：非相邻连杆胶囊两两之间（排除同一关节与其两根连杆）
   * 可达/限位/奇异：IK 解存在、关节角在软限位内、雅可比条件数不超限

姿态的**面级常量**（这是能不能上机的关键，不是可选项）：
   逐点独立求出的倾角会在相邻点之间跳（20°→45°），100mm 外的枪位会甩出几十毫米，
   真机跟不上；更要紧的是**工艺要求同一个面上喷角不变**。所以：
     ① 先求每点的**需求**倾角 need[i]（"这一点至少要多大的倾角才安全"）；
     ② 把相邻、法向一致、共面的行归成「面组」，一个面组只取一个 (倾角, 方位)；
     ③ 再按这个常量重新求解位姿；面内个别点若过不了硬约束，逐点退回本点最优并记
        reverted —— 那是物理上必须偏的点，不是决策抖动。
   早期做法是"只升不降的逐点限幅平滑"（每点最多变 3°），它能保证 |Δθ| ≤ 3°，
   但保不住"同一个面上倾角不变"：实测一笔 77 点里能从 30° 走到 55° 再落回 30°。

决策顺序（**最小修正原则**，改动量从小到大）：
   ① 基准倾角 ② 逐步加大 ③ 换倾斜方位（左右各 30°/60°）
每个点的选择都记录在报告里，能回答「这一点为什么要转 35°」。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import cloud as cloud_mod
from . import geom, robot, traj as traj_mod


@dataclass
class Params:
    standoff: float = 100.0            # 喷距（mm），master 的 standoff 参数
    safety: float = 20.0               # 间隙告警阈值（mm）：表面距离小于它 → P1 告警（软约束）
    hard_clearance: float = 2.0        # 绝对安全距离（mm）：平滑姿态低于它才回退（硬约束）；0=只要不穿透
    # 判定带 120mm：喷距本身就有 100mm，判定带比它小的话所有段都会"饱和"——
    # 报告里只剩一排"≥40"，等于什么都没测到。这是 DensoArmKine CloudCheck 的同一口径
    # （EffectiveSearchMm = max(120, safety*4)）。
    search: float = 120.0              # 判定带（mm）：超过它就报下界，不当实测值
    tilt_nonup_deg: float = 30.0       # α=90°（竖直侧壁）时的倾角；实际 = 它 × α/90°
    up_threshold_deg: float = 45.0     # 已不参与倾角决策（倾角改线性），保留以对齐 C++ ABI
    tilt_max_deg: float = 60.0         # 倾角上限
    tilt_step_deg: float = 5.0
    tilt_rate_deg_per_pt: float = 3.0  # 已不参与姿态决策（改为面级常量）：仅用作姿态跳变告警阈值
    max_spray_angle_deg: float = 30.0  # 工艺允许的最大喷角（超出报工艺违规，不是碰撞违规）
    tool_len: float = robot.VS087_TOOL_LEN
    include_base: bool = False         # 点云里若有地面/底座，打开它
    joint_margin_min: float = 1.0      # 关节限位余量下限（度）
    cond_max: float = 30.0             # 归一化雅可比条件数上限（奇异判定，只告警不阻断）
    link_stride: int = 1               # 连杆碰撞抽样步长（1=逐点全判，产品版默认）
    fix_transit: bool = True           # 过渡点（spray=0）碰撞时允许沿 +Z 抬升
    transit_lift_max: float = 400.0    # 过渡点最大抬升量（mm）
    az_fallback: bool = True           # 允许把「上翻方向」左右各偏 30/60° 作为兜底
    # 是否要求"这台机器人可达"。False = 只查几何（喷距/倾角/枪体），不解 IK、
    # 不建连杆碰撞体、不查自碰撞与关节限位。见 C++ 侧 RobotSafetyParams.requireReach。
    require_reach: bool = True


@dataclass
class PointResult:
    idx: int
    tilt_deg: float
    clearance: float
    min_part: str
    ik_ok: bool
    over_limit: bool
    margin_min: float
    cond: float
    theta: np.ndarray
    pos: np.ndarray
    gun: np.ndarray
    spray_angle_to_normal: float
    standoff_err: float
    violations: list[str] = field(default_factory=list)
    self_collision: float = math.inf    # 自碰撞最小间隙（负 = 干涉）
    lifted: float = 0.0                 # 过渡点抬升量
    az_deg: float = 0.0                 # 倾斜方位（相对"向上"方向绕法向的转角）
    reverted: bool = False              # 平滑位姿违反硬约束（低于绝对安全距离/碰撞/自碰撞/关节超限），已退回本点最优


@dataclass
class Result:
    points: list[PointResult]
    base: np.ndarray
    stats: dict
    violations: list[dict]
    params: Params
    cloud_info: dict


# ---------------------------------------------------------------- 基座位姿

AUTO_BASE_FWD = 0.50    # 基座沿工件外侧后退量（× reach）
AUTO_BASE_DOWN = 0.35   # 基座下沉量（× reach）


def auto_base(collider: cloud_mod.CloudCollider, traj: traj_mod.Traj, reach: float = robot.VS087_REACH) -> np.ndarray:
    """点云坐标系里「猜」一个基座位姿（T_base_in_cloud）。

    ⚠ 这是**假设**，不是标定值。真实的 T_cloud→base 必须来自手眼标定/工位图纸。
    猜法：基座沿工件外侧退 0.50·reach、下沉 0.35·reach，使工件质心落在 VS-087 的舒适区，+X 指向工件。
    报告里会显式标注「基座位置为默认假设」，并用可达率告诉用户猜得准不准 ——
    可达率低于 90% 基本就是猜错了，这时应当用 --base-pose / --base-json 给真值。
    """
    c = collider.centroid
    n_mean = geom.unit(traj.normal[traj.spray > 0].mean(axis=0)) if np.any(traj.spray > 0) else np.array([0.0, 0.0, 1.0])
    h = np.array([n_mean[0], n_mean[1], 0.0])
    if geom.norm(h) < 0.2:          # 面朝上/朝下：法向没有水平分量，用 +X
        h = np.array([1.0, 0.0, 0.0])
    h = geom.unit(h)
    # 摆位目标：让**工件质心**落在基座系的舒适区 (0.5·reach, 0, 0.35·reach)。
    # 这个"舒适区"不是我编的，是 VS-087 的连杆比例决定的：J2 轴在基座系 z=395，
    # 臂展 445+430。目标点若落在 z≈390、水平 450 附近，J2/J3/J5 都有大量余量。
    #
    # 走过的两个弯路（都留在这里，免得下次再踩）：
    #  · "基座比工件低 0.5·reach"：目标点跑到 z=-207（低于 J2 轴 600mm），
    #    所有种子解全部超限（J3=-284°），且连杆从工件底下钻上去穿板（surface=-21mm）；
    #  · "与工件齐平（低 0.1·reach）"：目标在 z=188，仍低于 J2 轴，依然所有解超限。
    origin = c - h * (AUTO_BASE_FWD * reach)
    origin[2] = c[2] - AUTO_BASE_DOWN * reach
    x = geom.unit(c - origin)
    z = np.array([0.0, 0.0, 1.0])
    y = geom.unit(np.cross(z, x))
    x = geom.unit(np.cross(y, z))
    T = np.eye(4)
    T[:3, 0], T[:3, 1], T[:3, 2] = x, y, z
    T[:3, 3] = origin
    return T


def base_from_euler(x: float, y: float, z: float, rx: float, ry: float, rz: float) -> np.ndarray:
    """基座位姿：平移 + ZYX 欧拉（度）。与 DENSO P 型数据的姿态约定一致（R = Rz·Ry·Rx）。"""
    T = np.eye(4)
    Rz = geom.rot_axis(np.array([0.0, 0.0, 1.0]), math.radians(rz))[:3, :3]
    Ry = geom.rot_axis(np.array([0.0, 1.0, 0.0]), math.radians(ry))[:3, :3]
    Rx = geom.rot_axis(np.array([1.0, 0.0, 0.0]), math.radians(rx))[:3, :3]
    T[:3, :3] = Rz @ Ry @ Rx
    T[:3, 3] = (x, y, z)
    return T


# ---------------------------------------------------------------- 评估器

class Evaluator:
    """把「一个枪位姿态」变成「一组结论」。碰撞/IK/限位的口径全在这里，别处不重算。"""

    def __init__(self, collider: cloud_mod.CloudCollider, base: np.ndarray, params: Params):
        self.cloud = collider
        self.base = base
        self.R = base[:3, :3]
        self.t = base[:3, 3]
        self.p = params
        self.tool = robot.Tool(length=params.tool_len)
        self.sizing = robot.ArmSizing.of()
        # 段对只依赖段数（include_base 固定 ⇒ 拓扑固定），一次算好
        segs0 = robot.arm_segments(robot.fk(np.zeros(6), self.tool), self.sizing,
                                   include_base=params.include_base)
        self._pairs = robot.self_collision_pairs(segs0)

    def to_base_p(self, p: np.ndarray) -> np.ndarray:
        return self.R.T @ (p - self.t)

    def to_base_R(self, Rc: np.ndarray) -> np.ndarray:
        return self.R.T @ Rc

    def to_base_v(self, v: np.ndarray) -> np.ndarray:
        """方向向量：点云系 → 基座系（只转不平移）。to_cloud_R 的反向。"""
        return self.R.T @ v

    def to_cloud_p(self, p: np.ndarray) -> np.ndarray:
        return self.R @ p + self.t

    def to_cloud_R(self, v: np.ndarray) -> np.ndarray:
        """方向向量：基座系 → 点云系（只转不平移）。ref_x 这类来自 FK 的轴必须过这一步。"""
        return self.R @ v

    def eval_pose(self, face: np.ndarray, d: np.ndarray, ref_x: np.ndarray | None,
                  seed: np.ndarray | None, tilt_deg: float, az_deg: float = 0.0,
                  check_links: bool = True):
        """给定枪身朝外方向 d（= 倾斜后的法向 n'），求 IK + 碰撞 + 自碰撞。返回 PointResult 或 None。

        ⚠ **工具 +Z 是 −d，不是 d**。这一条搞反过（代价：一整轮全部结论作废）：
          · `d`     = n' = 表面外法向，从工件指向枪；
          · 工具坐标系 +Z = 从法兰指向喷嘴 = **喷射方向** = −n' = −d（CSV 里那一列 IJK 就是它）；
          · 所以 TCP 在 `face + standoff·d`，而工具 +Z 必须取 `−d`。
        反过来的表现极具欺骗性：位置全对、喷距守恒、可读性也好，但 IK 为了把工具轴翻到
        反方向，只能选"反腕 + 大臂下探"的解（θ2=154°、超限 19°、小臂从工件底下穿过去 24.7mm）。
        当时差点把它当成"工位摆位不合理"去改 auto_base —— 记住这个症状，别再走一遍。
        """
        z_tool = -geom.unit(d)
        p_c = face + d * self.p.standoff

        # ---- 不要求可达：不解 IK，只查枪体与点云的间隙。
        # 枪体位置由「TCP + 枪轴」解析确定，与机器人无关；所以这不是"少查了一项"，
        # 而是把与型号无关的那部分单独算准。与 C++ eval() 的 requireReach==0 分支逐行对应。
        if not self.p.require_reach:
            # 碰撞段一律建在**基座系**：后面的 to_cloud_p 只负责 基座系→点云系。
            # 直接喂点云系的 TCP/枪轴、再被 to_cloud_p 变换一次，等于整根枪被平移了
            # 一个基座偏移，会算出"枪体撞了 1512 点、最深 −21.2mm"这种假结论。
            tcp_b = self.to_base_p(p_c)
            segs = robot.tool_segments(tcp_b - self.to_base_v(z_tool) * self.tool.length,
                                       tcp_b, self.sizing)
            worst = math.inf
            worst_part = ""
            for s in segs:
                hit = self.cloud.segment(self.to_cloud_p(s.a), self.to_cloud_p(s.b),
                                        s.r, self.p.search, s.part)
                if hit.surface_dist < worst:
                    worst = hit.surface_dist
                    worst_part = s.part
            return PointResult(
                idx=-1,
                tilt_deg=tilt_deg,
                clearance=worst,
                min_part=worst_part,
                ik_ok=False,          # 如实：没解过 IK
                over_limit=False,
                margin_min=math.inf,  # 没量到（不是 0°）
                cond=0.0,
                theta=np.zeros(6),
                pos=p_c,
                gun=-d,
                spray_angle_to_normal=0.0,
                standoff_err=abs(geom.norm(p_c - face) - self.p.standoff),
                self_collision=math.inf,
                az_deg=az_deg,
            )

        Rc = geom.pose_from_z(z_tool, 0.0)[:3, :3]
        if ref_x is not None:
            # 用上一点的工具 x 轴做参考（投影到新切平面）——保证姿态「最小旋转」。
            # 否则枪轴接近竖直时参考向量会在 +Z/+X 之间跳变，IK 会跳到另一个解支。
            # ⚠ ref_x 来自 fk()，是**基座系**；z_tool 是**点云系**。必须先转坐标系：
            #   混着用等于拿两个坐标系的向量做正交分解，工具 x 轴每步被额外偏转一份
            #   「基座偏航」并逐点累加 → 手腕绕枪轴持续自旋（complex 实测 J6 −17.6°/点，
            #   整行拧到限位；tray 基座偏航为 0 所以完全正常 —— 这就是混用的铁证）。
            rx = self.to_cloud_R(ref_x)
            x0 = rx - z_tool * float(np.dot(rx, z_tool))
            if geom.norm(x0) > 1e-6:
                x0 = geom.unit(x0)
                y0 = geom.unit(np.cross(z_tool, x0))
                Rc = np.column_stack([x0, y0, z_tool])

        sol = robot.ik(self.to_base_p(p_c), self.to_base_R(Rc), self.tool, seed=seed)
        if not sol.ok:
            return None

        fk = robot.fk(sol.theta, self.tool)
        segs = robot.arm_segments(fk, self.sizing, include_base=self.p.include_base)
        worst = math.inf
        worst_part = ""
        for s in segs:
            if not check_links and s.group != "tool":
                continue
            a = self.to_cloud_p(s.a)
            b = self.to_cloud_p(s.b)
            hit = self.cloud.segment(a, b, s.r, self.p.search, s.part)
            if hit.surface_dist < worst:
                worst = hit.surface_dist
                worst_part = s.part

        self_gap = math.inf
        for ia, ib in self._pairs:
            sa, sb = segs[ia], segs[ib]
            d_axis = geom.capsule_capsule_dist(sa.a, sa.b, sb.a, sb.b)
            self_gap = min(self_gap, d_axis - sa.r - sb.r)

        return PointResult(
            idx=-1,
            tilt_deg=tilt_deg,
            clearance=worst,
            min_part=worst_part,
            ik_ok=True,
            over_limit=bool(sol.over),
            margin_min=float(np.min(sol.margin)),
            cond=sol.cond,
            theta=sol.theta,
            pos=p_c,
            gun=-d,
            spray_angle_to_normal=0.0,   # 由调用方按原始法向填
            standoff_err=abs(geom.norm(p_c - face) - self.p.standoff),
            self_collision=self_gap,
            az_deg=az_deg,
        )


# ---------------------------------------------------------------- 倾角需求

def _dir_of(n: np.ndarray, u_up: np.ndarray, axis: np.ndarray, tilt_deg: float, az_deg: float) -> np.ndarray:
    """把法向 n 绕「倾斜轴」转 tilt_deg，方位可绕 n 偏 az_deg。返回单位向量 d（枪身朝外）。

    tilt=0 时严格返回 n（不做任何多余旋转），保证「不倾斜的点位姿一字不变」。
    """
    if tilt_deg <= 1e-9 or geom.norm(axis) < 1e-6:
        return geom.unit(n)
    u = u_up
    if abs(az_deg) > 1e-9:
        u = geom.rot_axis(n, math.radians(az_deg))[:3, :3] @ u_up
        u = geom.unit(u - n * float(np.dot(u, n)))
    return geom.unit(n * math.cos(math.radians(tilt_deg)) + u * math.sin(math.radians(tilt_deg)))


def _need_tilt(ev: Evaluator, traj: traj_mod.Traj, i: int, seed, ref_x, check_links: bool):
    """求第 i 点「至少需要多大倾角 / 要不要换方位」，返回 (tilt, az, PointResult|None)。

    决策顺序 = 最小修正：
      ① 线性基准倾角（α=法向与 +Z 夹角，tilt0 = tilt_nonup · α/90，α=0 → 0°）
      ② 逐步 +5° 到上限
      ③ 兜底：把「上翻方向」左右偏 30°/60° 再试（有些内壁往上扳会顶到对面的壁）

    两处刻意的加速（不改结论，只省时间，实测把单点最坏 2.8s 压到 0.3s）：
      * **早停**：连续 4 个候选的 clearance 变化 < 0.5mm —— 说明这个方向的倾角对碰撞毫无影响，
        再加也是白加（典型是"碰撞来自基座/连杆布局，不是枪体"，那是布局问题，倾角救不了）。
      * **退化直达**：法向与 +Z 平行时没有切向可扳，候选只剩 1 个，
        否则会拿同一个位姿算 13 遍（每一遍都 60~80ms，纯浪费）。
    """
    p = ev.p
    n = traj.normal[i]
    face = traj.pos[i] + traj.gun[i] * p.standoff
    up = np.array([0.0, 0.0, 1.0])

    # 倾角需求按法向**线性**给（不是"非朝上就一律 20°"的阶梯）：
    #   α = 法向与 +Z 的夹角（0°=正朝上，90°=竖直侧壁，180°=完全倒挂）
    #   tilt0 = tilt_nonup_deg · (α/90°)，封顶 tilt_max_deg
    # 90° 作满量程锚点是物理界限（"面立起来了"），更倒挂的面继续线性增长。
    # 阶梯在 45° 边界上会突跳 20°，落到轨迹上就是相邻点倾角跳变，实机过不了。
    alpha_deg = math.degrees(math.acos(float(np.clip(n[2], -1.0, 1.0))))
    tilt0 = min(p.tilt_nonup_deg * (alpha_deg / 90.0), p.tilt_max_deg)
    face_down = float(n[2]) < -0.7          # 朝下的悬垂面：往上扳会顶穿工件，不做线性外推
    if face_down:
        tilt0 = 0.0

    axis, u_up, _ = geom.tilt_axis_and_dir(n, up)
    degenerate = geom.norm(u_up) < 1e-6      # 夹角 < 5° 时 geom 已直接返回零向量
    if degenerate:
        axis = np.zeros(3)
        u_up = np.zeros(3)
        tilt0 = 0.0

    # 第一候选 = 线性基准原样（贴合需求，不无谓加大）；之后走**固定 5° 网格**。
    # 不能写成 tilt0 + k·step：tilt0 随法向连续变化，候选相位会跟着漂移，
    # 同一个物理点可能从"刚好过关的 40°"滑到"刚好撞上的 36.48°"，纯数值抖动。
    # 实测 curved idx 2990：网格 40° 自碰撞间隙 +1.18mm，漂移 36.48° 为 −1.16mm。
    cand: list[tuple[float, float]] = [(tilt0, 0.0)]
    if not degenerate:
        g = math.ceil(tilt0 / p.tilt_step_deg - 1e-9) * p.tilt_step_deg
        while g <= tilt0 + 1e-9:
            g += p.tilt_step_deg
        while g <= p.tilt_max_deg + 1e-9:
            cand.append((g, 0.0))
            g += p.tilt_step_deg

    best_any: PointResult | None = None
    best_any_key = (-math.inf, -math.inf)
    flat_run = 0
    prev: tuple[float, float, float] | None = None

    def _try(tilt_deg: float, az: float):
        nonlocal best_any, best_any_key, flat_run, prev
        d = _dir_of(n, u_up, axis, tilt_deg, az)
        pr = ev.eval_pose(face, d, ref_x, seed, tilt_deg, az, check_links=check_links)
        if pr is None:
            return None
        pr.spray_angle_to_normal = geom.angle_between(d, n)
        if (pr.clearance >= p.safety and pr.self_collision >= 0.0
                and not pr.over_limit and pr.margin_min >= p.joint_margin_min):
            return pr
        key = (pr.clearance if pr.self_collision >= 0 else -math.inf, pr.margin_min)
        if key > best_any_key:
            best_any_key = key
            best_any = pr
        # 早停计数：三个指标都不动才算"倾角对这条方向毫无影响"。
        # 只比 clearance 会误停——间隙常常早早饱和到 search 上限，而自碰撞/
        # 关节余量仍在随倾角变化，停了就漏掉后面能过关的档位。
        cur = (pr.clearance, pr.self_collision, pr.margin_min)
        if prev is not None and all(abs(a - b) < 0.5 for a, b in zip(cur, prev)):
            flat_run += 1
        else:
            flat_run = 0
        prev = cur
        return None

    for tilt_deg, az in cand:
        hit = _try(tilt_deg, az)
        if hit is not None:
            return tilt_deg, az, hit
        if flat_run >= 4:      # 倾角无效，别再徒劳加大
            break

    # ---- 兜底：换倾斜方位（只在"往上扳"整条路都无效时才展开）
    if p.az_fallback and not degenerate and flat_run < 4:
        grid0 = math.ceil(tilt0 / p.tilt_step_deg - 1e-9) * p.tilt_step_deg
        for az in (30.0, -30.0, 60.0, -60.0):
            for tt in (max(grid0, 10.0), max(grid0, 25.0), max(grid0, 40.0), max(grid0, 55.0)):
                if tt > p.tilt_max_deg + 1e-9:
                    continue
                hit = _try(tt, az)
                if hit is not None:
                    return tt, az, hit

    if best_any is not None:
        return best_any.tilt_deg, best_any.az_deg, best_any
    return tilt0, 0.0, None


def _face_groups(traj: traj_mod.Traj, rows: list[tuple[int, int]],
                 params: Params) -> list[tuple[int, int]]:
    """把相邻、法向一致、共面的行归成一个「面组」，返回面组的行区间 [r0, r1]。

    为什么需要它：工艺要求**同一个面上喷涂过程中喷角不变**。倾角如果逐点自由，
    就一定会随局部位姿抖动（实测一笔 77 点里能从 30° 走到 55° 再落回 30°，
    相邻点 TCP 沿球面来回摆几十毫米）。把姿态常量抬到「面」这一级，才有可能不变。

    判据（与 C++ 同口径）：
      · 行法向夹角 ≤ 2°      —— 同一张面（曲面会被自然切成多段，那是应该变的）
      · 沿法向偏移 ≤ 2mm     —— 排除"平行但不共面"的两张面（如壁厚两侧）
    """
    normals, mids = [], []
    for a, b in rows:
        n = traj.normal[a:b + 1].mean(axis=0)
        u = geom.norm(n)
        normals.append(geom.unit(n) if u > 1e-12 else np.zeros(3))
        m = (a + b) // 2
        mids.append(traj.pos[m] + traj.gun[m] * params.standoff)
    out: list[tuple[int, int]] = []
    cos_lim = math.cos(math.radians(2.0))
    for r in range(len(rows)):
        same = False
        if r > 0:
            same = (float(np.dot(normals[r], normals[r - 1])) >= cos_lim
                    and abs(float(np.dot(normals[r], mids[r] - mids[r - 1]))) <= 2.0)
        if same:
            out[-1] = (out[-1][0], r)
        else:
            out.append((r, r))
    return out


def _pose_key(pr: PointResult) -> tuple[float, float, float, float, float, float]:
    """硬约束二选一的排序键（越大越好）：IK 成功 > 自碰撞通过 > 关节不超限 >
    自间隙大 > 间隙大 > 关节余量大。

    碰撞（P0）永远压过姿态连续性（P2）；关节真正超限（over_limit）也单独判负，
    不能因为另一个位姿间隙大就选一个关节超限的位姿。
    """
    return (1.0 if pr.ik_ok else 0.0,
            1.0 if pr.self_collision >= 0.0 else 0.0,
            1.0 if not pr.over_limit else 0.0,
            pr.self_collision,
            pr.clearance,
            pr.margin_min)


def _hard_ok(pr: PointResult, params: Params) -> bool:
    """硬约束：IK 成功、不自碰撞、间隙不低于绝对安全距离、关节不超限。

    margin_min < joint_margin_min 只是 P1 软告警，不在这里判死；over_limit 才是硬约束。
    require_reach=False 时只剩「间隙不低于绝对安全距离」，而那时的 clearance 就是枪体间隙
    （没建连杆碰撞体），所以"枪不许戳进工件"这一条在两种口径下都成立。
    """
    if not (pr.clearance >= params.hard_clearance):
        return False
    if not params.require_reach:
        return True
    return (pr.ik_ok
            and pr.self_collision >= 0.0
            and not pr.over_limit)


# ---------------------------------------------------------------- L0 参数自检

def precheck(traj: traj_mod.Traj, collider: cloud_mod.CloudCollider, params: Params) -> list[dict]:
    """**L0：输入层自检**。这一层不通过，后面所有几何结论都是废的。

    两条判据都是硬的，来自几何而不是经验：

    ① **喷距 + 工具长必须留在臂展内**。法兰到面点的距离是 `standoff + tool_len`
       （工具 +Z 从法兰指向喷嘴=喷射方向 −n，所以法兰在面点**外侧** standoff+tool_len 处）。
       两者相加超过 reach 的 90% 时，末端姿态余量会急剧变小 —— 这是可达性问题，
       在后面的 IK 结果里表现为可达率下降或关节贴限位。

       ⚠ 这条曾经写反过：一度以为"tool_len > standoff ⇒ 法兰陷进工件"，并据此把
       VS-087 的 150mm 判成"参数不自洽"。真正的错误是工具 +Z 取成了 +n 而不是 −n
       （见 Evaluator.eval_pose 的注释）。工具越长，法兰离工件越远，**越安全**。

    ② **轨迹反推的喷距与参数一致**。轨迹点 = face + standoff·n；把每个点投到点云上量出来，
       就是生成这条轨迹时**实际用的** standoff。和参数差太多说明要么 standoff 填错、
       要么轨迹不是这套点云/这份参数生成的。
    """
    out: list[dict] = []

    # ① 喷距 + 工具长 vs 臂展
    #    只在「要求可达」时才有意义：不要求可达时臂展根本不参与判定。
    reach_need = params.standoff + params.tool_len
    limit = robot.VS087_REACH * 0.90
    if params.require_reach and reach_need > limit:
        out.append({
            "index": -1, "rule": "PARAM_REACH_MARGIN",
            "value": reach_need, "limit": limit, "severity": "P1",
            "note": (f"喷距 {params.standoff:.0f}mm + 工具长 {params.tool_len:.0f}mm = {reach_need:.0f}mm，"
                     f"超过臂展 90%（{limit:.0f}mm / reach {robot.VS087_REACH:.0f}mm）："
                     f"末端姿态余量会急剧变小，远端可能解不出逆解。"),
        })

    # ② 轨迹反推喷距（抽样 60 个喷涂点）
    #
    # ⚠ 这里必须用**沿枪轴的投影**，不能用最近点距离。
    #   最近点距离对斜面/曲面天然偏小（最近点根本不在法向上），实测 complex 的最近点
    #   中位数只有 72mm，于是把「喷距=100」误判成参数不一致；换成沿轴投影后
    #   全部样本都落在 98.5~100.5mm，说明 standoff=100 本来就是对的。
    #   几何定义也是沿轴投影：face = TCP + standoff·(I,J,K)。
    idx = np.where(traj.spray > 0)[0]
    if len(idx) > 0:
        sel = idx[np.linspace(0, len(idx) - 1, min(60, len(idx))).astype(int)]
        tube = max(3.0 * getattr(collider, "spacing", 1.0), 5.0)   # 管道半径：随点云间距自适应
        dists = []
        for i in sel:
            t = _ray_standoff(collider, traj.pos[i], traj.gun[i], tube)
            if t is not None:
                dists.append(t)
        if dists:
            med = float(np.median(dists))
            if abs(med - params.standoff) > max(8.0, 0.15 * params.standoff):
                out.append({
                    "index": -1, "rule": "PARAM_STANDOFF_MISMATCH",
                    "value": med, "limit": params.standoff, "severity": "P1",
                    "note": (f"轨迹点沿枪轴到点云的距离中位数 {med:.1f}mm，与参数 standoff={params.standoff:.0f}mm "
                             f"相差 {abs(med - params.standoff):.0f}mm —— 轨迹不是按这个喷距生成的，"
                             f"或点云与轨迹不是同一套坐标。"),
                })
    return out


def _ray_standoff(collider: cloud_mod.CloudCollider, p: np.ndarray, gun: np.ndarray,
                  tube: float, max_t: float = 300.0) -> float | None:
    """沿枪轴 gun 从轨迹点打到点云的距离（管道内点的投影中位数）。

    gun 指向工件（= −n）。取离轴不超过 tube 的点，求它们在 gun 方向上的正投影，
    取中位数 —— 比最近点稳健得多，且与 face/TCP 的几何定义一致。
    """
    rel = collider.pts - p
    t = rel @ gun
    perp2 = np.einsum("ij,ij->i", rel, rel) - t * t
    m = (t > 0) & (t < max_t) & (perp2 < tube * tube)
    if int(np.count_nonzero(m)) < 3:
        return None
    return float(np.median(t[m]))


# ---------------------------------------------------------------- 主流程

def plan(traj: traj_mod.Traj, collider: cloud_mod.CloudCollider, base: np.ndarray,
         params: Params, verbose: bool = True) -> Result:
    ev = Evaluator(collider, base, params)
    n_pts = len(traj)
    rows = traj.rows()

    l0 = precheck(traj, collider, params)
    if verbose:
        for v in l0:
            print(f"  [L0][{v['severity']}] {v['rule']}: {v['note']}")

    # ---- 阶段 1：每点的倾角需求（带 IK 种子连续性）
    need: list[float] = [0.0] * n_pts
    az_of: list[float] = [0.0] * n_pts
    cache: dict[int, PointResult] = {}
    prev_theta = None
    prev_x = None
    for i in range(n_pts):
        check_links = (i % max(1, params.link_stride) == 0)
        tilt, az, pr = _need_tilt(ev, traj, i, prev_theta, prev_x, check_links)
        need[i] = tilt
        az_of[i] = az
        if pr is not None:
            pr.idx = i
            cache[i] = pr
            prev_theta = pr.theta
            fk = robot.fk(pr.theta, ev.tool)
            prev_x = fk.tcp[:3, :3] @ np.array([1.0, 0.0, 0.0])
        if verbose and n_pts > 600 and (i + 1) % 600 == 0:
            print(f"  [需求] {i + 1}/{n_pts} ...", flush=True)

    # ---- 阶段 2：面级姿态常量（倾角 + 方位）
    # 同一个面组只有一个 (倾角, 方位)：姿态因此只可能在换面时改变。
    # 面内个别点在该常量下过不了硬约束，仍由阶段 3 退回本点最优（并记 reverted）——
    # 那是"物理上必须偏"的点，不是决策抖动。
    #
    # 倾角取「面内各点需求的最大值」：保守但零额外开销。
    # C++ 版（生产路径）在此之上做了抽样打分择优（同一面在别的倾角上碰撞余量可能更好），
    # 两份实现的选择规则不完全一致，跨实现对数时请以「姿态是否面级常量」为准，
    # 不要逐点比倾角数值。
    groups = _face_groups(traj, rows, params)
    final = list(need)
    az_row: dict[tuple[int, int], float] = {}
    up = np.array([0.0, 0.0, 1.0])
    for r0, r1 in groups:
        pts = [i for r in range(r0, r1 + 1) for i in range(rows[r][0], rows[r][1] + 1)
               if traj.spray[i] == 1]
        if not pts:
            continue
        n_mean = geom.unit(traj.normal[pts].mean(axis=0))
        alpha = math.degrees(math.acos(float(np.clip(n_mean[2], -1.0, 1.0))))
        tilt0 = min(params.tilt_nonup_deg * (alpha / 90.0), params.tilt_max_deg)
        if float(n_mean[2]) < -0.7:
            tilt0 = 0.0
        _, u_up, _ = geom.tilt_axis_and_dir(n_mean, up)
        if geom.norm(u_up) < 1e-6:
            tilt0 = 0.0
        tilt = max([tilt0] + [need[i] for i in pts])
        vals = [az_of[i] for i in pts if abs(az_of[i]) > 1e-9]
        az = 0.0
        if vals:
            uniq, cnt = np.unique(np.round(vals, 3), return_counts=True)
            az = float(uniq[int(np.argmax(cnt))])
        for r in range(r0, r1 + 1):
            a, b = rows[r]
            az_row[(a, b)] = az
            for i in range(a, b + 1):
                if traj.spray[i] == 1:
                    final[i] = tilt
    row_of: dict[int, tuple[int, int]] = {}
    for a, b in rows:
        for i in range(a, b + 1):
            row_of[i] = (a, b)

    # ---- 阶段 3：按 final 重新求解 + 归集违规
    results: list[PointResult] = []
    violations: list[dict] = []
    standoff_max_err = 0.0
    prev_theta = None
    prev_x = None
    recomputed = 0
    reverted = 0
    for i in range(n_pts):
        n = traj.normal[i]
        face = traj.pos[i] + traj.gun[i] * params.standoff
        tilt = final[i]
        rk = row_of.get(i)
        az = az_row.get(rk, 0.0) if rk else 0.0

        axis, u_up, _ = geom.tilt_axis_and_dir(n, np.array([0.0, 0.0, 1.0]))
        if geom.norm(u_up) < 1e-6:
            axis = np.zeros(3)
            u_up = np.zeros(3)

        cached = cache.get(i)
        reuse = (cached is not None
                 and abs(cached.tilt_deg - tilt) < 1e-6
                 and abs(cached.az_deg - az) < 1e-6)
        if reuse:
            pr = cached
            pr.idx = i
        else:
            d = _dir_of(n, u_up, axis, tilt, az)
            check_links = (i % max(1, params.link_stride) == 0)
            pr = ev.eval_pose(face, d, prev_x, prev_theta, tilt, az, check_links=check_links)
            recomputed += 1
            if pr is None:
                results.append(PointResult(
                    idx=i, tilt_deg=tilt, clearance=-math.inf, min_part="—", ik_ok=False,
                    over_limit=False, margin_min=0.0, cond=math.inf,
                    theta=np.full(6, np.nan), pos=traj.pos[i], gun=traj.gun[i],
                    spray_angle_to_normal=0.0, standoff_err=0.0, violations=["IK_FAIL"],
                ))
                violations.append({"index": i, "rule": "IK_FAIL", "value": None, "limit": None,
                                   "severity": "P0", "note": "所有候选倾角都解不出逆解（不可达）"})
                continue
            pr.idx = i
            # ---- 平滑回退：只回退硬约束，不再"差一点就回退"
            # 用户现场口径：有凸点没关系，只要不碰撞就以平滑为主。
            #   · 平滑位姿满足硬约束（IK / 自碰撞 / 间隙≥hard_clearance / 关节不超限）
            #     → 保留平滑，哪怕间隙比本点最优小几毫米；
            #   · 只有硬约束被破坏，且本点最优更安全时，才退回本点最优。
            # 实测（complex idx 3163）：本点最优 az 下间隙 -2.35mm，被整行统一成 az=30° 后
            # 变成 -5.81mm 真穿透 —— 这种硬约束被破坏时必须回退。
            if cached is not None and not _hard_ok(pr, params):
                if _hard_ok(cached, params) or _pose_key(pr) < _pose_key(cached):
                    pr = cached
                    pr.idx = i
                    pr.reverted = True
                    reverted += 1
        pr.spray_angle_to_normal = geom.angle_between(-pr.gun, n)
        standoff_max_err = max(standoff_max_err, pr.standoff_err)

        # ---- 违规归集
        # require_reach=False：SELF_COLLISION / JOINT_LIMIT 这两位不置。它们都是「这台机器人 +
        # 这条站位」的结论，点云本身没有对错。几何层（穿透 / 间隙 / 工艺喷角）照旧。
        _kin = params.require_reach
        if pr.clearance < 0:
            pr.violations.append("COLLISION")
            violations.append({"index": i, "rule": "COLLISION", "value": pr.clearance, "limit": 0.0,
                               "severity": "P0", "note": f"{pr.min_part} 穿透点云 {abs(pr.clearance):.1f}mm"})
        elif pr.clearance < params.safety:
            pr.violations.append("CLEARANCE")
            violations.append({"index": i, "rule": "CLEARANCE", "value": pr.clearance, "limit": params.safety,
                               "severity": "P1", "note": f"{pr.min_part} 表面距离 {pr.clearance:.1f}mm < {params.safety:.0f}mm"})
        if _kin and pr.self_collision < 0:
            pr.violations.append("SELF_COLLISION")
            violations.append({"index": i, "rule": "SELF_COLLISION", "value": pr.self_collision, "limit": 0.0,
                               "severity": "P0", "note": "连杆自碰撞"})
        if _kin and (pr.over_limit or pr.margin_min < params.joint_margin_min):
            pr.violations.append("JOINT_LIMIT")
            violations.append({"index": i, "rule": "JOINT_LIMIT", "value": pr.margin_min,
                               "limit": params.joint_margin_min, "severity": "P1",
                               "note": f"关节限位余量 {pr.margin_min:.1f}°"})
        if pr.spray_angle_to_normal > params.max_spray_angle_deg + 1e-9:
            pr.violations.append("SPRAY_ANGLE")
            violations.append({"index": i, "rule": "SPRAY_ANGLE", "value": pr.spray_angle_to_normal,
                               "limit": params.max_spray_angle_deg, "severity": "P1",
                               "note": f"喷角 {pr.spray_angle_to_normal:.1f}° 超工艺上限 {params.max_spray_angle_deg:.0f}°"})
        if pr.cond > params.cond_max:
            pr.violations.append("NEAR_SINGULAR")
            violations.append({"index": i, "rule": "NEAR_SINGULAR", "value": pr.cond, "limit": params.cond_max,
                               "severity": "P2", "note": "接近奇异位形"})

        results.append(pr)
        prev_theta = pr.theta
        fk = robot.fk(pr.theta, ev.tool)
        prev_x = fk.tcp[:3, :3] @ np.array([1.0, 0.0, 0.0])

    # ---- 过渡点抬升修正（spray=0，不喷涂 ⇒ 位置可改，只求不撞）
    lifted = 0
    if params.fix_transit:
        keep = []
        for v in violations:
            keep.append(v)
        for pr in results:
            if traj.spray[pr.idx] == 0 and np.isfinite(pr.clearance) and pr.clearance < params.safety:
                lift = min(params.transit_lift_max, params.safety - pr.clearance + 5.0)
                pr.pos = pr.pos + np.array([0.0, 0.0, lift])
                pr.lifted = lift
                pr.clearance += lift
                pr.violations = [v for v in pr.violations if v not in ("COLLISION", "CLEARANCE")]
                keep = [v for v in keep if not (v["index"] == pr.idx and v["rule"] in ("COLLISION", "CLEARANCE"))]
                lifted += 1
        violations = keep

    # ---- 姿态跳变：平滑回退的代价，必须报出来
    # 回退是为了不撞（P0），但它会让个别点偏离限幅、姿态突变。
    # 阈值取 3 倍限幅与 10° 的较大者：限幅内的正常波动不该刷告警，
    # 真正会让手腕"甩一下"的跳变才报。级别 P2 —— 它不是碰撞，是速度/加速度风险。
    if params.tilt_rate_deg_per_pt > 0 and len(results) == n_pts:
        thr = max(3.0 * params.tilt_rate_deg_per_pt, 10.0)
        for a, b in rows:
            for i in range(a + 1, b + 1):
                dlt = abs(results[i].spray_angle_to_normal - results[i - 1].spray_angle_to_normal)
                if dlt > thr:
                    results[i].violations.append("TILT_JUMP")
                    violations.append({"index": i, "rule": "TILT_JUMP", "value": dlt,
                                       "limit": params.tilt_rate_deg_per_pt, "severity": "P2",
                                       "note": f"倾角相对上一点跳变 {dlt:.1f}°（限幅 {params.tilt_rate_deg_per_pt:.1f}°/点）"})

    # L0 自检结论并进违规列表（P0 的 L0 问题一定让整体 REJECT）
    violations = l0 + violations

    angles = [p.spray_angle_to_normal for p in results]
    # 不要求可达时根本没解过 IK（ik_ok 恒 False）。按「默认全部通过」记成 n，
    # 但这只是记账口径 —— 报告/界面必须显著标成「未检查」，别把 100% 当成真测出来的。
    n_ok = n_pts if not params.require_reach else sum(1 for p in results if p.ik_ok)
    stats = {
        "points": n_pts,
        "spray_points": int(np.count_nonzero(traj.spray)),
        "reach_checked": bool(params.require_reach),
        "reachable_points": n_ok,
        "reach_fraction": (n_ok / n_pts) if n_pts else 0.0,
        "tilted_points": int(sum(1 for a in angles if a > 0.5)),
        "tilt_mean_deg": float(np.mean(angles)) if angles else 0.0,
        "tilt_p95_deg": float(np.percentile(angles, 95)) if angles else 0.0,
        "tilt_max_deg": float(np.max(angles)) if angles else 0.0,
        "standoff_max_err_mm": standoff_max_err,
        "min_clearance_mm": float(min((p.clearance for p in results), default=math.inf)),
        "transit_lifted": lifted,
        "recomputed_after_smoothing": recomputed,
        "pose_reverted": reverted,
        "violations": len(violations),
        "violations_by_rule": _count_rules(violations),
        "l0_issues": len(l0),
    }
    ci = {
        "points": len(collider.pts),
        "raw_points": getattr(collider, "raw_count", len(collider.pts)),
        "spacing_mm": collider.spacing,
        "aabb_min": collider.min.tolist(),
        "aabb_max": collider.max.tolist(),
        "centroid": collider.centroid.tolist(),
    }
    if verbose:
        reach_txt = (f"可达 {n_ok} ({stats['reach_fraction'] * 100:.1f}%)" if params.require_reach
                     else "可达 未检查（仅几何口径）")
        print(f"[plan] {n_pts} 点 · {reach_txt} · "
              f"倾斜 {stats['tilted_points']} 点 (均值 {stats['tilt_mean_deg']:.1f}°, 最大 {stats['tilt_max_deg']:.1f}°) · "
              f"违规 {len(violations)} · 最小间隙 {stats['min_clearance_mm']:.1f}mm")
        if standoff_max_err > 1e-3:
            print(f"  ⚠ 喷距保持失败：最大偏差 {standoff_max_err:.4f}mm")
    return Result(points=results, base=base, stats=stats, violations=violations, params=params, cloud_info=ci)


def _count_rules(vios: list[dict]) -> dict:
    out: dict[str, int] = {}
    for v in vios:
        out[v["rule"]] = out.get(v["rule"], 0) + 1
    return out


def apply_result(traj: traj_mod.Traj, res: Result) -> traj_mod.Traj:
    """把评估后的位姿写回轨迹（位置与 IJK 都换成倾斜后的）。"""
    pos = traj.pos.copy()
    gun = traj.gun.copy()
    for pr in res.points:
        if pr.idx < 0 or not np.all(np.isfinite(pr.theta)):
            continue
        pos[pr.idx] = pr.pos
        gun[pr.idx] = pr.gun
    return traj_mod.Traj(pos=pos, gun=gun, speed=traj.speed, spray=traj.spray,
                         extra=traj.extra, header=traj.header)
