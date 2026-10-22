# -*- coding: utf-8 -*-
"""关节跳变统计：读 spray_safety 输出的 checked CSV，统计相邻点单关节最大跳变。

用途：验证「IK 解支连续性」与「工具滚转不丢参考」两类改动，是否真的消除了
手腕换支 / 绕枪轴自旋造成的隐形跳变。只读文件，不改任何数据。

用法：
    python -m robot_safety.check_joints a_checked.csv [b_checked.csv ...]

口径：
  · 相邻 = CSV 里相邻两行（轨迹按输出顺序连续运动）
  · 每个相邻对取 6 个关节环上最短角差的最大值，记为该步的 jump
    （CSV 关节按 (-180,180] 输出，跨 ±180 的连续转动必须按环上角差算，
      否则 360° 的假跳变会淹没真实数据）
  · 按 spray 标志分成 spray→spray（喷涂段内，必须连续）与其余（过渡段，允许大幅重定向）

判据（实测标定，2026-10-22）：
  · 正常：喷涂段内均值 0.4~0.9°/点，p95 < 1.5°
  · 异常：均值 > 2° 或出现 180° 附近的整步跳变 → 查工具滚转参考（基座系/点云系混用）
"""
import sys
import numpy as np

JOINTS = ["J1", "J2", "J3", "J4", "J5", "J6"]


def load(path):
    with open(path, "r", encoding="utf-8-sig") as f:
        head = f.readline().strip().split(",")
    keys = [j + "_deg" for j in JOINTS]
    idx = {name: head.index(name) for name in keys + ["spray"]}
    cols = [idx[j + "_deg"] for j in JOINTS]
    # CSV 里有中文字符串列（worst_part/link_part），只能按列号取，不能整表 loadtxt
    th = np.loadtxt(path, delimiter=",", skiprows=1, usecols=cols, ndmin=2)
    spray = np.loadtxt(path, delimiter=",", skiprows=1, usecols=[idx["spray"]],
                       dtype=int, ndmin=1)
    return th, spray


def stats(path):
    th, spray = load(path)
    # 关节角按 (-180,180] 输出，跨 ±180 边界的连续转动会显示成 ~360° 跳变。
    # 物理上的跳变量取环上的最短角差：|((Δ+180) mod 360) - 180|
    raw = np.diff(th, axis=0)
    d = np.abs((raw + 180.0) % 360.0 - 180.0)
    step = d.max(axis=1)
    worst_joint = d.argmax(axis=1)
    inner = (spray[:-1] == 1) & (spray[1:] == 1)

    def rep(mask, label):
        s = step[mask]
        if s.size == 0:
            return f"  {label}: 无"
        return (f"  {label}: n={s.size} 均值={s.mean():.2f}° p95={np.percentile(s,95):.2f}° "
                f"p99={np.percentile(s,99):.2f}° 最大={s.max():.2f}° "
                f">5°={int((s>5).sum())} >10°={int((s>10).sum())} >20°={int((s>20).sum())}")

    print(f"\n{path}")
    print(f"  点数={len(th)}")
    print(rep(np.ones(len(step), bool), "全部相邻"))
    print(rep(inner, "喷涂段内"))

    # 喷涂段内最差的 8 步
    order = np.argsort(-np.where(inner, step, -1))[:8]
    bad = [i for i in order if inner[i] and step[i] > 5]
    if bad:
        print("  喷涂段内 >5° 的最差步骤：")
        for i in bad:
            print(f"    idx {i}->{i+1}: {JOINTS[worst_joint[i]]} 跳 {step[i]:.2f}°  "
                  + " ".join(f"{JOINTS[k]}={th[i+1,k]-th[i,k]:+.2f}" for k in range(6)
                             if abs(th[i+1, k] - th[i, k]) > 2))
    return step, inner


if __name__ == "__main__":
    for p in sys.argv[1:]:
        stats(p)
