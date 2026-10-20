"""机器人安全层（robot-safety 分支）。

只读 master 产出的轨迹 CSV，产出「校验后的轨迹 + 违规报告」，不改 master 的任何代码。
入口：`python -m robot_safety.run_safety --help`
"""

__all__ = ["geom", "robot", "cloud", "traj", "planner", "report"]
