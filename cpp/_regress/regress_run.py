# -*- coding: utf-8 -*-
"""喷涂轨迹算法全量回归

用法:
    python regress_run.py                 # 与基线比对，回归则非零退出
    python regress_run.py --update        # 刷新基线（只在确认改动符合预期时用）
    python regress_run.py --only box      # 只跑名字含 box 的用例
    python regress_run.py --no-clear      # 跳过安全距离复算（快）

每个用例测四项指标：
    faces   面级清单面数（几何理解是否变化）
    rows    轨迹行数（布行策略是否变化）
    points  轨迹点数（采样密度是否变化）
    minP    枪位到工件点云的最小距离（点级）
    minS    轨迹折线中段到工件的最小距离（穿模的直接指标）

后两项直接调用 traj_clear_scan 的测距内核，与人工体检口径完全一致。
"""
import os
import re
import sys
import json
import time
import subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from traj_clear_scan import read_ply, read_traj, Grid, pt_seg_dist   # noqa: E402

EXE = os.environ.get("SPRAY_EXE") or os.path.join(
    HERE, "..", "SprayTrajCore", "build", "bin", "Release", "spray_driver.exe")
BASE = os.path.join(HERE, "regress_baseline.json")
OUT = os.environ.get("REGRESS_OUT") or os.path.join(HERE, "regress_out")

# (名称, 点云, 参数, 安全半径mm, 说明)
CASES = [
    ("box",       "box.ply",       "-flat 40",                                              20.0, "开口盒 单层点云"),
    ("bin",       "bin.ply",       "-flat 40",                                              20.0, "开口槽 曲面+环"),
    ("cao",       "cao.ply",       "-flat 40",                                              20.0, "槽形件"),
    ("complex",   "complex.ply",   "-flat 40",                                              20.0, "高复杂度件"),
    ("curved",    "curved.ply",    "-flat 40",                                              20.0, "曲面件"),
    ("longstripe","longstripe.ply","-flat 40",                                              20.0, "长条面"),
    ("zhe",       "zhe.ply",       "-flat 40",                                              20.0, "折弯件"),
    ("tray",      "tray.ply",      "-flat 40",                                              20.0, "托盘"),
    ("big",       "big.ply",       "-flat 40",                                              20.0, "大件"),
    ("box_user",  "box.ply",       "-flat 40 -row 80 -pt 4 -so 100 -gc 20 -ee -1 -ea -1 -st 4 -sp 200",
                                                                                            20.0, "用户参数 行距80/离面100/双面4（穿模复现用例）"),
    ("rim",       "g_rim.ply",     "-flat 40",                                              20.0, "薄壁盒 带顶沿 单层点"),
    ("rim_st4",   "g_rim.ply",     "-flat 40 -st 4",                                        20.0, "薄壁盒 带顶沿 + 双面4"),
    ("norim",     "g_norim.ply",   "-flat 40",                                              20.0, "薄壁盒 无顶沿 真双层"),
]

# 允许波动
TOL_POINTS_PCT = 1.0      # 点数 ±1%
TOL_CLEAR_MM = 0.5        # 最小净空允许下降 0.5mm


def run_case(name, ply, args, clear):
    ply_path = os.path.join(HERE, ply)
    if not os.path.exists(ply_path):
        return {"skip": "点云不存在"}
    csv = os.path.join(OUT, name + ".csv")
    oply = os.path.join(OUT, name + ".ply")
    t0 = time.time()
    cp = subprocess.run([EXE, ply_path, csv, oply] + args.split(),
                        capture_output=True, text=True, errors="ignore", cwd=HERE)
    dt = time.time() - t0
    log = cp.stdout + cp.stderr

    if cp.returncode != 0:
        return {"fail": "退出码 %d" % cp.returncode, "log": log[-1500:], "sec": dt}
    m = re.search(r"OK: rows=(\d+) points=(\d+) rms=([\d.]+)", log)
    if not m:
        return {"fail": "无 OK 行", "log": log[-1500:], "sec": dt}
    rows, pts, rms = int(m.group(1)), int(m.group(2)), float(m.group(3))
    mf = re.search(r"窗口内 (\d+) 面", log)
    faces = int(mf.group(1)) if mf else -1
    n_lift = log.count("抬枪绕行")

    r = {"faces": faces, "rows": rows, "points": pts, "rms": round(rms, 3),
         "sec": round(dt, 1), "n_lift": n_lift}

    if clear > 0:
        Q = read_ply(ply_path)
        G, sp = read_traj(csv)
        grid = Grid(Q)
        dp = np.empty(len(G))
        for k in range(len(G)):
            c = grid.cand(G[k] - clear, G[k] + clear)
            dp[k] = 1e9 if c is None else float(np.linalg.norm(Q[c] - G[k], axis=1).min())
        ds = np.full(len(G) - 1, 1e9)
        for k in range(len(G) - 1):
            A, B = G[k], G[k + 1]
            c = grid.cand(np.minimum(A, B) - clear, np.maximum(A, B) + clear)
            if c is not None:
                ds[k] = float(pt_seg_dist(Q[c], A, B).min())
        r["minP"] = round(float(dp.min()), 2)
        r["minS"] = round(float(ds.min()), 2)
        r["badP"] = int((dp < clear).sum())
        r["badS"] = int((ds < clear).sum())
        # 查询半径外没有候选点 = 观测不到更近的东西，等价于"至少 clear"
        if r["minP"] >= 1e8: r["minP"] = clear
        if r["minS"] >= 1e8: r["minS"] = clear
    return r


def cmp(old, new):
    """返回 (是否回归, 说明列表)"""
    probs, notes = [], []
    if "fail" in new or "skip" in new:
        return True, [new.get("fail") or new.get("skip")]
    if "fail" in old:
        return False, ["基线为失败态，跳过比对"]
    for k in ("faces", "rows"):
        if old.get(k) != new.get(k):
            probs.append("%s %s→%s" % (k, old.get(k), new.get(k)))
    op, np_ = old.get("points", 0), new.get("points", 0)
    if op and abs(np_ - op) > op * TOL_POINTS_PCT / 100.0:
        probs.append("points %d→%d (%+.1f%%)" % (op, np_, 100.0 * (np_ - op) / op))
    elif op != np_:
        notes.append("points %d→%d" % (op, np_))
    for k, tol in (("minP", TOL_CLEAR_MM), ("minS", TOL_CLEAR_MM)):
        if k not in old or k not in new:
            continue
        d = new[k] - old[k]
        if d < -tol:
            probs.append("%s %.2f→%.2f (%.2f)" % (k, old[k], new[k], d))
        elif d > tol:
            notes.append("%s %.2f→%.2f (+%.2f)" % (k, old[k], new[k], d))
    return bool(probs), probs + notes


def main():
    argv = sys.argv[1:]
    update = "--update" in argv
    no_clear = "--no-clear" in argv
    only = None
    if "--only" in argv:
        only = argv[argv.index("--only") + 1]
    fast = "--fast" in argv          # 跳过 big.ply（单件 7 分钟，人工互检时才跑）
    os.makedirs(OUT, exist_ok=True)

    baseline = {}
    if os.path.exists(BASE):
        baseline = json.load(open(BASE, encoding="utf-8"))

    res, nreg = {}, 0
    print("=" * 100)
    print("%-11s %5s %5s %7s %7s %7s %7s %6s  %s"
          % ("用例", "面", "行", "点数", "minP", "minS", "耗时s", "绕行", "结论"))
    print("-" * 100)
    for name, ply, args, clear, desc in CASES:
        if only and only not in name:
            continue
        if fast and name == "big":
            continue
        r = run_case(name, ply, args, 0 if no_clear else clear)
        res[name] = r
        if update or name not in baseline:
            verdict = "基线"
            if "fail" in r:
                verdict = "崩溃"
        else:
            bad, msgs = cmp(baseline[name], r)
            verdict = ("REGRESS " + "; ".join(msgs[:2])) if bad else (
                "OK" if not msgs else "OK " + "; ".join(msgs[:2]))
            if bad:
                nreg += 1
        if "skip" in r:
            verdict = "跳过 " + r["skip"]
        print("%-11s %5s %5s %7s %7s %7s %7s %6s  %s"
              % (name,
                 r.get("faces", "-"), r.get("rows", "-"), r.get("points", "-"),
                 r.get("minP", "-"), r.get("minS", "-"), r.get("sec", "-"),
                 r.get("n_lift", "-"), verdict))
        if "fail" in r:
            print("  !! %s" % r["fail"])
            print("  " + (r.get("log", "") or "").replace("\n", "\n  ")[-900:])

    if update:
        merged = dict(baseline)
        merged.update(res)
        json.dump(merged, open(BASE, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print("\n基线已更新: %s (%d 用例)" % (BASE, len(merged)))
        return 0

    print("-" * 100)
    if nreg:
        print("结论: %d 个用例回归" % nreg)
        return 1
    print("结论: 全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
