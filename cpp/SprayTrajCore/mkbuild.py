# -*- coding: utf-8 -*-
"""Incremental build driver for SprayTrajCore.

Why this exists: build_release.bat calls cl.exe on spray_traj.cpp (8.6k lines)
and robot_safety.cpp (5.1k lines) with /O2 every single time. cl has no
incremental mode, so one changed line costs the full ~80s.

This script keeps the .obj files and only recompiles the translation units
whose .cpp / any .h is newer, then relinks only what depends on them.
Env (INCLUDE/LIB/PATH) must already be set by build_inc.bat.
"""
import os
import subprocess
import sys
import time

SRC = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(SRC, "build_rls")
OBJ = os.path.join(OUT, "obj")
CFLAGS = ["/nologo", "/std:c++17", "/EHsc", "/O2", "/MD", "/utf-8",
          "/DNDEBUG", "/DSPRAYTRAJ_EXPORTS", "/I" + SRC]

DLL = os.path.join(OUT, "SprayTraj.dll")
LIB = os.path.join(OUT, "SprayTraj.lib")


def cl_path():
    for d in os.environ.get("PATH", "").split(os.pathsep):
        p = os.path.join(d, "cl.exe")
        if os.path.isfile(p):
            return p
    return None


CL = cl_path()
if not CL:
    print("ERROR: cl.exe not on PATH - run through build_inc.bat")
    sys.exit(2)

os.makedirs(OBJ, exist_ok=True)
HEADERS = [os.path.join(SRC, f) for f in os.listdir(SRC) if f.endswith(".h")]


def stamp(p):
    try:
        return os.path.getmtime(p)
    except OSError:
        return -1.0


def is_stale(obj, cpp):
    """True when obj must be rebuilt."""
    t = stamp(obj)
    if t < 0:
        return True
    for dep in [cpp] + HEADERS:
        if stamp(dep) > t:
            return True
    return False


def run(args, what):
    t0 = time.time()
    # 不用 text=True：cl.exe 在中文 Windows 上输出 GBK，Python 默认按 UTF-8 解码会
    # 抛 UnicodeDecodeError，把真正的编译错误吞掉（表现为 [FAIL] 但没有任何 error 行）。
    r = subprocess.run(args, cwd=SRC, capture_output=True)
    out = (r.stdout or b"").decode("mbcs", "replace")
    err = (r.stderr or b"").decode("mbcs", "replace")
    dt = time.time() - t0
    if r.returncode != 0:
        print("[FAIL] %s (%.0fs)" % (what, dt))
        sys.stdout.write(out)
        sys.stderr.write(err)
        sys.exit(1)
    for line in out.splitlines():
        if "warning" in line.lower():
            print("  " + line)
    print("[ok] %s (%.1fs)" % (what, dt))
    return dt


def compile_one(name, always=False):
    cpp = os.path.join(SRC, name + ".cpp")
    obj = os.path.join(OBJ, name + ".obj")
    if not always and not is_stale(obj, cpp):
        print("[skip] %s.obj up to date" % name)
        return obj, False
    run([CL] + CFLAGS + ["/c", "/Fo" + obj, cpp], "compile " + name + ".cpp")
    return obj, True


t_all = time.time()
built = []

o_traj, c1 = compile_one("spray_traj")
o_safe, c2 = compile_one("robot_safety")
# stamp.cpp holds __DATE__/__TIME__ and must NEVER be cached, otherwise the
# DLL reports the time of some older build (see stamp.cpp).
o_stamp, _ = compile_one("stamp", always=True)
built += [x for x, c in (((o_traj, c1)), ((o_safe, c2))) if c]

need_dll = bool(built) or not os.path.isfile(DLL) or stamp(DLL) < max(stamp(o_traj), stamp(o_safe))
if need_dll:
    run([CL, "/nologo", "/LD", "/Fe:" + DLL, o_traj, o_safe, o_stamp,
         "/link", "/IMPLIB:" + LIB], "link SprayTraj.dll")

dll_t = stamp(DLL)

for exe, src in (("spray_safety", "safety_driver"), ("spray_driver", "driver")):
    o, changed = compile_one(src)
    exe_p = os.path.join(OUT, exe + ".exe")
    if changed or not os.path.isfile(exe_p) or stamp(exe_p) < dll_t:
        run([CL, "/nologo", "/Fe:" + exe_p, o,
             "/link", "/SUBSYSTEM:CONSOLE", "/LIBPATH:" + OUT,
             "SprayTraj.lib", "/STACK:16777216"], "link " + exe + ".exe")

print("BUILD OK: %s (%.1fs)" % (OUT, time.time() - t_all))
