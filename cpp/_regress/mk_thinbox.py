# 生成薄壁开口盒点云（400x300x300，壁厚4），用于验证"内壁顶部缺口"。
# 顶面开口，只有 4 条 4mm 宽的顶沿。采样步长 2mm，加 0.02mm 噪声。
import sys, random

W, D, H, T = 400.0, 300.0, 300.0, 4.0
STEP = float(__import__("os").environ.get("THINBOX_STEP","2.0"))
random.seed(7)


def grid(u0, u1, v0, v1, fn):
    pts = []
    nu = max(1, int(round((u1 - u0) / STEP)))
    nv = max(1, int(round((v1 - v0) / STEP)))
    for i in range(nu + 1):
        u = u0 + (u1 - u0) * i / nu
        for j in range(nv + 1):
            v = v0 + (v1 - v0) * j / nv
            pts.append(fn(u, v))
    return pts


P = []
# 外表面
P += grid(0, W, 0, D, lambda x, y: (x, y, 0.0))          # 外底 z=0
P += grid(0, W, 0, H, lambda x, z: (x, 0.0, z))          # 外壁 y=0
P += grid(0, W, 0, H, lambda x, z: (x, D, z))            # 外壁 y=D
P += grid(0, D, 0, H, lambda y, z: (0.0, y, z))          # 外壁 x=0
P += grid(0, D, 0, H, lambda y, z: (W, y, z))            # 外壁 x=W
# 内表面（底面 z=T，四壁从 z=T 到 z=H，顶面开口）
P += grid(T, W - T, T, D - T, lambda x, y: (x, y, T))    # 内底 z=T
P += grid(T, W - T, T, H, lambda x, z: (x, T, z))        # 内壁 y=T
P += grid(T, W - T, T, H, lambda x, z: (x, D - T, z))    # 内壁 y=D-T
P += grid(T, D - T, T, H, lambda y, z: (T, y, z))        # 内壁 x=T
P += grid(T, D - T, T, H, lambda y, z: (W - T, y, z))    # 内壁 x=W-T
# 顶沿（4mm 宽的环，z=H）
P += grid(0, W, 0, T, lambda x, y: (x, y, H))
P += grid(0, W, D - T, D, lambda x, y: (x, y, H))
P += grid(0, T, T, D - T, lambda x, y: (x, y, H))
P += grid(W - T, W, T, D - T, lambda x, y: (x, y, H))

out = sys.argv[1] if len(sys.argv) > 1 else "thinbox.ply"
# --inner-gap <mm>：把内壁顶部这一截（模拟扫描遮挡/边缘点被判给顶沿）从云里去掉，
# 用于复现"行端点比真实边沿短一截"造成的漏喷。
GAP = 0.0
if "--inner-gap" in sys.argv:
    GAP = float(sys.argv[sys.argv.index("--inner-gap") + 1])
    # 顶部 GAP 这一圈整体缺失（含顶沿环、内外壁顶部）：模拟扫描没扫到最上面一截
    P = [p for p in P if p[2] < H - GAP + 0.1]
with open(out, "w") as f:
    f.write("ply\nformat ascii 1.0\nelement vertex %d\n" % len(P))
    f.write("property float x\nproperty float y\nproperty float z\nend_header\n")
    for p in P:
        f.write("%.3f %.3f %.3f\n" % tuple(c + random.gauss(0, 0.02) for c in p))
print("wrote", out, len(P), "points")
