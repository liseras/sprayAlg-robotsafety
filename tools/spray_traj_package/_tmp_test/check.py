import csv, math, sys

holes = [(109.0, -61.0, 55.0), (-101.0, 49.0, 45.0)]

def check(path):
    rows = list(csv.DictReader(open(path)))
    P = [(float(r['x_mm']), float(r['y_mm']), float(r['z_mm'])) for r in rows]
    Z = [(float(r['toolz_x']), float(r['toolz_y']), float(r['toolz_z'])) for r in rows]
    maxstep = 0; maxang = 0; at = None
    for k in range(1, len(P)):
        d = math.dist(P[k], P[k-1])
        if d > maxstep: maxstep = d; at = k
        c = max(-1, min(1, sum(a*b for a, b in zip(Z[k], Z[k-1]))))
        maxang = max(maxang, math.degrees(math.acos(c)))
    # travel 点离孔心最近距离（判穿孔）
    minhole = 1e9; worst = None
    for k, r in enumerate(rows):
        if r['zone'] != 'travel': continue
        su, sv = float(r['su_mm']), float(r['sv_mm'])
        for (cu, cv, rr) in holes:
            d = math.hypot(su-cu, sv-cv)
            if d < minhole: minhole = d; worst = (k, su, sv)
    return dict(pts=len(rows), maxstep=round(maxstep,2), step_at=at,
                maxang=round(maxang,2), minhole=round(minhole,1), worst=worst)

for p in sys.argv[1:]:
    print(p, check(p))
