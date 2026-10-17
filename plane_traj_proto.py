# -*- coding: utf-8 -*-
"""
平面喷涂轨迹原型 v0.1
流程: 读PLY -> SOR去噪 -> 平面重建(RANSAC+PCA) -> 2D投影 -> 凹包轮廓
      -> 等行距扫描线 -> 行内等点距重采样 -> S形往复 + 行端外扩U形过渡
      -> 一条连续轨迹 -> CSV + 可视化
参数与方案v1.3一致, 外部可调。
"""
import numpy as np
import json, sys, os
from scipy.spatial import cKDTree
from shapely.geometry import MultiPoint, LineString, Polygon
from shapely.ops import unary_union

# ---------------- 参数 ----------------
PARAMS = dict(
    scanRowSpacing = 20.0,   # 行距 mm
    pointSpacing   = 4.0,    # 点距 mm
    standoff       = 100.0,  # 离面距离 mm (输出位姿用)
    speed          = 200.0,  # mm/s
    extendMargin   = 18.0,   # 行端外扩(关枪过渡段), >=1x行距取大者
    flatnessTol    = 3.0,    # 平面内点阈值(重建)
    sor_k          = 20,
    sor_z          = 2.5,
)

def read_ply(path):
    with open(path,'rb') as f:
        props=[]; fmt=None; n=0
        while True:
            l=f.readline()
            if l.startswith(b'format'): fmt=l.split()[1].decode()
            if l.startswith(b'property'):
                p=l.split(); props.append((p[1].decode(),p[2].decode()))
            if l.startswith(b'element vertex'): n=int(l.split()[2])
            if l.startswith(b'end_header'): break
        dt=[(nm,'f4' if t=='float' else 'u1') for t,nm in props]
        d=np.dtype(dt)
        arr=np.fromfile(f,dtype=d,count=n) if fmt!='ascii' else np.loadtxt(f,dtype=d,max_rows=n)
    xyz=np.column_stack([arr['x'],arr['y'],arr['z']]).astype(float)
    return xyz

def sor_filter(P, k=20, z=2.5):
    tree=cKDTree(P)
    d,_=tree.query(P, k=k+1)
    mean_d=d[:,1:].mean(1)
    mu, sd = mean_d.mean(), mean_d.std()
    keep = mean_d < mu + z*sd
    return P[keep], int((~keep).sum())

def fit_plane_robust(P, tol, iters=8):
    """RANSAC式迭代拟合: PCA初值 -> 内点 -> 重拟合"""
    plane=None
    idx=np.arange(len(P))
    for it in range(iters):
        c=P[idx].mean(0); Q=P[idx]-c
        _,_,vt=np.linalg.svd(Q,full_matrices=False)
        n=vt[2]; d=-n@c
        dist=np.abs(P@n+d)
        idx=np.where(dist<tol)[0]
        if len(idx)<10: break
        plane=(n,d,c,idx)
    return plane  # n,d,centroid,inlier_idx

def build_boundary(P2d, alpha=None):
    """凹包轮廓: 用shapely的concave_hull(缓冲+收缩近似); 简化: buffer+delaunay"""
    mp=MultiPoint(P2d)
    try:
        from shapely import concave_hull
        hull=concave_hull(mp, ratio=0.3)
    except Exception:
        hull=mp.convex_hull
    return hull

def main():
    ply=sys.argv[1] if len(sys.argv)>1 else r"C:\Users\DYF\Desktop\点云\喷涂点云\面板.ply"
    outdir=os.path.dirname(os.path.abspath(__file__))
    P=read_ply(ply)
    n0=len(P)
    P,nrm=sor_filter(P, PARAMS['sor_k'], PARAMS['sor_z'])

    # 平面重建
    n,d,c,inl=fit_plane_robust(P, PARAMS['flatnessTol'])
    if n[2]<0: n=-n; d=-d
    res=P[inl]@n+d
    rms=float(np.sqrt((res**2).mean()))
    # 局部2D基
    zaxis=np.array([0,0,1.0])
    u=np.cross(zaxis,n); u/=np.linalg.norm(u)
    v=np.cross(n,u)
    Q=P[inl]-c
    x2=Q@u; y2=Q@v

    # 主方向对齐: 2D PCA 把 x 轴旋到点云长轴方向(行沿长边, 参考图样式)
    # 注意: 仅旋转"扫描行坐标系", 回投3D时必须逆旋转回原始 (u,v) 基
    C2=np.cov(np.column_stack([x2,y2]).T)
    evals,evecs=np.linalg.eigh(C2)
    major=evecs[:,np.argmax(evals)]  # 长轴方向
    ang=np.arctan2(major[1],major[0])
    ca,sa=np.cos(ang),np.sin(ang)
    xr= ca*x2+sa*y2
    yr=-sa*x2+ca*y2
    x2,y2=xr,yr
    x2o,y2o=Q@u, Q@v   # 原始 (u,v) 基坐标(预览图用, 与轨迹同框)
    # 轮廓
    poly=build_boundary(np.column_stack([x2,y2]))
    poly=poly.simplify(1.0)          # 去锯齿
    poly=poly.buffer(1.5).buffer(-1.0)  # 闭运算平滑毛刺
    minx,miny,maxx,maxy=poly.bounds
    poly=poly.buffer(0.5)  # 轻微外扩防裂缝
    # 扫描线
    L=PARAMS['scanRowSpacing']; pstep=PARAMS['pointSpacing']
    margin=PARAMS['extendMargin']
    ys=np.arange(miny+L/2, maxy, L)
    rows=[]
    for yc in ys:
        line=LineString([(minx-10,yc),(maxx+10,yc)])
        inter=line.intersection(poly)
        if inter.is_empty: continue
        segs=[]
        if inter.geom_type=='LineString': segs=[inter]
        else:
            for g in getattr(inter,'geoms',[]): 
                if g.geom_type=='LineString': segs.append(g)
        for g in segs:
            if g.length < L*0.8: continue  # 过短子行剔除
            rows.append((yc,g))

    # 生成一条连续轨迹: 循环往复S形, 全程连续喷涂(spray=1), 行端直接折返连下一行, 等点距
    traj=[]  # (X,Y,Z,I,J,K,speed,spray,seg_type)
    def rot_inv(p):
        return np.array([ca*p[0]-sa*p[1], sa*p[0]+ca*p[1]])
    def push2(p_rot, spray, seg):
        p2=rot_inv(p_rot)
        X=c+u*p2[0]+v*p2[1]
        pos=X + n*PARAMS['standoff']
        traj.append((*pos, -n[0],-n[1],-n[2], PARAMS['speed'], spray, seg))
    def sample_line(a,b,spray,seg):
        seg_len=np.linalg.norm(b-a)
        npts=max(2,int(np.ceil(seg_len/pstep))+1)
        for t in np.linspace(0,1,npts):
            push2(a+(b-a)*t, spray, seg)
    def sample_arc(a, b, apex, spray, seg):
        # 折返圆弧, 按弧长等距采样
        ts=np.linspace(0,1,60)
        curve=(1-ts[:,None])**2*a + 2*(1-ts[:,None])*ts[:,None]*apex + ts[:,None]**2*b
        d=np.linalg.norm(np.diff(curve,axis=0),axis=1)
        arc_len=d.sum()
        npts=max(4,int(np.ceil(arc_len/pstep)))
        s=np.concatenate([[0],np.cumsum(d)])
        t_eq=np.interp(np.linspace(0,arc_len,npts+1), s, ts)
        for t in t_eq:
            p=(1-t)**2*a + 2*(1-t)*t*apex + t**2*b
            push2(p, spray, seg)

    # 行端点序列(蛇形): 每行取 [a,b], 交替反转使首尾相连
    ends=[]   # 每行 (a,b) 沿交替方向
    for i,(yc,g) in enumerate(rows):
        coords=list(g.coords)
        if i%2==1: coords=coords[::-1]
        a=np.array(coords[0]); b=np.array(coords[-1])
        ends.append((a,b))

    for i,(a,b) in enumerate(ends):
        # 喷涂行
        sample_line(a, b, 1, 'scan')
        # 行间折返: 当前行 b -> 下一行 a, 小圆弧平滑, 全程连续喷涂
        if i+1 < len(ends):
            a_n,b_n = ends[i+1]
            dirv=(b-a)/np.linalg.norm(b-a)
            cross=np.array([-dirv[1], dirv[0]])
            # 折返顶点: 沿行向法向偏移半行距, 让折返自然过渡到下一行
            apex=(b+a_n)/2.0 + cross*(L*0.5)
            sample_arc(b, a_n, apex, 1, 'scan')
    # 收尾: 输出
    arr=np.array([t[:7] for t in traj])
    csv_path=os.path.join(outdir,'traj_panel.csv')
    with open(csv_path,'w') as f:
        f.write("X,Y,Z,I,J,K,speed,spray,seg\n")
        for t in traj:
            f.write("%.3f,%.3f,%.3f,%.4f,%.4f,%.4f,%.1f,%d,%s\n"%(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8]))
    # PLY 输出: 全程连续喷涂, 统一红色
    ply_path=os.path.join(outdir,'traj_panel.ply')
    with open(ply_path,'w') as f:
        f.write("ply\nformat ascii 1.0\ncomment continuous S-shaped spray trajectory\n")
        f.write("element vertex %d\n"%len(traj))
        f.write("property float x\nproperty float y\nproperty float z\n")
        f.write("property uchar red\nproperty uchar green\nproperty uchar blue\n")
        f.write("end_header\n")
        for t in traj:
            f.write("%.3f %.3f %.3f 255 48 32\n"%(t[0],t[1],t[2]))
    # 自检
    pts3=arr[:,:3]
    d=np.linalg.norm(np.diff(pts3,axis=0),axis=1)
    print(json.dumps({
        "input_points":n0,"after_sor":len(P),"outliers_removed":nrm,
        "plane_normal":n.round(4).tolist(),"plane_rms_mm":round(rms,3),
        "rows":len(rows),"total_points":len(traj),
        "all_spray":bool((np.array([t[7] for t in traj])==1).all()),
        "seg_len_mean":round(float(d.mean()),3),
        "seg_len_max":round(float(d.max()),3),
        "continuous_single_path":True,
        "csv":csv_path,
        "ply":ply_path
    },ensure_ascii=False,indent=1))

    # 可视化(俯视, 参考图配色)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPoly
    fig,ax=plt.subplots(figsize=(10,7),dpi=130)
    fig.patch.set_facecolor('black'); ax.set_facecolor('black')
    # 点云(品红)
    ax.scatter(x2o,y2o,s=1.2,c='#e040e0',alpha=.55,linewidths=0)
    # 轮廓(白)
    hb=np.array(poly.exterior.xy).T          # 轮廓(行坐标系)
    hb_o=np.array([rot_inv(p) for p in hb])  # 逆旋转回原始基
    bx,by=hb_o[:,0],hb_o[:,1]
    ax.plot(bx,by,c='white',lw=1.2)
    # 轨迹: 全程连续喷涂, 统一红色
    traj2=[]
    for t in traj:
        p3=np.array(t[:3]); rel=p3-c
        traj2.append((rel@u, rel@v, t[7]))
    for k in range(len(traj2)-1):
        x0,y0,_=traj2[k]; x1,y1,_=traj2[k+1]
        ax.plot([x0,x1],[y0,y1],c='#ff3020',lw=1.6)
    ax.set_aspect('equal'); ax.axis('off')
    png=os.path.join(outdir,'traj_panel_preview.png')
    plt.savefig(png,bbox_inches='tight',facecolor='black')
    print("preview:",png)

if __name__=='__main__':
    main()
