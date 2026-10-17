// spray_traj.cpp —— 喷涂轨迹生成核心实现
// 流程：SOR去噪 -> RANSAC平面重建 -> 2D投影 + 长轴对齐 -> 凹包轮廓
//       -> 等行距扫描线 -> S形连续轨迹（行端圆弧折返，全程喷涂、等距重采样）
#include "spray_traj.h"
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <string>
#include <random>

// ---------------- 基础几何 ----------------
struct Vec3 { double x=0,y=0,z=0; };
static inline Vec3 operator+(Vec3 a, Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(Vec3 a, Vec3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(Vec3 a, double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dot(Vec3 a, Vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross(Vec3 a, Vec3 b){return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};}
static inline double len(Vec3 a){return std::sqrt(dot(a,a));}
static inline Vec3 norm(Vec3 a){double l=len(a); return l>1e-12? Vec3{a.x/l,a.y/l,a.z/l}: Vec3{0,0,0};}

struct Vec2 { double x=0,y=0; };
static inline Vec2 operator+(Vec2 a, Vec2 b){return {a.x+b.x,a.y+b.y};}
static inline Vec2 operator-(Vec2 a, Vec2 b){return {a.x-b.x,a.y-b.y};}
static inline Vec2 operator*(Vec2 a, double s){return {a.x*s,a.y*s};}
static inline Vec2 operator/(Vec2 a, double s){return {a.x/s,a.y/s};}
static inline double dot2(Vec2 a, Vec2 b){return a.x*b.x+a.y*b.y;}
static inline double cross2(Vec2 a, Vec2 b){return a.x*b.y-a.y*b.x;}
static inline double len2(Vec2 a){return std::sqrt(dot2(a,a));}
static inline Vec2 norm2(Vec2 a){double l=len2(a); return l>1e-12? Vec2{a.x/l,a.y/l}: Vec2{0,0};}

// ---------------- 简单 KD 树 ----------------
struct KdTree {
    struct Node { int idx; int left=-1, right=-1; };
    std::vector<Node> nodes;
    const std::vector<Vec3>* pts=nullptr;
    void build(std::vector<int> ids, int axis){
        if(ids.empty()) return;
        std::sort(ids.begin(), ids.end(), [&](int a, int b){
            double va = (axis==0)?(*pts)[a].x : (axis==1)?(*pts)[a].y:(*pts)[a].z;
            double vb = (axis==0)?(*pts)[b].x : (axis==1)?(*pts)[b].y:(*pts)[b].z;
            return va < vb;
        });
        int mid = ids.size()/2;
        int ni = (int)nodes.size();
        nodes.push_back({ids[mid], -1, -1});
        std::vector<int> l(ids.begin(), ids.begin()+mid);
        std::vector<int> r(ids.begin()+mid+1, ids.end());
        if(!l.empty()){ nodes[ni].left=(int)nodes.size(); build(l,(axis+1)%3);}
        if(!r.empty()){ nodes[ni].right=(int)nodes.size(); build(r,(axis+1)%3);}
    }
    void init(const std::vector<Vec3>& p){
        pts=&p; nodes.clear(); nodes.reserve(p.size());
        std::vector<int> ids(p.size()); std::iota(ids.begin(),ids.end(),0);
        build(ids,0);
    }
    // K 近邻距离和
    double knnSum(Vec3 q, int k) const {
        std::vector<double> best(k, 1e18);
        knnRec(0, q, k, 0, best);
        return std::accumulate(best.begin(), best.end(), 0.0);
    }
    void knnRec(int node, Vec3 q, int k, int axis, std::vector<double>& best) const {
        if(node<0) return;
        const Vec3& p = (*pts)[nodes[node].idx];
        double d = len(q-p);
        double worst = best[k-1];
        if(d < worst){
            best[k-1]=d;
            std::sort(best.begin(), best.end());
        }
        double diff = (axis==0)? q.x-p.x : (axis==1)? q.y-p.y : q.z-p.z;
        int near = diff<0 ? nodes[node].left : nodes[node].right;
        int far  = diff<0 ? nodes[node].right : nodes[node].left;
        knnRec(near, q, k, (axis+1)%3, best);
        if(std::fabs(diff) < best[k-1]) knnRec(far, q, k, (axis+1)%3, best);
    }
};

// ---------------- 主实现 ----------------
namespace {

struct TrajPoint { Vec3 pos; Vec3 dir; double speed; int spray; };

// SOR 去噪
std::vector<Vec3> sor_filter(const std::vector<Vec3>& P, int k, double z){
    KdTree tree; tree.init(P);
    std::vector<double> md(P.size());
    double sum=0;
    for(size_t i=0;i<P.size();++i){ md[i]=tree.knnSum(P[i], k)/k; sum+=md[i]; }
    double mu=sum/P.size();
    double var=0; for(double v:md) var+=(v-mu)*(v-mu); var/=P.size();
    double sd=std::sqrt(var);
    std::vector<Vec3> out;
    for(size_t i=0;i<P.size();++i) if(md[i] < mu + z*sd) out.push_back(P[i]);
    return out;
}

// 平面拟合：返回法向n、d、质心c、内点索引
struct Plane { Vec3 n; double d; Vec3 c; };

// 用一组点做最小二乘平面拟合（正规方程 + 逆迭代求最小特征向量）
static bool ls_plane(const std::vector<Vec3>& P, const std::vector<int>& idx,
                     Vec3& n, Vec3& c)
{
    if(idx.size()<3) return false;
    Vec3 m{0,0,0}; for(int i:idx){m.x+=P[i].x;m.y+=P[i].y;m.z+=P[i].z;}
    m=m*(1.0/(double)idx.size());
    double A00=0,A01=0,A02=0,A11=0,A12=0,A22=0;
    for(int i:idx){ double dx=P[i].x-m.x,dy=P[i].y-m.y,dz=P[i].z-m.z;
        A00+=dx*dx; A01+=dx*dy; A02+=dx*dz; A11+=dy*dy; A12+=dy*dz; A22+=dz*dz; }
    double invn=1.0/(double)idx.size();
    A00*=invn;A01*=invn;A02*=invn;A11*=invn;A12*=invn;A22*=invn;
    // 最小特征向量 = 法向（逆迭代）
    Vec3 x{0,0,1};
    for(int inner=0; inner<80; ++inner){
        double det = A00*(A11*A22-A12*A12) - A01*(A01*A22-A12*A02) + A02*(A01*A12-A11*A02);
        if(std::fabs(det)<1e-20) break;
        double i00=(A11*A22-A12*A12)/det;
        double i01=(A02*A12-A01*A22)/det;
        double i02=(A01*A12-A11*A02)/det;
        double i11=(A00*A22-A02*A02)/det;
        double i12=(A01*A02-A00*A12)/det;
        double i22=(A00*A11-A01*A01)/det;
        Vec3 y{ i00*x.x+i01*x.y+i02*x.z,
                i01*x.x+i11*x.y+i12*x.z,
                i02*x.x+i12*x.y+i22*x.z };
        double l=std::sqrt(dot(y,y));
        if(l<1e-12) break;
        x=Vec3{y.x/l,y.y/l,y.z/l};
    }
    n=x; c=m;
    return true;
}

// RANSAC 平面拟合：随机抽样 + 内点投票 + 内点最小二乘精化
// 只有到平面的距离 < tol 的点才算内点（凸起/凹陷超过 tol 会被"抹平"剔除）
bool fit_plane(const std::vector<Vec3>& P, double tol, Plane& pl, std::vector<int>& inl){
    int N=(int)P.size();
    if(N<10) return false;

    // 期望内点率下限 0.3，据此定迭代次数（95% 置信度至少抽中一次全内点样本）
    const double w = 0.3;
    int iters = (int)std::ceil(std::log(0.05) / std::log(1.0 - std::pow(w,3.0)));
    iters = std::max(200, std::min(800, iters));

    std::vector<int> best;
    Vec3 bn{0,0,1}, bc{0,0,0};

    std::mt19937 rng(12345);
    std::uniform_int_distribution<int> ud(0, N-1);

    for(int it=0; it<iters; ++it){
        // 抽 3 个不重复点
        int a=ud(rng), b=ud(rng), c2=ud(rng);
        if(a==b||b==c2||a==c2) continue;
        Vec3 n = norm(cross(P[b]-P[a], P[c2]-P[a]));
        if(len(n)<1e-6) continue;             // 三点共线，跳过
        double d = -dot(n, P[a]);
        // 统计内点
        std::vector<int> cur;
        for(int i=0;i<N;++i) if(std::fabs(dot(P[i],n)+d) < tol) cur.push_back(i);
        if(cur.size() > best.size()){
            best.swap(cur); bn=n;
        }
    }

    if(best.size() < std::max(3, (int)(N*0.1))) return false;

    // 内点最小二乘精化
    Vec3 nf, cf;
    if(!ls_plane(P, best, nf, cf)) return false;
    double df = -dot(nf, cf);

    // 精化后再按 tol 收一次内点，得到最终一致集
    std::vector<int> fin;
    for(int i=0;i<N;++i) if(std::fabs(dot(P[i],nf)+df) < tol) fin.push_back(i);
    if(fin.size() < std::max(3, (int)(N*0.1))) fin = best;

    pl.n=nf; pl.d=df; pl.c=cf;
    inl=fin;
    return inl.size() >= (size_t)(N*0.1);
}

// 顺序 RANSAC 多平面分割：反复拟合"最大平面"并剥离其内点，
// 直到剩余点不足。返回多个平面（每个含法向/质心/内点索引），
// 用于折型板（两面）、槽型盘（底盘+凸起侧边）这类多面工件分别重建。
//
// 参数语义（重要，修正过）：
//   tol (flatnessTol)  = 台阶抹平阈值：两共面片之间法向间距 < tol → 合并抹平；
//                        ≥ tol → 保持独立分割。它 *不是* 平面拟合的内点容差！
//                        若把大 tol 当拟合容差，凸边会被吸进底盘平面，
//                        最小二乘把平面拉斜（法向明显偏斜），重建面斜穿工件。
//   curvatureTol       = 曲面升级阈值：整片点云到主平面的残差 P99.5 <= 此值
//                        → 单一可喷面（曲率浅当平面）。
// 实现要点：
//   1) 平面拟合一律用紧容差 tightTol = min(tol, 3mm)，保证每个面片拟合准确；
//   2) 第一级：紧 tol 拟合主平面，全点到主平面残差 P99.5 <= curvatureTol
//      → 单面直接返回（不切分、不重叠）；
//   3) 第二级：顺序 RANSAC（紧 tol）切分出各面片；
//   4) 合并：法向夹角 < 15° 且两平面沿法向间距 < tol → 合并，
//      保留点数多的主段平面方程、吸收小段内点（"抹平到主平面"）。
//      RMS/P99.5 判据被弃用：RMS 对小面积高台阶不敏感（被底盘稀释），
//      曾导致凸边被错误合并、轨迹丢失。
struct PlaneSeg { Plane pl; std::vector<int> inl; };
std::vector<PlaneSeg> split_planes(const std::vector<Vec3>& P, double tol, double curvatureTol, int minPts){
    std::vector<PlaneSeg> out;
    const double tightTol = std::min(tol, 3.0);   // 拟合容差恒紧，与台阶阈值解耦

    // ---- 第一级：整体是否可当作"单一可喷面"？ ----
    {
        std::vector<Vec3> sub(P);   // 全点
        Plane pl; std::vector<int> inl;
        if(fit_plane(sub, tightTol, pl, inl)){
            if(dot(pl.n, Vec3{0,0,1}) < 0){ pl.n = pl.n*(-1.0); pl.d = -pl.d; }
            // 全点到主平面的残差 P99.5：捕捉少数但高的台阶/凸边点
            std::vector<double> res; res.reserve(P.size());
            for(size_t k=0;k<P.size();++k) res.push_back(std::fabs(dot(P[k],pl.n)+pl.d));
            double p995 = 0.0;
            {
                size_t kk = (size_t)((res.size()-1)*0.995);
                std::nth_element(res.begin(), res.begin()+kk, res.end());
                p995 = res[kk];
            }
            // 整体在曲面容差内（无真实台阶）→ 单面，直接返回（不切分、不重叠）
            if(p995 <= curvatureTol){
                PlaneSeg seg; seg.pl = pl; seg.inl = inl;
                fprintf(stderr, "[split] 整体单面: 内点=%zu/%zu P99.5=%.3f (<=curvatureTol=%.3f)\n",
                    inl.size(), P.size(), p995, curvatureTol);
                out.push_back(std::move(seg));
                return out;
            }
            fprintf(stderr, "[split] 整体 P99.5=%.3f > curvatureTol=%.3f，存在台阶/凸边，进入多面分割\n",
                p995, curvatureTol);
        }
    }

    // ---- 第二级：真实多面/台阶 → 顺序 RANSAC（紧容差）切分 ----
    std::vector<int> remain(P.size()); std::iota(remain.begin(), remain.end(), 0);
    // 最多 8 个平面，防止死循环
    for(int iter=0; iter<8 && (int)remain.size() >= minPts; ++iter){
        // 在剩余点集上跑 RANSAC 找最大平面
        Plane pl; std::vector<int> inl;
        std::vector<Vec3> sub; sub.reserve(remain.size());
        for(int i:remain) sub.push_back(P[i]);
        if(!fit_plane(sub, tightTol, pl, inl)) break;
        if((int)inl.size() < minPts) break;

        // inl 是 sub 的下标，映射回 P 的下标
        std::vector<int> globalInl; globalInl.reserve(inl.size());
        for(int k:inl) globalInl.push_back(remain[k]);

        PlaneSeg seg; seg.pl=pl; seg.inl=globalInl;
        fprintf(stderr, "[split] 段%d: 点数=%zu 法向=(%.3f,%.3f,%.3f) 质心=(%.1f,%.1f,%.1f)\n",
            iter, globalInl.size(), pl.n.x, pl.n.y, pl.n.z, pl.c.x, pl.c.y, pl.c.z);
        out.push_back(std::move(seg));

        // 从剩余集中剔除该平面的内点
        std::vector<bool> used(P.size(), false);
        for(int k:globalInl) used[k]=true;
        std::vector<int> nr; nr.reserve(remain.size()-globalInl.size());
        for(int k:remain) if(!used[k]) nr.push_back(k);
        remain.swap(nr);
    }

    // ---- 合并"台阶抹平"：法向一致 且 沿法向间距 < tol → 并入主段 ----
    // 保留点数多的段的平面方程（主平面），吸收小段内点。
    // 间距用两质心沿法向的投影差（平行平面间距），与 RMS 无关——
    // 语义就是"高差在抹平阈值内的台阶被抹平到主平面"。
    const double cosThresh = std::cos(15.0*3.141592653589793/180.0);
    bool merged = true;
    while(merged){
        merged = false;
        for(size_t i=0;i<out.size() && !merged;++i){
            for(size_t j=i+1;j<out.size() && !merged;++j){
                double cd = std::fabs(dot(out[i].pl.n, out[j].pl.n));
                if(cd < cosThresh) continue;       // 夹角大，不同朝向的面
                // 统一法向符号后算平行间距
                Vec3 nj = out[j].pl.n;
                if(dot(nj, out[i].pl.n) < 0) nj = nj*(-1.0);
                double dist = std::fabs(dot(out[j].pl.c - out[i].pl.c, out[i].pl.n));
                // 合并间距阈值：取台阶抹平阈值与浅曲面容差的较大者。
                // flatnessTol 覆盖"台阶高差在阈值内→抹平"；
                // curvatureTol 覆盖"缓曲底面被紧tol切成碎块、块间距略大→并回"。
                const double mergeDist = std::max(tol, curvatureTol*1.2);   // 20%余量：缓弯底面碎块间距常略超曲率容差
                if(dist >= mergeDist) continue;    // 台阶高差 >= 抹平阈值 → 保持独立

                // 小段并入大段（保留主段平面方程）
                size_t big = i, small = j;
                if(out[j].inl.size() > out[i].inl.size()){ big = j; small = i; }
                Vec3 nB = out[big].pl.n;
                if(dot(out[small].pl.n, nB) < 0){
                    out[small].pl.n = out[small].pl.n*(-1.0);
                    out[small].pl.d = -out[small].pl.d;
                }
                out[big].inl.insert(out[big].inl.end(),
                                    out[small].inl.begin(), out[small].inl.end());
                out.erase(out.begin() + small);   // erase 后重开 while 循环，下标重新取
                merged = true;
            }
        }
    }
    return out;
}

// 凹包（单调链凸包 + 边界膨胀；简单用凸包近似，平面板件足够）
std::vector<Vec2> convex_hull(std::vector<Vec2> pts){
    std::sort(pts.begin(),pts.end(),[](Vec2 a,Vec2 b){return a.x<b.x||(a.x==b.x&&a.y<b.y);});
    std::vector<Vec2> h;
    for(int pass=0; pass<2; ++pass){
        size_t base=h.size();
        for(auto&p:pts){
            while(h.size()>=base+2 && cross2(h[h.size()-1]-h[h.size()-2], p-h[h.size()-1])<=0) h.pop_back();
            h.push_back(p);
        }
        h.pop_back();
        std::reverse(pts.begin(),pts.end());
    }
    return h;
}

// 点在凸多边形内（含边界）
bool inside_poly(const std::vector<Vec2>& H, Vec2 p){
    int n=H.size();
    for(int i=0;i<n;++i){ Vec2 a=H[i],b=H[(i+1)%n];
        if(cross2(b-a, p-a) < -1e-9) return false; }
    return true;
}

// 线段与多边形求交的区间（水平扫描线 y=yc）
bool slice_interval(const std::vector<Vec2>& H, double yc, double& xs, double& xe){
    std::vector<double> hits;
    int n=H.size();
    for(int i=0;i<n;++i){ Vec2 a=H[i],b=H[(i+1)%n];
        if((a.y<=yc && b.y>yc)||(b.y<=yc&&a.y>yc)){
            double t=(yc-a.y)/(b.y-a.y);
            hits.push_back(a.x + t*(b.x-a.x));
        }
    }
    std::sort(hits.begin(),hits.end());
    if(hits.size()<2) return false;
    xs=hits[0]; xe=hits[hits.size()-1];
    return true;
}

} // namespace

// ---------------- 导出接口 ----------------
extern "C" {

SPRAYTRAJ_API SprayTrajParams spray_traj_default_params(void){
    SprayTrajParams p;
    p.scanRowSpacing = 12.0;
    p.pointSpacing   = 4.0;
    p.standoff       = 100.0;
    p.speed          = 200.0;
    p.turnRadiusFactor = 0.5;
    p.flatnessTol    = 3.0;
    p.curvatureTol   = 5.0;
    p.sorK           = 20;
    p.sorZ           = 2.5;
    p.normalFlip     = 0;
    return p;
}

SPRAYTRAJ_API SprayTrajResult spray_traj_generate(
    const double* xyz, int n, const SprayTrajParams* params)
{
    SprayTrajResult res;
    std::memset(&res, 0, sizeof(res));
    res.error = nullptr;

    std::string err;
    if(!xyz || n<10){ err="点云点数不足(<10)"; res.error=_strdup(err.c_str()); return res; }

    std::vector<Vec3> P(n);
    for(int i=0;i<n;++i) P[i]={xyz[i*3],xyz[i*3+1],xyz[i*3+2]};
    res.input_points=n;

    // 1. 去噪
    std::vector<Vec3> Pf = sor_filter(P, params->sorK, params->sorZ);
    res.removed_points = n - (int)Pf.size();
    if(Pf.size()<10){ err="去噪后点数不足"; res.error=_strdup(err.c_str()); return res; }

    // 2. 多平面分割重建：
    //   顺序 RANSAC 反复剥离最大平面。折型板 → 两个面；槽型盘 → 底盘 + 凸起侧边；
    //   高度差小的浅曲面/缓起伏会被第一个主平面吸收（flatnessTol 判定）。
    //   每个平面后续分别拟合、生成轨迹、沿各自法向抬升（避免撞击凸起）。
    std::vector<PlaneSeg> segs = split_planes(Pf, params->flatnessTol, params->curvatureTol, std::max(20, (int)(Pf.size()*0.02)));

    // 若分割结果为空（极少见），退回单平面
    if(segs.empty()){
        Plane pl; std::vector<int> inl;
        if(!fit_plane(Pf, params->flatnessTol, pl, inl)){ err="无法重建有效平面"; res.error=_strdup(err.c_str()); return res; }
        PlaneSeg seg; seg.pl=pl; seg.inl=inl; segs.push_back(std::move(seg));
    }

    // 法向定向：统一保证各面法向朝 +Z 半球（俯视相机外表面朝上），
    // 每个平面分别判向（折型板的两面法向可能相反，需各自翻转）。
    for(auto& seg : segs){
        Vec3& n = seg.pl.n;
        if(params->normalFlip){ n = n*(-1.0); seg.pl.d = -seg.pl.d; }
        else if(std::fabs(n.z) > 0.05){
            if(n.z < 0){ n = n*(-1.0); seg.pl.d = -seg.pl.d; }
        }
        // |n.z| 极小（近竖直面，如槽侧壁）：保持拟合方向，交由 normalFlip 纠正
    }

    // 报告：用面积/点数最大的主平面作为"重建面"统计信息
    {
        size_t bi=0; for(size_t i=1;i<segs.size();++i) if(segs[i].inl.size()>segs[bi].inl.size()) bi=i;
        Plane& pl = segs[bi].pl;
        const auto& inl = segs[bi].inl;
        double rss=0; for(int i:inl){ double e=dot(Pf[i],pl.n)+pl.d; rss+=e*e; }
        res.plane_rms=std::sqrt(rss/std::max(1,(int)inl.size()));
        res.normal_x=pl.n.x; res.normal_y=pl.n.y; res.normal_z=pl.n.z;
        res.plane_cx=pl.c.x; res.plane_cy=pl.c.y; res.plane_cz=pl.c.z;
        res.inlier_points=(int)inl.size();
    }

    // 输出所有分割平面（供 WPF 显示重建面）
    // 为每个面计算局部基(u,v)主轴 + 内点包围盒，WPF 直接画矩形（不做 PCA 近似）
    res.plane_count = (int)segs.size();
    for(size_t i=0;i<segs.size() && i<8;++i){
        Plane& P = segs[i].pl;
        res.planes[i].nx = P.n.x;
        res.planes[i].ny = P.n.y;
        res.planes[i].nz = P.n.z;
        res.planes[i].cx = P.c.x;
        res.planes[i].cy = P.c.y;
        res.planes[i].cz = P.c.z;
        res.planes[i].count = (int)segs[i].inl.size();

        // 局部基 u = cross(z, n)
        Vec3 zaxis{0,0,1};
        Vec3 uu = norm(cross(zaxis, P.n));
        if(len(uu)<1e-6) uu = Vec3{1,0,0};
        Vec3 vv = cross(P.n, uu);

        // 用该面内点做 2D PCA，得到主轴方向，避免矩形斜置
        double suu=0,svv=0,suv=0; int cn=0;
        for(int k : segs[i].inl){
            Vec3 d = Pf[k]-P.c;
            double cu=dot(d,uu), cv=dot(d,vv);
            suu+=cu*cu; svv+=cv*cv; suv+=cu*cv; cn++;
        }
        double ang = 0.0;
        if(cn>0){ suu/=cn; svv/=cn; suv/=cn; ang = 0.5*std::atan2(2*suv, suu-svv); }
        double ca=std::cos(ang), sa=std::sin(ang);
        Vec3 uf = uu*ca + vv*sa;
        Vec3 vf = vv*ca - uu*sa;
        uf = norm(uf); vf = norm(vf);

        // 内点在 (uf,vf) 下的包围盒：用 P0.5/P99.5 分位数（防长尾离群点撑大面片，
        // 与轨迹行端口径一致）
        std::vector<double> pu; pu.reserve(segs[i].inl.size());
        std::vector<double> pv; pv.reserve(segs[i].inl.size());
        for(int k : segs[i].inl){
            Vec3 d = Pf[k]-P.c;
            pu.push_back(dot(d,uf)); pv.push_back(dot(d,vf));
        }
        auto quant2=[&](std::vector<double>& v, double f)->double{
            if(v.empty()) return 0.0;
            size_t k=(size_t)((v.size()-1)*f);
            std::nth_element(v.begin(), v.begin()+k, v.end());
            return v[k];
        };
        double minU=quant2(pu,0.005), maxU=quant2(pu,0.995);
        double minV=quant2(pv,0.005), maxV=quant2(pv,0.995);

        res.planes[i].ux=uf.x; res.planes[i].uy=uf.y; res.planes[i].uz=uf.z;
        res.planes[i].vx=vf.x; res.planes[i].vy=vf.y; res.planes[i].vz=vf.z;
        res.planes[i].minU=minU; res.planes[i].maxU=maxU;
        res.planes[i].minV=minV; res.planes[i].maxV=maxV;
    }

    // 3. 对每个分割出的平面段，分别生成 S 形轨迹并拼接
    struct Row { double xs, xe, yc; int segIdx; };  // segIdx: 该行是扫描线与面交集的第几段（环形面会产生多段）
    double L=params->scanRowSpacing, pstep=params->pointSpacing;
    double turnR = params->turnRadiusFactor * L;

    // 主方向选择 + 行端点计算：对给定 2D 点集，穷举角度(每5°)选最平滑方向
    // 行端/行范围用 P0.5/P99.5 分位数包围盒（全长覆盖，防长尾离群点，也防稀疏端部被行过滤裁掉）
    auto quant=[&](std::vector<double> v, double f)->double{
        if(v.empty()) return 0.0;
        size_t k=(size_t)((v.size()-1)*f);
        std::nth_element(v.begin(), v.begin()+k, v.end());
        return v[k];
    };
    auto genRows=[&](const std::vector<Vec2>& pts, double ang, double L)->std::vector<Row>{
        double c=std::cos(ang), s=std::sin(ang);
        std::vector<Vec2> R(pts.size());
        std::vector<double> rx(pts.size()), ry(pts.size());
        for(size_t i=0;i<pts.size();++i){
            R[i]={c*pts[i].x+s*pts[i].y, -s*pts[i].x+c*pts[i].y};
            rx[i]=R[i].x; ry[i]=R[i].y;
        }
        double mnx=quant(rx,0.005), mxx=quant(rx,0.995);
        double mny=quant(ry,0.005), mxy=quant(ry,0.995);
        std::vector<Row> rs;
        // 支撑点阈值自适应：小平面（如窄侧边）点数少，用相对阈值；
        // 大平面用绝对阈值防噪声孤点。
        int supMin = (int)std::max(3, (int)(pts.size()*0.001));
        supMin = std::min(supMin, 20);
        // 首行从 mny+L*0.25 起、末行到 mxy+L*0.25，让轨迹两端外扩覆盖边缘，
        // 避免"首尾各损失半行距"导致窄侧面轨迹短一截。
        // 行端按"该行带内点的 x 聚类分段"：实心面 → 一整行（全长走刀）；
        // 环形/框形面（如槽型盘凸边环）→ 每行分成左右两段，中央空区不喷，
        // 否则环形面会生成满铺矩形轨迹、与底盘轨迹上下重叠。
        const double gapSeg = std::max(L*1.5, 20.0);   // x 间隙超过此值视为两段
        std::vector<Row> raw;
        for(double yc=mny+L*0.25; yc <= mxy+L*0.25+1e-9; yc+=L){
            // 收集该行带内的 x 值
            std::vector<double> bx; bx.reserve(R.size()/4);
            int sup=0;
            for(auto&p:R){
                if(std::fabs(p.y-yc)<=L*0.5){ sup++; bx.push_back(p.x); }
            }
            if(sup<supMin) continue;          // 跳过完全空白的带
            std::sort(bx.begin(), bx.end());
            // 按间隙聚类成连续段
            double segStart=bx[0], segPrev=bx[0];
            int segCnt=1, segIdx=0;
            auto flushSeg=[&](){
                // 段长 < 行距 的碎段丢弃：喷涂价值小、且会形成梳齿状碎轨迹
                if(segCnt>=supMin && (segPrev-segStart) >= L)
                    raw.push_back({segStart, segPrev, yc, segIdx});
                segIdx++;
            };
            for(size_t i=1;i<bx.size();++i){
                if(bx[i]-segPrev > gapSeg){ flushSeg(); segStart=bx[i]; segCnt=0; }
                segPrev=bx[i]; segCnt++;
            }
            flushSeg();
        }
        // 按段分组排序：先喷完所有第0段行（子蛇形），再第1段……
        // 避免左右段按行交替导致蛇形在环两端来回斜穿。
        std::stable_sort(raw.begin(), raw.end(), [](const Row&a, const Row&b){
            if(a.segIdx != b.segIdx) return a.segIdx < b.segIdx;
            return a.yc < b.yc;
        });
        rs = std::move(raw);
        if(rs.empty()) return rs;
        return rs;
    };

    std::vector<SprayTrajPoint> allPts;
    int totalRows = 0;

    for(auto& seg : segs){
        Plane& pl = seg.pl;
        const auto& inl = seg.inl;
        if(inl.size() < 20) continue;

        // 该平面的局部基
        Vec3 zaxis{0,0,1};
        Vec3 u = norm(cross(zaxis, pl.n));
        if(len(u)<1e-6) u = Vec3{1,0,0};
        Vec3 v = cross(pl.n, u);

        // 投影到 2D
        std::vector<Vec2> P2; P2.reserve(inl.size());
        for(int i:inl){ Vec3 q=Pf[i]-pl.c; P2.push_back({dot(q,u), dot(q,v)}); }

        // 主方向穷举
        double bestCost=1e18, theta=0; std::vector<Row> bestRows;
        for(int deg=0; deg<180; deg+=5){
            double ang=deg*3.141592653589793/180.0;
            auto rs=genRows(P2,ang,L);
            if(rs.size()<2) continue;
            double jag=0;
            for(size_t i=1;i<rs.size();++i)
                jag += std::fabs(rs[i].xs-rs[i-1].xs) + std::fabs(rs[i].xe-rs[i-1].xe);
            jag /= (rs.size()-1);
            double cost = jag + rs.size()*L*0.02;
            if(cost<bestCost){ bestCost=cost; theta=ang; bestRows=rs; }
        }
        if(bestRows.empty()){
            fprintf(stderr, "[traj] 段(点数=%zu, 法向=(%.2f,%.2f,%.2f)) 未生成有效行，跳过\n",
                inl.size(), pl.n.x, pl.n.y, pl.n.z);
            continue;
        }
        totalRows += (int)bestRows.size();
        fprintf(stderr, "[traj] 段(点数=%zu) 生成 %zu 行, theta=%.1f\n",
            inl.size(), bestRows.size(), theta*180/3.141592653589793);

        double ca=std::cos(theta), sa=std::sin(theta);
        auto rotInv=[&](Vec2 p)->Vec2{ return {ca*p.x - sa*p.y, sa*p.x + ca*p.y}; };

        // 生成该平面的 2D 轨迹（每个点带 spray 标记：1=喷涂，0=空移）
        std::vector<std::pair<Vec2,int>> traj2;
        auto sampleLine=[&](Vec2 a, Vec2 b, int spray){
            double LL=len2(b-a);
            int np=std::max(2,(int)std::ceil(LL/pstep)+1);
            for(int t=0;t<np;++t){ double f=(double)t/(np-1); traj2.push_back({a+(b-a)*f, spray}); }
        };
        auto sampleArc=[&](Vec2 a, Vec2 b, Vec2 apex){
            const int M=60; std::vector<Vec2> crv(M+1);
            for(int i=0;i<=M;++i){ double t=(double)i/M; crv[i]=a*((1-t)*(1-t)) + apex*(2*(1-t)*t) + b*(t*t);}
            std::vector<double> s(M+1,0);
            for(int i=1;i<=M;++i) s[i]=s[i-1]+len2(crv[i]-crv[i-1]);
            double arc=s[M];
            int np=std::max(4,(int)std::ceil(arc/pstep));
            for(int i=0;i<=np;++i){
                double target=(double)i/np*arc;
                int lo=0,hi=M; while(lo<hi){int mid=(lo+hi)/2; if(s[mid]>=target)hi=mid;else lo=mid+1;}
                traj2.push_back({crv[lo], 1});
            }
        };

        // 蛇形端点
        struct Seg2 { Vec2 a,b; };
        std::vector<Seg2> segs2;
        for(size_t i=0;i<bestRows.size();++i){
            Vec2 a{bestRows[i].xs, bestRows[i].yc}, b{bestRows[i].xe, bestRows[i].yc};
            if(i%2==1) std::swap(a,b);
            segs2.push_back({a,b});
        }
        for(size_t i=0;i<segs2.size();++i){
            sampleLine(segs2[i].a, segs2[i].b, 1);
            if(i+1<segs2.size()){
                double gap = len2(segs2[i+1].a - segs2[i].b);
                if(gap > turnR*3){
                    // 跨段大跳跃（如环形面的左右段切换、跨行大间隔）：
                    // 直线空移过渡，spray=0（不喷涂）。WPF 按 spray=0 断开显示。
                    sampleLine(segs2[i].b, segs2[i+1].a, 0);
                } else {
                    Vec2 dir=norm2(segs2[i].b - segs2[i].a);
                    Vec2 cr{-dir.y, dir.x};
                    Vec2 apex=(segs2[i].b + segs2[i+1].a)*0.5 + cr*turnR;
                    sampleArc(segs2[i].b, segs2[i+1].a, apex);
                }
            }
        }

        // 转 3D（用该平面自己的法向抬升，避免撞击凸起/侧壁）
        Vec3 gun = pl.n*(-1.0);
        for(auto& [p2, sp] : traj2){
            Vec2 q = rotInv(p2);
            Vec3 X = pl.c + u*q.x + v*q.y;
            Vec3 pos = X + pl.n*params->standoff;
            allPts.push_back({pos.x,pos.y,pos.z, gun.x,gun.y,gun.z, params->speed, sp});
        }
    }

    if(allPts.empty()){ err="无有效扫描行"; res.error=_strdup(err.c_str()); return res; }
    res.rows = totalRows;

    // 写回结果数组
    int cnt=(int)allPts.size();
    res.points=(SprayTrajPoint*)malloc(sizeof(SprayTrajPoint)*cnt);
    res.count=cnt;
    std::memcpy(res.points, allPts.data(), sizeof(SprayTrajPoint)*cnt);

    // 自检：点距（仅对连续段内；不同平面段之间有跳跃，跳过跨段检测）
    double sm=0,smax=0; int segCnt=0;
    for(int i=1;i<cnt;++i){
        Vec3 a{res.points[i].x,res.points[i].y,res.points[i].z};
        Vec3 b{res.points[i-1].x,res.points[i-1].y,res.points[i-1].z};
        double d=len(a-b);
        // 跨段跳跃（> 2×行距+点距）不计入等距统计
        if(d > L*2.0 + pstep) continue;
        sm+=d; smax=std::max(smax,d); segCnt++;
    }
    res.seg_len_mean = segCnt>0? sm/segCnt:0;
    res.seg_len_max = smax;
    return res;
}

SPRAYTRAJ_API void spray_traj_free(SprayTrajResult* res){
    if(!res) return;
    if(res->points) free(res->points);
    if(res->error)  free((void*)res->error);
    res->points=nullptr; res->error=nullptr; res->count=0;
}

} // extern "C"
