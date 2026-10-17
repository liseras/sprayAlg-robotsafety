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
// 关键设计（解决"底面被切成多块、范围重叠"）：
//   顺序 RANSAC 若用很紧的 tol（flatnessTol，毫米级）去切一个本身有轻微
//   起伏/弯曲、但应视为"单一可喷面"的底面，会把它沿弯曲方向切成 3~4 块，
//   各块用整块 P0.5/P99.5 分位数包围盒渲染后视觉互相重叠。
//   因此这里采用"先整体、后分割"两级判定。
//   —— 判据必须用高分位残差（P99.5）而非 RMS：RMS 是平方平均，槽型盘的
//      凸边、折型板的翻边只占少数点但高度差大（15~35mm），会被大面积平坦
//      底面"平均"掉，导致整体 RMS 很小（<5mm）却漏掉真实台阶、喷不全/撞枪。
//      P99.5 能稳定捕捉到这些"少数但高"的台阶点。
//   两级：
//     1) 整体平面拟合后，若残差 P99.5 <= curvatureTol，说明整片点云都在
//        曲面容差内（无真实台阶/凸边）→ 单一可喷面，直接返回，不切分；
//     2) 否则存在真实多面/台阶（折型板翻边、槽型盘凸边）→ 用紧 tol
//        （flatnessTol）顺序 RANSAC 切分，各段用各自内点的紧凑包围盒。
struct PlaneSeg { Plane pl; std::vector<int> inl; };
std::vector<PlaneSeg> split_planes(const std::vector<Vec3>& P, double tol, double curvatureTol, int minPts){
    std::vector<PlaneSeg> out;

    // ---- 第一级：整体是否可当作"单一可喷面"？用 P99.5 残差判定 ----
    {
        std::vector<int> all(P.size()); std::iota(all.begin(), all.end(), 0);
        Vec3 nf, cf;
        if(ls_plane(P, all, nf, cf)){
            if(dot(nf, Vec3{0,0,1}) < 0) nf = nf*(-1.0);
            double df = -dot(nf, cf);
            std::vector<double> res; res.reserve(all.size());
            double rss=0; for(int k:all){ double e=dot(P[k],nf)+df; res.push_back(std::fabs(e)); rss+=e*e; }
            double rms = std::sqrt(rss/std::max(1,(int)all.size()));
            // 高分位残差 P99.5：捕捉少数但高的台阶/凸边点
            double p995 = 0.0;
            {
                size_t kk = (size_t)((res.size()-1)*0.995);
                std::nth_element(res.begin(), res.begin()+kk, res.end());
                p995 = res[kk];
            }
            // 整体在曲面容差内（无真实台阶）→ 单面，直接返回（不切分、不重叠）
            if(p995 <= curvatureTol){
                PlaneSeg seg; seg.pl = {nf, df, cf}; seg.inl = all;
                fprintf(stderr, "[split] 整体单面: 点数=%zu RMS=%.3f P99.5=%.3f (<=curvatureTol=%.3f)\n",
                    all.size(), rms, p995, curvatureTol);
                out.push_back(std::move(seg));
                return out;
            }
            fprintf(stderr, "[split] 整体 P99.5=%.3f > curvatureTol=%.3f，存在台阶/凸边，进入多面分割\n",
                p995, curvatureTol);
        }
    }

    // ---- 第二级：真实多面/台阶 → 顺序 RANSAC 切分 ----
    std::vector<int> remain(P.size()); std::iota(remain.begin(), remain.end(), 0);
    // 最多 8 个平面，防止死循环
    for(int iter=0; iter<8 && (int)remain.size() >= minPts; ++iter){
        // 在剩余点集上跑 RANSAC 找最大平面
        Plane pl; std::vector<int> inl;
        std::vector<Vec3> sub; sub.reserve(remain.size());
        for(int i:remain) sub.push_back(P[i]);
        if(!fit_plane(sub, tol, pl, inl)) break;
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

    // 合并"法向一致且各自内点合并后 RMS 仍 <= curvatureTol"的相邻段：
    // 仅当两块是同一朝向、且合并后整体仍是浅曲面时才并，避免把真实台阶抹平。
    // 注意：不依赖质心距离（折型板两翻边质心可能相距很远但确是同一底面的延伸）。
    const double cosThresh = std::cos(15.0*3.141592653589793/180.0);
    bool merged = true;
    while(merged){
        merged = false;
        for(size_t i=0;i<out.size() && !merged;++i){
            for(size_t j=i+1;j<out.size() && !merged;++j){
                double cd = std::fabs(dot(out[i].pl.n, out[j].pl.n));
                if(cd < cosThresh) continue;       // 夹角大，不同朝向

                std::vector<int> uni = out[i].inl;
                uni.insert(uni.end(), out[j].inl.begin(), out[j].inl.end());
                Vec3 nf, cf;
                if(!ls_plane(P, uni, nf, cf)) continue;
                if(dot(nf, out[i].pl.n) < 0) nf = nf*(-1.0);
                double df = -dot(nf, cf);
                double rss=0; for(int k:uni){ double e=dot(P[k],nf)+df; rss+=e*e; }
                double rms = std::sqrt(rss/std::max(1,(int)uni.size()));
                if(rms <= curvatureTol){
                    out[i].pl.n=nf; out[i].pl.c=cf; out[i].pl.d=df;
                    out[i].inl=uni;
                    out.erase(out.begin()+j);
                    merged = true;
                }
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
    struct Row { double xs, xe, yc; };
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
        for(double yc=mny+L*0.25; yc <= mxy+L*0.25+1e-9; yc+=L){
            // 行端固定为分位数包围盒全长（分割出的面是规则平面，全长行走刀）
            int sup=0;
            for(auto&p:R){ if(std::fabs(p.y-yc)<=L*0.5) sup++; }
            if(sup<supMin) continue;          // 跳过完全空白的带
            if(mxx-mnx < L*0.2) continue;     // 面整体过窄
            rs.push_back({mnx,mxx,yc});
        }
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
        if(bestRows.empty()) continue;
        totalRows += (int)bestRows.size();

        double ca=std::cos(theta), sa=std::sin(theta);
        auto rotInv=[&](Vec2 p)->Vec2{ return {ca*p.x - sa*p.y, sa*p.x + ca*p.y}; };

        // 生成该平面的 2D 轨迹
        std::vector<Vec2> traj2;
        auto sampleLine=[&](Vec2 a, Vec2 b){
            double LL=len2(b-a);
            int np=std::max(2,(int)std::ceil(LL/pstep)+1);
            for(int t=0;t<np;++t){ double f=(double)t/(np-1); traj2.push_back(a+(b-a)*f); }
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
                traj2.push_back(crv[lo]);
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
            sampleLine(segs2[i].a, segs2[i].b);
            if(i+1<segs2.size()){
                Vec2 dir=norm2(segs2[i].b - segs2[i].a);
                Vec2 cr{-dir.y, dir.x};
                Vec2 apex=(segs2[i].b + segs2[i+1].a)*0.5 + cr*turnR;
                sampleArc(segs2[i].b, segs2[i+1].a, apex);
            }
        }

        // 转 3D（用该平面自己的法向抬升，避免撞击凸起/侧壁）
        Vec3 gun = pl.n*(-1.0);
        for(auto& p2 : traj2){
            Vec2 q = rotInv(p2);
            Vec3 X = pl.c + u*q.x + v*q.y;
            Vec3 pos = X + pl.n*params->standoff;
            allPts.push_back({pos.x,pos.y,pos.z, gun.x,gun.y,gun.z, params->speed, 1});
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
