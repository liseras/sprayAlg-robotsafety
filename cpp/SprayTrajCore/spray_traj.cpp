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
#include <functional>
#include <map>
#include <array>
#include <chrono>
#include <thread>
#include <atomic>

// ---------------- 简易并行 for ----------------
// 逐点阶段（SOR 邻域、法向估计、法向平滑）彼此完全独立，是最容易并行的部分。
// 按静态分块切给若干线程，每块内顺序执行 —— 结果与串行完全一致（每块写的
// 是自己的下标区间），不引入任何数值差异。
template<class F>
static void parallel_for(int n, F fn){
    if(n <= 0) return;
    unsigned hc = std::thread::hardware_concurrency();
    int T = (int)std::min<unsigned>(hc ? hc : 1u, 16u);
    // 小规模直接串行：线程创建/唤醒的开销比计算还大
    if(T <= 1 || n < 20000){ for(int i=0;i<n;++i) fn(i); return; }
    int chunk = (n + T - 1) / T;
    std::vector<std::thread> th;
    th.reserve(T);
    for(int t=0; t<T; ++t){
        int b = t*chunk, e = std::min(n, b+chunk);
        if(b >= e) break;
        th.emplace_back([&fn, b, e]{ for(int i=b;i<e;++i) fn(i); });
    }
    for(auto& x : th) x.join();
}

// ---------------- 分阶段计时 ----------------
// SPRAYTRAJ_PROF=1 时打印每个阶段的墙钟耗时，用于性能分析。
// 用高精度时钟 + 累加器：同一阶段被多次调用（如逐面重建）会合并统计。
struct ProfStage {
    const char* name;
    double* acc;
    std::chrono::steady_clock::time_point t0;
    explicit ProfStage(const char* n, double* a=nullptr) : name(n), acc(a), t0(std::chrono::steady_clock::now()) {}
    ~ProfStage(){
        double ms = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - t0).count();
        if(acc) *acc += ms;
        if(std::getenv("SPRAYTRAJ_PROF")) fprintf(stderr, "[prof] %-22s %8.1f ms\n", name, ms);
    }
};

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

// ---------------- KD 树（性能关键路径） ----------------
// 旧实现两个热点，在大点云(60万点)上是主要耗时来源：
//   1) 构建：每层 std::sort 全排 + 左右子数组各自拷贝 → O(n log²n) 且产生
//      2n 次 vector 分配。改成对同一份索引数组做 nth_element 原地中位切分，
//      O(n log n)、零分配。
//   2) 查询：每次查询 new 一个 vector 当结果桶，插入时对 k 个元素整体 sort，
//      距离用带 sqrt 的欧氏距离。
//
// 第二轮优化的两个点（查询是内存延迟瓶颈，不是算力瓶颈）：
//   a) 节点内联坐标：原来每访问一个节点要两次随机访存（节点表 → 点表），
//      现在坐标直接在节点里，一次访存就够；
//   b) 结果桶从"有序数组插入"改成定长最大堆：插入只需 O(log k) 次比较/交换，
//      而旧写法每次插入要 memmove 上百字节。邻居顺序对下游（PCA 拟合、
//      区域生长）无意义，因此不再排序；
//   c) 切分轴按"该子集跨度最大的轴"选，不再 x/y/z 循环。**这是查询慢的主因。**
//      线扫点云是薄壳（比如 1030×487×81mm），z 向跨度极小；循环选轴的写法
//      每三层就有一次拿 z 当切分轴，而 z 上几乎所有点挤在同一个薄片里，
//      这个节点等于没切开——子树仍然巨大，剪枝完全失效，每次查询要遍历
//      上千个节点。按跨度选轴后，切分总是发生在信息量最大的方向上。
//
// 注意：近邻集合与旧实现一致（等距并列时取哪个点可能不同，但 PCA 法向与
// SOR 判定对其不敏感，回归面数/行数不变）。
constexpr int KD_MAXK = 64;   // 单次查询支持的最大 K（栈上桶容量）

struct KnnBuf {
    double d2[KD_MAXK];   // 平方距离（最大堆，堆顶=当前最差）
    int    id[KD_MAXK];
    int    k = 0;         // 目标邻居数
    int    n = 0;         // 已收集数量
    void init(int kk){ k = std::min(kk, KD_MAXK); if(k<1) k=1; n = 0; }
    // 尚未收满时返回 +inf，让剪枝不生效；收满后返回堆顶
    inline double worst() const { return n < k ? 1e300 : d2[0]; }
    inline void push(double dd, int idx){
        if(n < k){
            int i = n++;
            while(i > 0){
                int p = (i-1)>>1;
                if(d2[p] >= dd) break;
                d2[i]=d2[p]; id[i]=id[p]; i=p;
            }
            d2[i]=dd; id[i]=idx;
        } else if(dd < d2[0]) {
            int i = 0;
            for(;;){
                int l = 2*i+1, r = l+1, m = i;
                double mx = dd;
                if(l < n && d2[l] > mx){ m=l; mx=d2[l]; }
                if(r < n && d2[r] > mx){ m=r; }
                if(m == i) break;
                d2[i]=d2[m]; id[i]=id[m]; i=m;
            }
            d2[i]=dd; id[i]=idx;
        }
    }
};

struct KdTree {
    // 坐标内联，避免查询时二次访存；axis 记在节点里（切分轴不再由深度推出来）
    struct Node { double x, y, z; int idx; int left=-1, right=-1; unsigned char axis=0; };
    std::vector<Node> nodes;        // 预分配到点数大小（一节点一点），支持并行构建
    std::vector<int>  ids;          // 构建期的索引工作区（原地切分，不拷贝子数组）
    const std::vector<Vec3>* pts=nullptr;
    std::atomic<int> nextNode{0};   // 节点槽位分配（并行构建时由各子树竞争）

    // lo..hi 为左闭右开区间；返回子树根节点下标
    int build_rec(int lo, int hi){
        if(lo >= hi) return -1;
        // 选该子集跨度最大的轴：避免在几乎不变形的方向上白切一刀
        double mn[3]={1e300,1e300,1e300}, mx[3]={-1e300,-1e300,-1e300};
        for(int t=lo;t<hi;++t){
            const Vec3& p=(*pts)[ids[t]];
            if(p.x<mn[0]) mn[0]=p.x; if(p.x>mx[0]) mx[0]=p.x;
            if(p.y<mn[1]) mn[1]=p.y; if(p.y>mx[1]) mx[1]=p.y;
            if(p.z<mn[2]) mn[2]=p.z; if(p.z>mx[2]) mx[2]=p.z;
        }
        double e0=mx[0]-mn[0], e1=mx[1]-mn[1], e2=mx[2]-mn[2];
        int axis = (e0>=e1 && e0>=e2) ? 0 : (e1>=e2 ? 1 : 2);

        int mid = lo + (hi - lo)/2;
        auto cmp = [&](int a, int b){
            double va = (axis==0)?(*pts)[a].x : (axis==1)?(*pts)[a].y:(*pts)[a].z;
            double vb = (axis==0)?(*pts)[b].x : (axis==1)?(*pts)[b].y:(*pts)[b].z;
            return va < vb;
        };
        std::nth_element(ids.begin()+lo, ids.begin()+mid, ids.begin()+hi, cmp);
        int ni = nextNode.fetch_add(1);
        int id = ids[mid];
        const Vec3& p = (*pts)[id];
        nodes[ni] = Node{p.x, p.y, p.z, id, -1, -1, (unsigned char)axis};
        // 左右子树规模都够大时并行构建：两侧互不相干（各自操作 ids 的不相交区间），
        // 树的结构与串行完全一致，只有节点槽位分配顺序不同（内部实现细节）。
        if(hi - lo > 60000){
            std::thread th([&]{ nodes[ni].left = build_rec(lo, mid); });
            nodes[ni].right = build_rec(mid+1, hi);
            th.join();
        } else {
            nodes[ni].left  = build_rec(lo,  mid);
            nodes[ni].right = build_rec(mid+1, hi);
        }
        return ni;
    }
    void init(const std::vector<Vec3>& p){
        pts=&p;
        // 节点数恒等于点数，预分配后并行构建就不必再扩容
        nodes.assign(p.size(), Node{});
        nextNode.store(0);
        ids.resize(p.size()); std::iota(ids.begin(),ids.end(),0);
        if(p.empty()){ std::vector<int>().swap(ids); return; }
        build_rec(0, (int)p.size());
        std::vector<int>().swap(ids);   // 查询阶段不再需要
    }

    void knn(Vec3 q, int k, KnnBuf& b) const {
        b.init(k);
        rec(0, q, b);
    }
    void rec(int node, const Vec3& q, KnnBuf& b) const {
        if(node < 0) return;
        const Node& nd = nodes[node];
        double dx=q.x-nd.x, dy=q.y-nd.y, dz=q.z-nd.z;
        b.push(dx*dx+dy*dy+dz*dz, nd.idx);
        int axis = nd.axis;
        double diff = (axis==0)? dx : (axis==1)? dy : dz;
        int near = diff<0 ? nd.left : nd.right;
        int far  = diff<0 ? nd.right : nd.left;
        rec(near, q, b);
        // 用平方距离剪枝：diff² < 当前最差平方距离才需要搜另一侧
        if(diff*diff < b.worst()) rec(far, q, b);
    }

    // K 近邻距离和（SOR 用，需要真实距离而非平方，故此处才开方）
    double knnSum(Vec3 q, int k) const {
        KnnBuf b; knn(q, k, b);
        double s=0;
        for(int i=0;i<b.n;++i) if(b.id[i]>=0) s += std::sqrt(b.d2[i]);
        return s;
    }
    // K 近邻索引（升序）
    void knnIdx(Vec3 q, int k, std::vector<int>& out) const {
        KnnBuf b; knn(q, k, b);
        std::vector<std::pair<double,int>> t;
        t.reserve(b.n);
        for(int i=0;i<b.n;++i) if(b.id[i]>=0) t.push_back({b.d2[i], b.id[i]});
        std::sort(t.begin(), t.end());
        out.clear(); out.reserve(t.size());
        for(auto& pr : t) out.push_back(pr.second);
    }
    // 零分配版本：直接写入调用方给定的 k 个 int 槽位（扁平数组用）。顺序未定义。
    int knnIdxInto(Vec3 q, int k, int* out) const {
        KnnBuf b; knn(q, k, b);
        int kk = std::min(k, KD_MAXK);
        for(int i=0;i<kk;++i) out[i] = -1;      // 邻居不足时留 -1，调用方跳过
        int n=0;
        for(int i=0;i<b.n;++i) if(b.id[i]>=0) out[n++]=b.id[i];
        return n;
    }
};

// ---------------- 主实现 ----------------
namespace {

struct TrajPoint { Vec3 pos; Vec3 dir; double speed; int spray; };

// SOR 去噪
std::vector<Vec3> sor_filter(const std::vector<Vec3>& P, int k, double z){
    // KnnBuf 是栈上定长桶，K 超过容量会被 init() 截断；这里显式夹一次并提示，
    // 避免"用户传了 sorK=80 却按 64 算"这种静默偏差。
    if(k > KD_MAXK){ fprintf(stderr, "[去噪] sorK=%d 超过上限，按 %d 计算\n", k, KD_MAXK); k = KD_MAXK; }
    if(k < 1) k = 1;
    KdTree tree;
    { ProfStage ps("1a.KD树构建"); tree.init(P); }
    std::vector<double> md(P.size());
    { ProfStage ps("1b.KNN邻域查询");
      parallel_for((int)P.size(), [&](int i){ md[i] = tree.knnSum(P[i], k)/(double)k; }); }
    double sum=0; for(double v:md) sum+=v;   // 串行求和，保证与单线程结果逐位一致
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
// 指针版本：供扁平 KNN 邻接表直接使用，避免每次查询都构造临时 vector
static bool ls_plane_p(const std::vector<Vec3>& P, const int* idx, int cnt,
                       Vec3& n, Vec3& c)
{
    if(cnt<3) return false;
    Vec3 m{0,0,0};
    for(int t=0;t<cnt;++t){ int i=idx[t]; if(i<0) continue; m.x+=P[i].x;m.y+=P[i].y;m.z+=P[i].z; }
    m=m*(1.0/(double)cnt);
    double A00=0,A01=0,A02=0,A11=0,A12=0,A22=0;
    for(int t=0;t<cnt;++t){ int i=idx[t]; if(i<0) continue;
        double dx=P[i].x-m.x,dy=P[i].y-m.y,dz=P[i].z-m.z;
        A00+=dx*dx; A01+=dx*dy; A02+=dx*dz; A11+=dy*dy; A12+=dy*dz; A22+=dz*dz; }
    double invn=1.0/(double)cnt;
    A00*=invn;A01*=invn;A02*=invn;A11*=invn;A12*=invn;A22*=invn;
    // 最小特征向量 = 法向（逆迭代）
    // 性能注意：A 在整个迭代过程中是不变的，所以逆矩阵只需算一次；
    // 旧写法把 det 和 6 个伴随矩阵元素放进循环里，并且用 `|det|<1e-20` 当退出条件
    // ——det 是循环不变量，这个条件永远为假，于是每次都老老实实跑满 80 轮。
    // 现在改为按特征向量收敛判据提前退出（对称阵逆迭代收敛很快，通常 4~8 轮）。
    double det = A00*(A11*A22-A12*A12) - A01*(A01*A22-A12*A02) + A02*(A01*A12-A11*A02);
    if(std::fabs(det)<1e-20){ n=Vec3{0,0,1}; c=m; return true; }   // 平凡/退化邻域
    const double i00=(A11*A22-A12*A12)/det;
    const double i01=(A02*A12-A01*A22)/det;
    const double i02=(A01*A12-A11*A02)/det;
    const double i11=(A00*A22-A02*A02)/det;
    const double i12=(A01*A02-A00*A12)/det;
    const double i22=(A00*A11-A01*A01)/det;
    Vec3 x{0,0,1};
    for(int inner=0; inner<80; ++inner){
        Vec3 y{ i00*x.x+i01*x.y+i02*x.z,
                i01*x.x+i11*x.y+i12*x.z,
                i02*x.x+i12*x.y+i22*x.z };
        double l=std::sqrt(dot(y,y));
        if(l<1e-12) break;
        Vec3 xn{y.x/l,y.y/l,y.z/l};
        // 收敛判据取"定点"：方向逐分量不再变化时，再迭代也只会得到同一个数。
        // 用 |dot|≈1 之类的相对阈值会停在 1e-13 量级，虽然肉眼无差别，
        // 但会在少数点的最后一位上留下漂移（回归时表现为个别坐标差 0.1mm 级）。
        bool fixed = (xn.x==x.x && xn.y==x.y && xn.z==x.z);
        x = xn;
        if(fixed) break;
    }
    n=x; c=m;
    return true;
}

// 容器版本（兼容原有调用点）
static bool ls_plane(const std::vector<Vec3>& P, const std::vector<int>& idx,
                     Vec3& n, Vec3& c)
{
    return ls_plane_p(P, idx.data(), (int)idx.size(), n, c);
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

    // 性能改造：先一次性生成全部假设，再并行打分，最后按原顺序取第一个最优。
    //
    // 原实现是"每轮 RANSAC 都全量扫一遍点云并 push_back 索引"，200 轮 × 60 万点
    // = 1.2 亿次 push_back + 反复扩容，是大点云上最大的一块单点开销。
    //
    // 等价性：原循环用 `cur.size() > best.size()` 严格大于才更新，所以最终选中的
    // 就是"按迭代顺序第一个达到全局最大内点数的假设"。这里并行算完所有假设的
    // 内点数后，同样按迭代顺序取第一个严格更优者 —— 结果完全相同。
    // RNG 抽样顺序也保持一致（无效样本仍照常消耗随机数），故抽样点集也相同。
    struct Hyp { Vec3 n; double d; };
    std::vector<Hyp> hyps;
    hyps.reserve(iters);
    for(int it=0; it<iters; ++it){
        int a=ud(rng), b=ud(rng), c2=ud(rng);
        if(a==b||b==c2||a==c2) continue;
        Vec3 n = norm(cross(P[b]-P[a], P[c2]-P[a]));
        if(len(n)<1e-6) continue;             // 三点共线，跳过
        hyps.push_back({n, -dot(n, P[a])});
    }

    if(!hyps.empty()){
        std::vector<int> cnt(hyps.size(), 0);
        parallel_for((int)hyps.size(), [&](int k){
            const Vec3& n = hyps[k].n; double d = hyps[k].d;
            int c=0;
            for(int i=0;i<N;++i) if(std::fabs(dot(P[i],n)+d) < tol) ++c;
            cnt[k]=c;
        });
        int bestK=-1, bestCnt=0;
        for(size_t k=0;k<hyps.size();++k)
            if(cnt[k] > bestCnt){ bestCnt=cnt[k]; bestK=(int)k; }
        if(bestK >= 0){
            bn = hyps[bestK].n;
            const Vec3& n = hyps[bestK].n; double d = hyps[bestK].d;
            best.reserve(bestCnt);
            for(int i=0;i<N;++i) if(std::fabs(dot(P[i],n)+d) < tol) best.push_back(i);
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

// 多平面分割：按"高度分层"替代顺序 RANSAC 切块。
// 顺序 RANSAC 会把凸边环（4边相连、同一高度）切成多块碎面（段1/2/3），
// 各自生成零散轨迹、面片重叠——视觉"乱七八糟"。
// 高度分层思路（通用、鲁棒）：
//   1) 先 RANSAC 拟合主平面得到法向 n（即"高度方向"）；
//   2) 全点投影到 n，得到高度分布，直方图聚类分成若干高度层；
//   3) 每层 = 该高度带内的所有点（无论实心还是环），各拟合一个平面、
//      各自生成轨迹。槽型盘 → 底盘层 + 凸边环层 = 2 面。
//
// 参数语义：
//   tol (flatnessTol)  = 高度层之间的最小间隔：相邻高度层峰值间距 < tol 视为同层
//                        （台阶抹平阈值，与"曲率浅当平面"一致）；
//   curvatureTol       = 整体高度极差 < 此值 → 单一可喷面（不切分）。
struct PlaneSeg {
    Plane pl;
    std::vector<int> inl;
    // 曲面块标记：true 时 pl 只是"平均平面"（用于定向/统计），真正的轨迹
    // 生成走曲面行追踪（gen_surface_rows），不是平面 2D 布行。
    bool curved = false;
};

// ---------------- 规则化轮廓重建（2D） ----------------
// 目的：把"面内点集"重建成边为直边的规则多边形（矩形优先），而不是逐 cell
// 打点的锯齿拼块。流程：细网格占用 → 形态学闭运算 → 按跨距判定填补孔洞 →
// 边界追踪 → RDP 直线简化 → 最小外接矩形吸附 → 三角化（含孔洞桥接）。
// 所有坐标为面局部基 (u,v) 下的 mm。

static double polyArea2(const std::vector<Vec2>& p){
    double a=0; size_t n=p.size();
    for(size_t i=0;i<n;++i){ const Vec2& q=p[i]; const Vec2& r=p[(i+1)%n]; a += q.x*r.y - r.x*q.y; }
    return a*0.5;
}

// 开链 Ramer-Douglas-Peucker：把抖动点列压成直边折线
static void rdpOpen(const std::vector<Vec2>& p, size_t i0, size_t i1, double eps, std::vector<char>& keep){
    if(i1<=i0+1) return;
    Vec2 a=p[i0], b=p[i1], ab=b-a;
    double L=len2(ab), dmax=-1; size_t idx=i0;
    for(size_t i=i0+1;i<i1;++i){
        double d = (L<1e-12)? len2(p[i]-a) : std::fabs(cross2(ab, p[i]-a))/L;
        if(d>dmax){ dmax=d; idx=i; }
    }
    if(dmax>eps){ keep[idx]=1; rdpOpen(p,i0,idx,eps,keep); rdpOpen(p,idx,i1,eps,keep); }
}

// 闭合轮廓简化：取最远两点切成两条开链分别 RDP
static std::vector<Vec2> simplifyClosed(const std::vector<Vec2>& p, double eps){
    size_t n=p.size();
    if(n<5) return p;
    size_t far=0; double dmax=-1;
    for(size_t i=1;i<n;++i){ double d=len2(p[i]-p[0]); if(d>dmax){ dmax=d; far=i; } }
    std::vector<char> keep(n,0); keep[0]=1; keep[far]=1;
    {
        std::vector<Vec2> c1(p.begin(), p.begin()+far+1);
        std::vector<char> k(c1.size(),0); k.front()=1; k.back()=1;
        rdpOpen(c1,0,c1.size()-1,eps,k);
        for(size_t i=0;i<c1.size();++i) if(k[i]) keep[i]=1;
    }
    {
        std::vector<Vec2> c2; c2.reserve(n-far+1);
        for(size_t i=far;i<n;++i) c2.push_back(p[i]);
        c2.push_back(p[0]);
        std::vector<char> k(c2.size(),0); k.front()=1; k.back()=1;
        rdpOpen(c2,0,c2.size()-1,eps,k);
        for(size_t i=0;i+1<c2.size();++i) if(k[i]) keep[far+i]=1;
    }
    std::vector<Vec2> out;
    for(size_t i=0;i<n;++i) if(keep[i]) out.push_back(p[i]);
    return out.size()>=3? out : p;
}

static std::vector<Vec2> convexHull(std::vector<Vec2> p){
    std::sort(p.begin(), p.end(), [](const Vec2&a, const Vec2&b){
        return a.x<b.x || (a.x==b.x && a.y<b.y); });
    p.erase(std::unique(p.begin(), p.end(), [](const Vec2&a, const Vec2&b){
        return std::fabs(a.x-b.x)<1e-9 && std::fabs(a.y-b.y)<1e-9; }), p.end());
    size_t n=p.size();
    if(n<3) return p;
    std::vector<Vec2> h(2*n);
    size_t k=0;
    for(size_t i=0;i<n;++i){
        while(k>=2 && cross2(h[k-1]-h[k-2], p[i]-h[k-2]) <= 1e-9) --k;
        h[k++]=p[i];
    }
    for(size_t i=n-1, t=k+1; i-->0; ){
        while(k>=t && cross2(h[k-1]-h[k-2], p[i]-h[k-2]) <= 1e-9) --k;
        h[k++]=p[i];
    }
    h.resize(k-1);
    return h;
}

// 最小外接矩形（旋转卡壳）：返回面积，rect 为 4 个角点（CCW）
static double minAreaRect(const std::vector<Vec2>& pts, std::vector<Vec2>& rect){
    std::vector<Vec2> h = convexHull(pts);
    if(h.size()<3) return -1.0;
    double best=1e18;
    for(size_t i=0;i<h.size();++i){
        Vec2 d = h[(i+1)%h.size()] - h[i];
        double L = len2(d); if(L<1e-9) continue;
        Vec2 e{d.x/L, d.y/L}, o{-d.y/L, d.x/L};
        double umin=1e18,umax=-1e18,vmin=1e18,vmax=-1e18;
        for(auto& q : h){
            double u=dot2(q,e), v=dot2(q,o);
            umin=std::min(umin,u); umax=std::max(umax,u);
            vmin=std::min(vmin,v); vmax=std::max(vmax,v);
        }
        double area=(umax-umin)*(vmax-vmin);
        if(area<best){
            best=area;
            rect = { e*umin+o*vmin, e*umax+o*vmin, e*umax+o*vmax, e*umin+o*vmax };
        }
    }
    return best;
}

static bool pointInTriStrict(Vec2 p, Vec2 a, Vec2 b, Vec2 c){
    double d1=cross2(b-a,p-a), d2=cross2(c-b,p-b), d3=cross2(a-c,p-c);
    bool neg = (d1<0)||(d2<0)||(d3<0);
    bool pos = (d1>0)||(d2>0)||(d3>0);
    return !(neg&&pos);
}

// 耳切三角化（凸多边形/凹多边形通用；带重复桥接点的退化情形按"共点不算遮挡"处理）
static bool earClip(const std::vector<Vec2>& poly, std::vector<int>& tris){
    int n=(int)poly.size();
    if(n<3) return false;
    std::vector<int> idx(n); std::iota(idx.begin(), idx.end(), 0);
    if(polyArea2(poly) < 0) std::reverse(idx.begin(), idx.end());
    int guard=0, limit=n*n+32;
    while((int)idx.size()>3 && guard<limit){
        bool clipped=false;
        int m=(int)idx.size();
        for(int i=0;i<m;++i){
            int i0=idx[(i+m-1)%m], i1=idx[i], i2=idx[(i+1)%m];
            Vec2 a=poly[i0], b=poly[i1], c=poly[i2];
            if(cross2(b-a,c-a) <= 1e-9) continue;      // 凹角/共线
            bool block=false;
            for(int k=0;k<m && !block;++k){
                int ki=idx[k];
                if(ki==i0||ki==i1||ki==i2) continue;
                if(len2(poly[ki]-a)<1e-9||len2(poly[ki]-b)<1e-9||len2(poly[ki]-c)<1e-9) continue;
                if(pointInTriStrict(poly[ki],a,b,c)) block=true;
            }
            if(block) continue;
            tris.push_back(i0); tris.push_back(i1); tris.push_back(i2);
            idx.erase(idx.begin()+i);
            clipped=true; ++guard; break;
        }
        if(!clipped) return false;
    }
    if(idx.size()==3){ tris.push_back(idx[0]); tris.push_back(idx[1]); tris.push_back(idx[2]); }
    return !tris.empty();
}

// 把孔洞桥接到外轮廓（孔的最右点向右投射，选最近且可见的外轮廓顶点）
static bool bridgeHole(std::vector<Vec2>& outer, const std::vector<Vec2>& hole){
    int hn=(int)hole.size(), on=(int)outer.size();
    if(hn<3||on<3) return false;
    int mi=0;
    for(int i=1;i<hn;++i) if(hole[i].x > hole[mi].x + 1e-12) mi=i;
    Vec2 M=hole[mi];
    auto crossT=[](Vec2 p, Vec2 q, Vec2 r){ return cross2(q-p, r-p); };
    auto segHit=[&](Vec2 a, Vec2 b, Vec2 c, Vec2 d)->bool{
        double d1=crossT(a,b,c), d2=crossT(a,b,d), d3=crossT(c,d,a), d4=crossT(c,d,b);
        const double e=1e-9;
        return ((d1>e&&d2<-e)||(d1<-e&&d2>e)) && ((d3>e&&d4<-e)||(d3<-e&&d4>e));
    };
    auto visible=[&](Vec2 V)->bool{
        for(int i=0;i<on;++i){
            Vec2 a=outer[i], b=outer[(i+1)%on];
            if(len2(a-V)<1e-9||len2(b-V)<1e-9) continue;
            if(segHit(M,V,a,b)) return false;
        }
        for(int i=0;i<hn;++i){
            Vec2 a=hole[i], b=hole[(i+1)%hn];
            if(len2(a-M)<1e-9||len2(b-M)<1e-9) continue;
            if(segHit(M,V,a,b)) return false;
        }
        return true;
    };
    int best=-1; double bestAng=1e18;
    for(int i=0;i<on;++i){
        Vec2 V=outer[i];
        if(V.x < M.x-1e-9) continue;                 // 优先右侧可见点
        if(!visible(V)) continue;
        double ang=std::atan2(std::fabs(V.y-M.y), std::max(1e-9, V.x-M.x));
        if(ang<bestAng){ bestAng=ang; best=i; }
    }
    if(best<0) return false;                         // 无右侧可见点：放弃该孔
    std::vector<Vec2> ins;                           // M ... M, V  （沿孔走一圈回到 M，再回 V）
    for(int k=0;k<=hn;++k) ins.push_back(hole[(mi+k)%hn]);
    ins.push_back(outer[best]);
    outer.insert(outer.begin()+best+1, ins.begin(), ins.end());
    return true;
}

static bool pointInPoly(Vec2 p, const std::vector<Vec2>& poly){
    bool in=false; size_t n=poly.size();
    for(size_t i=0,j=n-1;i<n;j=i++){
        Vec2 a=poly[i], b=poly[j];
        if(((a.y>p.y)!=(b.y>p.y)) && (p.x < (b.x-a.x)*(p.y-a.y)/((b.y-a.y)==0?1e-12:(b.y-a.y))+a.x)) in=!in;
    }
    return in;
}

// 诊断：重建网格的连通分量统计。
// 用途——识别"稀疏掠射回波被拟合成一个面"：真实可喷面重建出来必然是
// 一整块连续网格；边缘掠射点拟合出的伪面，网格会碎成若干互不相连的小块。
static void mesh_connect_stats(const std::vector<Vec2>& V, const std::vector<int>& T,
                               int& nComp, double& largestFrac, double& totalArea)
{
    nComp = 0; largestFrac = 0; totalArea = 0;
    int n = (int)V.size(), nt = (int)(T.size()/3);
    if(n <= 0 || nt <= 0) return;
    std::vector<int> par(n);
    for(int i=0;i<n;++i) par[i]=i;
    std::function<int(int)> find = [&](int x)->int{
        while(par[x]!=x){ par[x]=par[par[x]]; x=par[x]; }
        return x;
    };
    auto uni = [&](int a,int b){ a=find(a); b=find(b); if(a!=b) par[a]=b; };
    for(int t=0;t<nt;++t){ uni(T[t*3],T[t*3+1]); uni(T[t*3+1],T[t*3+2]); }
    std::map<int,double> area;
    for(int t=0;t<nt;++t){
        const Vec2& a=V[T[t*3]]; const Vec2& b=V[T[t*3+1]]; const Vec2& c=V[T[t*3+2]];
        double ar = 0.5*std::fabs((b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y));
        area[find(T[t*3])] += ar;
        totalArea += ar;
    }
    nComp = (int)area.size();
    double mx = 0;
    for(auto& kv : area) mx = std::max(mx, kv.second);
    largestFrac = totalArea > 1e-9 ? mx/totalArea : 0.0;
}

// 统计二值占用网格的 4-连通块数（诊断用）
static int count_occ_components(const std::vector<unsigned char>& occ, int nU, int nV)
{
    std::vector<char> seen(occ.size(), 0);
    int comp = 0;
    std::vector<int> st;
    const int du[4] = {1,-1,0,0}, dv[4] = {0,0,1,-1};
    for(int k0=0;k0<(int)occ.size();++k0){
        if(!occ[k0] || seen[k0]) continue;
        ++comp;
        seen[k0] = 1; st.push_back(k0);
        while(!st.empty()){
            int k = st.back(); st.pop_back();
            int iu = k % nU, iv = k / nU;
            for(int d=0; d<4; ++d){
                int ju = iu+du[d], jv = iv+dv[d];
                if(ju<0||jv<0||ju>=nU||jv>=nV) continue;
                int kk = jv*nU+ju;
                if(occ[kk] && !seen[kk]){ seen[kk]=1; st.push_back(kk); }
            }
        }
    }
    return comp;
}

// 环形窄面判定（(u,v) 平面上）：绕成一整圈、环宽窄、中间是孔。
// 两处要用同一个判据：
//   ① 轮廓重建的格距 —— 格距 ≥ 环宽时，开运算会把整圈削成几段弧。实测垃圾桶
//      顶环（r143~175、宽 32mm）在行距 30mm 下：格距 30 → 内点 6977 被剔除到
//      3653（只剩 250° 的环），轨迹只能盖住大半个环。
//   ② 轨迹布行 —— 平行扫描在环上要铺 11~12 条弦、3/4 是空走；同心圆 1~2 圈
//      就够（实测同心圆 1046mm vs 平行 3568mm）。
// 四条判据同时成立（缺一条都会误伤矩形面 / 实心圆盘 / 小垫片）：
//   ① 方位覆盖 ≥34/36 格；② 环宽 w=p97-p2 ≤ 2.5×行距；
//   ③ 内径 ≥ 0.45×外径（排除实心圆盘，底面内径≈0）；
//   ④ 外径 ≥ 3×环宽。
static bool detect_ring_uv(const std::vector<Vec2>& uv, double L,
                           Vec2& cOut, double& rin, double& rout,
                           double* wOut = nullptr, int* nsecOut = nullptr)
{
    if(uv.size() < 200) return false;
    Vec2 c{0,0};
    for(auto& q : uv){ c.x += q.x; c.y += q.y; }
    c.x /= (double)uv.size(); c.y /= (double)uv.size();
    std::vector<double> rr; rr.reserve(uv.size());
    for(auto& q : uv) rr.push_back(len2(q-c));
    std::sort(rr.begin(), rr.end());
    const double r2  = rr[(size_t)(0.02*(double)(rr.size()-1))];
    const double r97 = rr[(size_t)(0.97*(double)(rr.size()-1))];
    const double w   = r97 - r2;
    const int NS = 36; std::vector<unsigned char> sec(NS,0);
    for(auto& q : uv){
        Vec2 d{q.x-c.x, q.y-c.y};
        if(len2(d) < 1e-9) continue;
        double th = std::atan2(d.y, d.x);
        int k = (int)((th + 3.14159265358979)/(2.0*3.14159265358979)*(double)NS);
        if(k<0) k=0; if(k>=NS) k=NS-1;
        sec[k]=1;
    }
    int nsec=0; for(int k=0;k<NS;++k) if(sec[k]) ++nsec;
    if(wOut) *wOut = w;
    if(nsecOut) *nsecOut = nsec;
    cOut = c; rin = r2; rout = r97;          // 诊断用：判否时也给出实测值
    const bool ringish = (nsec >= 34) && (w > 1e-6) && (w <= 2.5*L)
                      && (r2 >= 0.45*r97) && (r97 >= 3.0*w);
    return ringish;
}

// 由面内点集重建规则三角网格（局部 uv 坐标，mm）。
//   rowGap   扫描线间距（跨线空隙要用闭运算桥接，否则每行独立成条）
//   fillHole 孔洞填补阈值：连通孔洞的外接跨距 <= 此值 → 填实（覆盖重建）
static bool build_plane_mesh(const std::vector<Vec2>& uv, double rowGap, double fillHole,
                             std::vector<Vec2>& outVerts, std::vector<int>& outTris,
                             int maxVerts, int maxTris,
                             std::vector<char>* keepMask = nullptr,
                             std::vector<Vec2>* fillRectOut = nullptr)
{
    const int n = (int)uv.size();
    if(n < 20) return false;
    double umin=1e18, umax=-1e18, vmin=1e18, vmax=-1e18;
    for(const auto& q : uv){
        umin=std::min(umin,q.x); umax=std::max(umax,q.x);
        vmin=std::min(vmin,q.y); vmax=std::max(vmax,q.y);
    }
    double wU=umax-umin, wV=vmax-vmin;
    if(wU<2.0 || wV<2.0) return false;

    // 1) 细网格占用（目标 ~8000 cell，边长夹在 [2,6]mm，四周留 1 cell 空边）
    double cs = std::sqrt(wU*wV/8000.0);
    cs = std::min(6.0, std::max(2.0, cs));
    int nU=(int)std::ceil(wU/cs)+3, nV=(int)std::ceil(wV/cs)+3;
    std::vector<unsigned char> occ((size_t)nU*nV, 0);
    for(const auto& q : uv){
        int iu=(int)((q.x-umin)/cs)+1, iv=(int)((q.y-vmin)/cs)+1;
        iu=std::min(std::max(iu,0),nU-1); iv=std::min(std::max(iv,0),nV-1);
        occ[(size_t)iv*nU+iu]=1;
    }

    // 2) 形态学闭运算：先膨胀再腐蚀，核半径按扫描线间距取
    //    （线扫点云沿行进方向很密、跨行很疏，不做闭运算每行都是一条独立细条）
    {
        int rr=(int)std::ceil(std::max(rowGap,1.0)/(2.0*cs))+1;
        rr=std::min(rr,6);
        auto dilate=[&](std::vector<unsigned char>& g, int r){
            std::vector<unsigned char> t(g.size(),0), o(g.size(),0);
            for(int iv=0;iv<nV;++iv) for(int iu=0;iu<nU;++iu){
                if(!g[(size_t)iv*nU+iu]) continue;
                int a=std::max(0,iu-r), b=std::min(nU-1,iu+r);
                for(int k=a;k<=b;++k) t[(size_t)iv*nU+k]=1;
            }
            for(int iv=0;iv<nV;++iv){
                int a=std::max(0,iv-r), b=std::min(nV-1,iv+r);
                for(int iu=0;iu<nU;++iu){
                    if(!t[(size_t)iv*nU+iu]) continue;
                    for(int k=a;k<=b;++k) o[(size_t)k*nU+iu]=1;
                }
            }
            g.swap(o);
        };
        dilate(occ, rr);
        std::vector<unsigned char> inv(occ.size(),0);
        for(size_t k=0;k<occ.size();++k) inv[k]=occ[k]?0:1;
        dilate(inv, rr);
        for(size_t k=0;k<occ.size();++k) occ[k]= inv[k]?0:1;   // 腐蚀
        if(const char* d=std::getenv("SPRAYTRAJ_DUMP")) if(d[0]=='1')
            fprintf(stderr, "[mesh] 闭运算后连通块=%d (rr=%d cs=%.2f)\n",
                    count_occ_components(occ,nU,nV), rr, cs);
        // 开运算（先腐蚀后膨胀）：去掉"细尾巴"——相邻面在折棱处混入的
        // 一两条扫描线会拖出细长条，不做开运算会让轮廓长出一根刺。
        {
            std::vector<unsigned char> inv2(occ.size(),0);
            for(size_t k=0;k<occ.size();++k) inv2[k]=occ[k]?0:1;
            dilate(inv2, rr);
            for(size_t k=0;k<occ.size();++k) occ[k]= inv2[k]?0:1;   // 腐蚀
            dilate(occ, rr);                                        // 膨胀回来
        }
        if(const char* d=std::getenv("SPRAYTRAJ_DUMP")) if(d[0]=='1')
            fprintf(stderr, "[mesh] 开运算后连通块=%d\n", count_occ_components(occ,nU,nV));
    }

    // 3) 消除"对角夹点"（两格斜对角相连、正交两侧空）：保证边界是简单闭合环
    for(int pass=0; pass<4; ++pass){
        bool changed=false;
        for(int iv=0;iv+1<nV;++iv) for(int iu=0;iu+1<nU;++iu){
            unsigned char a=occ[(size_t)iv*nU+iu], b=occ[(size_t)iv*nU+iu+1];
            unsigned char c=occ[(size_t)(iv+1)*nU+iu], d=occ[(size_t)(iv+1)*nU+iu+1];
            if(a&&d&&!b&&!c){ occ[(size_t)iv*nU+iu+1]=1; changed=true; }
            else if(b&&c&&!a&&!d){ occ[(size_t)iv*nU+iu]=1; changed=true; }
        }
        if(!changed) break;
    }

    // 4) 孔洞：外部泛洪，未被外部到达的空区即内部孔洞；仅当跨距 <= fillHole 才填
    {
        std::vector<unsigned char> out((size_t)nU*nV,0);
        std::vector<int> st;
        auto push=[&](int iu,int iv){
            if(iu<0||iv<0||iu>=nU||iv>=nV) return;
            size_t k=(size_t)iv*nU+iu;
            if(occ[k]||out[k]) return;
            out[k]=1; st.push_back((int)k);
        };
        for(int iu=0;iu<nU;++iu){ push(iu,0); push(iu,nV-1); }
        for(int iv=0;iv<nV;++iv){ push(0,iv); push(nU-1,iv); }
        while(!st.empty()){
            int k=st.back(); st.pop_back();
            int iu=k%nU, iv=k/nU;
            push(iu+1,iv); push(iu-1,iv); push(iu,iv+1); push(iu,iv-1);
        }
        // 内部空区按连通域统计跨距
        std::vector<char> seen((size_t)nU*nV,0);
        for(size_t k0=0;k0<occ.size();++k0){
            if(occ[k0]||out[k0]||seen[k0]) continue;
            std::vector<int> comp; comp.push_back((int)k0); seen[k0]=1;
            int lo_u=nU,hi_u=0,lo_v=nV,hi_v=0;
            for(size_t t=0;t<comp.size();++t){
                int k=comp[t], iu=k%nU, iv=k/nU;
                lo_u=std::min(lo_u,iu); hi_u=std::max(hi_u,iu);
                lo_v=std::min(lo_v,iv); hi_v=std::max(hi_v,iv);
                const int du[4]={1,-1,0,0}, dv[4]={0,0,1,-1};
                for(int d=0;d<4;++d){
                    int ju=iu+du[d], jv=iv+dv[d];
                    if(ju<0||jv<0||ju>=nU||jv>=nV) continue;
                    size_t kk=(size_t)jv*nU+ju;
                    if(occ[kk]||out[kk]||seen[kk]) continue;
                    seen[kk]=1; comp.push_back((int)kk);
                }
            }
            double spanU=(hi_u-lo_u+1)*cs, spanV=(hi_v-lo_v+1)*cs;
            if(std::max(spanU,spanV) <= fillHole)
                for(int k : comp) occ[(size_t)k]=1;   // 覆盖重建：小洞填实
        }
    }

    // 4.5) 回填"哪些点属于这个重建面"：开运算已经去掉了稀疏细尾巴，
    // 落在被删掉区域里的点（相邻面混入的孤立扫描线）不应继续参与轨迹生成，
    // 否则轨迹会在工件外面拉出梳齿。
    if(keepMask){
        keepMask->assign(uv.size(), 0);
        for(size_t t=0;t<uv.size();++t){
            int iu=(int)((uv[t].x-umin)/cs)+1, iv=(int)((uv[t].y-vmin)/cs)+1;
            if(iu<0||iv<0||iu>=nU||iv>=nV) continue;
            if(occ[(size_t)iv*nU+iu]) (*keepMask)[t]=1;
        }
    }

    // 5) 边界追踪（材料的左边界），外环 CCW、内环 CW
    const int NW=nU+1;
    std::vector<int> nxt((size_t)NW*(nV+1), -1);
    auto addEdge=[&](int iu,int iv,int ju,int jv){
        int a=iv*NW+iu, b=jv*NW+ju;
        if(nxt[(size_t)a]<0) nxt[(size_t)a]=b;
    };
    for(int iv=0;iv<nV;++iv) for(int iu=0;iu<nU;++iu){
        if(!occ[(size_t)iv*nU+iu]) continue;
        bool bEmpty = (iv==0)   || !occ[(size_t)(iv-1)*nU+iu];
        bool rEmpty = (iu==nU-1)|| !occ[(size_t)iv*nU+iu+1];
        bool tEmpty = (iv==nV-1)|| !occ[(size_t)(iv+1)*nU+iu];
        bool lEmpty = (iu==0)   || !occ[(size_t)iv*nU+iu-1];
        if(bEmpty) addEdge(iu,   iv,   iu+1, iv  );   // 下边：+u
        if(rEmpty) addEdge(iu+1, iv,   iu+1, iv+1);   // 右边：+v
        if(tEmpty) addEdge(iu+1, iv+1, iu,   iv+1);   // 上边：-u
        if(lEmpty) addEdge(iu,   iv+1, iu,   iv  );   // 左边：-v
    }
    std::vector<char> usedEdge(nxt.size(),0);
    std::vector<std::vector<Vec2>> loops;
    for(size_t s=0;s<nxt.size();++s){
        if(nxt[s]<0||usedEdge[s]) continue;
        std::vector<Vec2> lp; int cur=(int)s; int guard=0;
        while(cur>=0 && !usedEdge[(size_t)cur] && guard++ < (int)nxt.size()+8){
            usedEdge[(size_t)cur]=1;
            int iu=cur%NW, iv=cur/NW;
            lp.push_back(Vec2{ umin+(iu-1)*cs, vmin+(iv-1)*cs });
            cur=nxt[(size_t)cur];
        }
        if(lp.size()>=4) loops.push_back(std::move(lp));
    }

    // 6) 简化 + 凸包规整 + 矩形吸附
    std::vector<std::vector<Vec2>> outers, holes;
    bool dbg = false;
    if(const char* d=std::getenv("SPRAYTRAJ_DUMP")) dbg = (d[0]=='1');
    if(dbg) fprintf(stderr, "[mesh] 网格 %dx%d cs=%.2f 环数=%zu\n", nU, nV, cs, loops.size());
    if(dbg){
        char nm[64]; snprintf(nm,sizeof(nm),"mask_%dx%d.pgm", nU, nV);
        if(FILE* fm=fopen(nm,"wb")){
            fprintf(fm,"P5\n%d %d\n255\n", nU, nV);
            for(int iv=nV-1; iv>=0; --iv)
                for(int iu=0; iu<nU; ++iu)
                    fputc(occ[(size_t)iv*nU+iu]?255:0, fm);
            fclose(fm);
        }
    }
    // 最大环面积：用于剔除"同一个面上孤立的小碎块"。
    // 现场件边缘是掠射扫描，点云沿边缘断续，形态学闭运算也接不上，
    // 边界追踪就会得到几个互不相连的小环——渲染上表现为一个面上叠了几块补丁。
    // 面积不足主环 10% 的环直接丢弃（这点面积也生成不出有效行）。
    double maxLoopArea = 0.0;
    for(auto& lp : loops) maxLoopArea = std::max(maxLoopArea, std::fabs(polyArea2(lp)));

    for(auto& lp : loops){
        double a0=polyArea2(lp);
        if(dbg) fprintf(stderr, "[mesh]   环: 点=%zu 面积=%.0f\n", lp.size(), a0);
        if(std::fabs(a0) < 9.0*cs*cs) { if(dbg) fprintf(stderr, "[mesh]     -> 太小丢弃\n"); continue; }
        if(loops.size()>1 && std::fabs(a0) < 0.10*maxLoopArea){
            if(dbg) fprintf(stderr, "[mesh]     -> 碎块丢弃(占主环 %.0f%%)\n", std::fabs(a0)/maxLoopArea*100.0);
            continue;
        }
        // 直线简化的容差取"扫描线间距量级"：线扫点云边界天然是"梳齿"状
        // （每行末端位置不同，齿深约等于行距），容差小于行距就会把梳齿保留
        // 下来，看起来像锯子。取 max(2×cell, 1.5×行距) 把梳齿抹成直边。
        std::vector<Vec2> s = simplifyClosed(lp, std::max(2.0*cs, 1.5*rowGap));
        double a=polyArea2(s);
        if(std::fabs(a)<1e-6) continue;

        // 凸包规整：相邻面在折棱处会互相混入几个点，在轮廓上咬出浅缺口。
        // 只有当"轮廓相对凸包的最大凹陷深度 <= max(孔洞填补阈值, 2.5×行距)"时
        // 才用凸包拉平——即只抹平噪声咬边；真实缺口/凹槽（如长板宽板的中间
        // 残缺，深达 100mm）保留原轮廓，不会被"填成矩形"。
        {
            std::vector<Vec2> hull = convexHull(s);
            if(hull.size()>=3){
                double ha = polyArea2(hull);
                double maxDev = 0.0;
                for(const auto& q : s){
                    double dmin=1e18;
                    for(size_t e=0;e<hull.size();++e){
                        Vec2 A=hull[e], B=hull[(e+1)%hull.size()], AB=B-A;
                        double L2=dot2(AB,AB);
                        double t = L2>1e-12? dot2(q-A,AB)/L2 : 0.0;
                        t = std::min(1.0, std::max(0.0, t));
                        Vec2 proj = A + AB*t;
                        dmin = std::min(dmin, len2(q-proj));
                    }
                    maxDev = std::max(maxDev, dmin);
                }
                // 容差取 max(5×行距, 孔洞填补阈值)：行距决定了"扫描线间混入"的
                // 量级，边界的浅缺口都在这个尺度内；真实大缺口（长板宽板的中间
                // 残缺，深 100mm）超出容差，仍然保留原轮廓。
                double devTol = std::max(5.0*rowGap, fillHole);
                if(std::fabs(ha)>1e-9 && maxDev <= devTol && std::fabs(a)/std::fabs(ha) >= 0.80){
                    if((a>0)!=(ha>0)) std::reverse(hull.begin(), hull.end());
                    s = hull;
                    a = polyArea2(s);
                }
            }
        }

        // 矩形吸附：轮廓与其最小外接矩形面积比很高 → 认定是矩形面/矩形孔，
        // 直接用 4 个角点，边界彻底变成直边（这是"规则重建"的关键一步）
        std::vector<Vec2> rect;
        double ra=minAreaRect(s, rect);
        double ratio = (ra>1e-9)? std::fabs(a)/ra : 0.0;
        double thr = (a0>0)? 0.85 : 0.80;
        // 只剩 3 个点的轮廓（小面被简化过头）→ 直接用最小外接矩形，
        // 三角形看起来像"缺口"，矩形才是规整的重建面。
        if(s.size() < 4 && ra > 1e-9){
            if((a>0)!=(polyArea2(rect)>0)) std::reverse(rect.begin(), rect.end());
            s = rect; ratio = 1.0;
        }
        if(dbg) fprintf(stderr, "[mesh]    简化后 点=%zu 面积=%.0f 矩形度=%.3f (阈值%.2f)\n",
                        s.size(), a, ratio, thr);
        if(ra>1e-9 && ratio>=thr){
            if((a>0)!=(polyArea2(rect)>0)) std::reverse(rect.begin(), rect.end());
            s = rect;
        }
        if(a0>0) outers.push_back(std::move(s));
        else     holes.push_back(std::move(s));
    }
    if(outers.empty()) return false;

    // 6.5) 窄条面矩形补全
    // 场景：折型板的折边带只有 ~43mm 宽、700mm 长，激光对它是掠射，点云沿长向
    // 断续、短向只在棱线附近密。形态学接不上断口（闭运算半径已开到上限），
    // 边界追踪出来是一条 4~8mm 的细带 + 若干孤立点，重建覆盖只有外接矩形的 40%。
    // 但工件本身是连续的（只是扫描缺陷），工程上必须按条带补齐，否则那段没有面、
    // 也就没有轨迹。判据三条同时满足才补：
    //   a) 短向跨度 <= 4×行距（确实是"窄条"，不是大平面；大平面的大缺口必须保留）
    //   b) 重建覆盖 < 65%（确实缺得厉害）
    //   c) 长向跨度 > 2×短向（是条带而非小碎块）
    {
        // 条带范围必须用"回筛后保留下来的点"来定，否则个别边缘离群点会把
        // 矩形拉长到没有数据的地方（实测折边带因此多出 237mm 的空矩形，
        // 面比轨迹长一截）。
        double u0=umin, u1=umax, v0=vmin, v1=vmax;
        if(keepMask && keepMask->size()==uv.size()){
            double a0=1e18, a1=-1e18, b0=1e18, b1=-1e18; int cnt=0;
            for(size_t t=0;t<uv.size();++t){
                if(!(*keepMask)[t]) continue;
                a0=std::min(a0,uv[t].x); a1=std::max(a1,uv[t].x);
                b0=std::min(b0,uv[t].y); b1=std::max(b1,uv[t].y);
                ++cnt;
            }
            if(cnt>=20){ u0=a0; u1=a1; v0=b0; v1=b1; }
        }
        double wU=u1-u0, wV=v1-v0;
        double wShort=std::min(wU,wV), wLong=std::max(wU,wV);
        double recon=0;
        for(auto& o : outers) recon += std::fabs(polyArea2(o));
        for(auto& h : holes)  recon -= std::fabs(polyArea2(h));
        double box = wU*wV;
        if(wShort > 4.0 && wShort <= 4.0*rowGap && box > 1.0 && wLong > 2.0*wShort
           && recon/box < 0.65){
            std::vector<Vec2> rect = {
                {u0,v0},{u1,v0},{u1,v1},{u0,v1}
            };
            outers.clear(); holes.clear();
            outers.push_back(rect);
            if(fillRectOut) *fillRectOut = rect;      // 供轨迹阶段按同一矩形拉通行
            if(dbg) fprintf(stderr, "[mesh] 窄条补全: 短向=%.1f 长向=%.1f 覆盖%.0f%% -> 按矩形补齐\n",
                            wShort, wLong, recon/box*100.0);
        }
    }

    // 7) 三角化（外环 + 内含孔洞桥接）
    std::sort(outers.begin(), outers.end(), [](const std::vector<Vec2>&a, const std::vector<Vec2>&b){
        return std::fabs(polyArea2(a)) > std::fabs(polyArea2(b)); });
    std::vector<char> holeUsed(holes.size(),0);
    for(auto& outer : outers){
        std::vector<int> hs;
        for(size_t h=0;h<holes.size();++h){
            if(holeUsed[h]) continue;
            if(pointInPoly(holes[h][0], outer)){ hs.push_back((int)h); holeUsed[h]=1; }
        }
        std::sort(hs.begin(), hs.end(), [&](int a,int b){
            return std::fabs(polyArea2(holes[a])) > std::fabs(polyArea2(holes[b])); });
        std::vector<Vec2> poly = outer;
        for(int h : hs) bridgeHole(poly, holes[h]);            // 失败则该孔不挖（实心显示）

        std::vector<int> tris;
        bool ok = earClip(poly, tris);
        if(!ok){                                                // 兜底：凸包扇形
            std::vector<Vec2> hull = convexHull(poly);
            if(hull.size()<3) continue;
            poly = hull; tris.clear();
            for(size_t k=1;k+1<hull.size();++k){
                tris.push_back(0); tris.push_back((int)k); tris.push_back((int)k+1);
            }
        }
        int base=(int)outVerts.size();
        if(base + (int)poly.size() > maxVerts) break;
        if((int)outTris.size()/3 + (int)tris.size()/3 > maxTris) break;
        for(auto& q : poly) outVerts.push_back(q);
        for(size_t k=0;k+2<(size_t)tris.size();k+=3){
            outTris.push_back(base+tris[k]);
            outTris.push_back(base+tris[k+1]);
            outTris.push_back(base+tris[k+2]);
        }
    }
    return !outTris.empty();
}

// 逐点法向（KD 树 + KNN PCA + 2 轮保边平滑）。
// 原来是 split_planes 内部的一个局部块。抽出来是因为"曲面块判定"也要用同一份
// 法向——否则要么重复算一遍（60 万点约 0.5s，白扔），要么两处口径不一致。
// K 取 25：线扫点云沿线密、线间疏，K 过小(10)时邻域落在同一条扫描线上近乎
// 共线，PCA 法向退化（区域生长碎片化的根源）。K=25 保证邻域跨越多条扫描线。
struct CloudNormals {
    static const int K = 25;
    KdTree tree;
    std::vector<Vec3> nor;
    std::vector<int>  knn;      // N*K 扁平邻接表（-1 表示不足）

    void build(const std::vector<Vec3>& P){
        const int n = (int)P.size();
        { ProfStage ps("2a.KD树构建"); tree.init(P); }
        nor.assign(n, Vec3{0,0,1});
        knn.assign((size_t)n*(size_t)K, -1);
        {
            ProfStage ps("2b.KNN法向");
            parallel_for(n, [&](int i){
                int* row = &knn[(size_t)i*(size_t)K];
                tree.knnIdxInto(P[i], K, row);
                Vec3 nf, cf;
                if(ls_plane_p(P, row, K, nf, cf) && len(nf)>1e-9) nor[i]=nf;
                else nor[i]=Vec3{0,0,1};
            });
        }
        // 法向图平滑 2 轮：nor[i] = 归一化(自身 + 邻居法向均值)，消除逐点 PCA 噪声。
        // 关键：必须"保边"——只平均与自身法向接近(<30°)的邻居。
        // 否则薄壁/陡壁点(其 KNN 邻域绝大多数是地面点)的法向会被拉平，
        // 侧壁在区域生长里根本分不出来（棕色槽型盘的近垂直侧壁就是这样丢的）。
        const double cosKeep = std::cos(30.0*3.141592653589793/180.0);
        for(int pass=0; pass<2; ++pass){
            std::vector<Vec3> sm(n);
            parallel_for(n, [&](int i){
                Vec3 s = nor[i];
                const int* row = &knn[(size_t)i*(size_t)K];
                for(int t=0;t<K;++t){
                    int j = row[t]; if(j<0) continue;
                    if(std::fabs(dot(nor[j], nor[i])) < cosKeep) continue;   // 跨棱不平均
                    s = s + nor[j];
                }
                sm[i] = norm(s);
                if(len(sm[i])<1e-9) sm[i]=nor[i];
            });
            nor.swap(sm);
        }
    }
};

// 多平面分割：法向区域生长（技术方案 2.4 节设计）。
// 思路：
//   1) 主平面 RANSAC 确定高度方向；全点到主平面高度极差 < curvatureTol
//      → 单一可喷面（曲率浅当平面），直接返回；
//   2) 否则对全部点估计法向（KNN PCA），按"法向与块种子法向夹角 <
//      normalGrowDeg"做 BFS 区域生长，把点云分成若干"法向一致"的面块；
//   3) 每块最小二乘拟合平面；近共面小块（夹角<3° 且偏移<flatnessTol）
//      向大块合并；碎块（< minPts）丢弃。
// 为什么不用"平面剥离/高度分层"：
//   - 顺序 RANSAC 会把环形凸边切成碎面；高度分层对连续斜坡折边无峰可分；
//   - 剥离出的"主平面外点"里混着主面自身的过渡带碎带，二次拟合拟出
//     散布多处的伪面。区域生长按局部法向归属，折型板两折面各成一整块。
std::vector<PlaneSeg> split_planes(const std::vector<Vec3>& P, const CloudNormals& cn,
                                   double tol, double curvatureTol, int minPts){
    std::vector<PlaneSeg> out;

    // ---- 第一步：RANSAC 拟合主平面，确定高度方向 n ----
    Plane mainPl; std::vector<int> mainInl;
    bool ransacOk;
    { ProfStage ps("2c.RANSAC主平面"); ransacOk = fit_plane(P, std::min(tol,3.0), mainPl, mainInl); }
    if(!ransacOk){
        // 拟合失败：退化为整体一个面
        std::vector<int> all(P.size()); std::iota(all.begin(),all.end(),0);
        Vec3 nf, cf;
        if(ls_plane(P, all, nf, cf)){
            if(nf.z<0){ nf=nf*(-1.0); }
            PlaneSeg seg; seg.pl={nf, -dot(nf,cf), cf}; seg.inl=all;
            out.push_back(std::move(seg));
        }
        return out;
    }
    if(mainPl.n.z < 0){ mainPl.n = mainPl.n*(-1.0); mainPl.d = -mainPl.d; }
    Vec3 n = mainPl.n;

    // ---- 高度极差判定：极差 < curvatureTol → 单一可喷面 ----
    {
        double hmin=1e18, hmax=-1e18;
        for(const auto& p:P){ double h=dot(p,n)+mainPl.d; hmin=std::min(hmin,h); hmax=std::max(hmax,h); }
        double span = hmax-hmin;
        if(span <= curvatureTol){
            PlaneSeg seg; seg.pl = mainPl; seg.inl = mainInl;
            fprintf(stderr, "[split] 台阶高度极差 %.2f <= curvatureTol %.2f，单一可喷面\n", span, curvatureTol);
            out.push_back(std::move(seg));
            return out;
        }
        fprintf(stderr, "[split] 台阶高度范围[%.1f,%.1f] span=%.1f > curvatureTol，启用法向区域生长分割\n", hmin, hmax, span);
    }

    // ---- 法向：由调用方统一算好传进来（CloudNormals） ----
    // 以前在这里现算。抽出去是因为"曲面块判定"要用同一份法向。
    const std::vector<Vec3>& nor = cn.nor;
    const std::vector<int>&  knnFlat = cn.knn;
    const int K = CloudNormals::K;

    // ---- 区域生长：种子法向固定，夹角阈值 normalGrowDeg ----
    const double cosGrow = std::cos(5.0*3.141592653589793/180.0);  // 与块种子法向夹角 < 5°
    ProfStage psGrow("2d.区域生长");
    std::vector<int> regionId(P.size(), -1);
    std::vector<Vec3> seedN;
    std::vector<std::vector<int>> regions;
    std::vector<int> stack;
    for(size_t s=0; s<P.size(); ++s){
        if(regionId[s] >= 0) continue;
        int id = (int)regions.size();
        Vec3 sn = nor[s];
        regions.push_back({});
        std::vector<int>& reg = regions.back();
        regionId[s] = id; reg.push_back((int)s);
        stack.clear(); stack.push_back((int)s);
        while(!stack.empty()){
            int i = stack.back(); stack.pop_back();
            const int* row = &knnFlat[(size_t)i*(size_t)K];
            for(int t=0;t<K;++t){
                int j = row[t]; if(j<0) continue;
                if(regionId[j] >= 0) continue;
                if(std::fabs(dot(nor[j], sn)) >= cosGrow){
                    regionId[j] = id; reg.push_back(j);
                    stack.push_back(j);
                }
            }
        }
        seedN.push_back(sn);
    }
    fprintf(stderr, "[split] 区域生长: 原始块数=%zu\n", regions.size());
    if(const char* d=std::getenv("SPRAYTRAJ_DUMP")){
        if(d[0]=='1'){
            std::vector<std::pair<int,double>> sz;   // (点数, 与主面夹角)
            for(auto& reg : regions){
                Vec3 nf, cf;
                if(!ls_plane(P, reg, nf, cf)) continue;
                double a = std::acos(std::min(1.0,std::fabs(dot(nf,n))))*180/3.141592653589793;
                sz.push_back({(int)reg.size(), a});
            }
            std::stable_sort(sz.begin(), sz.end(), [](auto&a,auto&b){return a.first>b.first;});
            fprintf(stderr, "[dbg] 块大小/夹角 Top20:");
            for(size_t i=0;i<sz.size() && i<20;++i)
                fprintf(stderr, " %d(%.0f°)", sz[i].first, sz[i].second);
            fprintf(stderr, "\n[dbg] 夹角>40° 的点数合计=%d, 其中最大块=%d\n",
                [&]{ int s=0; for(auto&e:sz) if(e.second>40) s+=e.first; return s; }(),
                [&]{ int m=0; for(auto&e:sz) if(e.second>40) m=std::max(m,e.first); return m; }());
        }
    }

    // ---- 每块拟合平面（先不过滤！碎块必须参与合并）----
    // 为什么不能先按 minPts 过滤：区域生长会把同一片缓翘曲面切成几十块，
    // 其中大多数小于 minPts。若先过滤，这些点就再也回不到主面上——
    // 重建面缺角、轨迹只覆盖半张板（棕色槽型盘地面只拿到 45% 的点）。
    // 正确顺序是"先按平面一致性合并成簇，再按簇的总点数过滤"。
    struct RawReg { Plane pl; std::vector<int> inl; double ang; };   // ang: 与主法向夹角(°)
    const int minSmall = std::max(10, (int)(P.size()*0.0015));
    std::vector<RawReg> raw;
    for(auto& reg : regions){
        if((int)reg.size() < minSmall) continue;
        Vec3 nf, cf;
        if(!ls_plane(P, reg, nf, cf)) continue;
        Plane pl; pl.n=nf; pl.c=cf; pl.d=-dot(nf,cf);
        double a = std::acos(std::min(1.0,std::fabs(dot(nf,n))))*180.0/3.141592653589793;
        raw.push_back({pl, std::move(reg), a});
    }
    // 按点数降序：大块先当种子，小块往大块并
    std::stable_sort(raw.begin(), raw.end(), [](const RawReg&a, const RawReg&b){
        return a.inl.size() > b.inl.size(); });

    // ---- 陡壁（侧壁）的定义 ----
    // 实测棕色槽型盘外沿有一圈 ~35mm 高、面接近垂直的裙边，与地面法向夹角 70~90°。
    // 这类面在区域生长里会被切成 20~180 点的碎弧段（法向沿圆弧连续变化，
    // 固定种子法向的 BFS 必然分段），不能靠在块层面合并来恢复，
    // 而是单独用 RANSAC 在"陡壁点集合"上逐面剥离（见下面的侧壁剥离）。
    const double wallDeg = 40.0;
    // 陡壁块不在这里成面（它们会被区域生长切成碎弧段），只把点留给后面的
    // "侧壁 RANSAC 剥离"统一处理，避免同一面被两条路径重复建出来。
    std::vector<RawReg> flatRaw;
    for(auto& r : raw)
        if(r.ang <= wallDeg) flatRaw.push_back(std::move(r));

    // 合并簇：法向差 <= angDeg 且 块内点到簇平面的最大残差 <= residCap 才并。
    // 单趟、按点数降序（大块先成簇），不做反复重扫以控制 O(n^2) 成本。
    auto mergeClusters = [&](std::vector<RawReg>& blocks, double angDeg, double residCap){
        std::vector<RawReg> res2;
        if(blocks.empty()) return res2;
        const double cs = std::cos(angDeg*3.141592653589793/180.0);
        std::vector<char> used(blocks.size(), 0);
        for(size_t i=0;i<blocks.size();++i){
            if(used[i]) continue;
            used[i] = 1;
            RawReg cur = blocks[i];
            const Vec3 seedN = blocks[i].pl.n;   // 固定种子法向：门限不随合并漂移
            for(size_t j=i+1;j<blocks.size();++j){
                if(used[j]) continue;
                // 用"对种子法向"的夹角判定，而不是对合并后重新拟合的法向判定：
                // 否则 13° 的斜坡种子可以吸收 0.5° 的地面碎块，重新拟合后法向
                // 又"平均"回地面方向，簇的身份完全乱了。
                if(std::fabs(dot(blocks[j].pl.n, seedN)) < cs) continue;
                double maxRes = 0.0;
                for(int k : blocks[j].inl)
                    maxRes = std::max(maxRes, std::fabs(dot(P[k]-cur.pl.c, cur.pl.n)));
                if(const char* d=std::getenv("SPRAYTRAJ_DUMP")){
                    if(d[0]=='1' && maxRes > residCap)
                        fprintf(stderr, "[dbg] 合并拒绝: 块%zu(%zu点, 夹角%.1f°) 最大残差%.2f > %.2f\n",
                            j, blocks[j].inl.size(),
                            std::acos(std::min(1.0,std::fabs(dot(blocks[j].pl.n,cur.pl.n))))*180.0/3.141592653589793,
                            maxRes, residCap);
                }
                if(maxRes > residCap) continue;
                cur.inl.insert(cur.inl.end(), blocks[j].inl.begin(), blocks[j].inl.end());
                Vec3 nf, cf;
                if(ls_plane(P, cur.inl, nf, cf)){
                    if(dot(nf, cur.pl.n) < 0) nf = nf*(-1.0);
                    cur.pl = {nf, -dot(nf,cf), cf};
                }
                used[j] = 1;
            }
            res2.push_back(std::move(cur));
        }
        return res2;
    };

    // 平缓面：角度硬隔离 8°（实测盘底碎块 0.3~1.0°、折型板折角 13.2°，
    // 8° 把"同面翘曲碎块"和"真实折边"分开），距离由 flatnessTol 决定。
    // 平缓面：两级合并。
    // 第一级在"块"上做，用固定种子法向——保证簇的身份是真一致的（不会出现
    // 13° 的斜坡种子吸收 0.5° 地面碎块后重新拟合又"平均"回地面的怪簇）。
    // 缺点是这样只能合并到种子 ±8°，曲面上会切出一堆碎片；所以第二级在"簇"上
    // 再合并一次：此时每个簇都有了明确的法向，跨簇合并不会再出现法向漂移问题。
    std::vector<RawReg> flatMerged = mergeClusters(flatRaw, 8.0, std::max(tol, 2.0));
    flatMerged = mergeClusters(flatMerged, 8.0, std::max(tol, 2.0));
    // 侧壁的最小点数与平面度门槛：太小的弧段（圆角过渡带）不算独立可喷面，
    // 否则一堆碎面各自生成轨迹，视觉上又是"乱七八糟"。
    const int    minWallPts = std::max(60, (int)(P.size()*0.008));
    const double wallRmsCap = std::max(tol*0.5, 2.5);

    for(auto& m : flatMerged){
        if((int)m.inl.size() < minPts) continue;
        // 平面性检验：一堆方向各异的斜壁段拼在一起时，最小二乘平面会"平均"回
        // 地面法向，看着像个平面其实不是（棕色槽型盘盘沿的斜坡环就这样）。
        double rss=0; for(int k : m.inl){ double e=dot(P[k],m.pl.n)+m.pl.d; rss+=e*e; }
        double rms = std::sqrt(rss/(double)std::max<size_t>(1,m.inl.size()));
        if(const char* d=std::getenv("SPRAYTRAJ_DUMP"))
            if(d[0]=='1')
                fprintf(stderr, "[dbg] 平面簇 %zu点 角度%.1f° RMS=%.2f (上限%.2f)\n",
                    m.inl.size(), m.ang, rms, std::max(tol, 3.0));
        if(rms > std::max(tol, 3.0)){
            fprintf(stderr, "[split] 丢弃非平面簇: 点数=%zu RMS=%.2f (>%.2f) 夹角%.1f°\n",
                m.inl.size(), rms, std::max(tol,3.0), m.ang);
            continue;
        }
        PlaneSeg seg; seg.pl = m.pl; seg.inl = std::move(m.inl);
        fprintf(stderr, "[split] 面块: 点数=%zu 法向=(%.3f,%.3f,%.3f) 种子夹角=%.1f° RMS=%.2f\n",
            seg.inl.size(), seg.pl.n.x, seg.pl.n.y, seg.pl.n.z, m.ang, rms);
        out.push_back(std::move(seg));
    }

    // ---- 侧壁（近垂直面）RANSAC 逐面剥离 ----
    // 为什么不用"区域生长块合并"：侧壁在俯视扫描里是一条很窄的带子，
    // 点稀疏、法向沿圆弧连续变化，区域生长会把它切成几十个 20~180 点的弧段，
    // 再按块合并也拼不回"一个面"（棕色槽型盘 4 条侧壁只能捞到 1 条）。
    // 改用直接在"陡壁点集合"上做 RANSAC：每次剥出最大的一个共面集 = 一个侧壁面，
    // 直到剩下的点凑不出一个面为止。这样方盘的 4 条侧壁各自成面。
    {
        std::vector<int> spIdx;
        spIdx.reserve(P.size()/4);
        for(size_t i=0;i<P.size();++i){
            double a = std::acos(std::min(1.0,std::fabs(dot(nor[i],n))))*180.0/3.141592653589793;
            if(a > wallDeg) spIdx.push_back((int)i);
        }
        std::vector<char> gone(spIdx.size(), 0);
        // 自适应容差：先按 1mm 找"真平"的侧壁；找不到再逐级放宽到 1.6/2.2/3.0mm。
        // 俯视扫描对侧向面是掠射，点又稀又抖，固定紧容差只能捞到正对相机的那一面。
        const double wallTols[4] = {1.0, 1.6, 2.2, 3.0};
        for(int it=0; it<8; ++it){
            std::vector<Vec3> cur;
            std::vector<int> curSp;                        // cur[k] 对应 spIdx[curSp[k]]
            for(size_t k=0;k<spIdx.size();++k)
                if(!gone[k]){ cur.push_back(P[spIdx[k]]); curSp.push_back((int)k); }
            if((int)cur.size() < minWallPts) break;

            bool accepted = false;
            for(int ti=0; ti<4 && !accepted; ++ti){
                Plane pl; std::vector<int> inl;
                if(!fit_plane(cur, wallTols[ti], pl, inl)) continue;
                if((int)inl.size() < minWallPts) continue;
                // 必须真的是陡壁，而不是混进陡壁集合里的地面点
                double angF = std::acos(std::min(1.0,std::fabs(dot(pl.n,n))))*180.0/3.141592653589793;
                if(angF <= wallDeg) continue;
                double rss=0; for(int k : inl){ double e=dot(cur[k],pl.n)+pl.d; rss+=e*e; }
                double rms = std::sqrt(rss/(double)inl.size());
                if(rms > wallRmsCap) continue;
                if(const char* d=std::getenv("SPRAYTRAJ_DUMP")) if(d[0]=='1')
                    fprintf(stderr, "[dbg] 侧壁候选(tol=%.1f): %zu点 夹角%.1f° RMS=%.2f 剩余%zu\n",
                        wallTols[ti], inl.size(), angF, rms, cur.size());

                // 与已接受的侧壁重复（同向且同位置）：把点并进去，别丢
                int dupAt = -1;
                for(size_t si=0; si<out.size(); ++si){
                    const auto& s = out[si];
                    if(std::fabs(dot(s.pl.n, pl.n)) < std::cos(8.0*3.141592653589793/180.0)) continue;
                    if(std::fabs(s.pl.d - pl.d) < std::max(tol*0.5, 2.0)) dupAt = (int)si;
                }
                PlaneSeg seg; seg.pl = pl;
                for(int k : inl){ seg.inl.push_back(spIdx[curSp[k]]); gone[curSp[k]] = 1; }
                accepted = true;                            // 这些点已消费
                if(dupAt >= 0){
                    auto& dst = out[dupAt].inl;
                    dst.insert(dst.end(), seg.inl.begin(), seg.inl.end());
                    if(const char* d=std::getenv("SPRAYTRAJ_DUMP")) if(d[0]=='1')
                        fprintf(stderr, "[dbg] 侧壁并入已有面%d: +%zu 点\n", dupAt, seg.inl.size());
                    break;
                }
                fprintf(stderr, "[split] 侧壁面: 点数=%zu 法向=(%.3f,%.3f,%.3f) RMS=%.2f 与地面夹角=%.1f°\n",
                    seg.inl.size(), seg.pl.n.x, seg.pl.n.y, seg.pl.n.z, rms, angF);
                out.push_back(std::move(seg));
            }
            if(!accepted) break;
        }
    }

    // 主面（点数最多）排到最前，便于报告与显示
    std::stable_sort(out.begin(), out.end(), [](const PlaneSeg&a, const PlaneSeg&b){
        return a.inl.size() > b.inl.size(); });

    // ---- 塌缩判定：缓翘曲/浅曲面 → 单一可喷面 ----
    // 区域生长在"法向连续渐变"的缓翘曲板上（现场大板常见波浪变形）会切出
    // 多块近共面小面，轨迹呈拼块状。按"曲率浅一律当平面"原则，若满足：
    //   1) 所有面块法向与主平面法向夹角 < 8°（无真实折边/侧壁）；
    //   2) 整体对主平面的 RMS <= max(2×flatnessTol, 6mm)（起伏平缓）；
    //   3) 高度直方图无第二峰（排除槽型盘"同法向不同高度"的台阶结构）；
    // 则塌缩为单一可喷面（主平面），全部点作为内点。
    // 注意：即使只分割出 1 个面也要做这一步——缓翘曲的大板在做区域生长时
    // 会被切成上千个法向略有差异的小块，绝大多数块小于 minPts 被丢弃，
    // 只剩下一块"主区"，重建面就只覆盖板子的一部分（缺角）。塌缩判定成立时
    // 把整片点云都归给这一个面，重建轮廓才是完整板面。
    if(!out.empty()){
        const double cosNear = std::cos(10.0*3.141592653589793/180.0);
        // 近共面判定（允许少数派）：偏离主面 >10° 的面块点数占比 <=8% 时
        // 视为翘曲边缘的少数派碎块，不影响"整体近共面"判定。
        // 注意：不能放宽太多——槽型盘的侧壁（法向偏离 10~18°，面本身接近垂直）
        // 占比约 10%，必须保留为独立面；折型板折面占 46%，也必须保留。
        double outlierPts = 0, totalPts = 0;
        for(auto& s : out){
            totalPts += (double)s.inl.size();
            if(std::fabs(dot(s.pl.n, n)) < cosNear) outlierPts += (double)s.inl.size();
        }
        bool allNear = (outlierPts <= totalPts*0.08);
        // 已识别出侧壁面时绝不塌缩——塌缩会把侧壁整片倒进地面，
        // 侧壁（近垂直）恰恰是需要单独喷的面。
        bool hasWall = false;
        for(auto& s : out)
            if(std::acos(std::min(1.0,std::fabs(dot(s.pl.n,n))))*180.0/3.141592653589793 > wallDeg)
                hasWall = true;
        if(hasWall) allNear = false;
        double gRss=0; for(const auto& p:P){ double e=dot(p,n)+mainPl.d; gRss+=e*e; }
        double gRms = std::sqrt(gRss/P.size());
        // 高度直方图峰值检查：>=2 个相距明显的峰说明存在"同法向不同高度"
        // 的台阶结构（如槽型盘底盘+凸边环顶），绝不能塌缩。
        double hmin=1e18,hmax=-1e18;
        for(const auto& p:P){ double h=dot(p,n)+mainPl.d; hmin=std::min(hmin,h); hmax=std::max(hmax,h); }
        double binW2=2.0; int nb=(int)std::ceil((hmax-hmin)/binW2)+1;
        std::vector<double> hist2(nb,0.0);
        for(const auto& p:P){ double h=dot(p,n)+mainPl.d; int b=(int)((h-hmin)/binW2); if(b>=nb)b=nb-1; if(b<0)b=0; hist2[b]+=1.0; }
        std::vector<int> peaks2;
        for(int b=0;b<nb;++b){
            bool isMax=(b==0||hist2[b]>=hist2[b-1])&&(b==nb-1||hist2[b]>=hist2[b+1]);
            if(isMax && hist2[b] >= P.size()*0.02) peaks2.push_back(b);
        }
        // 合并间距 < tol 的相邻峰
        {
            std::vector<int> mg;
            for(int b:peaks2){
                if(mg.empty() || (b-mg.back())*binW2 >= tol) mg.push_back(b);
                else if(hist2[b]>hist2[mg.back()]) mg.back()=b;
            }
            peaks2=mg;
        }
        double rmsCap = std::max(2.0*tol, 6.0);
        fprintf(stderr, "[split] 塌缩检查: 主面法向=(%.3f,%.3f,%.3f) 面块=%zu 全部近共面(10°)=%d 全局RMS=%.2f(上限%.1f) 高度峰数=%zu\n",
            n.x, n.y, n.z, out.size(), (int)allNear, gRms, rmsCap, peaks2.size());
        for(auto& s : out)
            fprintf(stderr, "[split]   面块法向=(%.3f,%.3f,%.3f) 与主面夹角=%.1f°\n",
                s.pl.n.x, s.pl.n.y, s.pl.n.z,
                std::acos(std::min(1.0,std::fabs(dot(s.pl.n,n))))*180/3.141592653589793);
        if(allNear && gRms <= rmsCap && peaks2.size() <= 1){
            fprintf(stderr, "[split] 缓翘曲/浅曲面 → 塌缩为单一可喷面\n");
            out.clear();
            PlaneSeg seg; seg.pl = mainPl;
            seg.inl.resize(P.size()); std::iota(seg.inl.begin(), seg.inl.end(), 0);
            out.push_back(std::move(seg));
        }
    }

    // ---- 同向同高度面并入主面 ----
    // 盘沿/圆角的过渡带会被区域生长切成若干"弧段簇"，各自拟合出来方向接近地面、
    // 平面度也不差（RMS≈1mm），但它们是同一个面的一部分，不是独立的可喷面。
    // 判据：与主面夹角 <30° 且 全部点落在主面容差内（<= max(tol,2mm)）→ 并回主面。
    // 真实台阶（槽型盘凸边、折型板折面）要么有高度差、要么夹角大，不会被并。
    if(out.size() > 1){
        size_t mi = 0;
        for(size_t i=1;i<out.size();++i) if(out[i].inl.size() > out[mi].inl.size()) mi = i;
        const Plane& mp = out[mi].pl;
        const double cosSame = std::cos(30.0*3.141592653589793/180.0);
        const double sameD = std::max(tol, 2.0);
        std::vector<int> absorbed;
        for(size_t i=0;i<out.size();++i){
            if(i == mi) continue;
            if(std::fabs(dot(out[i].pl.n, mp.n)) < cosSame) continue;
            double maxD = 0.0;
            for(int k : out[i].inl) maxD = std::max(maxD, std::fabs(dot(P[k],mp.n)+mp.d));
            if(maxD > sameD) continue;
            fprintf(stderr, "[split] 同向同高度并入主面: %zu点 夹角%.1f° 最大高差%.2fmm\n",
                out[i].inl.size(),
                std::acos(std::min(1.0,std::fabs(dot(out[i].pl.n,mp.n))))*180.0/3.141592653589793, maxD);
            absorbed.push_back((int)i);
        }
        if(!absorbed.empty()){
            std::vector<PlaneSeg> keep;
            std::vector<int> addPts;
            for(size_t i=0;i<out.size();++i){
                if(i == mi) continue;              // 主面最后单独移动，避免被 move 掉
                bool ab = false;
                for(int a : absorbed) if((size_t)a == i) ab = true;
                if(ab) addPts.insert(addPts.end(), out[i].inl.begin(), out[i].inl.end());
                else   keep.push_back(std::move(out[i]));
            }
            out[mi].inl.insert(out[mi].inl.end(), addPts.begin(), addPts.end());
            keep.push_back(std::move(out[mi]));
            out.swap(keep);
        }
    }

    if(out.empty()){
        PlaneSeg seg; seg.pl = mainPl; seg.inl = mainInl;
        out.push_back(std::move(seg));
    }

    // ---- 残点归并：没分到任何面的点，并到"距离最近且在容差内"的面上 ----
    // 为什么要做：区域生长切出的碎弧段（圆角过渡、斜坡环）够不上平面门槛会被
    // 丢掉，这些点如果不归并，重建面的轮廓就停在它们内侧，轨迹喷不到边缘。
    // 归并只看点到面的距离（<= max(flatnessTol, curvatureTol)），
    // 让主面轮廓扩到真实边缘，同时不至于把远处的结构拉进来。
    {
        std::vector<char> owned(P.size(), 0);
        for(auto& s : out) for(int k : s.inl) owned[k] = 1;
        const double mergeD = std::max(tol, curvatureTol);
        std::vector<std::vector<int>> extra(out.size());
        int moved = 0;
        for(size_t i=0;i<P.size();++i){
            if(owned[i]) continue;
            int best = -1; double bd = 1e18;
            for(size_t j=0;j<out.size();++j){
                double d = std::fabs(dot(P[i], out[j].pl.n)+out[j].pl.d);
                if(d < bd){ bd = d; best = (int)j; }
            }
            if(best >= 0 && bd <= mergeD){ extra[best].push_back((int)i); ++moved; }
        }
        for(size_t j=0;j<out.size();++j){
            if(extra[j].empty()) continue;
            out[j].inl.insert(out[j].inl.end(), extra[j].begin(), extra[j].end());
        }
        if(moved) fprintf(stderr, "[split] 残点归并: %d 个点并入最近的平面（容差%.1fmm）\n", moved, mergeD);
    }

    // ---- 连续曲面合并：相邻 + 同向 → 必须是同一个面 ----
    // 现象（现场长板/翘曲板）：整块板翘曲 20~50mm，flatnessTol=3 时区域生长会把它
    // 切成 5~6 块"法向只差 2~6°、空间上首尾相连"的平行片。它们各有自己的行方向，
    // 一片一套行，拼起来就是"面很多 + 轨迹乱"。从喷涂角度看，这些片是同一张脸上
    // 的同一个朝向，拆开毫无意义。
    // 判据（用几何连续性，不用残差阈值——残差无法区分"缓翘"和"台阶"）：
    //   ① 与"本组种子面"法向夹角 <= 6°（固定种子，禁止链式漂移，折边不会被吞）；
    //   ② 两面在点云里直接相连（存在一条 KNN 边横跨两个面）；
    //      真实台阶中间隔着竖壁点，两侧水平面并不直接相连；
    //   ③ 该面点到种子面的高度中位数 <= warpCap(=max(3*tol,20mm))，
    //      防止把上下两层的翻边/凸台并进来。
    // 副作用（正向）：棕色槽型盘在 flat=3 下盘底被切碎的问题同样被这条规则治好，
    // 因为碎块本来就相邻同向。
    if(out.size() > 1 && out.size() <= 24){
        const size_t F = out.size();
        std::vector<int> lab(P.size(), -1);
        for(size_t i=0;i<F;++i) for(int k : out[i].inl) lab[k] = (int)i;
        // 邻接矩阵：KNN 图里横跨两个面的边 → 这两个面空间相连
        std::vector<char> adj(F*F, 0);
        for(size_t i=0;i<P.size();++i){
            int a = lab[i]; if(a < 0) continue;
            const int* row = &knnFlat[i*(size_t)K];
            for(int t=0;t<K;++t){
                int j = row[t]; if(j<0) continue;
                int b = lab[j];
                if(b < 0 || b == a) continue;
                adj[(size_t)a*F+b] = adj[(size_t)b*F+a] = 1;
            }
        }
        const double cosCon   = std::cos(8.0*3.141592653589793/180.0);
        const double warpCap  = std::max(3.0*tol, 20.0);
        std::vector<int> grp(F, -1);
        std::vector<std::vector<int>> groupMembers;
        for(size_t i=0;i<F;++i){
            if(grp[i] >= 0) continue;
            int gid = (int)groupMembers.size();
            grp[i] = gid;
            groupMembers.push_back({(int)i});
            const Vec3 seedN = out[i].pl.n;        // 固定种子：比较永远对种子做
            const Vec3 seedC = out[i].pl.c;
            bool grew = true;
            while(grew){
                grew = false;
                for(size_t j=0;j<F;++j){
                    if(grp[j] >= 0) continue;
                    if(std::fabs(dot(out[j].pl.n, seedN)) < cosCon) continue;
                    bool touch = false;
                    for(int m : groupMembers[gid])
                        if(adj[(size_t)j*F+(size_t)m]){ touch = true; break; }
                    if(!touch){
                        if(const char* d=std::getenv("SPRAYTRAJ_DUMP")) if(d[0]=='1')
                            fprintf(stderr, "[dbg] 连续合并跳过: 面%zu 与组%d 不相邻 (夹角%.1f°)\n",
                                j, gid, std::acos(std::min(1.0,std::fabs(dot(out[j].pl.n,seedN))))*180.0/3.141592653589793);
                        continue;
                    }
                    std::vector<double> hh; hh.reserve(out[j].inl.size());
                    for(int k : out[j].inl) hh.push_back(std::fabs(dot(P[k]-seedC, seedN)));
                    size_t mid = hh.size()/2;
                    std::nth_element(hh.begin(), hh.begin()+mid, hh.end());
                    if(hh[mid] > warpCap){
                        if(const char* d=std::getenv("SPRAYTRAJ_DUMP")) if(d[0]=='1')
                            fprintf(stderr, "[dbg] 连续合并跳过: 面%zu 高度中位%.1f > %.1f\n", j, hh[mid], warpCap);
                        continue;
                    }
                    grp[j] = gid;
                    groupMembers[gid].push_back((int)j);
                    grew = true;
                }
            }
        }
        if((int)groupMembers.size() < (int)F){
            std::vector<PlaneSeg> merged;
            for(auto& mem : groupMembers){
                PlaneSeg seg;
                seg.pl = out[mem[0]].pl;
                for(int m : mem){
                    seg.inl.insert(seg.inl.end(), out[m].inl.begin(), out[m].inl.end());
                    if(m != mem[0])
                        fprintf(stderr, "[split] 连续曲面合并: 面%zu(%zu点) 并入 面%d(%zu点) 夹角%.1f°\n",
                            (size_t)m, out[m].inl.size(), mem[0], out[mem[0]].inl.size(),
                            std::acos(std::min(1.0,std::fabs(dot(out[m].pl.n, out[mem[0]].pl.n))))*180.0/3.141592653589793);
                }
                std::sort(seg.inl.begin(), seg.inl.end());
                seg.inl.erase(std::unique(seg.inl.begin(), seg.inl.end()), seg.inl.end());
                Vec3 nf, cf;
                if(ls_plane(P, seg.inl, nf, cf)){
                    if(dot(nf, seg.pl.n) < 0) nf = nf*(-1.0);
                    seg.pl = {nf, -dot(nf,cf), cf};
                }
                merged.push_back(std::move(seg));
            }
            out.swap(merged);
        }
    }

    // ---- 终局同面合并：法向接近 + 距离接近 → 必须并成同一张脸 ----
    // 现象（现场点云 2026-08-21-11-33-13 的板边倒角）：同一条倒角带被拆成两条平行窄带
    // （36.0°/28mm 宽 与 47.9°/18mm 宽），法向只差 11.9°、平面偏移只差 1.1mm、
    // 点云直接相连。两块各有一套行方向、各布一套行，中间还多一次空移——
    // 肉眼看就是"同一个面上两套乱轨迹"。
    // 判据（必须同时满足，缺一条都会把真折边/台阶吞掉）：
    //   ① 同向：法向夹角 <= 15°（法向接近）；
    //   ② 距离接近：小块点到大块平面的高度 |h| 的 p95 <= 6mm（真折边过不了这条）；
    //   ③ 相邻：KNN 图里有边横跨两面（同一次扫描的连续数据，排除同平面孤岛）；
    //   ④ 并后仍是平面：重拟合 RMS <= max(flatnessTol, 8mm)。
    // 与上一条"连续曲面合并"的分工：那条是 8° 硬门限（防缓翘板链式漂移），
    // 这条用"距离极近"换角度余量，专治扫描条带把同一张窄面切成多条。
    if(out.size() > 1 && out.size() <= 24){
        const size_t F0 = out.size();
        std::vector<int> lab(P.size(), -1);
        for(size_t i=0;i<F0;++i) for(int k : out[i].inl) lab[k] = (int)i;
        std::vector<char> adj(F0*F0, 0);
        for(size_t i=0;i<P.size();++i){
            int a = lab[i]; if(a < 0) continue;
            const int* row = &knnFlat[i*(size_t)K];
            for(int t=0;t<K;++t){
                int j = row[t]; if(j<0) continue;
                int b = lab[j];
                if(b >= 0 && b != a) adj[(size_t)a*F0+b] = adj[(size_t)b*F0+a] = 1;
            }
        }
        double angG = 15.0, hCap = 6.0;   // 现场实测标定；环境变量仅供标定用
        if(const char* e = std::getenv("SPRAYTRAJ_MERGE_ANG")){ double v = atof(e); if(v > 0) angG = v; }
        if(const char* e = std::getenv("SPRAYTRAJ_MERGE_H"))  { double v = atof(e); if(v > 0) hCap = v; }
        const char* dumpEnv = std::getenv("SPRAYTRAJ_DUMP");
        const bool dbg = dumpEnv && dumpEnv[0] == '1';
        const double cosG   = std::cos(angG*3.141592653589793/180.0);
        const double rmsCap = std::max(tol, 8.0);
        // 法向定向（与出函数后的统一定向口径一致）：先 +Z 半球，近竖直按背离质心。
        // 必须用有向夹角：薄板双面扫描时上下两面法向相反，夹角=0 但绝不是同一张脸。
        Vec3 partC{0,0,0};
        for(const auto& p : P) partC = partC + p;
        partC = partC*(1.0/(double)std::max<size_t>(1,P.size()));
        auto orient = [&](Vec3 n, const Vec3& c){
            if(std::fabs(n.z) > 0.05){ if(n.z < 0) n = n*(-1.0); }
            else if(dot(n, c - partC) < 0) n = n*(-1.0);
            return n;
        };
        std::vector<std::vector<int>> members(F0);          // out[i] 由哪些原始面组成
        for(size_t i=0;i<F0;++i) members[i] = {(int)i};
        bool did = true;
        while(did && out.size() > 1){
            did = false;
            for(size_t a=0; a<out.size() && !did; ++a){
                for(size_t b=a+1; b<out.size() && !did; ++b){
                    double ang = std::acos(std::max(-1.0, std::min(1.0, dot(orient(out[a].pl.n, out[a].pl.c),
                                                                         orient(out[b].pl.n, out[b].pl.c)))))*180.0/3.141592653589793;
                    if(ang > angG) continue;
                    bool touch = false;
                    for(int x : members[a]){ for(int y : members[b]) if(adj[(size_t)x*F0+y]){ touch = true; break; } if(touch) break; }
                    if(dbg) fprintf(stderr, "[dbg] 同面候选: 面%zu-%zu 夹角%.1f° 相邻=%d\n", a, b, ang, (int)touch);
                    if(!touch) continue;
                    const PlaneSeg& big = out[a].inl.size() >= out[b].inl.size() ? out[a] : out[b];
                    const PlaneSeg& sml = out[a].inl.size() >= out[b].inl.size() ? out[b] : out[a];
                    if(sml.inl.empty()) continue;
                    std::vector<double> hh; hh.reserve(sml.inl.size());
                    for(int k : sml.inl) hh.push_back(std::fabs(dot(P[k]-big.pl.c, big.pl.n)));
                    size_t q = (size_t)(0.95*(double)(hh.size()-1));
                    std::nth_element(hh.begin(), hh.begin()+q, hh.end());
                    double h95 = hh[q];
                    if(dbg) fprintf(stderr, "[dbg] 同面判距: 面%zu-%zu 高度p95=%.2fmm\n", a, b, h95);
                    if(h95 > hCap) continue;
                    std::vector<int> all = out[a].inl;
                    all.insert(all.end(), out[b].inl.begin(), out[b].inl.end());
                    std::sort(all.begin(), all.end());
                    all.erase(std::unique(all.begin(), all.end()), all.end());
                    Vec3 nf, cf;
                    if(!ls_plane(P, all, nf, cf)) continue;
                    double rms = 0;
                    for(int k : all){ double d = dot(P[k]-cf, nf); rms += d*d; }
                    rms = std::sqrt(rms/(double)std::max<size_t>(1,all.size()));
                    if(rms > rmsCap){
                        if(dbg) fprintf(stderr, "[dbg] 同面合并放弃: 面%zu-%zu 重拟合RMS=%.2f > %.2f\n", a, b, rms, rmsCap);
                        continue;
                    }
                    fprintf(stderr, "[split] 同面合并(法向%.1f° 距离p95=%.2fmm): 面%zu(%zu点) + 面%zu(%zu点) -> %zu点 RMS=%.2f(门限%.1f)\n",
                            ang, h95, a, out[a].inl.size(), b, out[b].inl.size(), all.size(), rms, rmsCap);
                    if(dot(nf, out[a].pl.n) < 0) nf = nf*(-1.0);
                    out[a].pl = {nf, -dot(nf,cf), cf};
                    out[a].inl = std::move(all);
                    out.erase(out.begin()+(long)b);
                    members[a].insert(members[a].end(), members[b].begin(), members[b].end());
                    members.erase(members.begin()+(long)b);
                    did = true;
                }
            }
        }
    }

    // 最终按点数降序：保证 plane0 永远是主面（报告/显示/配色都按此约定）
    std::stable_sort(out.begin(), out.end(), [](const PlaneSeg&a, const PlaneSeg&b){
        return a.inl.size() > b.inl.size(); });

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

// ==================== 曲面块（自由曲面）行追踪 ====================
// 为什么不做"参数化展开"：
//   一开始想按柱面展开（轴 + 剖面脊线 -> (弧长, 轴向) 二维域），再复用平面
//   那套 2D 布行。实测不行——这块"曲型板"沿任意候选轴展开后，剖面残差
//   p90 都在 20mm 以上（最好的轴 res90=22.7mm，p99=47.6mm），说明它既不是
//   柱面也不是回转面，是双曲率自由曲面。任何"投影/展开成平面"的做法都会
//   在陡峭段产生空间上的挤压拉伸，行距失控。
//
// 所以直接做**曲面上的行追踪**（和工业 CAM 的等距偏置走刀同一个思路）：
//   1. 一圈行：从曲面一侧边缘起步，按点距 pstep 步进；
//   2. 每一步把"目标点"投影到局部拟合平面上（用该点邻域 30 个点做 PCA 平面），
//      投影位移过大 = 已经出了工件，本行结束；
//   3. 切向做平行传输（投影掉法向分量），所以行是贴着曲面走的；
//   4. 行尾沿 t×n 侧移行距 L 起下一行，方向取反 —— 天生就是蛇形，
//      行与行相邻、不需要额外的换行过渡；
//   5. 两遍定距：先由投影包络估出横向跨度 D，N=round(D/L)，实际行距 D/N，
//      保证首末两行都贴边（否则最后一条缝会差半行距）。
// 好处：行距是曲面上的真实弦长，法向是局部拟合法向，与曲率无关。

struct SurfPatch {
    std::vector<Vec3> P;        // 曲面块点（世界坐标）
    std::vector<Vec3> N;        // 对应法向（已定向朝外）
    KdTree tree;                // 建在 P 上
    Vec3 c0{0,0,0};             // 平均平面上的参考点
    Vec3 n0{0,0,1};             // 平均法向（朝外）
    Vec3 eL{1,0,0}, eS{0,1,0};  // 平均平面内正交基（长/短）
    double spacing = 2.0;       // 点云中位点间距（决定局部拟合半径）
    double rFit = 6.0;          // 局部拟合查询半径
    // 行追踪结果：所有行按蛇形顺序首尾相接
    std::vector<Vec3> rowP;
    std::vector<Vec3> rowN;
    std::vector<int>  rowBeg;   // 每行起始下标，末尾再加一个哨兵
    double rowGap = 0.0;        // 实际行距
    double cover = 0.0;         // 覆盖率 0..1
    double fitRms = 0.0;        // 点到局部拟合面的残差 RMS
    int    nRow = 0;
    // ---- 重建曲面模型（与显示网格共用同一个 (u,v) 参数域网格）----
    // 曲面块原先只有"点云 + 逐点法向"：显示网格另算一套格内中位高度，
    // 轨迹逐点投影、法向逐点 PCA —— 同一个工件的"曲面"存在两套互不相干
    // 的表达，而且都不平滑。这里把几何与姿态统一到同一个网格上：
    //   几何 mH  ：格内中位高度 -> 归一化高斯（σ≈1.5 格），曲率连续
    //   姿态 mNm ：格内逐点法向取中位（抗坏点）-> 3x3 归一化平滑
    // 只用于显示与姿态，不改行点位置，因此覆盖率/贴合残差/行数完全不变。
    int    mNu=0, mNv=0;
    double mCell=0, mU0=0, mV0=0;
    std::vector<double> mH;             // 平滑高度（沿 n0），1e300=无效
    std::vector<Vec3>   mNm;            // 格中位法向（已定向朝外）
    std::vector<unsigned char> mOcc;    // 1=该格有实测点（插值才可信）
    bool   modelOk=false;
};

// 把 q 投到曲面块上：返回曲面点 pos 与局部法向 nrm。
// 失败表示 q 已经不在工件上（邻域点太少 / 离点云太远 / 邻域不共面）。
//
// ⚠ 邻域必须"跨扫描线"：线扫点云沿线 0.3mm、线间 1.5mm，半径给小了邻域
// 退化成一条线，ls_plane_p 会返回 (0,0,1) 这种垃圾法向（它的 det≈0 分支），
// 于是投影沿 Z 跑偏几毫米、行追踪第一步就被判"出界"。
// 所以这里用固定 K=36 个近邻，并对拟合结果做残差校验，不共面就判失败。
static const char* g_surfWhy = "";
static double g_surfLocal = 1.0;   // 最近一次投影查询的局部采样尺度
static bool surf_project(const SurfPatch& sp, Vec3 q, Vec3& pos, Vec3& nrm)
{
    g_surfWhy = "ok";
    KnnBuf b; sp.tree.knn(q, 36, b);
    if(b.n < 6){ g_surfWhy="邻域不足"; return false; }
    // 局部采样尺度：取最近的 8 个邻居的平均距离。
    // 线扫点云有"横向密、纵向疏"的条纹，个别位置 4mm 内只有 1 个点，
    // 用全局中位点间距当门槛会把正常曲面点误判成"出界"（实测行追踪第 2 步就断）。
    double loc8 = 0; int n8 = 0;
    {
        double best[8];
        for(int i=0;i<8;++i) best[i] = 1e300;
        for(int i=0;i<b.n;++i){
            double d2v = b.d2[i];
            if(d2v >= best[7]) continue;
            int k = 7;
            while(k > 0 && best[k-1] > d2v){ best[k] = best[k-1]; --k; }
            best[k] = d2v;
        }
        for(int i=0;i<8;++i) if(best[i] < 1e299){ loc8 += std::sqrt(best[i]); ++n8; }
        if(n8 > 0) loc8 /= (double)n8;
    }
    if(n8 < 3){ g_surfWhy="近邻不足8"; return false; }
    // 局部尺度必须封顶：不封顶会形成正反馈（稀疏->门槛放松->走得更偏->更稀疏），
    // 实测第一行直接飘到工件外 1 米处，第二行起点变成 (‑49.8,‑1398.8,722.2)。
    const double rCap = std::max(2.0, 8.0*sp.spacing);
    if(loc8 > 2.5*rCap){ g_surfWhy="远离点云"; return false; }
    const double rLocal = std::min(std::max(0.5, loc8), rCap);
    g_surfLocal = rLocal;

    int idx[36], cnt = 0;
    for(int i=0;i<b.n;++i) if(b.id[i] >= 0) idx[cnt++] = b.id[i];
    if(cnt < 6) return false;
    Vec3 nf, cf;
    if(!ls_plane_p(sp.P, idx, cnt, nf, cf)) return false;
    nf = norm(nf); if(len(nf) < 0.5) return false;
    // 残差校验：邻域必须真的共面，否则这个法向不可信
    {
        double rss=0;
        for(int i=0;i<cnt;++i){ double e=dot(sp.P[idx[i]]-cf, nf); rss+=e*e; }
        double rms = std::sqrt(rss/(double)cnt);
        if(rms > std::max(0.8, 0.9*rLocal)){
            static char buf[96]; snprintf(buf,sizeof(buf),"邻域不共面 rms=%.2f > %.2f", rms, std::max(0.8,0.9*rLocal));
            g_surfWhy = buf; return false;
        }
    }
    // 平面法向对齐到最近点的法向（ls_plane 的符号是随机的）
    int nearest = idx[0]; double nd = 1e300;
    for(int i=0;i<cnt;++i){ double d=len(sp.P[idx[i]]-q); if(d<nd){ nd=d; nearest=idx[i]; } }
    if(dot(nf, sp.N[nearest]) < 0) nf = nf*(-1.0);
    // 投影：把 q 沿 nf 投到过 cf 的平面上
    pos = q - nf*dot(q-cf, nf);
    nrm = nf;
    // 投完必须仍贴着点云（否则是从工件边缘外掉的笛卡尔点）
    KnnBuf b2; sp.tree.knn(pos, 1, b2);
    if(b2.n < 1) return false;
    double dmin = std::sqrt(b2.d2[0]);
    const double accept = std::max(2.0, 1.8*rLocal);
    if(dmin > accept){
        static char buf2[96]; snprintf(buf2,sizeof(buf2),"离点云 %.2f > %.2f", dmin, accept); g_surfWhy = buf2;
        if(std::getenv("SPRAYTRAJ_DUMP"))
            fprintf(stderr,"[曲面dbg]   投影后离点云 %.2fmm > %.2fmm (rLocal=%.2f cnt=%d off=%.2f)\n",
                            dmin, accept, rLocal, cnt, dot(q-cf, nf));
        return false;
    }
    return true;
}

// 从 start 沿横向等值线追踪一条行。返回行点/法向（追加到 out）。
//
// 行进方向不是"把 eRow 投到切平面"，而是每步重算的方向场 t = n × eIdx：
//   eIdx 是横向（步进）轴，行是"横向坐标 s(p)=p·eIdx 的等值线"。
//   等值线的面内梯度是 ∇f = eIdx - n(eIdx·n)，其正交方向正是 n × eIdx，
//   且 (n×eIdx)·eIdx ≡ 0 —— 沿它走横向坐标严格不变。
//   早期用"平行传输 eRow"会因 eRow 只与**平均**法向 n0 正交而长出横向分量：
//   实测一行 584mm 横漂 123mm（行距才 11.8mm），行成了斜穿全宽的乱线，
//   覆盖率只有 31%。
//
// bridge：数据空洞的跨越上限。线扫点云在反光/遮挡处会缺一小段，直接判"到边"
// 会把长行截断成两截。缺口服 fillHole 语义处理：能接上就补直线点。
static void surf_march(const SurfPatch& sp, Vec3 start, Vec3 eIdx, Vec3 t0, double step,
                       double bridge, int maxSteps, std::vector<Vec3>& outP, std::vector<Vec3>& outN,
                       bool isoH = false)
{
    Vec3 p, n;
    if(!surf_project(sp, start, p, n)) return;
    // ---- 上坡方向：决定"行"是什么曲线 ----
    // isoH=false（老路，可展柱面）：上坡方向固定为 eIdx，行 = 过轴平面截线。
    //   曲型板沿长边可展，截线就是一条沿长边的直线，正确。
    // isoH=true（回转/双曲率面）：上坡方向取**高度场的梯度**
    //   h(p)=dot(p-c0,n0) → ∇h = n0·(n·n0) - n（法向的面内分量取反）。
    //   行 = 等高线，即绕轴的闭合圈。
    //   梯度≈0 处是平面（底面/盖板），等值线无定义，直接判定行到此为止 ——
    //   这正好把"平面区"从行走向里摘出去，不会像老路那样直线横穿整个曲面。
    auto upOf = [&](const Vec3& nn)->Vec3{
        if(!isoH) return eIdx;
        Vec3 g = sp.n0*dot(nn, sp.n0) - nn;
        double l = len(g);
        if(l < 0.15) return Vec3{0,0,0};
        return g*(1.0/l);
    };
    // 起点切向用方向场给（t0 只用来定符号，保证朝期望的一侧走）
    Vec3 u0 = upOf(n);
    if(len(u0) < 0.2) return;
    Vec3 t = cross(n, u0);
    if(len(t) < 0.2) return;
    t = norm(t);
    if(dot(t, t0) < 0) t = t*(-1.0);
    outP.push_back(p); outN.push_back(n);
    // ---- 闭合判据（仅等高线模式）----
    // 等高线是**闭合圈**。没有闭合判据时，追踪走到起点附近不会停，会继续
    // 绕第二圈、第三圈……直到 maxSteps 用完。实测垃圾桶（锥面 Rmax≈140，
    // 一圈 383mm≈96 步）每行被追成 2956mm≈741 步（7.7 圈），7 层叠起来就是
    // 一根螺距 2.5mm 的螺旋弹簧 —— 界面上看到"一条轨迹绕那么多圈、设 30mm
    // 行距还是很密"，就是这个。
    // 判据：走过的弧长够长（防种子附近原地打转）且回到起点附近 → 闭合收工。
    const Vec3  pStart   = p;
    const double closeTol = 1.8*step;
    double sRun = 0.0;
    // 单步投影位移上限：按局部采样尺度自适应。稀疏区（线扫条纹间隙）局部平面
    // 拟合本身就有毫米级偏差，用固定值会把正常行进误判成"跳层"。
    double maxOff = std::max(2.0, 1.2*std::max(2.0,8.0*sp.spacing));
    const int nBridge = (int)std::floor(bridge/step);
    const bool dbg = (std::getenv("SPRAYTRAJ_DUMP") != nullptr);
    if(dbg) fprintf(stderr, "[曲面dbg] 行起点=(%.1f,%.1f,%.1f) t=(%.3f,%.3f,%.3f) step=%.2f maxOff=%.2f bridge=%d步\n",
                    p.x,p.y,p.z, t.x,t.y,t.z, step, maxOff, nBridge);
    for(int s=0; s<maxSteps; ++s){
        Vec3 q = p + t*step;
        Vec3 qp, qn;
        bool projOk = surf_project(sp, q, qp, qn);
        double offQ = projOk ? len(qp-q) : -1.0;
        if(projOk) maxOff = std::max(maxOff, 1.2*g_surfLocal);
        bool ok = projOk && (offQ <= maxOff);
        if(dbg && s<6) fprintf(stderr, "[曲面dbg]  step%d: proj=%d off=%.2f why=%s\n", s, (int)projOk, offQ, g_surfWhy);
        if(!ok){
            // ---- 数据空洞：沿切向直线探一段，接上就补点继续 ----
            int got = -1; Vec3 bp, bn;
            for(int j=2; j<=nBridge; ++j){
                Vec3 c = p + t*(step*j);
                bool pj = surf_project(sp, c, bp, bn);
                if(pj && len(bp-c) <= maxOff){ got = j; break; }
            }
            if(got < 0){
                if(dbg && s<10) fprintf(stderr, "[曲面dbg]  step%d: 出界（空探 %d 步无果）原因=%s\n", s, nBridge, g_surfWhy);
                break;
            }
            if(dbg && s<10) fprintf(stderr, "[曲面dbg]  step%d: 跨空洞 %d 步\n", s, got);
            // 补点法向：从洞前法向 slerp 到落点法向。
            // 早期用"常数 n"（而且 n 还是种子法向，见上面的更新修复）——
            // 洞两端各留一个姿态突跳。实测 row04 过洞处 26.8mm 枪位跳跃，
            // slerp 把它摊到每个补点上（每步 ~4mm 级）。
            {
                double cdot = std::max(-1.0, std::min(1.0, dot(n, bn)));
                double ang  = (cdot > 1.0-1e-12) ? 0.0 : std::acos(cdot);
                for(int j=1;j<got;++j){
                    double f = (double)j/(double)got;
                    Vec3 nj = n;
                    if(ang > 1e-6){
                        double s2 = std::sin(ang);
                        nj = norm(n*(std::sin((1.0-f)*ang)/s2) + bn*(std::sin(f*ang)/s2));
                    }
                    outP.push_back(p + t*(step*j));
                    outN.push_back(nj);
                }
            }
            p = bp; n = bn;
            Vec3 uu = upOf(n);
            Vec3 tn = (len(uu) > 0.2) ? cross(n, uu) : Vec3{0,0,0};
            if(len(tn) > 0.2){ tn = norm(tn); if(dot(tn,t) < 0) tn = tn*(-1.0); t = tn; }
            outP.push_back(p); outN.push_back(n);
            sRun += step*(double)got;
            if(isoH && sRun > 8.0*step && len(p-pStart) < closeTol){
                if(dbg) fprintf(stderr,"[曲面dbg] 闭合(跨洞后) 弧长=%.0fmm\n", sRun);
                break;
            }
            continue;
        }
        // 方向场给新切向，符号延续当前走向
        Vec3 uu = upOf(qn);
        if(len(uu) < 0.2) break;                                  // 平面区：等高线无定义，行到此为止
        Vec3 tn = cross(qn, uu);
        if(len(tn) < 0.2) break;                                  // 法向与横向轴近平行，行在此退化
        tn = norm(tn);
        if(dot(tn, t) < 0) tn = tn*(-1.0);
        if(dot(tn, t) < 0.342){                                   // 拐角 > 70°，视为到边
            if(dbg && s<10) fprintf(stderr, "[曲面dbg] step%d: 拐角过大 dot=%.3f\n", s, dot(tn,t));
            break;
        }
        t = tn;
        p = qp;
        // 法向必须随步进更新：过洞补点（bridge）会拿当前 n 当补点法向，
        // 之前这里漏了赋值，n 一直是种子点法向 —— 实测补点的法向与真实
        // 局部法向差 ~120°，枪位（X+100·n）被甩出 100+mm 的假跳跃。
        n = qn;
        outP.push_back(p); outN.push_back(qn);
        sRun += step;
        if(isoH && sRun > 8.0*step && len(p-pStart) < closeTol){
            if(dbg) fprintf(stderr,"[曲面dbg] 闭合 弧长=%.0fmm 步数=%d\n", sRun, s+1);
            break;
        }
    }
}

// 行内法向平滑：低通（行内 [1,2,3,2,1]） + 双向角速率限幅。位置一律不动。
// 抽成函数是因为它要跑两遍：行追踪后一遍（压逐点 PCA 噪声），
// 跨行姿态统一后再一遍（统一后的场在格边界处仍有跳变，必须再过一次限幅）。
static void smooth_row_normals(std::vector<Vec3>& rn, int lpPasses, double slewDeg)
{
    if(rn.size() < 5) return;
    static const double w[5] = {1.0, 2.0, 3.0, 2.0, 1.0};
    std::vector<Vec3> sm(rn.size());
    for(int pass=0; pass<lpPasses; ++pass){
        for(int i=0;i<(int)rn.size();++i){
            Vec3 acc{0,0,0};
            for(int d=-2; d<=2; ++d){
                int j = i+d;
                if(j<0 || j>=(int)rn.size()) continue;
                acc = acc + rn[j]*w[d+2];
            }
            acc = norm(acc);
            if(len(acc) > 0.5){
                if(dot(acc, rn[i]) < 0) acc = acc*(-1.0);   // 保持原朝向（防拟合符号随机）
                sm[i] = acc;
            } else sm[i] = rn[i];
        }
        rn.swap(sm);
    }
    if(slewDeg <= 0.0) return;
    const double slewRad = slewDeg*3.141592653589793/180.0;
    const double cs = std::cos(slewRad);
    for(int pass=0; pass<2; ++pass){                  // 0=正向 1=反向
        const bool fwd = (pass==0);
        const int stop = fwd ? (int)rn.size() : -1;
        const int step = fwd ? 1 : -1;
        for(int i = fwd ? 1 : (int)rn.size()-2; i!=stop; i+=step){
            const Vec3 a = rn[i-step], b = rn[i];
            const double c = dot(a,b);
            if(c >= cs) continue;                     // 没超限
            const double ang = std::acos(std::max(-1.0,std::min(1.0,c)));
            Vec3 perp = b - a*c;                      // b 在"垂直于 a"方向上的分量
            const double pl = len(perp);
            if(pl < 1e-9 || ang < 1e-6) continue;
            perp = perp*(1.0/pl);
            Vec3 nb = norm(a*std::cos(slewRad) + perp*std::sin(slewRad));
            if(len(nb) > 0.5 && dot(nb, b) > 0) rn[i] = nb;
        }
    }
}

// ==================== 重建曲面模型 ====================
// 把曲面块表达成一个**可求值、可求导、曲率连续**的模型：以平均平面 (eL,eS)
// 为参数域建规则网格，每格同时给出几何与姿态。
//
// 为什么需要它（实测依据，曲型板 u08 / curved）：
//   1) 显示网格只用"格内中位高度 + 空洞扩散"，从未做正则化 ——
//      相邻三角面二面角 p50=4.0° p95=26.0° max=84.9°，相邻顶点法向
//      p95=36.4°。真实拱面在 18mm 格距下法向变化应只有 ~2.6°。
//      所以用户说"重建曲面不规则、曲率不一致"是**成立的**。
//   2) 轨迹用的逐点 PCA 法向是重尾分布：p50 误差 0.65°、p90 6.7°、
//      p95 17.2°、max 47.7°（坏点集中在扫描条纹间隙与边缘外推处）。
//      经 standoff=100mm 放大成 30~74mm 的枪位横移，而实际行距只有
//      11.76mm —— 相邻行枪位因此大面积互相穿越（反序率 14%）。
//   3) 现有平滑只在**行内**做（[1,2,3,2,1]x4 + 6°/步限幅），管不了跨行
//      不一致；限幅还会把单点异常摊成连续 6°/步斜坡，形成可见"甩尾"。
//
// 数值策略为什么"几何用平滑高度、姿态用中位法向"：
//   直接从平滑高度场求梯度会把残余噪声再放大一次（离线实测法向改动
//   p95 38°、枪位反而更差），而"格内法向取中位再平滑"实测把枪位偏离
//   从 17.6mm 压到 7.3mm、反序率 14%->1%。二者共用同一网格，逻辑上仍是
//   一个模型。
static void sm_blur_scalar(std::vector<double>& A, int nu, int nv, double sigma)
{
    if(sigma <= 0.0) return;
    const int hw = std::max(1, (int)std::ceil(3.0*sigma));
    std::vector<double> k(2*hw+1); double ks=0;
    for(int i=-hw;i<=hw;++i){ double w=std::exp(-0.5*(i/sigma)*(i/sigma)); k[i+hw]=w; ks+=w; }
    for(auto& v : k) v/=ks;
    std::vector<double> tmp(A.size(), 1e300);
    for(int iv2=0;iv2<nv;++iv2) for(int iu=0;iu<nu;++iu){
        double s=0, ws=0;
        for(int d=-hw; d<=hw; ++d){
            int j=iu+d; if(j<0||j>=nu) continue;
            double v=A[(size_t)j*nv+iv2];
            if(v < 1e299){ s+=v*k[d+hw]; ws+=k[d+hw]; }
        }
        if(ws > 1e-12) tmp[(size_t)iu*nv+iv2] = s/ws;
    }
    for(int iv2=0;iv2<nv;++iv2) for(int iu=0;iu<nu;++iu){
        double s=0, ws=0;
        for(int d=-hw; d<=hw; ++d){
            int j=iv2+d; if(j<0||j>=nv) continue;
            double v=tmp[(size_t)iu*nv+j];
            if(v < 1e299){ s+=v*k[d+hw]; ws+=k[d+hw]; }
        }
        if(ws > 1e-12) A[(size_t)iu*nv+iv2] = s/ws;
    }
}

static void sm_blur_vec(std::vector<Vec3>& A, const std::vector<unsigned char>& occ,
                        int nu, int nv, int passes)
{
    std::vector<Vec3> tmp(A.size());
    for(int pass=0; pass<passes; ++pass){
        for(int iv2=0;iv2<nv;++iv2) for(int iu=0;iu<nu;++iu){
            size_t c=(size_t)iu*nv+iv2;
            if(!occ[c]) continue;
            Vec3 s{0,0,0}; int n=0;
            for(int a=-1;a<=1;++a) for(int b=-1;b<=1;++b){
                int ju=iu+a, jv=iv2+b;
                if(ju<0||jv<0||ju>=nu||jv>=nv) continue;
                size_t cc=(size_t)ju*nv+jv;
                if(!occ[cc]) continue;
                s = s + A[cc]; ++n;
            }
            if(n>0) tmp[c] = norm(s);
        }
        for(int iv2=0;iv2<nv;++iv2) for(int iu=0;iu<nu;++iu){
            size_t c=(size_t)iu*nv+iv2;
            if(occ[c] && len(tmp[c]) > 0.5) A[c] = tmp[c];
        }
    }
}

// 构建模型。cellHint<=0 时按点云跨度自适应（6~12mm）。
static bool build_surf_model(SurfPatch& sp, double cellHint)
{
    const int m = (int)sp.P.size();
    if(m < 100) return false;
    std::vector<Vec2> uv(m);
    double minU=1e18,maxU=-1e18,minV=1e18,maxV=-1e18;
    for(int i=0;i<m;++i){
        Vec3 d = sp.P[i]-sp.c0;
        uv[i] = Vec2{ dot(d, sp.eL), dot(d, sp.eS) };
        minU=std::min(minU,uv[i].x); maxU=std::max(maxU,uv[i].x);
        minV=std::min(minV,uv[i].y); maxV=std::max(maxV,uv[i].y);
    }
    const double wU=maxU-minU, wV=maxV-minV;
    if(wU < 2.0 || wV < 2.0) return false;

    double cell = (cellHint > 0.0) ? cellHint
                                   : std::max(6.0, std::min(12.0, std::max(wU,wV)/64.0));
    const int nu = (int)std::ceil(wU/cell)+1;
    const int nv = (int)std::ceil(wV/cell)+1;
    if(nu < 4 || nv < 4) return false;
    const size_t NG = (size_t)nu*(size_t)nv;
    if(NG > 400000) return false;

    std::vector<std::vector<int>> bk(NG);
    for(int i=0;i<m;++i){
        int iu=(int)std::floor((uv[i].x-minU)/cell);
        int iv2=(int)std::floor((uv[i].y-minV)/cell);
        iu=std::min(nu-1,std::max(0,iu)); iv2=std::min(nv-1,std::max(0,iv2));
        bk[(size_t)iu*nv+iv2].push_back(i);
    }

    std::vector<double> h(NG, 1e300);
    std::vector<Vec3> nn(NG, Vec3{0,0,0});
    std::vector<unsigned char> occ(NG, 0);
    int multi=0, empty=0;
    for(size_t c=0;c<NG;++c){
        const std::vector<int>& b = bk[c];
        if((int)b.size() < 3){ ++empty; continue; }
        std::vector<double> hs; hs.reserve(b.size());
        for(int i : b) hs.push_back(dot(sp.P[i]-sp.c0, sp.n0));
        // 多值格：同一 (u,v) 上存在两个表面（薄板背面/翻边/侧面），
        // 中位高度会取到二者之间，梯度完全错误 —— 直接丢弃，交给邻域扩散。
        {
            std::vector<double> t = hs;
            size_t k1=(size_t)((t.size()-1)*0.25), k3=(size_t)((t.size()-1)*0.75);
            std::nth_element(t.begin(), t.begin()+k1, t.end()); double q1=t[k1];
            std::nth_element(t.begin(), t.begin()+k3, t.end()); double q3=t[k3];
            if(q3 - q1 > 8.0){ ++multi; continue; }
        }
        {
            std::vector<double> t = hs;
            size_t k=(t.size()-1)/2;
            std::nth_element(t.begin(), t.begin()+k, t.end());
            h[c] = t[k];
        }
        // 中位法向（逐分量中位，抗坏点；与逐点法向重尾正好互补）
        {
            std::vector<double> a,bb,cc;
            a.reserve(b.size()); bb.reserve(b.size()); cc.reserve(b.size());
            for(int i : b){ a.push_back(sp.N[i].x); bb.push_back(sp.N[i].y); cc.push_back(sp.N[i].z); }
            size_t k=(a.size()-1)/2;
            std::nth_element(a.begin(), a.begin()+k, a.end()); double mx=a[k];
            std::nth_element(bb.begin(), bb.begin()+k, bb.end()); double my=bb[k];
            std::nth_element(cc.begin(), cc.begin()+k, cc.end()); double mz=cc[k];
            Vec3 v = norm(Vec3{mx,my,mz});
            if(len(v) < 0.5) v = sp.n0;
            if(dot(v, sp.n0) < 0) v = v*(-1.0);
            nn[c] = v;
        }
        occ[c] = 1;
    }
    int occN=0; for(size_t c=0;c<NG;++c) if(occ[c]) ++occN;
    if(occN < 20) return false;

    // 空格：邻域已填值扩散（保持网格连续；这些格 mOcc=0，插值时不采用）
    for(int it=0; it<8; ++it){
        bool any=false;
        std::vector<double> nh = h;
        std::vector<Vec3> nv2 = nn;
        for(int iv2=0;iv2<nv;++iv2) for(int iu=0;iu<nu;++iu){
            size_t c=(size_t)iu*nv+iv2;
            if(h[c] < 1e299) continue;
            double s=0; int k=0; Vec3 vs{0,0,0};
            for(int a=-1;a<=1;++a) for(int b=-1;b<=1;++b){
                int ju=iu+a, jv=iv2+b;
                if(ju<0||jv<0||ju>=nu||jv>=nv) continue;
                size_t cc=(size_t)ju*nv+jv;
                if(h[cc] >= 1e299) continue;
                s += h[cc]; vs = vs + nn[cc]; ++k;
            }
            if(k){ nh[c]=s/(double)k; Vec3 w=norm(vs); if(len(w)>0.5) nv2[c]=w; any=true; }
        }
        h.swap(nh); nn.swap(nv2);
        if(!any) break;
    }

    // 几何：归一化高斯（σ=1.5 格），把"格内中位"的台阶磨成连续曲率
    sm_blur_scalar(h, nu, nv, 1.5);
    // 姿态：3x3 归一化平滑一遍（只在实测格上，避免空洞带偏）
    sm_blur_vec(nn, occ, nu, nv, 1);

    sp.mNu=nu; sp.mNv=nv; sp.mCell=cell; sp.mU0=minU; sp.mV0=minV;
    sp.mH.swap(h); sp.mNm.swap(nn); sp.mOcc.swap(occ);
    sp.modelOk = true;
    fprintf(stderr, "[曲面] 曲面模型: %dx%d 格=%.1fmm 占用=%d/%d 多值格=%d 空格=%d\n",
            nu, nv, cell, occN, (int)NG, multi, empty);
    return true;
}

// 只查**几何**（平滑高度）：允许用空格扩散出来的值 —— 高度场本身连续，
// 扩散值只差亚毫米，用它保住网格连续比留个洞好得多。
static bool surf_model_height(const SurfPatch& sp, const Vec3& p, double& hout)
{
    if(!sp.modelOk) return false;
    Vec3 d = p - sp.c0;
    double fu = (dot(d, sp.eL)-sp.mU0)/sp.mCell;
    double fv = (dot(d, sp.eS)-sp.mV0)/sp.mCell;
    int i0=(int)std::floor(fu), j0=(int)std::floor(fv);
    if(i0 < 0 || j0 < 0 || i0+1 >= sp.mNu || j0+1 >= sp.mNv) return false;
    size_t c[4] = { (size_t)i0*sp.mNv+j0,       (size_t)(i0+1)*sp.mNv+j0,
                    (size_t)i0*sp.mNv+j0+1,     (size_t)(i0+1)*sp.mNv+j0+1 };
    for(int q=0;q<4;++q) if(sp.mH[c[q]] > 1e299) return false;
    const double a = fu-i0, b = fv-j0;
    const double w[4] = { (1-a)*(1-b), a*(1-b), (1-a)*b, a*b };
    double hs=0;
    for(int q=0;q<4;++q) hs += w[q]*sp.mH[c[q]];
    hout = hs;
    return true;
}

// 查**姿态**：四邻格必须都是实测格，否则返回 false（调用方回退到逐点法向）。
static bool surf_model_query(const SurfPatch& sp, const Vec3& p, double& hout, Vec3& nout)
{
    if(!sp.modelOk) return false;
    Vec3 d = p - sp.c0;
    double fu = (dot(d, sp.eL)-sp.mU0)/sp.mCell;
    double fv = (dot(d, sp.eS)-sp.mV0)/sp.mCell;
    int i0=(int)std::floor(fu), j0=(int)std::floor(fv);
    if(i0 < 0 || j0 < 0 || i0+1 >= sp.mNu || j0+1 >= sp.mNv) return false;
    size_t c[4] = { (size_t)i0*sp.mNv+j0,       (size_t)(i0+1)*sp.mNv+j0,
                    (size_t)i0*sp.mNv+j0+1,     (size_t)(i0+1)*sp.mNv+j0+1 };
    for(int q=0;q<4;++q) if(!sp.mOcc[c[q]]) return false;
    const double a = fu-i0, b = fv-j0;
    const double w[4] = { (1-a)*(1-b), a*(1-b), (1-a)*b, a*b };
    double hs=0; Vec3 ns{0,0,0};
    for(int q=0;q<4;++q){
        if(sp.mH[c[q]] > 1e299) return false;
        hs += w[q]*sp.mH[c[q]];
        ns = ns + sp.mNm[c[q]]*w[q];
    }
    Vec3 v = norm(ns);
    if(len(v) < 0.5) return false;
    if(dot(v, sp.n0) < 0) v = v*(-1.0);
    hout = hs; nout = v;
    return true;
}

// 在曲面块上生成蛇形行（就地写回 sp.rowP/rowN/rowBeg）
static bool gen_surface_rows(SurfPatch& sp, double L, double pstep, double fillHole,
                             int dirMode, double speed)
{
    const int m = (int)sp.P.size();
    if(m < 50 || L <= 0 || pstep <= 0) return false;
    if(!sp.modelOk) build_surf_model(sp, 0.0);   // 显示网格已建则复用，避免重复计算
    // ---- 姿态限幅上限：按末端允许角速度反算 ----
    //   Δmax[°/点] = ω_max[°/s] · pstep[mm] / speed[mm/s]
    // 老版本固定 6°/点：点距 4mm、速度 200mm/s 时等于 **300°/s** 的末端角速度，
    // 多数六轴机器人（典型 100~200°/s）跟不上，实际表现是到位滞后/超调，
    // 视觉上就是"轨迹在拐角处甩出去"。默认 ω_max=150°/s → 3.0°/点。
    // 夹在 [0.5, 6.0]：点距很大或速度很慢时不至于把真实折边也抹平。
    const double wMax = 150.0;   // °/s
    const double slewDeg = std::min(6.0, std::max(0.5, wMax*pstep/std::max(1.0, speed)));
    fprintf(stderr, "[曲面] 姿态限幅 %.2f°/点 = %.0f°/s (点距%.1fmm 速度%.0fmm/s)\n",
            slewDeg, slewDeg*std::max(1.0,speed)/std::max(1e-9,pstep), pstep, speed);
    // 空洞跨越上限：跟平面路径的 fillHole 语义一致，但曲面行不宜跨太远
    const double bridge = std::min(fillHole, 3.0*L);

    // ---- 面内主轴：把曲面点投到平均平面，取 2D PCA ----
    std::vector<Vec2> uv(m);
    double minU=1e18,maxU=-1e18,minV=1e18,maxV=-1e18;
    for(int i=0;i<m;++i){
        Vec3 d = sp.P[i]-sp.c0;
        uv[i] = Vec2{ dot(d, sp.eL), dot(d, sp.eS) };
        minU=std::min(minU,uv[i].x); maxU=std::max(maxU,uv[i].x);
        minV=std::min(minV,uv[i].y); maxV=std::max(maxV,uv[i].y);
    }
    double sxx=0,syy=0,sxy=0;
    for(int i=0;i<m;++i){ sxx+=uv[i].x*uv[i].x; syy+=uv[i].y*uv[i].y; sxy+=uv[i].x*uv[i].y; }
    double ang = 0.5*std::atan2(2*sxy, sxx-syy);
    Vec3 a1 = norm(sp.eL*std::cos(ang) + sp.eS*std::sin(ang));
    Vec3 a2 = norm(cross(sp.n0, a1));
    if(len(a1) < 0.5 || len(a2) < 0.5){ a1 = sp.eL; a2 = sp.eS; }

    // 行方向 = 较长那一维；扫描推进方向 = 另一维
    double ext1=0, ext2=0;
    {
        double lo1=1e18,hi1=-1e18,lo2=1e18,hi2=-1e18;
        for(int i=0;i<m;++i){
            double t1 = dot(sp.P[i], a1), t2 = dot(sp.P[i], a2);
            lo1=std::min(lo1,t1); hi1=std::max(hi1,t1);
            lo2=std::min(lo2,t2); hi2=std::max(hi2,t2);
        }
        ext1 = hi1-lo1; ext2 = hi2-lo2;
    }
    Vec3 eRow = a1, eIdx = a2;
    if(ext1 < ext2) std::swap(eRow, eIdx);          // 保证 eRow 是长的那一维
    if(dirMode == 2) std::swap(eRow, eIdx);         // 沿短边喷：行平行短边

    // ---- 行走向自检：可展柱面 vs 回转/双曲率面 ----
    // "过轴平面截线"（沿 eRow 走的等 eIdx 曲线）只有在曲面沿 eRow **可展**
    // 时才是一条像样的行。回转体（圆锥/圆柱/槽盘）沿 eRow 走一圈方位角要变
    // 几十度，截线会退化成沿母线的一束双曲线：外圈挤、内圈疏，端点在轴心汇聚。
    // 实测垃圾桶(锥台+底面)：27 行全部从上沿直插到底面中心，总长 12.7m、
    // 母线间距 20.4mm(外)~7.0mm(内)；等高线圈只要 17 圈 / 11.4m，间距均匀。
    //
    // 判据用"沿长轴的高程起伏"而不是"法向离散"：后者在双曲率板上也会很大，
    // 实测曲型板 p90=50.8° 被误判成回转体，行数从 18 涨到 34、rms 8.9→18.9。
    // 高程起伏ratio = 箱内 h 的极差 / 箱内 eRow 跨度：
    //   柱面的母线是沿 eRow 的直线 → h 沿 eRow 不变 → ratio ~0.02（只剩箱宽引起的那点）
    //   圆锥的截线是双曲线      → h 沿 eRow 剧变   → ratio ~0.3
    // 门限 0.10，两边各有 3~4 倍余量。
    bool isoMode = false;
    double isoG = 1.0;          // 高度方向走 1mm 时曲面上实际走多少 mm（= 1/|∇h|）
    {
        const int NB = 12;
        double rlo=1e18, rhi=-1e18, ilo=1e18, ihi=-1e18;
        for(int i=0;i<m;++i){
            double tr=dot(sp.P[i], eRow), ti=dot(sp.P[i], eIdx);
            rlo=std::min(rlo,tr); rhi=std::max(rhi,tr);
            ilo=std::min(ilo,ti); ihi=std::max(ihi,ti);
        }
        if(ihi-ilo > 1e-6 && rhi-rlo > 1e-6){
            std::vector<double> rat;
            for(int b=0;b<NB;++b){
                double a0 = ilo + (ihi-ilo)*(double)b/(double)NB;
                double a1 = ilo + (ihi-ilo)*(double)(b+1)/(double)NB;
                double tlo=1e18, thi=-1e18, hlo=1e18, hhi=-1e18;
                int cnt=0;
                for(int i=0;i<m;++i){
                    double ti = dot(sp.P[i], eIdx);
                    if(ti < a0 || ti >= a1) continue;
                    double tr = dot(sp.P[i], eRow);
                    double h  = dot(sp.P[i], sp.n0);
                    tlo=std::min(tlo,tr); thi=std::max(thi,tr);
                    hlo=std::min(hlo,h);  hhi=std::max(hhi,h);
                    ++cnt;
                }
                if(cnt < 40 || thi-tlo < 1e-6) continue;
                rat.push_back((hhi-hlo)/(thi-tlo));
            }
            double ratio = 0.0;
            if(rat.size() >= 4){
                std::nth_element(rat.begin(), rat.begin()+rat.size()/2, rat.end());
                ratio = rat[rat.size()/2];
            }
            // 第二个条件：投影是否绕满一圈。回转体的 (u,v) 投影必然覆盖 360° 方位；
            // 一片板（哪怕弯成 60° 的弧板）投影出来是个矩形/细长条，方位只覆盖
            // 100°~150°。实测曲型板 ratio=0.372（它确实是双曲率的），只靠起伏
            // 判据会被误判成回转体，行数从 18 涨到 30、覆盖率 92%→55%。
            const int NS = 36;
            std::vector<unsigned char> sec(NS, 0);
            int npj = 0;
            for(int i=0;i<m;++i){
                Vec3 d = sp.P[i]-sp.c0;
                double aa = dot(d, sp.eL), bb = dot(d, sp.eS);
                if(std::fabs(aa) < 1e-9 && std::fabs(bb) < 1e-9) continue;
                double th = std::atan2(bb, aa);
                int k = (int)((th + 3.14159265358979)/(2.0*3.14159265358979)*(double)NS);
                if(k<0) k=0; if(k>=NS) k=NS-1;
                sec[k]=1; ++npj;
            }
            int nsec=0; for(int k=0;k<NS;++k) if(sec[k]) ++nsec;
            const double cover = npj>0 ? (double)nsec/(double)NS : 0.0;
            // 第三个条件：外接框不能太狭长。回转体的投影接近正方形/圆；
            // 一片长条板（曲型板实测 548×176，宽高比 3.1）方位也能覆盖满一圈，
            // 但它显然不是回转体 —— 少了这条会把它带进等高线模式，行被切成
            // 6 个点的碎片（覆盖率 92%→55%）。
            const double aspect = (ext1 > ext2) ? ext1/std::max(1e-9,ext2)
                                                : ext2/std::max(1e-9,ext1);
            isoMode = (ratio > 0.10) && (cover >= 5.0/6.0) && (aspect < 1.8);
            fprintf(stderr, "[曲面] 行走向自检: 高程起伏=%.3f 方位覆盖=%.0f%% 宽高比=%.2f → %s\n",
                    ratio, cover*100.0, aspect,
                    isoMode ? "回转/双曲率，按等高线圈布行" : "可展/开放面，按平面截线布行");
        }
    }

    // ---- 定距：两条路分别是"横向等距锚点"和"等高线逐层等弧长推进" ----
    std::vector<double> iv(m);
    std::vector<unsigned char> isoOk(m, 1);
    std::vector<double> anchors;          // 每一行的锚（平面=横向坐标，iso=高度）
    auto qt=[&](std::vector<double>& v, double f)->double{
        if(v.empty()) return 0.0;
        size_t k=(size_t)((v.size()-1)*f);
        std::nth_element(v.begin(), v.begin()+k, v.end());
        return v[k];
    };
    double lo=0.0, hi=0.0, D=0.0, gap=L;
    if(isoMode){
        // 锚点沿**高度**分层。只有"法向明显偏离 n0"的点参与 —— 平面区
        // （底面/盖板）的高度是常数，既会把范围压扁，也生不出等高线。
        double hlo=1e18, hhi=-1e18;
        for(int i=0;i<m;++i){
            iv[i]=dot(sp.P[i]-sp.c0, sp.n0);
            Vec3 g = sp.n0*dot(sp.N[i], sp.n0) - sp.N[i];
            if(len(g) < 0.15){ isoOk[i]=0; continue; }
            hlo=std::min(hlo,iv[i]); hhi=std::max(hhi,iv[i]);
        }
        if(hhi-hlo < 1e-6) return false;
        lo=hlo; hi=hhi;
        // |∇h| 随高度的查找表（1mm 箱取中位）。**不能用全局平均值**：
        // 曲型板（碗形）的 |∇h| 从 0.1 到 1.0，用平均 0.6 会让底部行距 72mm、
        // 陡壁 7mm，实测覆盖率只有 54.7%、贴合残差 18.9mm。
        const int NH = (int)std::floor(hhi-hlo)+1;
        std::vector<std::vector<double>> gb(NH);
        double gAll=0; int gAllN=0;
        for(int i=0;i<m;++i){
            if(!isoOk[i]) continue;
            Vec3 g = sp.n0*dot(sp.N[i], sp.n0) - sp.N[i];
            double gl = len(g);
            int b = (int)std::floor(iv[i]-hlo); if(b<0) b=0; if(b>=NH) b=NH-1;
            gb[b].push_back(gl); gAll += gl; ++gAllN;
        }
        const double gGlob = gAllN>0 ? std::max(0.10, gAll/(double)gAllN) : 0.5;
        std::vector<double> gm(NH, gGlob);
        for(int b=0;b<NH;++b){
            if(gb[b].size() < 5) continue;
            std::nth_element(gb[b].begin(), gb[b].begin()+gb[b].size()/2, gb[b].end());
            gm[b] = std::max(0.08, gb[b][gb[b].size()/2]);
        }
        // 逐层推进：每层的间隔 = L·|∇h|(该层)，等价于"沿曲面量每层都走 L"
        double hcur = hlo;
        for(int it=0; it<2000 && hcur <= hhi; ++it){
            anchors.push_back(hcur);
            int b = (int)std::floor(hcur-hlo); if(b<0) b=0; if(b>=NH) b=NH-1;
            hcur += L*gm[b];
        }
        // 收尾贴边：逐层推进必然在 hhi 前停下，留下最多一整个层距的空带。
        // 实测垃圾桶（hlo=3.7 hhi=187）停在 161.3，顶部 30mm 没有行 —— 覆盖率
        // 100%→97.2%，贴合残差 6.5→10.2mm。这里把末层推到 hhi，或差得多时补一层。
        if(anchors.size() >= 2){
            const double lastStep = anchors.back() - anchors[anchors.size()-2];
            const double rest = hhi - anchors.back();
            if(rest > 0.55*lastStep)            anchors.push_back(hhi);
            else if(rest > 1e-6)                anchors.back() = hhi;
        }
        if(anchors.size() < 2) return false;
        gap = L;                 // 逐层推进的间隔就是 L，实际行距也就是 L
        fprintf(stderr, "[曲面] 等高线分层: 高差=%.1fmm %zu 层 行距=%.1fmm\n",
                hhi-hlo, anchors.size(), gap);
    }else{
        for(int i=0;i<m;++i) iv[i]=dot(sp.P[i], eIdx);
        std::vector<double> ivs = iv;
        lo = qt(ivs, 0.005); hi = qt(ivs, 0.995);
        D = hi - lo;
        if(D < L*0.5) return false;
        int N = (int)std::floor(D/L + 0.5);
        if(N < 1) N = 1;
        if(N > 2000) N = 2000;
        gap = D/(double)N;
        anchors.reserve(N);
        for(int k=0;k<N;++k) anchors.push_back(lo + (k+0.5)*gap);
    }
    sp.rowGap = gap;                 // 沿曲面的实际行距
    if(std::getenv("SPRAYTRAJ_DUMP"))
        fprintf(stderr, "[曲面dbg] ext(a1)=%.1f ext(a2)=%.1f eRow=(%.3f,%.3f,%.3f) eIdx=(%.3f,%.3f,%.3f) "
                        "lo=%.1f hi=%.1f D=%.1f N=%zu gap=%.2f iso=%d\n",
                ext1, ext2, eRow.x,eRow.y,eRow.z, eIdx.x,eIdx.y,eIdx.z, lo, hi, D,
                anchors.size(), gap, (int)isoMode);

    // ---- 逐行独立起步 ----
    // 早期版本是"上一行末点沿 t×n 侧移一格 gap 起下一行"，结果行会随曲面法向的
    // 微小偏差逐行错位累积：实测行距看似 11.8mm，但行与行互相穿插，覆盖率只有
    // 31%。改成每行都以"横向坐标 = lo+(k+0.5)*gap"为锚点从点云里独立取种子，
    // 侧移只用来决定行的走向（不负责定位），错位不会传递到下一行。
    //
    // 单行从种子向**两侧**分别追踪再拼接：从边缘起步时行会被截断成半行，
    // 从中间起步则能走满整条弦。
    sp.rowP.clear(); sp.rowN.clear(); sp.rowBeg.clear();
    sp.rowBeg.push_back(0);          // 首行起点哨兵（发射段按 rowBeg[r]..rowBeg[r+1] 切行）
    // 单行步数上限。平面截线是"一串从边到边的开放曲线"，2×(外接尺寸) 够用；
    // 等高线圈是**闭合**的，最长一圈 = 2π·Rmax，实测曲型板（碗形，Rmax=334）
    // 一圈要 525 步，而按外接尺寸只给 373 步 —— 行被腰斩成 340mm 的碎片，
    // 覆盖率只有 54.8%、残差 15.8mm。
    int maxSteps = (int)std::ceil((ext1+ext2)*2.0/pstep) + 32;
    if(isoMode){
        double rmax = 0.0;
        for(int i=0;i<m;++i) rmax = std::max(rmax, len(sp.P[i]-sp.c0));
        maxSteps = (int)std::ceil(2.0*3.14159265358979*rmax*1.3/pstep) + 32;
    }
    const bool dbg = (std::getenv("SPRAYTRAJ_DUMP") != nullptr);
    int rowsMade = 0;

    // 一条完整行：iso 模式追一圈闭合圈；平面截线模式两侧各追一遍再拼接
    auto marchBoth = [&](Vec3 seed, std::vector<Vec3>& oP, std::vector<Vec3>& oN){
        oP.clear(); oN.clear();
        if(isoMode){
            // 等高线是闭合圈：从种子朝一侧追满一圈就回到种子，另一侧不需要再追。
            // 以前两侧都追、且没有闭合判据，两侧各自绕到步数上限（各 3.85 圈），
            // 拼起来 7.7 圈 —— 行距 30mm 下锥面 7 层变成一根密集螺旋。
            surf_march(sp, seed, eIdx, eRow, pstep, bridge, maxSteps, oP, oN, true);
            return;
        }
        std::vector<Vec3> fP, fN, bP, bN;
        surf_march(sp, seed, eIdx, eRow*(-1.0), pstep, bridge, maxSteps, bP, bN, isoMode);
        surf_march(sp, seed, eIdx, eRow,        pstep, bridge, maxSteps, fP, fN, isoMode);
        oP.clear(); oN.clear();
        oP.reserve(bP.size()+fP.size()); oN.reserve(bP.size()+fN.size());
        // 反向前半是"seed, seed-1步, seed-2步…"，倒序后末尾那个点就是 seed
        for(size_t i=bP.size(); i-->1;){ oP.push_back(bP[i]); oN.push_back(bN[i]); }
        for(size_t i=0;i<fP.size();++i){ oP.push_back(fP[i]); oN.push_back(fN[i]); }
    };

    for(size_t k=0; k<anchors.size(); ++k){
        const double tk = anchors[k];
        int seed = -1; double bestd = 1e18;
        for(int i=0;i<m;++i){
            if(!isoOk[i]) continue;              // 平面点不参与等高线取种
            double d = std::fabs(iv[i]-tk);
            if(d<bestd){ bestd=d; seed=i; }
        }
        if(seed < 0) continue;
        std::vector<Vec3> rp, rn;
        marchBoth(sp.P[seed], rp, rn);
        if(rp.size() < 2) continue;
        // 闭合等高线：向两侧各追一遍会绕满一圈并在种子附近相遇，尾部是重复段。
        // 不去掉的话蛇形串接会在圈末原地折返 180°（枪位图上是个尖刺）。
        if(isoMode && rp.size() > 8){
            while(rp.size() > 8 && len(rp.back()-rp.front()) < 1.5*pstep){
                rp.pop_back(); rn.pop_back();
            }
        }
        // ---- 姿态来自重建曲面模型（跨行一致化）----
        // 逐点 PCA 法向的重尾误差（p95 17°、max 48°）经 100mm 枪位放大成
        // 30~74mm 横移，而行距只有 11.76mm，相邻行枪位就大面积互相穿越。
        // 行内低通只在**行内**做，跨行不一致它管不到。这里把模型网格上的
        // 中位法向插值回来：位置一个点都不动（覆盖率/残留不变），只换姿态。
        // 查不到（点云覆盖不到的外推区）就保留原法向，不做猜测。
        if(sp.modelOk){
            int hit=0;
            for(size_t i=0;i<rp.size();++i){
                double hq; Vec3 nq;
                if(surf_model_query(sp, rp[i], hq, nq)){ rn[i]=nq; ++hit; }
            }
            if(dbg) fprintf(stderr,"[曲面dbg] 行%02d 模型法向命中 %d/%d\n",
                            rowsMade, hit, (int)rp.size());
        }
        // ---- 行内法向平滑（枪姿连续性）----
        // 逐点局部拟合的法向有 1~2° 级噪声（线扫条纹间隙里更大），加上换墙/
        // 卷边处的几度~几十度摆动，4mm 步进的枪位（X+100·n）会冒出 8~50mm 的
        // "姿态跳跃"，渲染出来就是"轨迹乱七八糟"。枪头旋转本身有惯性/限速，
        // 这里分两段处理（位置 rp 一律不动，覆盖率/贴合残差完全不受影响）：
        //   1) 行内 [1,2,3,2,1] 低通 ×4 遍（等效高斯 σ≈6.5mm）：压逐点法向噪声；
        //   2) 双向角速率限幅：换墙/卷边的折角摊成缓转。
        //      折角处一两个点转几十度，枪位（100mm 外）横甩 ±50mm，
        //      是"轨迹乱"的最大来源；限幅后枪头像真机一样匀速转过去。
        smooth_row_normals(rn, 4, slewDeg);
        if((rowsMade & 1) != 0){                     // 蛇形：奇数行反向存放
            std::reverse(rp.begin(), rp.end());
            std::reverse(rn.begin(), rn.end());
        }
        if(dbg){
            double s0 = dot(rp.front(), eIdx), s1 = dot(rp.back(), eIdx);
            fprintf(stderr,"[曲面dbg] 行%02d 锚=%.1f 实到=%.1f..%.1f(漂%.1f) 点数=%zu\n",
                    rowsMade, tk, std::min(s0,s1), std::max(s0,s1), std::fabs(s0-s1), rp.size());
        }
        sp.rowP.insert(sp.rowP.end(), rp.begin(), rp.end());
        sp.rowN.insert(sp.rowN.end(), rn.begin(), rn.end());
        sp.rowBeg.push_back((int)sp.rowP.size());
        ++rowsMade;
    }
    sp.nRow = rowsMade;

    // ---- 跨行姿态一致化：同一位置、不同行的法向必须一致 ----
    // 行内低通（上面每一步）只在**沿行方向**平滑，管不了"相邻行在同一 (u,v)
    // 处法向不一致"。而相邻行法向差 6° 就足以在 100mm 枪位处横移 10.5mm >>
    // 行距 11.76mm，视觉上就是两行互相穿越。
    // 做法：把全部行的法向按 (u,v) 网格（与几何模型同一套网格，8.5mm）取中位，
    // 再叠 4 遍 3x3 归一化平滑（有效尺度 ~50mm，跨 3~4 行），最后插值回每个行点。
    // 位置 rowP 一个点都不动 —— 覆盖率/贴合残差不受影响，只统一姿态。
    if(sp.modelOk && !sp.rowP.empty()){
        // 姿态场的格宽必须**大于行距**：几何模型用 8.5mm 格，比行距 11.76mm 还小，
        // 相邻两行会落在不同格里，"跨行取中位"根本无从谈起（实测命中率只有 53%）。
        // 这里独立用 1.3×行距（≈15mm）的分辨率，仍在同一个 (u,v) 参数域上。
        const double cell2 = std::max(sp.mCell, 1.3*L);
        const double spanU = (double)(sp.mNu-1)*sp.mCell;
        const double spanV = (double)(sp.mNv-1)*sp.mCell;
        const int nu=(int)std::ceil(spanU/cell2)+1;
        const int nv=(int)std::ceil(spanV/cell2)+1;
        const size_t NG=(size_t)nu*(size_t)nv;
        std::vector<std::vector<int>> bk(NG);
        for(size_t k=0;k<sp.rowP.size();++k){
            Vec3 d=sp.rowP[k]-sp.c0;
            int iu=(int)std::floor((dot(d,sp.eL)-sp.mU0)/cell2);
            int iv2=(int)std::floor((dot(d,sp.eS)-sp.mV0)/cell2);
            if(iu<0||iv2<0||iu>=nu||iv2>=nv) continue;
            bk[(size_t)iu*nv+iv2].push_back((int)k);
        }
        std::vector<Vec3> fn(NG, sp.n0);
        std::vector<unsigned char> foc(NG,0);
        for(size_t c=0;c<NG;++c){
            if(bk[c].empty()) continue;
            const size_t n=bk[c].size();
            std::vector<double> a(n),b(n),e(n);
            for(size_t q=0;q<n;++q){
                const Vec3& v=sp.rowN[bk[c][q]];
                a[q]=v.x; b[q]=v.y; e[q]=v.z;
            }
            size_t k2=(n-1)/2;
            std::nth_element(a.begin(),a.begin()+k2,a.end()); double mx=a[k2];
            std::nth_element(b.begin(),b.begin()+k2,b.end()); double my=b[k2];
            std::nth_element(e.begin(),e.begin()+k2,e.end()); double mz=e[k2];
            Vec3 v=norm(Vec3{mx,my,mz});
            if(len(v)<0.5) v=sp.n0;
            if(dot(v,sp.n0)<0) v=v*(-1.0);
            fn[c]=v; foc[c]=1;
        }
        sm_blur_vec(fn, foc, nu, nv, 4);   // 只在"确有行点"的格上平滑，别把点云法向混进来
        int rep=0;
        for(size_t k=0;k<sp.rowP.size();++k){
            Vec3 d=sp.rowP[k]-sp.c0;
            double fu=(dot(d,sp.eL)-sp.mU0)/cell2;
            double fv=(dot(d,sp.eS)-sp.mV0)/cell2;
            int i0=(int)std::floor(fu), j0=(int)std::floor(fv);
            if(i0<0||j0<0||i0+1>=nu||j0+1>=nv) continue;
            const double aa=fu-i0, bb=fv-j0;
            const double w4[4]={ (1-aa)*(1-bb), aa*(1-bb), (1-aa)*bb, aa*bb };
            const size_t cc[4]={ (size_t)i0*nv+j0, (size_t)(i0+1)*nv+j0,
                                 (size_t)i0*nv+j0+1, (size_t)(i0+1)*nv+j0+1 };
            Vec3 v{0,0,0};
            if(foc[cc[0]]&&foc[cc[1]]&&foc[cc[2]]&&foc[cc[3]]){
                Vec3 s{0,0,0};
                for(int q=0;q<4;++q) s = s + fn[cc[q]]*w4[q];
                v = norm(s);
            }else{
                // 四邻不全有行点（行端 / 点云边缘）：取 3x3 邻域里最近的有效格。
                // 不能保留原来的噪声法向 —— 那正是枪位图上两处"尖峰"的来源
                // （实测 max 26.5mm，够真机甩一下）。
                double bd=1e18;
                for(int di=-1;di<=1;++di) for(int dj=-1;dj<=1;++dj){
                    int ii=i0+di, jj=j0+dj;
                    if(ii<0||jj<0||ii>=nu||jj>=nv) continue;
                    size_t cx=(size_t)ii*nv+jj;
                    if(!foc[cx] || len(fn[cx]) < 0.5) continue;
                    double dd=(double)(di*di+dj*dj);
                    if(dd < bd){ bd=dd; v=fn[cx]; }
                }
            }
            if(len(v)<0.5) continue;
            if(dot(v, sp.rowN[k])<0) v=v*(-1.0);
            sp.rowN[k]=v; ++rep;
        }
        // 统一后的场在格边界上仍可能有几度跳变，按行再过一遍低通+限幅：
        // 否则 100mm 枪位会把这几度放大成二十几毫米的"枪位跳跃"（seg_max 实测
        // 从 12.9mm 涨到 26mm 就是这个原因）。
        for(int r=0;r<sp.nRow;++r){
            const int b=sp.rowBeg[r], e=sp.rowBeg[r+1];
            if(e-b < 5) continue;
            std::vector<Vec3> seg(sp.rowN.begin()+b, sp.rowN.begin()+e);
            smooth_row_normals(seg, 2, slewDeg);
            std::copy(seg.begin(), seg.end(), sp.rowN.begin()+b);
        }
        if(std::getenv("SPRAYTRAJ_DUMP"))
            fprintf(stderr,"[曲面dbg] 跨行姿态一致化: 行点 %d/%d 被统一\n", rep, (int)sp.rowP.size());
    }
    if(rowsMade < 1) return false;

    // 调试：SPRAYTRAJ_DUMP=1 时导出行点（含横向坐标 s = p·eIdx）与行边界。
    // 曲面件的 3D 截图分不清"行内抖动"和"端部卷边"，离线按 s 看最直观。
    if(dbg){
        fprintf(stderr, "[曲面dbg] eIdx=(%.4f,%.4f,%.4f) eRow=(%.4f,%.4f,%.4f) n0=(%.4f,%.4f,%.4f) gap=%.3f\n",
                eIdx.x,eIdx.y,eIdx.z, eRow.x,eRow.y,eRow.z, sp.n0.x,sp.n0.y,sp.n0.z, gap);
        if(FILE* fr = fopen("surf_rows.csv","w")){
            fprintf(fr,"row,k,x,y,z,s,nx,ny,nz\n");
            for(int r=0;r<sp.nRow;++r)
                for(int k=sp.rowBeg[r]; k<sp.rowBeg[r+1]; ++k)
                    fprintf(fr,"%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", r, k,
                            sp.rowP[k].x, sp.rowP[k].y, sp.rowP[k].z, dot(sp.rowP[k], eIdx),
                            sp.rowN[k].x, sp.rowN[k].y, sp.rowN[k].z);
            fclose(fr);
        }
    }

    // ---- 覆盖率：每个曲面点到最近行点的距离 ----
    {
        KdTree rt; rt.init(sp.rowP);
        int cov=0;
        const double tol = std::max(gap*0.75, L*0.75);
        double rss=0; int rn=0;
        for(int i=0;i<m;++i){
            KnnBuf b; rt.knn(sp.P[i], 1, b);
            double d = (b.n>0)? std::sqrt(b.d2[0]) : 1e9;
            if(d <= tol) ++cov;
            rss += d*d; ++rn;
        }
        sp.cover = (double)cov/(double)m;
        sp.fitRms = std::sqrt(rss/std::max(1,rn));
    }
    fprintf(stderr, "[曲面] 行追踪: %d 行 行距=%.2fmm 点距=%.2fmm 覆盖率=%.1f%% 残差RMS=%.2fmm\n",
            sp.nRow, sp.rowGap, pstep, sp.cover*100.0, sp.fitRms);
    return true;
}

// 曲面块重建网格：**独立的高度场结构化网格**，不走 build_plane_mesh。
// 为什么不能复用：build_plane_mesh 末尾有"最小外接矩形吸附"——矩形度 ≥0.85 时
// 直接把整个面替换成 4 个角点，曲面用它做网格会被压成一块平板（WPF 里看不出曲面）。
// 做法：
//   1) 点投到平均平面得 (u,v)；网格分辨率在 ABI 预算内取最细：
//      for K=64..8，cell=max(wU,wV)/K，取第一个满足
//      nu*nv<=maxVerts 且 2(nu-1)(nv-1)<=maxTris 的 K；
//   2) 每个网格点取"最近格"内点云的中位高度（防单点噪声顶尖刺），
//      空格用最多 8 轮邻域平均扩散填；
//   3) 闭运算（膨胀 1 格再腐蚀 1 格）桥接扫描线间隙，否则网格中间漏缝；
//   4) 只对"四角都占用且四角高度有效"的格发 2 个三角（CCW，法向朝外，
//      与平面路径的绕向约定一致）；顶点/三角数在循环里逐项对照 ABI 预算保护。
static bool build_surface_mesh(const SurfPatch& sp, double rowGap, double fillHole,
                               std::vector<Vec3>& verts, std::vector<int>& tris,
                               int maxVerts, int maxTris)
{
    (void)rowGap; (void)fillHole;   // 网格分辨率改由 ABI 预算决定，不再用行距推
    const int m = (int)sp.P.size();
    if(m < 20) return false;
    std::vector<Vec2> uv(m);
    double minU=1e18,maxU=-1e18,minV=1e18,maxV=-1e18;
    for(int i=0;i<m;++i){
        Vec3 d = sp.P[i]-sp.c0;
        uv[i] = Vec2{ dot(d, sp.eL), dot(d, sp.eS) };
        minU=std::min(minU,uv[i].x); maxU=std::max(maxU,uv[i].x);
        minV=std::min(minV,uv[i].y); maxV=std::max(maxV,uv[i].y);
    }
    const double wU=maxU-minU, wV=maxV-minV;
    if(wU<2.0 || wV<2.0) return false;

    // ---- 1) 网格分辨率：ABI 预算内取最细 ----
    int nU=0,nV=0; double cell=0;
    for(int K=64; K>=8; --K){
        double c = std::max(wU,wV)/(double)K;
        int a=(int)std::ceil(wU/c)+1, b=(int)std::ceil(wV/c)+1;
        if((long long)a*b <= (long long)maxVerts &&
           2LL*(a-1)*(b-1) <= (long long)maxTris){ nU=a; nV=b; cell=c; break; }
    }
    if(nU<=0){   // 兜底：预算小到 K=8 都放不下时逐维收缩
        double c = std::max(wU,wV)/8.0;
        nU=(int)std::ceil(wU/c)+1; nV=(int)std::ceil(wV/c)+1;
        while((nU>2 || nV>2) &&
              ((long long)nU*nV > (long long)maxVerts ||
               2LL*(nU-1)*(nV-1) > (long long)maxTris)){
            if(nU>=nV){ if(nU>2) --nU; else --nV; }
            else      { if(nV>2) --nV; else --nU; }
        }
        cell=c;
    }
    // 顶点格均匀铺满包围盒（最后一列/行贴边，避免浮点累加越界）
    const double du = wU/(double)(nU-1), dv = wV/(double)(nV-1);

    // ---- 2) 每个网格点：最近格内点云的中位高度 ----
    const size_t NG = (size_t)nU*(size_t)nV;
    std::vector<std::vector<double>> bucket(NG);
    for(int i=0;i<m;++i){
        int iu=(int)std::lround((uv[i].x-minU)/du);
        int iv2=(int)std::lround((uv[i].y-minV)/dv);
        iu=std::min(nU-1,std::max(0,iu)); iv2=std::min(nV-1,std::max(0,iv2));
        bucket[(size_t)iu*nV+iv2].push_back(dot(sp.P[i]-sp.c0, sp.n0));
    }
    std::vector<double> h(NG, 1e300);
    std::vector<unsigned char> occ(NG, 0);
    for(size_t i=0;i<NG;++i){
        std::vector<double>& bk=bucket[i];
        if(bk.empty()) continue;
        std::nth_element(bk.begin(), bk.begin()+bk.size()/2, bk.end());
        h[i]=bk[bk.size()/2];       // 中位高：单点噪声顶不出尖刺
        occ[i]=1;
    }
    // 空格：最多 8 轮"邻域已填值平均"的扩散填空
    for(int it=0; it<8; ++it){
        bool any=false;
        std::vector<double> nh = h;
        for(int iv2=0;iv2<nV;++iv2) for(int iu=0;iu<nU;++iu){
            size_t k=(size_t)iu*nV+iv2;
            if(h[k] < 1e299) continue;
            double s=0; int c=0;
            for(int a=-1;a<=1;++a) for(int b=-1;b<=1;++b){
                int ju=iu+a, jv=iv2+b;
                if(ju<0||jv<0||ju>=nU||jv>=nV) continue;
                double v=h[(size_t)ju*nV+jv];
                if(v < 1e299){ s+=v; ++c; }
            }
            if(c){ nh[k]=s/c; any=true; }
        }
        h.swap(nh);
        if(!any) break;
    }
    // ---- 3) 闭运算（膨胀 1 格再腐蚀 1 格，8 邻域）：桥接扫描线间隙 ----
    // 网格外一律当作"已占用"参与腐蚀：否则贴着包围盒边界的整圈会被腐蚀掉，
    // 网格比点云范围小一圈（视觉上整个面缩水半格）。
    {
        std::vector<unsigned char> dil(NG,0), ero(NG,0);
        for(int iv2=0;iv2<nV;++iv2) for(int iu=0;iu<nU;++iu){
            if(!occ[(size_t)iu*nV+iv2]) continue;
            for(int a=-1;a<=1;++a) for(int b=-1;b<=1;++b){
                int ju=iu+a, jv=iv2+b;
                if(ju<0||jv<0||ju>=nU||jv>=nV) continue;
                dil[(size_t)ju*nV+jv]=1;
            }
        }
        for(int iv2=0;iv2<nV;++iv2) for(int iu=0;iu<nU;++iu){
            if(!dil[(size_t)iu*nV+iv2]) continue;
            bool full=true;
            for(int a=-1;a<=1 && full;++a) for(int b=-1;b<=1 && full;++b){
                int ju=iu+a, jv=iv2+b;
                if(ju<0||jv<0||ju>=nU||jv>=nV) continue;   // 网格外算已占用
                if(!dil[(size_t)ju*nV+jv]) full=false;
            }
            if(full) ero[(size_t)iu*nV+iv2]=1;
        }
        occ.swap(ero);
    }
    // ---- 3.5) 有曲面模型时：格高换成模型的平滑高度，再在网格自身尺度上补一次高斯 ----
    // 模型格(8.5mm)比显示网格格(18.2mm)细，直接插值仍会留下格间台阶；
    // 在网格自身尺度补一次 σ=1.2 格（≈22mm）的正则化，顶点二面角才真正落下来。
    if(sp.modelOk){
        int rep=0;
        for(int iv2=0;iv2<nV;++iv2) for(int iu=0;iu<nU;++iu){
            size_t k=(size_t)iu*nV+iv2;
            if(!occ[k] || h[k] > 1e299) continue;
            double u = minU + wU*(double)iu/(double)(nU-1);
            double vv= minV + wV*(double)iv2/(double)(nV-1);
            double hq;
            if(surf_model_height(sp, sp.c0 + sp.eL*u + sp.eS*vv + sp.n0*h[k], hq)){
                h[k] = hq; ++rep;
            }
        }
        if(rep > 0) sm_blur_scalar(h, nU, nV, 1.2);
        fprintf(stderr, "[曲面] 网格高度正则化: %d 个格取自曲面模型 + σ=1.2格高斯\n", rep);
    }
    // ---- 4) 三角化：只对"四角都占用且高度有效"的格发 2 个三角 ----
    verts.clear(); tris.clear();
    std::vector<int> vid(NG, -1);           // 顶点惰性创建，按 nU*nV 预分配
    bool budgetHit=false;
    for(int iv2=0; iv2+1<nV && !budgetHit; ++iv2) for(int iu=0; iu+1<nU; ++iu){
        const size_t cc[4] = { (size_t)iu*nV+iv2,       (size_t)(iu+1)*nV+iv2,
                               (size_t)(iu+1)*nV+iv2+1, (size_t)iu*nV+iv2+1 };
        if(!occ[cc[0]]||!occ[cc[1]]||!occ[cc[2]]||!occ[cc[3]]) continue;
        if(h[cc[0]]>1e299||h[cc[1]]>1e299||h[cc[2]]>1e299||h[cc[3]]>1e299) continue;
        int idx4[4];
        for(int q=0;q<4;++q){
            if(vid[cc[q]] < 0){
                if((int)verts.size() >= maxVerts){ budgetHit=true; break; }   // ABI 保护
                double u = minU + wU*(double)(cc[q]/(size_t)nV)/(double)(nU-1);
                double vv= minV + wV*(double)(cc[q]%(size_t)nV)/(double)(nV-1);
                // 顶点高度优先取"重建曲面模型"的平滑高度：模型是格内中位高度
                // 再做归一化高斯正则化得到的，曲率连续；直接用本函数的格内
                // 中位高度会让相邻三角面二面角 p95 到 26°（用户看到的就是
                // "曲面不规则、曲率不一致"）。模型查不到时回退原高度。
                double hh = h[cc[q]];
                verts.push_back(sp.c0 + sp.eL*u + sp.eS*vv + sp.n0*hh);
                vid[cc[q]]=(int)verts.size()-1;
            }
            idx4[q]=vid[cc[q]];
        }
        if(budgetHit) break;
        if((int)(tris.size()/3)+2 > maxTris) break;      // 三角预算保护
        tris.push_back(idx4[0]); tris.push_back(idx4[1]); tris.push_back(idx4[2]);
        tris.push_back(idx4[0]); tris.push_back(idx4[2]); tris.push_back(idx4[3]);
    }
    int occN=0; for(size_t i=0;i<NG;++i) if(occ[i]) ++occN;
    fprintf(stderr, "[曲面] 曲面网格: %dx%d 顶点=%d 三角=%d (cell=%.1fmm 占用格=%d/%d)\n",
            nU, nV, (int)verts.size(), (int)(tris.size()/3), cell, occN, (int)NG);
    return !verts.empty() && !tris.empty();
}

// 多平面件排除（第十四轮）。
// 现象：开口盒子这种"底面 + 4 个侧壁"的纯多面体，被当成曲面块走行追踪，
//       覆盖率只有 57.3%、残差 59.6mm，轨迹完全不对。
// 原因：连通判据是"相邻点法向夹角 < 20°"。棱边两侧的点，KNN 邻域同时包含两个面，
//       法向落在中间过渡带；这一串过渡点把 0° 的底面和 90° 的侧壁一段段接起来，
//       于是 91.2% 的点连成一个"光滑块"，整体平面拟合 p95 自然高达 167.7mm，
//       超过"曲面"门限（max(4*curvatureTol, 2*行距) = 24mm）而被判成曲面。
//
// 判据：真正的曲面没法被少数几个平面**精确**覆盖；多面体可以。
//       用紧公差顺序 RANSAC 剥离，若 <=nmax 个平面就能覆盖 >=coverReq 的点 → 多面体。
// 实测（tol=2mm, nmax=6）：
//       开口盒子  5 个平面覆盖 100.0%（底面 23.2 / 侧壁 22.3+22.3 / 端壁 16.0+16.0）
//       曲型板    6 个平面只覆盖 43.5%（曲面是连续的，每个 2mm 薄层只吃一条带）
// 余量约 46 个百分点，两边都不靠边。
static bool patch_is_polyhedral(const std::vector<Vec3>& P, const std::vector<int>& blk,
                                int nmax, double tolPlane, double coverReq,
                                int& nPlaneOut, double& covOut)
{
    nPlaneOut = 0; covOut = 0.0;
    const int N = (int)blk.size();
    if(N < 200) return false;

    // 抽样到 <=20000 点：判据是"覆盖比例"，抽样不改变结论，但把 RANSAC 成本压到常数级
    const int cap = 20000;
    const int step = std::max(1, N/cap);
    std::vector<Vec3> Q;
    Q.reserve((size_t)(N/step) + 1);
    for(int i=0;i<N;i+=step) Q.push_back(P[blk[i]]);
    const int M = (int)Q.size();
    if(M < 200) return false;

    const int iters = 200;
    const int minPts = std::max(50, (int)(M*0.02));
    std::mt19937 rng(20771);
    std::vector<int> cur(M), rest;
    std::iota(cur.begin(), cur.end(), 0);
    int covered = 0;
    for(int it=0; it<nmax; ++it){
        const int m = (int)cur.size();
        if(m < minPts) break;
        int bestCnt = -1;
        Vec3 bn{0,0,1}; double bd = 0;
        std::uniform_int_distribution<int> ud(0, m-1);
        for(int k=0;k<iters;++k){
            int a=ud(rng), b=ud(rng), c=ud(rng);
            if(a==b || b==c || a==c) continue;
            Vec3 nn = cross(Q[cur[b]]-Q[cur[a]], Q[cur[c]]-Q[cur[a]]);
            double ln = len(nn);
            if(ln < 1e-9) continue;
            nn = nn*(1.0/ln);
            double d = -dot(nn, Q[cur[a]]);
            int cnt = 0;
            for(int t=0;t<m;++t) if(std::fabs(dot(Q[cur[t]], nn)+d) < tolPlane) ++cnt;
            if(cnt > bestCnt){ bestCnt = cnt; bn = nn; bd = d; }
        }
        if(bestCnt < minPts) break;
        rest.clear(); rest.reserve((size_t)(m-bestCnt));
        int got = 0;
        for(int t=0;t<m;++t){
            if(std::fabs(dot(Q[cur[t]], bn)+bd) < tolPlane) ++got;
            else rest.push_back(cur[t]);
        }
        covered += got;
        ++nPlaneOut;
        cur.swap(rest);
        covOut = (double)covered/(double)M;
        if(covOut >= coverReq) return true;
    }
    return false;
}

// 在给定点集（按 blk 索引）上跑一轮 RANSAC，返回最大紧公差平面的内点（全局下标）。
// 与 patch_is_polyhedral 同源，但只取第一个平面 —— 用于"块内大平盘"剥离。
static bool ransac_one_plane(const std::vector<Vec3>& P, const std::vector<int>& blk,
                             double tol, int iters, Vec3& nOut, double& dOut,
                             std::vector<int>& inlOut)
{
    inlOut.clear();
    const int N = (int)blk.size();
    if(N < 100) return false;
    std::mt19937 rng(20771);
    std::uniform_int_distribution<int> ud(0, N-1);
    int bestCnt = 0; Vec3 bn{0,0,1}; double bd = 0;
    for(int k=0;k<iters;++k){
        int a=ud(rng), b=ud(rng), c=ud(rng);
        if(a==b || b==c || a==c) continue;
        Vec3 nn = cross(P[blk[b]]-P[blk[a]], P[blk[c]]-P[blk[a]]);
        double ln = len(nn);
        if(ln < 1e-9) continue;
        nn = nn*(1.0/ln);
        double d = -dot(nn, P[blk[a]]);
        int cnt = 0;
        for(int t=0;t<N;++t) if(std::fabs(dot(P[blk[t]], nn)+d) < tol) ++cnt;
        if(cnt > bestCnt){ bestCnt = cnt; bn = nn; bd = d; }
    }
    if(bestCnt < 100) return false;
    nOut = bn; dOut = bd;
    inlOut.reserve(bestCnt);
    for(int t=0;t<N;++t) if(std::fabs(dot(P[blk[t]], bn)+bd) < tol) inlOut.push_back(blk[t]);
    return true;
}

// 曲面块检测：找"法向局部连续"的最大连通块，再判断它平不平。
//   * 连通条件用**相邻点**法向夹角 < 20°，不是"与种子夹角"——否则弯过 60°
//     的曲面会被种子法向切成几块（这正是原来被硬切成 8 个平面的原因）；
//   * 平整判据用整体平面拟合的 p95 残差。实测 10 个样本里只有曲型板
//     p95=72.9mm，其余最大 16.5mm（complex），所以门限取
//     max(4*curvatureTol, 2*行距) 有 4 倍以上余量，不会误伤其他工件。
static bool detect_curved_patch(const std::vector<Vec3>& P, const CloudNormals& cn,
                               const Vec3& partC, double curvatureTol, double rowSpacing,
                               SurfPatch& sp, std::vector<int>& outInl)
{
    const int n = (int)P.size();
    if(n < 500) return false;
    const int K = CloudNormals::K;

    // ---- 并查集：相邻点法向夹角 < 20° 连通 ----
    std::vector<int> par(n);
    for(int i=0;i<n;++i) par[i]=i;
    // 路径压缩的非递归 find（不用 std::function，14M 次调用差一个间接跳转）
    auto findRoot = [&](int x)->int{
        int r = x;
        while(par[r] != r) r = par[r];
        while(par[x] != r){ int nx = par[x]; par[x] = r; x = nx; }
        return r;
    };
    const double cosGrow = std::cos(20.0*3.141592653589793/180.0);
    for(int i=0;i<n;++i){
        const int* row = &cn.knn[(size_t)i*(size_t)K];
        for(int t=1;t<K;++t){
            int j = row[t]; if(j<0) continue;
            if(std::fabs(dot(cn.nor[i], cn.nor[j])) < cosGrow) continue;
            int a=findRoot(i), b=findRoot(j);
            if(a!=b) par[a]=b;
        }
    }
    std::vector<int> cnt(n,0);
    for(int i=0;i<n;++i) cnt[findRoot(i)]++;
    int bestRoot=-1, bestCnt=0;
    for(int i=0;i<n;++i) if(cnt[i] > bestCnt){ bestCnt=cnt[i]; bestRoot=i; }
    if(bestRoot < 0 || bestCnt < 50) return false;
    double share = (double)bestCnt/(double)n;

    outInl.clear(); outInl.reserve(bestCnt);
    for(int i=0;i<n;++i) if(findRoot(i)==bestRoot) outInl.push_back(i);
    if((int)outInl.size() < 50) return false;

    // ---- 平不平：整体平面拟合的 p95 残差 ----
    Vec3 nf, cf;
    if(!ls_plane(P, outInl, nf, cf)) return false;
    nf = norm(nf);
    const double thr = std::max(4.0*curvatureTol, 2.0*rowSpacing);
    {
        std::vector<double> ds; ds.reserve(outInl.size());
        for(int i : outInl) ds.push_back(std::fabs(dot(P[i]-cf, nf)));
        std::vector<double> t = ds;
        size_t k = (size_t)((t.size()-1)*0.95);
        std::nth_element(t.begin(), t.begin()+k, t.end());
        double p95 = t[k];
        fprintf(stderr, "[曲面] 最大光滑块: %d/%d 点(%.1f%%) 平面残差 p95=%.1fmm (门限=%.1fmm)\n",
                bestCnt, n, share*100.0, p95, thr);
        if(share < 0.5 || p95 < thr) return false;
    }

    // ---- 多平面件排除：能拆成几个精确平面 → 不是曲面，退回平面分割 ----
    {
        // 公差取 0.2×行距（11.76mm → 2.4mm）：既容得下现场点云的墙面噪声（~0.5mm），
        // 又远小于"浅曲面被切成薄带"的厚度，避免把 sag 30mm 级的缓曲面误判成多面体
        // （sag 30mm 的板在 2.4mm 带宽下要 12 个平面才铺满，nmax=8 只够盖 ~75%）。
        const double tolPl = std::min(3.0, std::max(1.5, 0.20*rowSpacing));
        int np = 0; double cov = 0.0;
        bool poly = patch_is_polyhedral(P, outInl, 8, tolPl, 0.93, np, cov);
        fprintf(stderr, "[曲面] 多平面校验: 紧公差(%.1fmm)剥离 %d 个平面覆盖 %.1f%%%s\n",
                tolPl, np, cov*100.0,
                poly ? " → 多面体，退回平面分割" : "（未达 93%，仍按曲面处理）");
        if(poly) return false;
    }

    // ---- 剥离块内的"平行于 n0 的大平盘"（底面 / 盖板 / 环形唇边）----
    // 曲面块靠"相邻点法向差 < 20°"连通，棱边处的过渡点会把平盘和曲面串成一块。
    // 平盘对曲面行追踪是致命的：沿 n0 的高度场在盘上是常数，等值线无定义，
    // 追踪会退化成"沿一条直线横穿整个曲面"。
    // 实测垃圾桶（顶部环形唇边 + 锥台 + 底面圆盘）三者连成 100% 的块，27 行
    // 全部从上沿直插到底面中心，母线间距外圈 20.4mm / 内圈 7.0mm，端点在轴心抱团。
    //
    // 找法**不能用 RANSAC**：锥面本身在 2.4mm 公差下就能凑出"伪平面窄带"
    // （实测 3043 点、法向 (0.85,-0.25,0.46)，把 2726 点的底面比下去）。
    // 改用**沿 nf 的 1mm 高度直方图**：平盘在高度场上必然是一个尖峰 ——
    //   垃圾桶：底面 2700 点/箱 vs 锥面 ~36 点/箱；曲型板：~150 点/箱、占 1.7%。
    // 两条判据同时满足才剥：
    //   ① 峰值箱 ±tol 的点占块 ≥ 5%；
    //   ② 这些点在垂直于 nf 的平面上最小宽度 ≥ 3×行距 —— 排除弯曲板被 2.4mm
    //      薄层切出的窄条（曲型板实测宽度只有 4~5mm）。
    // 循环最多 3 轮：一个件可能同时有顶部唇边和底面两个平盘。
    for(int pass=0; pass<3; ++pass){
        const double tolPl = std::min(3.0, std::max(1.5, 0.20*rowSpacing));
        {   // 每轮用"剩下这些点"重拟合法向，第二个盘才不会被第一个盘带偏
            Vec3 nf2, cf2;
            if(outInl.size() >= 200 && ls_plane(P, outInl, nf2, cf2)){
                nf2 = norm(nf2);
                if(dot(nf2, nf) < 0) nf2 = nf2*(-1.0);
                nf = nf2;
            }
        }
        double hlo=1e18, hhi=-1e18;
        for(int i : outInl){ double h=dot(P[i],nf); hlo=std::min(hlo,h); hhi=std::max(hhi,h); }
        if(hhi-hlo < 2.0*tolPl) break;
        const int nb = (int)std::floor(hhi-hlo)+1;
        std::vector<int> hc(nb, 0);
        for(int i : outInl){
            int b = (int)std::floor(dot(P[i],nf)-hlo);
            if(b<0) b=0; if(b>=nb) b=nb-1;
            hc[b]++;
        }
        int bb=-1, bc=0;
        for(int b=0;b<nb;++b) if(hc[b]>bc){ bc=hc[b]; bb=b; }
        if(bb<0 || bc < 100) break;
        const double hc0 = hlo + (double)bb + 0.5;
        std::vector<int> pin; pin.reserve(outInl.size());
        for(int i : outInl) if(std::fabs(dot(P[i],nf)-hc0) <= tolPl) pin.push_back(i);
        if(pin.size() < 100) break;
        const double share2 = (double)pin.size()/(double)outInl.size();
        // 平面内 2D 最小宽度（4σ）
        double wmin = 0.0;
        {
            Vec3 u = cross(nf, Vec3{0,0,1});
            if(len(u) < 0.5) u = cross(nf, Vec3{1,0,0});
            u = norm(u);
            Vec3 v = norm(cross(nf, u));
            double mu=0, mv=0;
            for(int i : pin){ mu += dot(P[i],u); mv += dot(P[i],v); }
            mu /= (double)pin.size(); mv /= (double)pin.size();
            double suu=0,svv=0,suv=0;
            for(int i : pin){
                double a = dot(P[i],u)-mu, b = dot(P[i],v)-mv;
                suu += a*a; svv += b*b; suv += a*b;
            }
            suu /= (double)pin.size(); svv /= (double)pin.size(); suv /= (double)pin.size();
            double tr = suu+svv, det = suu*svv - suv*suv;
            double lm = tr*0.5 - std::sqrt(std::max(0.0, tr*tr*0.25 - det));
            wmin = 4.0*std::sqrt(std::max(0.0, lm));
        }
        const bool big  = (share2 >= 0.05);
        const bool wide = (wmin >= 3.0*rowSpacing);
        if(big && wide){
            fprintf(stderr, "[曲面] 剥离块内平盘[%d]: %zu 点(%.1f%%) 沿nf高度=%.1fmm 宽%.1fmm\n",
                    pass, pin.size(), share2*100.0, hc0, wmin);
            std::vector<char> drop(n, 0);
            for(int i : pin) drop[i] = 1;
            std::vector<int> keep; keep.reserve(outInl.size());
            for(int i : outInl) if(!drop[i]) keep.push_back(i);
            if((int)keep.size() < 200) break;
            outInl.swap(keep);
            continue;                    // 可能还有第二个平盘
        }
        if(pass == 0)
            fprintf(stderr, "[曲面] 块内平盘检查: 最大高度峰 %zu 点(%.1f%%) 宽%.1fmm → 不剥离\n",
                    pin.size(), share2*100.0, wmin);
        break;
    }

    // ---- 组装 SurfPatch ----
    sp.P.clear(); sp.N.clear();
    sp.P.reserve(outInl.size()); sp.N.reserve(outInl.size());
    Vec3 nsum{0,0,0};
    for(int i : outInl){
        sp.P.push_back(P[i]);
        Vec3 ni = norm(cn.nor[i]);
        if(len(ni) < 0.5) ni = nf;
        sp.N.push_back(ni);
        nsum = nsum + ni;
    }
    Vec3 n0 = norm(nsum);
    if(len(n0) < 0.5) n0 = nf;
    // 定向：和平面路径一致 —— 朝 +Z 半球；近竖直面则背离工件质心
    if(std::fabs(n0.z) > 0.05){
        if(n0.z < 0) n0 = n0*(-1.0);
    } else if(dot(n0, cf-partC) < 0){
        n0 = n0*(-1.0);
    }
    for(auto& v : sp.N) if(dot(v, n0) < 0) v = v*(-1.0);
    sp.n0 = n0;
    sp.c0 = cf;
    sp.eL = norm(cross(Vec3{0,0,1}, n0));
    if(len(sp.eL) < 1e-6) sp.eL = Vec3{1,0,0};
    sp.eS = norm(cross(n0, sp.eL));

    sp.tree.init(sp.P);

    // 点间距：抽 200 个点取中位最近邻距离
    {
        std::vector<double> ds;
        int step = std::max(1, (int)sp.P.size()/200);
        for(size_t i=0;i<sp.P.size(); i+=step){
            KnnBuf b; sp.tree.knn(sp.P[i], 2, b);
            if(b.n>=2) ds.push_back(std::sqrt(b.d2[0] > b.d2[1] ? b.d2[0] : b.d2[1]));
        }
        // 注：KnnBuf 是最大堆，d2[0] 是"最差"的那个，所以取两者较大值即第 2 近
        if(ds.empty()) return false;
        size_t k = ds.size()/2;
        std::nth_element(ds.begin(), ds.begin()+k, ds.end());
        sp.spacing = ds[k];
    }
    sp.rFit = std::min(30.0, std::max(3.0, 4.0*sp.spacing));
    fprintf(stderr, "[曲面] 判定为曲面块: 点间距=%.2fmm 拟合半径=%.1fmm 平均法向=(%.3f,%.3f,%.3f)\n",
            sp.spacing, sp.rFit, n0.x, n0.y, n0.z);
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
    p.flatnessTol    = 40.0;
    p.curvatureTol   = 5.0;
    p.sorK           = 20;
    p.sorZ           = 2.5;
    p.normalFlip     = 0;
    p.fillHoleSize   = 300.0;
    p.scanDirMode    = 0;
    return p;
}

SPRAYTRAJ_API SprayTrajResult spray_traj_generate(
    const double* xyz, int n, const SprayTrajParams* params)
{
    SprayTrajResult res;
    std::memset(&res, 0, sizeof(res));
    res.error = nullptr;
    ProfStage psTotal("0.总耗时");   // 任意返回路径都会打印

    std::string err;
    if(!xyz || n<10){ err="点云点数不足(<10)"; res.error=_strdup(err.c_str()); return res; }

    std::vector<Vec3> P(n);
    for(int i=0;i<n;++i) P[i]={xyz[i*3],xyz[i*3+1],xyz[i*3+2]};
    res.input_points=n;

    double tSor=0, tSplit=0, tMesh=0, tTraj=0;
    // 1. 去噪
    std::vector<Vec3> Pf;
    {
        ProfStage ps("1.SOR去噪", &tSor);
        Pf = sor_filter(P, params->sorK, params->sorZ);
    }
    res.removed_points = n - (int)Pf.size();
    if(Pf.size()<10){ err="去噪后点数不足"; res.error=_strdup(err.c_str()); return res; }
    {
        double zmin=Pf[0].z, zmax=Pf[0].z;
        for(auto&p:Pf){ zmin=std::min(zmin,p.z); zmax=std::max(zmax,p.z); }
        fprintf(stderr, "[去噪] 点数 %zu->%zu, Z范围[%.1f,%.1f]\n", P.size(), Pf.size(), zmin, zmax);
    }

    // 2. 多平面分割重建：
    //   顺序 RANSAC 反复剥离最大平面。折型板 → 两个面；槽型盘 → 底盘 + 凸起侧边；
    //   高度差小的浅曲面/缓起伏会被第一个主平面吸收（flatnessTol 判定）。
    //   每个平面后续分别拟合、生成轨迹、沿各自法向抬升（避免撞击凸起）。
    //
    //   曲面件例外：先试曲面块判定。命中则**最大光滑块整体走曲面行追踪**，
    //   剩下的点（凸台/翻边/耳朵等）仍交给 split_planes 按平面处理。
    //   没命中时 split_planes 的输入与以前完全一致，输出逐字节不变。
    std::vector<PlaneSeg> segs;
    SurfPatch surf;
    bool isCurved = false;
    {
        ProfStage ps("2.平面分割", &tSplit);
        CloudNormals cn;
        cn.build(Pf);

        Vec3 partC{0,0,0};
        for(const auto& p : Pf) partC = partC + p;
        partC = partC*(1.0/(double)std::max<size_t>(1,Pf.size()));

        std::vector<int> curvedInl;
        isCurved = detect_curved_patch(Pf, cn, partC, params->curvatureTol,
                                       params->scanRowSpacing, surf, curvedInl);
        if(isCurved){
            // 剩余点单独分割（索引要映射回 Pf）
            std::vector<char> isCurvedPt(Pf.size(), 0);
            for(int i : curvedInl) isCurvedPt[i] = 1;
            std::vector<Vec3> Pres; std::vector<int> back;
            Pres.reserve(Pf.size() - curvedInl.size());
            back.reserve(Pf.size() - curvedInl.size());
            for(size_t i=0;i<Pf.size();++i) if(!isCurvedPt[i]){ Pres.push_back(Pf[i]); back.push_back((int)i); }
            segs.reserve(8);
            if(Pres.size() >= 100){
                CloudNormals cnRes;
                cnRes.build(Pres);
                auto sub = split_planes(Pres, cnRes, params->flatnessTol, params->curvatureTol,
                                        std::max(20, (int)(Pf.size()*0.02)));
                for(auto& s : sub){
                    for(int& k : s.inl) k = back[k];
                    segs.push_back(std::move(s));
                }
            }
            // 曲面块排在最前：它覆盖绝大部分点，先喷
            PlaneSeg cs;
            double rss=0;
            for(int i : curvedInl){ double e=dot(Pf[i]-surf.c0, surf.n0); rss+=e*e; }
            cs.pl.n = surf.n0;
            cs.pl.d = -dot(surf.n0, surf.c0);
            cs.pl.c = surf.c0;
            cs.inl = curvedInl;
            cs.curved = true;
            fprintf(stderr, "[曲面] 剩余 %zu 点走平面分割，得到 %zu 个面\n", Pres.size(), segs.size());
            segs.insert(segs.begin(), std::move(cs));
            (void)rss;
        } else {
            segs = split_planes(Pf, cn, params->flatnessTol, params->curvatureTol, std::max(20, (int)(Pf.size()*0.02)));
        }
    }

    // 若分割结果为空（极少见），退回单平面
    if(segs.empty()){
        Plane pl; std::vector<int> inl;
        if(!fit_plane(Pf, params->flatnessTol, pl, inl)){ err="无法重建有效平面"; res.error=_strdup(err.c_str()); return res; }
        PlaneSeg seg; seg.pl=pl; seg.inl=inl; segs.push_back(std::move(seg));
    }

    // 过滤零碎小平面：主面之外的面点数占比 < 3% 视为碎面丢弃。
    // 否则顺序 RANSAC 会切出几百点的碎面，各自生成零散轨迹，
    // 视觉上"乱七八糟"。折型板翻边/槽型盘凸边/深槽侧壁都是占比>3%的实体结构，
    // 不会误删；几百点的小碎面（噪声、槽底碎片）占比<3%被正确丢弃。
    //
    // ⚠ 但"点数占比"单独不足以判碎面：它会把【朝向独立的真实结构】误杀。
    // 实测（cao 槽型盘）：66.1° 的侧壁 578 点，占比 2.75% < 3%，被当成碎面删掉；
    // 而 flat=30 时同一个侧壁因残点归并多吸了点（750 点）就活了 —— 表现为
    // "阈值放大反而面更多"，与用户"阈值越大越能合并"的语义正好相反。
    // 判据修正：与主面夹角 > 40°（陡壁/翻边/裙边）视为独立结构，只要求绝对点数下限，
    // 不参与占比竞赛。近共面（<40°）的碎片仍然按 3% 丢弃——那才是区域生长产生的碎屑。
    if(segs.size() > 1){
        size_t maxIdx=0; for(size_t i=1;i<segs.size();++i)
            if(segs[i].inl.size() > segs[maxIdx].inl.size()) maxIdx=i;
        double total = (double)Pf.size();
        const Plane& mainPl = segs[maxIdx].pl;
        const size_t absFloor = std::max<size_t>(80, (size_t)(total*0.005));
        std::vector<PlaneSeg> keep;
        for(auto& s : segs){
            double ratio = (double)s.inl.size()/total;
            double angMain = std::acos(std::min(1.0, std::fabs(dot(s.pl.n, mainPl.n))))*180.0/3.141592653589793;
            bool big      = ratio >= 0.03;
            bool distinct = (angMain > 40.0) && (s.inl.size() >= absFloor);
            if(big || distinct){
                if(!big) fprintf(stderr, "[split] 保留独立朝向小面: 点数=%zu (占比%.1f%%, 与主面夹角%.1f° > 40°)\n",
                                 s.inl.size(), ratio*100.0, angMain);
                keep.push_back(std::move(s));
            }
            else fprintf(stderr, "[split] 丢弃碎面: 点数=%zu (占比%.1f%% < 3%%, 与主面夹角%.1f°)\n",
                          s.inl.size(), ratio*100.0, angMain);
        }
        if(!keep.empty()) segs = std::move(keep);
    }

    // 法向定向：统一保证各面法向朝 +Z 半球（俯视相机外表面朝上），
    // 每个平面分别判向（折型板的两面法向可能相反，需各自翻转）。
    // 近竖直面（侧壁）另按"背离工件质心"定向：PCA 拟合法向符号是随机的，
    // 一旦朝内，抬升 standoff 会把喷头送进工件内部。
    Vec3 partC{0,0,0};
    for(const auto& p : Pf) partC = partC + p;
    partC = partC*(1.0/(double)std::max<size_t>(1,Pf.size()));
    for(auto& seg : segs){
        Vec3& n = seg.pl.n;
        bool flip = false;
        if(params->normalFlip){ flip = true; }
        else if(std::fabs(n.z) > 0.05){ if(n.z < 0) flip = true; }
        else if(dot(n, seg.pl.c - partC) < 0){ flip = true; }   // 近竖直：背离质心
        if(!flip) continue;
        n = n*(-1.0);
        seg.pl.d = -seg.pl.d;
        // 曲面块：法向场要一起翻，否则计算抬升时喷头会钻到工件里
        if(seg.curved && isCurved){
            surf.n0 = surf.n0*(-1.0);
            for(auto& v : surf.N) v = v*(-1.0);
        }
    }

    // 报告：用面积/点数最大的主平面作为"重建面"统计信息
    // 曲面块：平面残差天然很大（它不是平面），但曲面贴合残差要等行追踪跑完才知道，
    // 所以这里先记下标记，函数末尾再用 surf.fitRms 覆盖（否则对外报的是 57mm 的假值）。
    bool reportCurved = false;
    {
        size_t bi=0; for(size_t i=1;i<segs.size();++i) if(segs[i].inl.size()>segs[bi].inl.size()) bi=i;
        Plane& pl = segs[bi].pl;
        const auto& inl = segs[bi].inl;
        double rss=0; for(int i:inl){ double e=dot(Pf[i],pl.n)+pl.d; rss+=e*e; }
        res.plane_rms=std::sqrt(rss/std::max(1,(int)inl.size()));
        res.normal_x=pl.n.x; res.normal_y=pl.n.y; res.normal_z=pl.n.z;
        res.plane_cx=pl.c.x; res.plane_cy=pl.c.y; res.plane_cz=pl.c.z;
        res.inlier_points=(int)inl.size();
        reportCurved = (segs[bi].curved && isCurved);
    }

    // 输出所有分割平面（供 WPF 显示重建面）
    // 为每个面计算局部基(u,v)主轴 + 内点包围盒，WPF 直接画矩形（不做 PCA 近似）
    res.plane_count = 0;
    size_t wi = 0;                     // 只统计"网格重建成功"的面，避免空面被前端画成包围盒
    std::vector<char> meshOk(segs.size(), 0);   // 没重建出面片的面不生成轨迹
    // 窄条补全记录：补全面要让轨迹也铺满整个矩形，不能只跟着点云走
    std::vector<char> planeFilled(segs.size(), 0);
    std::vector<std::array<Vec3,4>> planeFillRect(segs.size());
    for(size_t i=0;i<segs.size() && wi<8;++i){
        Plane& P = segs[i].pl;
        res.planes[wi].nx = P.n.x;
        res.planes[wi].ny = P.n.y;
        res.planes[wi].nz = P.n.z;
        res.planes[wi].cx = P.c.x;
        res.planes[wi].cy = P.c.y;
        res.planes[wi].cz = P.c.z;
        res.planes[wi].count = (int)segs[i].inl.size();

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

        res.planes[wi].ux=uf.x; res.planes[wi].uy=uf.y; res.planes[wi].uz=uf.z;
        res.planes[wi].vx=vf.x; res.planes[wi].vy=vf.y; res.planes[wi].vz=vf.z;
        res.planes[wi].minU=minU; res.planes[wi].maxU=maxU;
        res.planes[wi].minV=minV; res.planes[wi].maxV=maxV;

        // ---- 规则化重建面片：轮廓 → 直边简化 → 矩形吸附 → 三角化 ----
        // 旧做法是"逐 cell 打小方片"，边界必然是阶梯锯齿（面看起来像拼图）。
        // 这里改成矢量化轮廓：矩形面直接输出 4 个角点，孔洞用直边多边形/矩形
        // 表达，三角化后交给 WPF 渲染，重建出来就是干净的规则平面。
        if(segs[i].curved && isCurved){
            // 曲面块：独立的高度场结构化网格（不复用平面那套矢量化轮廓，
            // 矩形吸附会把曲面压成 4 角平板），顶点按格内中位高抬回曲面，
            // 所以 WPF 看到的是贴着工件的曲面片。
            std::vector<Vec3> verts3; std::vector<int> tris;
            bool okMesh = false;
            {
                ProfStage ps("3.曲面重建", &tMesh);
                if(!surf.modelOk) build_surf_model(surf, 0.0);
                okMesh = build_surface_mesh(surf, params->scanRowSpacing, params->fillHoleSize,
                                            verts3, tris,
                                            SPRAYTRAJ_PLANE_MAX_POLY_VERTS,
                                            SPRAYTRAJ_PLANE_MAX_POLY_TRIS);
            }
            res.planes[wi].polyVertCount = 0;
            res.planes[wi].polyTriCount  = 0;
            if(okMesh){
                int nv = std::min((int)verts3.size(), SPRAYTRAJ_PLANE_MAX_POLY_VERTS);
                for(int v=0; v<nv; ++v){
                    res.planes[wi].polyVerts[v*3+0]=verts3[v].x;
                    res.planes[wi].polyVerts[v*3+1]=verts3[v].y;
                    res.planes[wi].polyVerts[v*3+2]=verts3[v].z;
                }
                int nt = std::min((int)(tris.size()/3), SPRAYTRAJ_PLANE_MAX_POLY_TRIS);
                for(int t=0; t<nt; ++t){
                    res.planes[wi].polyTris[t*3+0]=tris[t*3+0];
                    res.planes[wi].polyTris[t*3+1]=tris[t*3+1];
                    res.planes[wi].polyTris[t*3+2]=tris[t*3+2];
                }
                res.planes[wi].polyVertCount=nv;
                res.planes[wi].polyTriCount=nt;
                fprintf(stderr, "[重建] 面片%zu: 顶点=%d 三角=%d（曲面块）内点=%zu\n",
                        wi, nv, nt, segs[i].inl.size());
                meshOk[i] = 1;
                ++wi;
            } else {
                fprintf(stderr, "[重建] 面片%zu: 曲面网格重建失败，不输出该面\n", i);
            }
            continue;
        }
        {
            std::vector<Vec2> uv; uv.reserve(segs[i].inl.size());
            for(int k : segs[i].inl){
                Vec3 d = Pf[k]-P.c;
                uv.push_back(Vec2{ dot(d,uf), dot(d,vf) });
            }
            std::vector<Vec2> verts; std::vector<int> tris;
            // 调试：SPRAYTRAJ_DUMP=1 时导出该面内点的局部 (u,v) 坐标，供离线复现
            if(const char* dump = std::getenv("SPRAYTRAJ_DUMP")){
                if(dump[0]=='1'){
                    char nm[128]; snprintf(nm,sizeof(nm),"uv_plane%d.csv",(int)i);
                    if(FILE* fu=fopen(nm,"w")){
                        for(auto& q : uv) fprintf(fu,"%.3f,%.3f\n", q.x, q.y);
                        fclose(fu);
                    }
                }
            }
            std::vector<char> keepMask;
            std::vector<Vec2> fillRect;
            bool okMesh = false;
            {
                ProfStage ps("3.轮廓重建", &tMesh);
                // 环形窄面：重建格距必须能分辨环宽。默认格距 = 行距，行距调到 30mm
                // 时格距 30 > 环宽 32 的一半，开运算把整圈削成几段弧，内点被剔掉
                // 一半（6977→3653），轨迹只能盖住 250° 的环。
                double cellMesh = params->scanRowSpacing;
                {
                    Vec2 rc; double rin=0, rout=0, rw=0; int ns=0;
                    if(detect_ring_uv(uv, params->scanRowSpacing, rc, rin, rout, &rw, &ns)){
                        cellMesh = std::min(params->scanRowSpacing,
                                            std::max(1.5, rw/6.0));
                        fprintf(stderr, "[重建] 面片%zu 环形(内径%.1f 外径%.1f 宽%.1f 方位%d/36): 重建格距 %.1fmm (原 %.1f)\n",
                                i, rin, rout, rw, ns, cellMesh, params->scanRowSpacing);
                    }
                }
                okMesh = build_plane_mesh(uv, cellMesh, params->fillHoleSize,
                                          verts, tris,
                                          SPRAYTRAJ_PLANE_MAX_POLY_VERTS,
                                          SPRAYTRAJ_PLANE_MAX_POLY_TRIS,
                                          &keepMask, &fillRect);
            }
            // 补全面：记下矩形四角的世界坐标，轨迹阶段据此在补全区布行
            if(!fillRect.empty() && fillRect.size() >= 4){
                planeFilled[i] = 1;
                for(int k=0;k<4;++k)
                    planeFillRect[i][k] = P.c + uf*fillRect[k].x + vf*fillRect[k].y;
            }
            res.planes[wi].polyVertCount = 0;
            res.planes[wi].polyTriCount  = 0;
            // 用重建面掩膜回筛内点：稀疏混入的孤立扫描线已被开运算剔除，
            // 这些点不能再参与后面按面生成轨迹，否则轨迹会拉到工件外面。
            if(okMesh && keepMask.size()==segs[i].inl.size()){
                std::vector<int> kept; kept.reserve(segs[i].inl.size());
                for(size_t t=0;t<keepMask.size();++t) if(keepMask[t]) kept.push_back(segs[i].inl[t]);
                if(kept.size() >= segs[i].inl.size()/2 && kept.size() >= 20){
                    if(kept.size() != segs[i].inl.size())
                        fprintf(stderr, "[重建] 面片%zu: 剔除离群内点 %zu -> %zu\n",
                                i, segs[i].inl.size(), kept.size());
                    segs[i].inl.swap(kept);
                    res.planes[wi].count = (int)segs[i].inl.size();
                }
            }
            if(okMesh){
                int nv=(int)verts.size();
                if(nv > SPRAYTRAJ_PLANE_MAX_POLY_VERTS) nv = SPRAYTRAJ_PLANE_MAX_POLY_VERTS;
                for(int v=0; v<nv; ++v){
                    Vec3 pos = P.c + uf*verts[v].x + vf*verts[v].y;
                    res.planes[wi].polyVerts[v*3+0]=pos.x;
                    res.planes[wi].polyVerts[v*3+1]=pos.y;
                    res.planes[wi].polyVerts[v*3+2]=pos.z;
                }
                int nt=(int)(tris.size()/3);
                if(nt > SPRAYTRAJ_PLANE_MAX_POLY_TRIS) nt = SPRAYTRAJ_PLANE_MAX_POLY_TRIS;
                for(int t=0; t<nt; ++t){
                    res.planes[wi].polyTris[t*3+0]=tris[t*3+0];
                    res.planes[wi].polyTris[t*3+1]=tris[t*3+1];
                    res.planes[wi].polyTris[t*3+2]=tris[t*3+2];
                }
                res.planes[wi].polyVertCount=nv;
                res.planes[wi].polyTriCount=nt;
                {
                    int nc=0; double lf=0, ar=0;
                    mesh_connect_stats(verts, tris, nc, lf, ar);
                    fprintf(stderr, "[重建] 面片%zu: 顶点=%d 三角=%d（规则轮廓）连通块=%d 最大块占比=%.0f%% 面积=%.0fmm2 内点=%zu\n",
                            wi, nv, nt, nc, lf*100.0, ar, segs[i].inl.size());
                }
                // 诊断：非主面内点到主面平面的垂直距离分布。
                // 极端案例——折型板左侧那片掠射点被拟合成 63° 的"面"，
                // 但每个点到主板面的距离都在 flatnessTol 之内，属于同一块板的起伏。
                if(i>0 && !segs[0].inl.empty()){
                    std::vector<double> ds; ds.reserve(segs[i].inl.size());
                    for(int k : segs[i].inl)
                        ds.push_back(std::fabs(dot(Pf[k], segs[0].pl.n) + segs[0].pl.d));
                    std::vector<double> t=ds;
                    std::nth_element(t.begin(), t.begin()+t.size()/2, t.end());
                    double med=t[t.size()/2];
                    std::nth_element(t.begin(), t.begin()+t.size()*95/100, t.end());
                    double p95=t[t.size()*95/100];
                    double mx=*std::max_element(ds.begin(),ds.end());
                    fprintf(stderr, "[diag] 面%zu 到主面距离: 中位=%.1f p95=%.1f max=%.1f mm (flatTol=%.0f)\n",
                            i, med, p95, mx, params->flatnessTol);
                }
                meshOk[i] = 1;
                ++wi;
            } else {
                fprintf(stderr, "[重建] 面片%zu: 轮廓提取失败（点太少/全散点），不输出该面\n", i);
            }
        }
    }
    res.plane_count = (int)wi;

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
    auto genRows=[&](const std::vector<Vec2>& pts, double ang, double L, double fillHole)->std::vector<Row>{
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
        // 环形/框形面（如槽型盘凸边环）→ 每行分成左右两段，中央空区不喷。
        // 点云缺陷覆盖重建：间隙 < fillHole 的小洞/缺口视为可填补——直接
        // 跨越生成连续行（重建面已把该洞补齐，行应贴重建面连续无断点）；
        // 间隙 >= fillHole 的真实孔/槽才断开为两段。
        const double gapSeg = std::max(fillHole, L*0.5);   // x 间隙超过此值视为两段；以下则跨越填补
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
        // 覆盖重建：整条行带内完全没有点的区段（凹陷/孔洞）在 fillHole 之内时，
        // 按上下两行线性插值补行——重建面已经把这块补齐，轨迹也应该连续覆盖，
        // 否则中间会空一块、看起来"缺行、不规整"。
        if(fillHole > L*1.5 && raw.size() >= 2){
            std::vector<Row> filled; filled.reserve(raw.size()*2);
            for(size_t i=0;i<raw.size();++i){
                filled.push_back(raw[i]);
                if(i+1>=raw.size() || raw[i+1].segIdx != raw[i].segIdx) continue;
                double dy = raw[i+1].yc - raw[i].yc;
                if(dy <= L*1.5 || dy > fillHole) continue;
                int k = (int)std::floor(dy/L);
                for(int t=1;t<k;++t){
                    double f = (double)t*L/dy;
                    Row r;
                    r.yc = raw[i].yc + dy*f;
                    r.xs = raw[i].xs + (raw[i+1].xs-raw[i].xs)*f;
                    r.xe = raw[i].xe + (raw[i+1].xe-raw[i].xe)*f;
                    r.segIdx = raw[i].segIdx;
                    if(r.xe - r.xs >= L) filled.push_back(r);
                }
            }
            raw.swap(filled);
        }
        rs = std::move(raw);
        if(rs.empty()) return rs;
        return rs;
    };

    std::vector<SprayTrajPoint> allPts;
    int totalRows = 0;
    ProfStage psTraj("4.轨迹生成", &tTraj);

    // ---- 面间过渡的安全高度 ----
    // 面间转移是空移，但控制器会把相邻两个点用直线插值走完。原来只补一个
    // spray=0 哨兵点（位置不变），下一点直接是另一面的首点 —— 实际运动就是一
    // 根穿过工件凸台/槽壁的斜线。工业上必须"抬枪绕行"，否则撞件。
    // 处理：沿上一面法向退出 clear → 抬到全局安全高度 zSafe → 平面内平移 → 下降接下一面首点。
    //   clear 取 standoff：喷头本来就离面 100mm，再退 100mm 即彻底离开该面。
    //   zSafe 取工件最高点 + 2×standoff：最高面上的工作点正好等于 zSafe，更低的面则向上抬。
    double zPartMax = -1e18;
    for(const auto& p : Pf){ zPartMax = std::max(zPartMax, p.z); }
    const double clear = params->standoff;
    const double zSafe = zPartMax + 2.0*params->standoff;
    Vec3 prevNormal{0,0,1};
    int transCount = 0;

    for(size_t si=0; si<segs.size(); ++si){
        if(!meshOk[si]) continue;          // 轮廓都重建不出来的面，不喷
        auto& seg = segs[si];
        Plane& pl = seg.pl;
        const auto& inl = seg.inl;
        if(inl.size() < 20) continue;

        // ================= 曲面块：曲面上的行追踪 =================
        // 不走平面的 2D 布行——曲面投影成 2D 会在陡峭段挤压，行距失控。
        // 行点/法向已由 gen_surface_rows 在曲面上直接算好，这里只做抬升与串接。
        if(seg.curved && isCurved){
            if(surf.nRow <= 0){
                ProfStage ps("4b.曲面行追踪", &tTraj);
                if(!gen_surface_rows(surf, L, pstep, params->fillHoleSize,
                                     params->scanDirMode, params->speed)){
                    fprintf(stderr, "[曲面] 行追踪失败，跳过该块\n");
                    continue;
                }
            }
            const int nr = surf.nRow;
            if(nr < 1 || surf.rowBeg.size() < (size_t)nr+1) continue;

            auto emitPt=[&](const Vec3& X, const Vec3& nrm, int spray){
                Vec3 pos = X + nrm*params->standoff;
                allPts.push_back({pos.x,pos.y,pos.z, -nrm.x,-nrm.y,-nrm.z, params->speed, spray});
            };
            // 行间连接：两行端点相距约一个行距，按点距插值，短距保持喷涂不关枪
            auto emitGap=[&](const Vec3& a, const Vec3& na, const Vec3& b, const Vec3& nb){
                double d = len(b-a);
                if(d <= pstep*0.5) return;
                int np = std::max(2, (int)std::ceil(d/pstep));
                int spray = (d <= 3.0*L) ? 1 : 0;
                for(int t=1;t<np;++t){
                    double f = (double)t/(double)np;
                    emitPt(a + (b-a)*f, na + (nb-na)*f, spray);
                }
            };

            // 面间过渡：从上一面终点抬枪绕行到曲面块首点（与平面路径同一套）
            {
                Vec3 firstX = surf.rowP[surf.rowBeg[0]];
                Vec3 firstN = surf.rowN[surf.rowBeg[0]];
                Vec3 firstPos = firstX + firstN*params->standoff;
                if(!allPts.empty()){
                    Vec3 pEnd{allPts.back().x, allPts.back().y, allPts.back().z};
                    Vec3 L1 = pEnd + prevNormal*clear;
                    allPts.push_back({L1.x,L1.y,L1.z, prevNormal.x,prevNormal.y,prevNormal.z, params->speed, 0});
                    allPts.push_back({L1.x,L1.y,zSafe, prevNormal.x,prevNormal.y,prevNormal.z, params->speed, 0});
                    allPts.push_back({firstPos.x,firstPos.y,zSafe, -firstN.x,-firstN.y,-firstN.z, params->speed, 0});
                    ++transCount;
                }
                for(int r=0; r<nr; ++r){
                    int b = surf.rowBeg[r], e = surf.rowBeg[r+1];
                    if(e-b < 2) continue;
                    if(r > 0){
                        emitGap(surf.rowP[b-1], surf.rowN[b-1], surf.rowP[b], surf.rowN[b]);
                    }
                    for(int k=b;k<e;++k) emitPt(surf.rowP[k], surf.rowN[k], 1);
                    totalRows++;
                }
                prevNormal = surf.rowN[surf.rowBeg[nr]-1];
            }
            {
                double angN = std::acos(std::min(1.0,std::fabs(dot(surf.n0,Vec3{0,0,1}))))*180.0/3.141592653589793;
                fprintf(stderr, "[traj] 面%zu(曲面): 内点=%zu 行数=%d 与Z夹角=%.1f° 行距=%.2fmm 覆盖率=%.1f%%\n",
                        si, inl.size(), nr, angN, surf.rowGap, surf.cover*100.0);
            }
            continue;
        }

        // 该平面的局部基
        Vec3 zaxis{0,0,1};
        Vec3 u = norm(cross(zaxis, pl.n));
        if(len(u)<1e-6) u = Vec3{1,0,0};
        Vec3 v = cross(pl.n, u);

        // 投影到 2D
        std::vector<Vec2> P2; P2.reserve(inl.size());
        for(int i:inl){ Vec3 q=Pf[i]-pl.c; P2.push_back({dot(q,u), dot(q,v)}); }

        // 补全面：把补全矩形作为"虚拟支撑点"并入点集，让 genRows 在补全区也布行。
        // 否则会出现"面补齐了、轨迹还是只走有数据的那一段"。
        //
        // 关键：虚拟点必须**沿矩形的两条边**铺，不能铺成轴对齐包围盒（AABB）。
        // 矩形在轨迹系里与行方向存在几度偏差时，990×28.5 的长条取 AABB 会变成
        // 990×107（4.5° 偏差就能撑出 78mm），于是窄边被排出 9 行、轨迹张成扇形。
        // 同时把该面的行方向**强制取矩形长边**，让矩形在旋转后严格轴对齐。
        double forcedTheta = -1.0;
        if(planeFilled[si]){
            Vec3 c00 = planeFillRect[si][0];
            Vec3 ex  = planeFillRect[si][1] - c00;      // 矩形长边（build 的 u 向）
            Vec3 ey  = planeFillRect[si][3] - c00;      // 矩形短边（build 的 v 向）
            double lx = len(ex), ly = len(ey);
            if(lx > 1e-6 && ly > 1e-6){
                double exu = dot(ex,u), exv = dot(ex,v);
                forcedTheta = std::atan2(exv, exu);
                double step = std::max(2.0, params->scanRowSpacing*0.5);
                int nx = std::max(1, (int)std::ceil(lx/step));
                int ny = std::max(1, (int)std::ceil(ly/step));
                for(int ix=0; ix<=nx; ++ix)
                    for(int iy=0; iy<=ny; ++iy){
                        Vec3 pos = c00 + ex*((double)ix/nx) + ey*((double)iy/ny);
                        Vec3 q = pos - pl.c;
                        P2.push_back({dot(q,u), dot(q,v)});
                    }
                fprintf(stderr, "[traj] 面%zu: 补全面 %.0f x %.0f mm（含 %.0fmm 补全区），行方向锁定长边 %.1f°\n",
                        si, lx, ly, ly, forcedTheta*180.0/3.141592653589793);
            }
        }

        // 主方向选择：
        //   scanDirMode=0 自动：穷举角度(每5°)，按"行端锯齿 + 行数代价"选最平滑方向
        //   scanDirMode=1 沿长边喷：行平行长边（换行次数最少），行方向 = 长轴
        //   scanDirMode=2 沿短边喷：行平行短边（换行次数多、行长短），行方向 = 短轴
        // 行端/行范围用 P0.5/P99.5 分位数包围盒（全长覆盖，防长尾离群点，
        // 也防稀疏端部被行过滤裁掉）
        double theta = 0; std::vector<Row> bestRows;
        // 行方向必须贴合"重建矩形的长边"：重建面已经规整成矩形，行若与矩形边
        // 斜交就会看起来歪斜。原来按 5° 网格穷举搜索，可能落在离边几度到十几度
        // 的方向上；这里直接用最小外接矩形的边作为基准方向。
        bool haveRectAng = false; double thLong = 0.0;
        {
            std::vector<Vec2> rect;
            double ra = minAreaRect(P2, rect);
            if(ra > 1e-9 && rect.size() == 4){
                Vec2 e0 = rect[1]-rect[0], e1 = rect[2]-rect[1];
                Vec2 e = (len2(e0) >= len2(e1)) ? e0 : e1;
                thLong = std::atan2(e.y, e.x);
                haveRectAng = true;
            }
        }
        auto normAng=[&](double a){ while(a<0) a+=3.141592653589793; while(a>=3.141592653589793) a-=3.141592653589793; return a; };
        // 兜底校验：行方向必须落在"面内更长的那一维"上。
        // minAreaRect 的边序对窄带点集不可靠（现场长板端部 1012×30 的棱条上，
        // thLong 被判成了 30mm 那一边），于是行沿 30mm 往返、沿 1012mm 前行，
        // 生成 83 条极短行——视觉上就是那条"密密麻麻的小锯齿"。
        auto extentAlong=[&](double a, double& e1, double& e2){
            double c=std::cos(a), s=std::sin(a);
            double lo=1e18,hi=-1e18,lo2=1e18,hi2=-1e18;
            for(auto& q : P2){
                double t=c*q.x+s*q.y, u=-s*q.x+c*q.y;
                lo=std::min(lo,t); hi=std::max(hi,t);
                lo2=std::min(lo2,u); hi2=std::max(hi2,u);
            }
            e1=hi-lo; e2=hi2-lo2;
        };
        if(forcedTheta >= 0){
            // 补全面：行方向由矩形长边唯一确定，不需要 minAreaRect 也不需要穷举。
            // scanDirMode=2（沿短边喷）时转 90°。
            theta = normAng((params->scanDirMode == 2)
                            ? forcedTheta + 1.5707963267948966 : forcedTheta);
            bestRows = genRows(P2, theta, L, params->fillHoleSize);
        }
        else if(haveRectAng){
            // mode1/0 = 沿长边；mode2 = 沿短边（长边方向转 90°）
            theta = normAng((params->scanDirMode == 2) ? thLong + 1.5707963267948966 : thLong);
            double e1=0,e2=0; extentAlong(theta,e1,e2);
            bool wantShort = (params->scanDirMode == 2);
            if((!wantShort && e1 < e2) || (wantShort && e1 > e2)){
                theta = normAng(theta + 1.5707963267948966);
                std::swap(e1,e2);
            }
            bestRows = genRows(P2, theta, L, params->fillHoleSize);
            // 只有"面足够宽、却只排出不到 3 行"才说明方向选错了、需要重选；
            // 窄面（垂直行向的宽度撑不满 3 行）本来就只能铺 1~2 行，必须保留，
            // 否则会掉进下面的穷举搜索，选到沿窄边的方向 → 上百条 20mm 的碎行。
            if(bestRows.size() < 3 && e2 >= 3.0*L) bestRows.clear();
        }
        if(bestRows.empty()){
            double bestCost=1e18;
            for(int deg=0; deg<180; deg+=5){
                double ang=deg*3.141592653589793/180.0;
                auto rs=genRows(P2,ang,L,params->fillHoleSize);
                if(rs.size()<2) continue;
                double jag=0;
                for(size_t i=1;i<rs.size();++i)
                    jag += std::fabs(rs[i].xs-rs[i-1].xs) + std::fabs(rs[i].xe-rs[i-1].xe);
                jag /= (rs.size()-1);
                double cost = jag + rs.size()*L*0.02;
                if(cost<bestCost){ bestCost=cost; theta=ang; bestRows=rs; }
            }
        }
        if(bestRows.empty()){
            fprintf(stderr, "[traj] 段(点数=%zu) 没有可用行，跳过\n", inl.size());
            continue;
        }
        // ---- 环形窄面：平行扫描极度浪费，改走同心圆 ----
        // 顶环唇边 r143~175、宽 32mm、外径 350。按行距 30mm 平行扫描要铺 12 条弦、
        // 总长 3.6m，而环面积只有 0.032m² —— 74% 的行程在环内空洞上空走，界面上
        // 就是"一堆密集平行线"。走同心圆只要 1 圈 ~1.0m（实测 3568→1046mm）。
        // 四条判据同时成立才改（缺一条都会误伤矩形面 / 实心圆盘）：
        //   ① 方位覆盖 ≥34/36 格 —— 确实绕成一圈；
        //   ② 环宽 w = p97-p2 ≤ 2.5×行距 —— 窄带；宽了本来就该正常铺行；
        //   ③ 内径 ≥ 0.45×外径 —— 排除实心圆盘（底面 r≤57，内径≈0）；
        //   ④ 外径 ≥ 3×环宽 —— 排除小圆垫片（一圈周长还没行距长）。
        bool ringMode = false;
        Vec2 ringC{0,0}; double ringIn=0, ringOut=0; int nRing=0;
        if(forcedTheta < 0){
            double rw=0; int ns=0;
            const bool ringish = detect_ring_uv(P2, L, ringC, ringIn, ringOut, &rw, &ns);
            fprintf(stderr, "[traj] 面%zu 环形判据: 内径%.1f 外径%.1f 环宽%.1f 方位%d/36 → %s\n",
                    si, ringIn, ringOut, rw, ns, ringish ? "同心圆" : "平行扫描");
            if(ringish){
                ringMode = true;
                nRing = std::max(1, (int)std::floor(rw/L + 0.5));
            }
        }
        totalRows += ringMode ? nRing : (int)bestRows.size();
        {
            double angN = std::acos(std::min(1.0,std::fabs(dot(pl.n,Vec3{0,0,1}))))*180.0/3.141592653589793;
            if(ringMode)
                fprintf(stderr, "[traj] 面%zu: 内点=%zu 行数=%d(同心圆) 与Z夹角=%.1f°\n",
                    si, inl.size(), nRing, angN);
            else
                fprintf(stderr, "[traj] 面%zu: 内点=%zu 行数=%zu 与Z夹角=%.1f° 行方向=%.1f°\n",
                    si, inl.size(), bestRows.size(), angN, theta*180.0/3.141592653589793);
        }

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
        if(ringMode){
            // 同心圆：每圈首尾闭合，相邻圈反向下刀，圈间只在角度 0 处横移 w/nRing。
            // 圈间那一小段也要按点距插点 —— 实测 2 圈时跳 12mm（pstep 的 3 倍），
            // 会让 seg_max 从 4.0 跳到 12.0，看着像点距失控。
            auto toRow = [&](Vec2 q)->Vec2{ return {ca*q.x + sa*q.y, -sa*q.x + ca*q.y}; };
            Vec2 prevEnd{0,0}; bool hasPrev = false;
            for(int k=0;k<nRing;++k){
                const double rad  = ringIn + (ringOut-ringIn)*((double)k+0.5)/(double)nRing;
                const double circ = 2.0*3.14159265358979*rad;
                int np = std::max(8, (int)std::ceil(circ/pstep));
                const bool rev = ((k & 1) != 0);
                {   // 圈间过渡（角度 0 处，沿半径方向 w/nRing）
                    Vec2 q0{ ringC.x + rad, ringC.y };
                    Vec2 first = toRow(q0);
                    if(hasPrev){
                        double g = len2(first - prevEnd);
                        int nI = std::max(1, (int)std::ceil(g/pstep));
                        for(int t=1;t<nI;++t)
                            traj2.push_back({prevEnd + (first-prevEnd)*((double)t/(double)nI), 1});
                    }
                }
                for(int t=0;t<=np;++t){
                    int ti = rev ? (np - t) : t;
                    double aa = 2.0*3.14159265358979*(double)ti/(double)np;
                    Vec2 q{ ringC.x + rad*std::cos(aa), ringC.y + rad*std::sin(aa) };
                    traj2.push_back({toRow(q), 1});
                }
                prevEnd = traj2.back().first; hasPrev = true;
            }
        }
        if(!ringMode)
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

        // 面间过渡：从上一面终点出发，抬枪绕行到本面首点。
        // 全部为 spray=0 空移点，WPF 与上位机据此断开喷涂段。
        if(!allPts.empty() && !traj2.empty()){
            Vec3 pEnd{ allPts.back().x, allPts.back().y, allPts.back().z };
            Vec2 q0 = rotInv(traj2[0].first);
            Vec3 X0 = pl.c + u*q0.x + v*q0.y;
            Vec3 pStart = X0 + pl.n*params->standoff;

            auto pushMove=[&](const Vec3& p){
                allPts.push_back({p.x,p.y,p.z, gun.x,gun.y,gun.z, params->speed, 0});
            };
            // ① 沿上一面法向退出（离开该面，不擦面）
            Vec3 L1 = pEnd + prevNormal*clear;
            pushMove(L1);
            // ② 抬到全局安全高度（最高面上二者等高，不会出现"先冲高再回落"）
            pushMove({L1.x, L1.y, zSafe});
            // ③ 在安全高度上平移到本面首点正上方
            pushMove({pStart.x, pStart.y, zSafe});
            ++transCount;
        }
        prevNormal = pl.n;

        for(auto& [p2, sp] : traj2){
            Vec2 q = rotInv(p2);
            Vec3 X = pl.c + u*q.x + v*q.y;
            Vec3 pos = X + pl.n*params->standoff;
            allPts.push_back({pos.x,pos.y,pos.z, gun.x,gun.y,gun.z, params->speed, sp});
        }
    }
    if(transCount)
        fprintf(stderr, "[traj] 面间过渡: %d 处，安全高度 z=%.1f（工件最高 %.1f）\n",
                transCount, zSafe, zPartMax);

    if(allPts.empty()){ err="无有效扫描行"; res.error=_strdup(err.c_str()); return res; }
    res.rows = totalRows;

    // 曲面件：把"平面残差"换成"行到曲面的贴合残差"（行追踪跑完才有的值）
    if(reportCurved && surf.nRow > 0){
        fprintf(stderr, "[曲面] 平面残差 RMS=%.2fmm（曲面件正常），改报曲面贴合残差=%.2fmm\n",
                res.plane_rms, surf.fitRms);
        res.plane_rms = surf.fitRms;
    }

    // 写回结果数组
    int cnt=(int)allPts.size();
    res.points=(SprayTrajPoint*)malloc(sizeof(SprayTrajPoint)*cnt);
    res.count=cnt;
    std::memcpy(res.points, allPts.data(), sizeof(SprayTrajPoint)*cnt);

    // 自检：点距（仅统计"喷涂段内"的相邻点）。
    // 空移点（面间抬枪/跨段跳转）必须排除：它们的间距由避障路径决定（本例 5.8mm、
    // 100mm 不等），混进来会把 seg_len_max 抬高，看起来像"点距失控"。
    double sm=0,smax=0; int segCnt=0;
    for(int i=1;i<cnt;++i){
        if(res.points[i].spray==0 || res.points[i-1].spray==0) continue;
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
