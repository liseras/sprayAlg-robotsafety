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
#include <queue>
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
    // mask = 子树包含的面标签位图（面数 <=64 时有效）；用于"跳过本面点"查询：
    // 本面点都在枪位下方 standoff 处，不可能进入安全半径，整棵子树可直接剪掉。
    struct Node { double x, y, z; int idx; int left=-1, right=-1; unsigned char axis=0; unsigned long long mask=~0ULL; };
    std::vector<Node> nodes;        // 预分配到点数大小（一节点一点），支持并行构建
    std::vector<int>  ids;          // 构建期的索引工作区（原地切分，不拷贝子数组）
    const std::vector<Vec3>* pts=nullptr;
    std::atomic<int> nextNode{0};   // 节点槽位分配（并行构建时由各子树竞争）

    // lo..hi 为左闭右开区间；返回子树根节点下标
    int build_rec(int lo, int hi, const std::vector<int>* lab){
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
        nodes[ni] = Node{p.x, p.y, p.z, id, -1, -1, (unsigned char)axis, ~0ULL};
        unsigned long long own = ~0ULL;
        if(lab){
            const int f = (*lab)[id];
            if(f >= 0 && f < 64) own = 1ULL << f;
        }
        nodes[ni].mask = own;
        // 左右子树规模都够大时并行构建：两侧互不相干（各自操作 ids 的不相交区间），
        // 树的结构与串行完全一致，只有节点槽位分配顺序不同（内部实现细节）。
        if(hi - lo > 60000){
            std::thread th([&]{ nodes[ni].left = build_rec(lo, mid, lab); });
            nodes[ni].right = build_rec(mid+1, hi, lab);
            th.join();
        } else {
            nodes[ni].left  = build_rec(lo,  mid, lab);
            nodes[ni].right = build_rec(mid+1, hi, lab);
        }
        unsigned long long lm = (nodes[ni].left >=0) ? nodes[nodes[ni].left ].mask : 0ULL;
        unsigned long long rm = (nodes[ni].right>=0) ? nodes[nodes[ni].right].mask : 0ULL;
        nodes[ni].mask = own | lm | rm;
        return ni;
    }
    void init(const std::vector<Vec3>& p, const std::vector<int>* lab=nullptr){
        pts=&p;
        // 节点数恒等于点数，预分配后并行构建就不必再扩容
        nodes.assign(p.size(), Node{});
        nextNode.store(0);
        ids.resize(p.size()); std::iota(ids.begin(),ids.end(),0);
        if(p.empty()){ std::vector<int>().swap(ids); return; }
        build_rec(0, (int)p.size(), lab);
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

    // 半径上界最近邻：只关心"有没有点比 r 更近"时用。
    // 枪位在工件外 standoff(100mm) 处，最近点本身就在 100mm 外；普通 knn(k=1)
    // 会把 best 初始化到 100mm，于是所有横向 100mm 内的分支都要展开，6 万点云上
    // 单次查询可达 0.6ms。这里把剪枝上界直接压到安全阈值 r（20~25mm），
    // 超出 r 的分支立即剪掉；返回 min(真实最近距离, r)。
    // 所有调用点都只做 "< r" / ">= r" 判定，因此结果与精确 KNN 逐位一致。
    double nearestWithin(Vec3 q, double r) const {
        double best2 = r*r;
        recWithin(0, q, best2);
        return std::sqrt(best2);
    }
    void recWithin(int node, const Vec3& q, double& best2) const {
        if(node < 0) return;
        const Node& nd = nodes[node];
        double dx=q.x-nd.x, dy=q.y-nd.y, dz=q.z-nd.z;
        double d2=dx*dx+dy*dy+dz*dz;
        if(d2 < best2) best2 = d2;
        int axis = nd.axis;
        double diff = (axis==0)? dx : (axis==1)? dy : dz;
        int near = diff<0 ? nd.left : nd.right;
        int far  = diff<0 ? nd.right : nd.left;
        recWithin(near, q, best2);
        if(diff*diff < best2) recWithin(far, q, best2);
    }

    // 同上，但跳过 faceMask 对应的面：仅当该面所有点都远离查询点（> r）时使用。
    double nearestWithinExcept(Vec3 q, double r, unsigned long long faceMask) const {
        if(faceMask == 0ULL) return nearestWithin(q, r);
        double best2 = r*r;
        recWithinExcept(0, q, best2, faceMask);
        return std::sqrt(best2);
    }
    void recWithinExcept(int node, const Vec3& q, double& best2, unsigned long long faceMask) const {
        if(node < 0) return;
        const Node& nd = nodes[node];
        if((nd.mask & ~faceMask) == 0ULL) return;   // 子树只剩被跳过的面
        double dx=q.x-nd.x, dy=q.y-nd.y, dz=q.z-nd.z;
        double d2=dx*dx+dy*dy+dz*dz;
        if(d2 < best2) best2 = d2;
        int axis = nd.axis;
        double diff = (axis==0)? dx : (axis==1)? dy : dz;
        int near = diff<0 ? nd.left : nd.right;
        int far  = diff<0 ? nd.right : nd.left;
        recWithinExcept(near, q, best2, faceMask);
        if(diff*diff < best2) recWithinExcept(far, q, best2, faceMask);
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

    // ---- 半径查询（射线自由空间判定用）----
    // 需要"我在射线上前进时，前方有没有点"，而 KNN 只能给最近邻居：
    // 同面邻点永远是最邻近的，会把真正的远处命中挡住。所以按球收集候选，
    // 再由调用方用"沿射线的投影 + 横向偏移"筛。
    // 沿枪轴找圆柱内第一个实测点，保留点号供覆盖归属判断。
    int firstRayPoint(Vec3 origin, Vec3 dir, double maxDist, double radius) const {
        double best = maxDist;
        int hit = -1;
        const double radius2 = radius*radius;
        auto visit = [&](auto&& self, int id) -> void {
            if(id < 0) return;
            const Node& node = nodes[id];
            const Vec3 off = Vec3{node.x, node.y, node.z} - origin;
            const double dist = dot(off, dir);
            if(dist >= 0.0 && dist <= best && dot(off, off) - dist*dist <= radius2 &&
               (dist < best || hit < 0 || node.idx < hit)){
                best = dist;
                hit = node.idx;
            }
            const double coord = node.axis == 0 ? off.x : (node.axis == 1 ? off.y : off.z);
            const double slope = node.axis == 0 ? dir.x : (node.axis == 1 ? dir.y : dir.z);
            if(coord >= std::min(0.0, best*slope) - radius) self(self, node.left);
            if(coord <= std::max(0.0, best*slope) + radius) self(self, node.right);
        };
        if(!nodes.empty()) visit(visit, 0);
        return hit;
    }

    int radiusCollect(Vec3 q, double r, int* out, int maxOut) const {
        int n = 0;
        recCollect(0, q, r*r, out, n, maxOut);
        return n;
    }
    void recCollect(int node, const Vec3& q, double r2, int* out, int& n, int maxOut) const {
        if(node < 0 || n >= maxOut) return;
        const Node& nd = nodes[node];
        double dx=q.x-nd.x, dy=q.y-nd.y, dz=q.z-nd.z;
        double d2=dx*dx+dy*dy+dz*dz;
        if(d2 <= r2) out[n++] = nd.idx;
        int axis = nd.axis;
        double diff = (axis==0)? dx : (axis==1)? dy : dz;
        int near = diff<0 ? nd.left : nd.right;
        int far  = diff<0 ? nd.right : nd.left;
        recCollect(near, q, r2, out, n, maxOut);
        if(diff*diff <= r2) recCollect(far, q, r2, out, n, maxOut);
    }
};

// ---------------- 主实现 ----------------
namespace {

struct TrajPoint { Vec3 pos; Vec3 dir; double speed; int spray; };

// ---------------- 体素降采样（大点云的第一道闸） ----------------
// 为什么需要：算法里所有"逐点"阶段（SOR 邻域、法向 PCA、区域生长、枪位净空过滤）
// 成本都是 O(N × 局部密度)。线扫点云到 100 万点时点间距已到 0.3~0.5mm，而工艺
// 只要点距 4mm、行距 12mm —— 超出工艺需求的采样对轨迹没有任何贡献，只是把每一段
// 都拖慢十几倍。实测 212 万点整条流程 128s，其中 116s 花在「内点枪位过滤」这一个
// 逐点 KD 查询上。降采样后同一块板 15 万点，全流程进 2s 量级，且**点数不再随
// 扫描仪分辨率增长**（分辨率翻倍只是把体素填得更满，输出点数不变）。
//
// 体素代表点：把点云按 leaf 立方格分块，每块取"离格心最近的**实测点**"。
// 不取格内均值，理由有两条：
//   · 不产生点云里不存在的坐标，棱边/薄壁不会被均值拉平（薄壁件最怕这个）；
//   · 抽稀后的点集是原始点的**子集**，几何关系只会变稀、不会变形。
//
// 分工口径（用户 2026-11-20 15:46）：抽稀**只服务布行与安全判据**，
// 分割/重建/碎面判定/轮廓一律走全分辨率。实测依据：p047 现场件 226114 点，
// 全分辨率分割出 2 个面；2mm 体素后出 2 面 + 1 线状 + 1 碎面（还得靠环级丢弃
// 才拦住）。那个"碎面"是体素把主面撕下来的一块，不是真碎面 —— 碎面过滤在
// 相当程度上是在给体素擦屁股。全分辨率下它被主面正常吸收。
//
// 分辨率口径见 auto_voxel_leaf()：锚定重建尺度 kReconGap/6，下限守点距/2
// （采样定理：点距 p 的点云几何信息只到 p/2），上限 10mm。与工艺行距无关。
//
// 体素代表点的**索引**版本：不复制点，只挑出"每格最靠近格心"的那个原始下标。
// 返回下标而不是复制点：分割结果 seg.inl 用的是原始点下标，两套口径共用
// 同一个索引空间才能互相对照（抽稀点只是其中一个子集）。
static std::vector<int> voxel_downsample_idx(const std::vector<Vec3>& P, double leaf)
{
    std::vector<int> out;
    if(leaf <= 1e-9 || P.size() < 20000){ out.resize(P.size()); for(size_t i=0;i<P.size();++i) out[i]=(int)i; return out; }
    const double inv = 1.0/leaf;
    struct KV { unsigned long long key; int idx; };
    std::vector<KV> kv(P.size());
    const long long OFF = 1LL<<20;
    for(size_t i=0;i<P.size();++i){
        long long ix = (long long)std::floor(P[i].x*inv) + OFF;
        long long iy = (long long)std::floor(P[i].y*inv) + OFF;
        long long iz = (long long)std::floor(P[i].z*inv) + OFF;
        ix = std::max(0LL, std::min(0x1FFFFFLL, ix));
        iy = std::max(0LL, std::min(0x1FFFFFLL, iy));
        iz = std::max(0LL, std::min(0x1FFFFFLL, iz));
        kv[i].key = ((unsigned long long)ix<<42) | ((unsigned long long)iy<<21) | (unsigned long long)iz;
        kv[i].idx = (int)i;
    }
    std::sort(kv.begin(), kv.end(), [](const KV&a, const KV&b){ return a.key < b.key; });
    out.reserve(P.size()/4+16);
    size_t i=0;
    while(i<kv.size()){
        const unsigned long long key = kv[i].key;
        const long long gx = (long long)((key>>42)&0x1FFFFFULL) - OFF;
        const long long gy = (long long)((key>>21)&0x1FFFFFULL) - OFF;
        const long long gz = (long long)(key & 0x1FFFFFULL) - OFF;
        const Vec3 c{ ((double)gx+0.5)*leaf, ((double)gy+0.5)*leaf, ((double)gz+0.5)*leaf };
        int bestIdx = -1; double bestD = 1e300;
        for(; i<kv.size() && kv[i].key==key; ++i){
            const Vec3& p = P[kv[i].idx];
            const double dx=p.x-c.x, dy=p.y-c.y, dz=p.z-c.z;
            const double d = dx*dx+dy*dy+dz*dz;
            if(d < bestD){ bestD = d; bestIdx = kv[i].idx; }
        }
        if(bestIdx >= 0) out.push_back(bestIdx);
    }
    return out;
}

// 重建/分割专用的几何尺度（mm），与工艺参数 scanRowSpacing 完全解耦。
// 为什么不能用行距：行距是"布行"参数，同一份点云可以填 12/20/60，重建结果
// 必须完全一致。现场点云实测线间距仅 1~8mm（KNN 估计），12mm 是历史默认
// 重建尺度，对稀疏边缘仍有余量，同时尽量保持既有回归基线。
static constexpr double kReconGap = 12.0;

// 自动档体素边长：**只由重建尺度定，绝不用工艺行距**。
// 历史口径是 clamp(点距/2, 行距/8, 10)：行距 12 → 2mm、行距 120 → 10mm，
// 于是同一份点云填不同行距会得到**不同的降采样点云 → 不同的平面重建**，
// 用户 2026-11-20 15:28 实测 p047 行距 12 出 2 个面、行距 120 出 3 个面，
// 就是这个口径漏出来的（重建与工艺必须是两条独立的链）。
// 现在：分辨率锚定重建尺度 kReconGap/6 = 2mm（"重建格距的 1/6"，比一个重
// 建格子细一档，足够分辨轮廓与窄筋），下限守点距/2（采样定理），上限 10mm
// （实测 10mm 在 8° 斜面上出现 15mm 台阶）。
static double auto_voxel_leaf(const SprayTrajParams* p){
    const double byGeo   = kReconGap/6.0;
    const double byPoint = p->pointSpacing*0.5;
    return std::min(10.0, std::max(byPoint, byGeo));
}
// 兜底重试用的临时档位（>0 时覆盖自动档；-1 = 不覆盖）。
// 为什么需要：双层薄壳在"体素边长大致等于壁厚"的窄带里，两层会被切成两个面、
// 法向朝壳内，枪位验收判"两侧都伸不进枪"→ 整件不出轨迹（实测 4.0mm 踩中，3/5/7mm 都正常）。
// 这不是算法错（拒喷是安全侧），但用户没必要为此重跑，退一档重新算即可。
static double g_voxelLeafOverride = -1.0;

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
    // 阈值就是经典的 mean + z·σ。
    // 曾试过加下限"中位数×1.6"去保住棱边点（边界点邻域只有半边支撑，meanKNN
    // 比内部高 28%，干净点云 σ 小时会被整圈删掉）。方向是对的，但改在这里
    // 代价太大：保住的边界点/重复点改变了后续法向聚类与拆层，实测薄壁盒
    // g_norim 从 10 面塌成 1 面、曲面件 curved 从 2 面变 3 面并冒出 0.05mm
    // 的贴件轨迹。棱边看不见的问题改在"枪位净空"里解决——那里有解析平面，
    // 不依赖点云是否把棱边采到（见 buildBlkFaces / gunClearAt）。
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

// 多平面分割：法向区域生长（技术方案 2.4 节），当前实现见 split_planes。
//
// 参数语义（2026-11-19 重定，替代旧的"高度分层"口径）：
//   tol (flatnessTol) = **吸收/覆盖阈值**：离这张面多远以内的点，仍然算这张面覆盖到
//                       （进入 seg.inl 参与轮廓、覆盖率、残点归并），但**不参与解平面方程**。
//                       它不再是"平面内点容差"，也不再被静默截断成 3mm；真正的上限由
//                       净空预算给出：absorbTol = min(tol, 0.6×standoff)。
//   拟合尺度           = 固定 3mm（解方程的内点尺度，不承担工艺语义）。噪声 σ 只用于
//                       诊断上报：σ 明显偏大说明点云噪声/边缘多，但不自动放宽拟合尺度——
//                       实测薄壁盒放宽到 4mm，主面就会从开口面翻到侧壁。
//   curvatureTol      = 整体高度极差 < 此值 → 单一可喷面（不切分）。
//
// 为什么必须拆成两个尺度：RANSAC 投票容差一旦跟着 tol 放大，每个假设的内点数都≈全点，
// 打分失去区分度，选出的假设等于随机抽样——平面会被拉斜、并跨台阶。而"吸收"如果参与
// 重拟合，吸收多少就歪多少。两者是互斥需求，只能各用一个数。
// ---------------- 面级法向诊断（内壁 / 外壁） ----------------
// 背景：PCA 法向本身无符号（符号由数值求解决定），原来的定向只有两条规则：
//   |n.z|>0.05 → 归 +Z 半球；否则"背离整件质心"。
// 这条规则的隐含假设是"工件凸、只扫到外表面"，对开口盒的**内壁**恰好全反——
// 内壁的可喷侧朝腔体，正是"朝向质心"的一侧。
//
// 这里改用射线自由空间判材料侧：沿 ±n 各打若干条射线，
//   某一侧在 nearWall 内命中点云 → 该侧是材料（壁的另一层 / 工件的实体侧）
//   另一侧没有任何命中       → 该侧是可喷侧 = 真实法向应指向的方向
//   两侧都自由（单层壳：只扫到内表面、外表面没建出来）→ 射线无信息，
//   保留原定向并标 src=2，界面提示人工确认。这种情况几何上确实二义。
struct FaceOrient {
    Vec3   n{0,0,1};      // 真实法向（可喷侧）
    double angZ    = 0;   // 与 +Z 夹角 0~180°
    double tilt    = 0;   // 相对水平面仰角 -90~90°（= 90 - angZ）
    int    side    = 1;   // +1 正向 / -1 反向 / 0 竖直
    int    isInner = 0;   // 1 = 内壁
    int    src     = 2;   // 1=射线判定 2=质心兜底 3=强制翻转
    double freePlus  = -1;// 沿 +n 自由距离 mm（-1 = 量程内无命中）
    double freeMinus = -1;// 沿 -n 自由距离 mm
};

struct PlaneSeg {
    Plane pl;
    std::vector<int> inl;
    // 曲面块标记：true 时 pl 只是"平均平面"（用于定向/统计），真正的轨迹
    // 生成走曲面行追踪（gen_surface_rows），不是平面 2D 布行。
    bool curved = false;
    // 通过缓翘曲塌缩检查，保留平面布行，仅豁免后续平面性拒喷。
    bool collapsed = false;
    // 载体拟合或拆层时的实测支撑点；吸收扩展覆盖后仍用它们核对当前平面。
    std::vector<int> fitSupport;
    bool preservedCarrier = false;   // 合并重拟合失效，沿通过校验的外侧载体布行。
    FaceOrient fo;        // 面级朝向诊断；angZ 参与轨迹生成（朝下超窗口不喷），其余只上报
    // 单层壳 / 薄壁双层：1=单层 2=双层 0=未判定；layerThick 是沿法向的分位跨度
    int    layerCnt = 0;
    double layerThick = 0.0;
    // 平面冻结：true = 这张面的平面方程由"吸收"一次性定死，后续任何阶段的
    // 合并/重拟合都不得改写它（只能往 inl 里加点）。语义见 split_planes 里的吸收段：
    // 用户设的"面吸收阈值"是工艺口径，吸收完这张面就定了，不该被下游重算。
    bool   fixed = false;
    // 枪位验收专用核心点：最终合并闸门把侧面折边并进主面时，先在这里存一份
    // **合并前的本面点**。原因（2026-11-23 实测）：merge 把 83~88° 的侧壁并进主面后，
    // s.inl 从 9142 涨到 17386，枪位验收按 s.inl 均匀取 6 点，一旦取到侧壁根部的点，
    // 沿**主面法向**打射线 1.1mm 就撞到侧壁本体 → 整张主面判"枪伸不进去"而拒喷
    // （p001 主面 8433 点、p009 主面 15631 点都是这么掉的）。
    // 枪位验收问的是"从这张面自己的面上打出去空不空"，侧壁不回答这个问题。
    // 空 = 没发生过合并，验收与布行都用 inl（与改动前逐位相同）。
    std::vector<int> core;
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

// 面内最小外接矩形的短边（mm）。用于"最小可喷宽度"门限：
//   点云近似共线时，任意包含这条线的平面残差都≈0，flatnessTol 拦不住；
//   但它的二维凸包只有一条细带，最小外接矩形短边≈重建格距/点距。
// 返回值 <0 表示点太少、凸包退化；长边写入 longSide（可为空）。
static double minRectShortSide(const std::vector<Vec2>& pts, double* longSide = nullptr){
    std::vector<Vec2> rect;
    if(minAreaRect(pts, rect) < 0.0){
        if(longSide) *longSide = -1.0;
        return -1.0;
    }
    const double a = len2(rect[1] - rect[0]);
    const double b = len2(rect[2] - rect[1]);
    if(longSide) *longSide = std::max(a, b);
    return std::min(a, b);
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

// ---------------- 环带射线模型：环形 / 类环形 ----------------
// detect_ring_uv 只认"正圆 + 与形心同心"：它拿点到形心的半径 p2/p97 当环宽。
// 方框（矩形焊框）在角上的半径是边中点的 1.41 倍，p97-p2 立刻从 30mm 涨到 300mm，
// 于是所有**非圆**的环带都被判成"平行扫描"：一条连续环带被切成几十条短行，
// 每条还要被中央空洞断成两段，换行/抬枪比有效行程还多。
//
// 这里改成"从环内一点向 N 个方位发射线"的带宽模型：
//   射线上被材料占据的那一段 = 环带在该方位的径向跨度 [rIn, rOut]。
//   对方位不做任何形状假设 → 圆环、椭圆环、方框、跑道形、圆角矩形、
//   甚至开口环（射线落在缺口上就记 rIn<0）都能统一描述，
//   轨迹就是"同一方位按等分取点、沿方位角连成的闭合/开口曲线"。
//
// 判据（缺一条都会误伤实心盘 / 带小孔的板 / 实心条）：
//   ① 有材料的方位 >= 55%（开口环允许缺口）；
//   ② 单条射线上出现多段材料的方位 <= 30%（凹陷形状/中间有凸台时才会多段）；
//   ③ 环宽中位 w <= 3×行距（更宽的带子平行扫描本来就正常）；
//   ④ 内径中位 >= 0.35×外径中位 —— 中间确实是空的（实心盘 rIn≈0）；
//   ⑤ 外径中位 >= 3×环宽 —— 一圈周长显著大于带宽（排除小垫片）。
struct BandUV {
    Vec2  c{0,0};                      // 环内参考点（射线原点）
    int   N = 0;                       // 方位采样数
    std::vector<float> rIn, rOut;      // 每方位内/外半径；rIn<0 = 该方位无材料
    double wMed=0, wP90=0;             // 射线截长（中位 / p90）——只作诊断
    double wPerp=0;                    // 垂直带宽：环面积 / 中线周长（判据用它）
    double rInMed=0, rOutMed=0;        // 内径 / 外径中位
    double perim=0;                    // 中线长度（覆盖方位上的积分）
    double closure=0;                  // 有材料的方位比例
    double multi=0;                    // 射线上命中多段材料的方位比例
    bool   ok=false;
};

// 内部：在已膨胀的占据栅格上做一次射线扫描
// 注意：这里**不接受工艺行距**。判"是不是一条环带"是纯几何问题，同一块面
// 的行距填 12 还是 120 都必须给出同一个答案（用户 2026-11-20 15:28 反馈：
// 行距 120 时小碎面反而不被过滤了 —— 就是这里曾用 6×行距 当带宽上限，
// 行距 120 等于把上限放到 720mm，几乎任何面都能冒充环带而被豁免）。
static BandUV band_eval_seed(const std::vector<unsigned char>& dil, int nx, int ny,
                             double x0, double y0, double cell,
                             Vec2 c)
{
    BandUV o; o.c = c;
    const int N = std::min(512, std::max(96, (int)std::lround(
                    3.14159265358979*std::sqrt((nx*cell)*(nx*cell)+(ny*cell)*(ny*cell))/25.0)));
    o.N = N;
    o.rIn.assign((size_t)N, -1.0f);
    o.rOut.assign((size_t)N, -1.0f);
    double tMax = 0.0;
    for(int sx=0;sx<2;++sx) for(int sy=0;sy<2;++sy){
        Vec2 cor{x0 + (double)(nx-1)*cell*sx, y0 + (double)(ny-1)*cell*sy};
        tMax = std::max(tMax, len2(cor - c));
    }
    tMax += 2.0*cell;
    const double step = cell*0.5;
    int nHit=0, nMulti=0;
    std::vector<double> ws; ws.reserve((size_t)N);
    std::vector<double> ri, ro; ri.reserve((size_t)N); ro.reserve((size_t)N);
    for(int a=0;a<N;++a){
        const double th = 2.0*3.14159265358979*(double)a/(double)N;
        const double dx = std::cos(th), dy = std::sin(th);
        int nRun = 0;
        double best0 = -1, best1 = -1;             // 最长的一段材料
        double runS = -1, runE = -1;
        for(double t=0;t<=tMax;t+=step){
            const double px = c.x + dx*t, py = c.y + dy*t;
            const int ix = (int)std::floor((px - x0)/cell);
            const int iy = (int)std::floor((py - y0)/cell);
            const bool h = (ix>=0 && ix<nx && iy>=0 && iy<ny) && dil[(size_t)iy*nx + ix];
            if(h && runS < 0) runS = t;
            if(!h && runS >= 0){
                runE = t; ++nRun;
                if(best0 < 0 || (runE-runS) > (best1-best0)){ best0=runS; best1=runE; }
                runS = -1;
            }
        }
        if(runS >= 0){ ++nRun; runE = tMax;
            if(best0 < 0 || (runE-runS) > (best1-best0)){ best0=runS; best1=runE; } }
        if(nRun == 1 && best1 - best0 > 1e-6){
            o.rIn[a]  = (float)best0;
            o.rOut[a] = (float)best1;
            ++nHit;
            ws.push_back(best1 - best0);
            ri.push_back(best0); ro.push_back(best1);
        } else if(nRun >= 2){
            // 一条射线穿过两处材料：中间有凸台/两条臂。这里既不能当环带，
            // 也不能用其中一段冒充 —— 记为缺口，让轨迹在那里断开，
            // 由换行连接段跨过去（换行是正常连续喷涂，kind=0，不标色）。
            ++nMulti;
        }
    }
    o.closure = (double)nHit/(double)N;
    o.multi   = (double)nMulti/(double)N;
    if(ws.empty()) return o;
    std::sort(ws.begin(), ws.end());
    std::sort(ri.begin(), ri.end());
    std::sort(ro.begin(), ro.end());
    o.wMed   = ws[ws.size()/2];
    o.wP90   = ws[(size_t)((double)(ws.size()-1)*0.9)];
    o.rInMed = ri[ri.size()/2];
    o.rOutMed= ro[ro.size()/2];
    // ---- 带宽必须用"垂直宽度"，不能用射线截长 ----
    // 方框的角上射线是斜着穿过环带材料的：712×616 外框 + 650×554 内孔，
    // 45° 射线截长 111mm，而垂直带宽只有 31mm。用中位截长当带宽，会把所有
    // 矩形框/跑道环都判成"太宽，走平行扫描"——正是这次要修的那类误判。
    // 做法：用中线（(rIn+rOut)/2 的轨迹）的局部切向，把截长折算成垂直宽度
    //   w⊥ = run · |sin∠(射线, 中线切向)|
    // 取中位：角上少数几个方位偏大，压不过整条边上的方位。
    {
        std::vector<Vec2> mid((size_t)N);
        std::vector<char> has((size_t)N, 0);
        for(int a=0;a<N;++a){
            if(o.rIn[(size_t)a] < 0.0f) continue;
            const double th = 2.0*3.14159265358979*(double)a/(double)N;
            const double rm = 0.5*((double)o.rIn[(size_t)a] + (double)o.rOut[(size_t)a]);
            mid[(size_t)a] = Vec2{ rm*std::cos(th), rm*std::sin(th) };
            has[(size_t)a] = 1;
        }
        std::vector<double> wp; wp.reserve((size_t)N);
        for(int a=0;a<N;++a){
            if(!has[(size_t)a]) continue;
            const int ap = (a+N-1)%N, an = (a+1)%N;
            if(!has[(size_t)ap] || !has[(size_t)an]) continue;
            Vec2 T = norm2(mid[(size_t)an] - mid[(size_t)ap]);
            if(len2(T) < 1e-9) continue;
            const double th = 2.0*3.14159265358979*(double)a/(double)N;
            const Vec2 D{ std::cos(th), std::sin(th) };
            const double run = (double)o.rOut[(size_t)a] - (double)o.rIn[(size_t)a];
            wp.push_back(run*std::fabs(cross2(D, T)));
        }
        if(!wp.empty()){
            std::sort(wp.begin(), wp.end());
            o.wPerp = wp[wp.size()/2];
        } else o.wPerp = o.wMed;
    }
    // 中线长度（覆盖方位上的中位半径 × 角度），用于"一圈是否显著长于带宽"
    double perim = 0.0;
    {
        const double dth = 2.0*3.14159265358979/(double)N;
        for(int a=0;a<N;++a){
            if(o.rIn[(size_t)a] < 0.0f) continue;
            perim += 0.5*((double)o.rIn[(size_t)a] + (double)o.rOut[(size_t)a])*dth;
        }
    }
    o.perim = perim;
    // 带宽上限：不再用"6×行距"这种工艺量，改成纯几何的 rOut >= 3·wPerp
    // （见下面 o.ok 里的这一条）。半带宽撑不到外径 1/3 的"环"其实是整块板
    // 或厚壁筒，平行扫描本来就对；这一条已足够，且与行距无关。
    // 历史口径 6L 的副作用：行距 120 → 上限 720mm → 小碎面全被当成环带豁免。
    o.ok = (o.closure >= 0.55) && (o.multi <= 0.30)
        && (o.wPerp > 1.0)
        && (o.rInMed >= 0.35*o.rOutMed)
        && (o.rOutMed >= 3.0*o.wPerp)
        && (perim >= 3.0*o.wPerp);
    return o;
}

// uv：面内点（任意正交基）。判据里没有任何工艺量，重建/分割/布行三处同源。
static bool detect_band_uv(const std::vector<Vec2>& uv, BandUV& out)
{
    out = BandUV();
    if(uv.size() < 200) return false;
    double x0=1e18,x1=-1e18,y0=1e18,y1=-1e18;
    for(const auto& q : uv){
        x0=std::min(x0,q.x); x1=std::max(x1,q.x);
        y0=std::min(y0,q.y); y1=std::max(y1,q.y);
    }
    const double ext = std::max(x1-x0, y1-y0);
    if(ext < 1e-6) return false;
    const double cell = std::min(6.0, std::max(2.5, ext/150.0));
    const int nx = (int)std::floor((x1-x0)/cell) + 1;
    const int ny = (int)std::floor((y1-y0)/cell) + 1;
    if(nx < 4 || ny < 4 || (long long)nx*ny > 4000000LL) return false;
    std::vector<unsigned char> occ((size_t)nx*ny, 0);
    for(const auto& q : uv){
        int ix=(int)std::floor((q.x-x0)/cell), iy=(int)std::floor((q.y-y0)/cell);
        if(ix<0) ix=0; if(ix>=nx) ix=nx-1;
        if(iy<0) iy=0; if(iy>=ny) iy=ny-1;
        occ[(size_t)iy*nx + ix] = 1;
    }
    // 1 格膨胀：闭合扫描线间距（现场 1~8mm），否则一条带子会被射线采成
    // "带空隙的几段"，环宽被低估、方位覆盖率也被低估。
    std::vector<unsigned char> dil(occ.size(), 0);
    for(int iy=0;iy<ny;++iy)
        for(int ix=0;ix<nx;++ix){
            if(!occ[(size_t)iy*nx + ix]) continue;
            for(int dy=-1;dy<=1;++dy)
                for(int dx=-1;dx<=1;++dx){
                    const int jx=ix+dx, jy=iy+dy;
                    if(jx<0||jx>=nx||jy<0||jy>=ny) continue;
                    dil[(size_t)jy*nx + jx] = 1;
                }
        }
    // 距离场（到最近材料的格数，chamfer 近似）：用来给射线挑一个"最空的种子点"。
    std::vector<float> dt(dil.size(), 1e9f);
    for(size_t k=0;k<dil.size();++k) if(dil[k]) dt[k]=0.0f;
    auto dget=[&](int ix,int iy)->float{
        if(ix<0||ix>=nx||iy<0||iy>=ny) return 1e9f;
        return dt[(size_t)iy*nx + ix];
    };
    for(int iy=0;iy<ny;++iy)
        for(int ix=0;ix<nx;++ix){
            float& d = dt[(size_t)iy*nx + ix];
            if(d == 0.0f) continue;
            float b = d;
            b = std::min(b, dget(ix-1,iy)+1.0f);
            b = std::min(b, dget(ix,iy-1)+1.0f);
            b = std::min(b, dget(ix-1,iy-1)+1.4142f);
            b = std::min(b, dget(ix+1,iy-1)+1.4142f);
            d = b;
        }
    for(int iy=ny-1;iy>=0;--iy)
        for(int ix=nx-1;ix>=0;--ix){
            float& d = dt[(size_t)iy*nx + ix];
            if(d == 0.0f) continue;
            float b = d;
            b = std::min(b, dget(ix+1,iy)+1.0f);
            b = std::min(b, dget(ix,iy+1)+1.0f);
            b = std::min(b, dget(ix+1,iy+1)+1.4142f);
            b = std::min(b, dget(ix-1,iy+1)+1.4142f);
            d = b;
        }
    // 内部空洞 = 从边界进不去的空格。种子优先取它里面离材料最远的一格：
    // 这样开口环（缺口连到外部）也能退化到"形心/包围盒中心"再试一次。
    std::vector<unsigned char> outside(dil.size(), 0);
    std::vector<int> stk; stk.reserve(1024);
    auto pushOut=[&](int ix,int iy){
        if(ix<0||ix>=nx||iy<0||iy>=ny) return;
        const size_t k=(size_t)iy*nx + ix;
        if(outside[k] || dil[k]) return;
        outside[k]=1; stk.push_back((int)k);
    };
    for(int ix=0;ix<nx;++ix){ pushOut(ix,0); pushOut(ix,ny-1); }
    for(int iy=0;iy<ny;++iy){ pushOut(0,iy); pushOut(nx-1,iy); }
    while(!stk.empty()){
        const int k = stk.back(); stk.pop_back();
        const int ix = k%nx, iy = k/nx;
        pushOut(ix-1,iy); pushOut(ix+1,iy); pushOut(ix,iy-1); pushOut(ix,iy+1);
    }
    std::vector<Vec2> seeds;
    {
        int bi=-1; float bd=-1.0f;
        for(int iy=0;iy<ny;++iy)
            for(int ix=0;ix<nx;++ix){
                const size_t k=(size_t)iy*nx + ix;
                if(dil[k] || outside[k]) continue;
                if(dt[k] > bd){ bd = dt[k]; bi = (int)k; }
            }
        if(bi >= 0 && bd >= 1.5f)
            seeds.push_back(Vec2{x0 + ((double)(bi%nx)+0.5)*cell,
                                 y0 + ((double)(bi/nx)+0.5)*cell});
    }
    {
        // 缺口/漏边会让洪水填充从孔里漏到外面（outside 全 1），此时上面的
        // "内部空洞"找不到；再补一个"离材料最远、且不在 2 格边界带内"的种子。
        int bi=-1; float bd=-1.0f;
        for(int iy=2;iy<ny-2;++iy)
            for(int ix=2;ix<nx-2;++ix){
                const size_t k=(size_t)iy*nx + ix;
                if(dil[k]) continue;
                if(dt[k] > bd){ bd = dt[k]; bi = (int)k; }
            }
        if(bi >= 0 && bd >= 1.5f)
            seeds.push_back(Vec2{x0 + ((double)(bi%nx)+0.5)*cell,
                                 y0 + ((double)(bi/nx)+0.5)*cell});
    }
    {
        Vec2 c{0,0};
        for(const auto& q : uv){ c.x += q.x; c.y += q.y; }
        c.x /= (double)uv.size(); c.y /= (double)uv.size();
        seeds.push_back(c);
        seeds.push_back(Vec2{0.5*(x0+x1), 0.5*(y0+y1)});
    }
    // 选种子不能只看"方位覆盖"：形心落在材料上时覆盖率更高（射线从材料内部出发），
    // 但 rIn≈0，那是实心件的特征。真正的环带要从**空腔里的种子**看：
    // 覆盖率 × (1-多段率) × (0.25+0.75·内/外径) —— 既要求绕得满，也要求中间是空的。
    auto quality=[&](const BandUV& e)->double{
        if(e.closure < 0.30) return -1.0;
        const double hole = (e.rOutMed > 1e-9) ? std::min(1.0, e.rInMed/e.rOutMed) : 0.0;
        return e.closure*(1.0 - e.multi)*(0.25 + 0.75*hole);
    };
    BandUV best; bool haveBest=false; double bestQ=-2.0;
    std::vector<Vec2> tried;
    const bool dbg = (std::getenv("SPRAYTRAJ_BAND_DBG") != nullptr);
    for(const Vec2& s : seeds){
        bool dup = false;
        for(const Vec2& t0 : tried) if(len2(t0 - s) < 2.0*cell){ dup = true; break; }
        if(dup) continue;
        tried.push_back(s);
        BandUV e = band_eval_seed(dil, nx, ny, x0, y0, cell, s);
        const double q = quality(e);
        if(dbg)
            fprintf(stderr, "  [band]   种子(%.1f,%.1f): 方位%.0f%% 多段%.0f%% 带宽%.1f(内%.1f 外%.1f) q=%.2f %s\n",
                    s.x, s.y, e.closure*100.0, e.multi*100.0, e.wPerp, e.rInMed, e.rOutMed,
                    q, e.ok ? "OK" : "");
        if(q > bestQ){ best = e; bestQ = q; haveBest = true; }
    }
    out = haveBest ? best : BandUV();
    return out.ok;
}

// kReconGap 在文件前部定义（auto_voxel_leaf 之前），此处不再重复。

// 由面内点集重建规则三角网格（局部 uv 坐标，mm）。
//   reconGap 重建几何尺度（不是工艺行距；跨线空隙用闭运算桥接，否则每行独立成条）
//   fillHole 孔洞填补阈值：连通孔洞的外接跨距 <= 此值 → 填实（覆盖重建）
// 值 2 标出阈值主动填补的格，供布行区分补面与无支撑区域。
static int fill_reconstruction_holes(std::vector<unsigned char>& occ, int width, int height,
                                     double dx, double dy, double maxSpan,
                                     const std::vector<unsigned char>* rimSupport = nullptr)
{
    if(maxSpan <= 0.0) return 0;
    const double mouthLimit = maxSpan/3.0;
    const int rx = (int)std::min((double)width, std::floor(mouthLimit/(2.0*dx)));
    const int ry = (int)std::min((double)height, std::floor(mouthLimit/(2.0*dy)));
    const int w = width+2*rx+2, h = height+2*ry+2;
    const int x0 = rx+1, y0 = ry+1;
    std::vector<unsigned char> original((size_t)w*h, 0);
    for(int y=0; y<height; ++y)
        for(int x=0; x<width; ++x)
            original[(size_t)(y+y0)*w+x+x0] = occ[(size_t)y*width+x] != 0;

    // 只在副本上临时封住窄口，再查哪些空区由外部变成了内孔。
    // 外轮廓不使用闭运算结果，避免把开阔凹口或不同面之间的外部空区填掉。
    auto closeMask = [&](int radX, int radY){
        auto filter = [&](const std::vector<unsigned char>& src, int radius, bool vertical, bool dilate){
            auto dst = src;
            if(radius == 0) return dst;
            const int lines = vertical ? w : h, length = vertical ? h : w;
            std::vector<int> sums(length+1, 0);
            for(int j=0; j<lines; ++j){
                auto index = [&](int i){ return vertical ? (size_t)i*w+j : (size_t)j*w+i; };
                for(int i=0; i<length; ++i) sums[i+1] = sums[i]+src[index(i)];
                for(int i=0; i<length; ++i){
                    const int sum = sums[std::min(length,i+radius+1)]-sums[std::max(0,i-radius)];
                    dst[index(i)] = dilate ? sum>0 : sum==2*radius+1;
                }
            }
            return dst;
        };
        auto mask = filter(original, radX, false, true);
        mask = filter(mask, radY, true, true);
        mask = filter(mask, radX, false, false);
        return filter(mask, radY, true, false);
    };
    auto visit = [&](int p, auto&& fn){
        const int x=p%w, y=p/w;
        if(x>0) fn(p-1);
        if(x+1<w) fn(p+1);
        if(y>0) fn(p-w);
        if(y+1<h) fn(p+w);
    };
    int filled = 0;
    std::vector<unsigned char> originalExterior;
    const std::array<std::pair<int,int>,3> radii{{{0,0}, {(rx+1)/2,(ry+1)/2}, {rx,ry}}};
    std::pair<int,int> previous{-1,-1};
    for(const auto& radius : radii){
        if(radius == previous) continue;
        previous = radius;
        const auto mask = closeMask(radius.first,radius.second);
        std::vector<unsigned char> seen(mask.size(),0);
        std::vector<int> cells{0};
        seen[0]=1;
        auto addEmpty = [&](int p){
            if(!mask[p] && !seen[p]){ seen[p]=1; cells.push_back(p); }
        };
        for(size_t i=0; i<cells.size(); ++i) visit(cells[i],addEmpty);
        if(radius.first==0 && radius.second==0) originalExterior = seen;
        for(int y=y0; y<y0+height; ++y) for(int x=x0; x<x0+width; ++x){
            const int start = y*w+x;
            if(mask[start] || seen[start]) continue;
            cells.assign(1,start);
            seen[start]=1;
            int left=x, right=x, bottom=y, top=y;
            for(size_t i=0; i<cells.size(); ++i){
                const int p=cells[i];
                left=std::min(left,p%w); right=std::max(right,p%w);
                bottom=std::min(bottom,p/w); top=std::max(top,p/w);
                visit(p,addEmpty);
            }
            const bool open = radius.first || radius.second;
            if(open && !originalExterior[start]) continue;
            if(open){
                // 补上临时封口新增的格；只从孔内向窄口长，且严格限定在孔周围。
                const int xl=left-2*radius.first, xr=right+2*radius.first;
                const int yl=bottom-2*radius.second, yr=top+2*radius.second;
                for(size_t i=0; i<cells.size(); ++i) visit(cells[i],[&](int p){
                    const int px=p%w, py=p/w;
                    if(seen[p] || !mask[p] || original[p] || !originalExterior[p] ||
                       px<xl || px>xr || py<yl || py>yr) return;
                    seen[p]=1; cells.push_back(p);
                    left=std::min(left,px); right=std::max(right,px);
                    bottom=std::min(bottom,py); top=std::max(top,py);
                });
            }
            const double spanX = (right-left+1)*dx, spanY = (top-bottom+1)*dy;
            if(std::max(spanX,spanY) > maxSpan) continue;
            if(rimSupport){
                // 曲面投影里的空区也可能是被剥离的腔底。只有边缘有足够朝向一致的
                // 实测表面支持，才在该基面内补洞。
                int observed=0, facing=0;
                for(int cell : cells) visit(cell,[&](int p){
                    const int px=p%w-x0, py=p/w-y0;
                    if(!original[p] || px<0 || px>=width || py<0 || py>=height) return;
                    const auto v=(*rimSupport)[(size_t)py*width+px];
                    if(v==0) return;
                    ++observed;
                    if(v==2) ++facing;
                });
                if(observed==0 || 2*facing<observed) continue;
            }
            const int keep = filled;
            for(int p : cells){
                const int px=p%w-x0, py=p/w-y0;
                if(px<0 || px>=width || py<0 || py>=height) continue;
                auto& cell = occ[(size_t)py*width+px];
                if(cell==0){ cell=2; ++filled; }
            }
            if(filled != keep)
                fprintf(stderr, "[mesh] %s孔洞填补: %.1f×%.1fmm，阈值 %.1fmm，填补 %d 格\n",
                        open ? "开口" : "闭合",
                        spanX, spanY, maxSpan, filled-keep);
        }
    }
    return filled;
}

static bool build_plane_mesh(const std::vector<Vec2>& uv, double reconGap, double fillHole,
                             double minFaceWidth,
                             std::vector<Vec2>& outVerts, std::vector<int>& outTris,
                             int maxVerts, int maxTris,
                             std::vector<char>* keepMask = nullptr,
                             std::vector<Vec2>* fillRectOut = nullptr,
                             std::vector<std::vector<Vec2>>* loopsOut = nullptr)
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
        int rr=(int)std::ceil(std::max(reconGap,1.0)/(2.0*cs))+1;
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

    fill_reconstruction_holes(occ, nU, nV, cs, cs, fillHole);

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
    // 只比较外环面积。小而紧凑的孤岛用更严格的面积门槛，细长折边仍按原门槛。
    double maxLoopArea = 0.0;
    for(const auto& lp : loops) maxLoopArea = std::max(maxLoopArea, polyArea2(lp));
    std::vector<const std::vector<Vec2>*> droppedIslands;

    for(auto& lp : loops){
        double a0=polyArea2(lp);
        if(dbg) fprintf(stderr, "[mesh]   环: 点=%zu 面积=%.0f\n", lp.size(), a0);
        if(a0 > 0.0){
            bool drop = a0 < 9.0*cs*cs || a0 < 0.10*maxLoopArea;
            if(!drop && minFaceWidth > 0.0 && a0 < 0.25*maxLoopArea){
                double longSide = 0.0;
                minRectShortSide(lp, &longSide);
                drop = longSide < 4.0*minFaceWidth;
            }
            if(drop){
                droppedIslands.push_back(&lp);
                fprintf(stderr, "[mesh] 剔除孤立碎块: %.0fmm²（占主块 %.1f%%）\n",
                        a0, 100.0*a0/maxLoopArea);
                continue;
            }
        }
        // 直线简化的容差取"重建几何尺度"：线扫点云边界天然是"梳齿"状
        // （每行末端位置不同，齿深约等于扫描线间距），容差小于齿深就会把梳齿保留
        // 下来，看起来像锯子。取 max(2×cell, 1.5×重建尺度) 把梳齿抹成直边。
        std::vector<Vec2> s = simplifyClosed(lp, std::max(2.0*cs, 1.5*reconGap));
        double a=polyArea2(s);
        if(std::fabs(a)<1e-6) continue;

        // 凸包规整：相邻面在折棱处会互相混入几个点，在轮廓上咬出浅缺口。
        // 只有当"轮廓相对凸包的最大凹陷深度 <= max(孔洞填补阈值, 5×重建尺度)"时
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
                // 容差取 max(5×重建尺度, 孔洞填补阈值)：重建尺度决定了"扫描线间混入"的
                // 量级，边界的浅缺口都在这个尺度内；真实大缺口（长板宽板的中间
                // 残缺，深 100mm）超出容差，仍然保留原轮廓。
                double devTol = std::max(5.0*reconGap, fillHole);
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
    if(keepMask && !droppedIslands.empty()){
        for(size_t k=0; k<uv.size(); ++k){
            if(!(*keepMask)[k]) continue;
            for(const auto* island : droppedIslands){
                if(pointInPoly(uv[k], *island)){ (*keepMask)[k]=0; break; }
            }
        }
    }

    // 6.5) 窄条面矩形补全
    // 场景：折型板的折边带只有 ~43mm 宽、700mm 长，激光对它是掠射，点云沿长向
    // 断续、短向只在棱线附近密。形态学接不上断口（闭运算半径已开到上限），
    // 边界追踪出来是一条 4~8mm 的细带 + 若干孤立点，重建覆盖只有外接矩形的 40%。
    // 但工件本身是连续的（只是扫描缺陷），工程上必须按条带补齐，否则那段没有面、
    // 也就没有轨迹。判据三条同时满足才补：
    //   a) 短向跨度 <= 4×重建尺度（确实是"窄条"，不是大平面；大平面的大缺口必须保留）
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
        if(wShort > 4.0 && wShort <= 4.0*reconGap && box > 1.0 && wLong > 2.0*wShort
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
        // 输出"重建轮廓多边形"（已是世界无关的 2D 局部坐标），供轨迹阶段把行裁剪
        // 在重建面内——这是抑制"轨迹出框"最直接的一环：面的轮廓已经算出来了，
        // 布行却仍按 P0.5/P99.5 包围盒走，等于算完不用。
        if(loopsOut) loopsOut->push_back(poly);
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

// 点云局部噪声尺度 σ（mm）：对每个抽样点，用其 K 邻域解一个局部平面，把"邻点到该
// 平面的平均残差"换算成该点的局部 σ（高斯下 平均|r| = 0.7979·s，s = σ·sqrt(1−3/K)，
// K=25 时合计 0.7484σ）；再取全体局部 σ 的 **25% 分位** 作为整片点云的噪声尺度。
//
// 为什么用 25% 分位而不是中位数：薄壁件上"绝大多数点都挨着棱边"（壁厚 2mm、高 300mm，
// 几乎每个点的邻域都跨折角），中位数被折角残差污染——实测 g_rim/g_norim 中位数给
// 2.4mm，而它们的壁面拟合 RMS 只有 0.02mm。25% 分位落在真正的平坦区，才代表噪声。
//
// 为什么拟合尺度不能跟着用户的大阈值走：
//   RANSAC 的投票是"到假设平面距离 < tol 就记一票"。容差放大到几十上百毫米后，
//   任何一个假设的内点数都≈全部点，打分失去区分度，选出来的假设等价于随机抽的
//   三点 —— 这不是"拟合更宽松"，是"RANSAC 坏掉"。
//   最终实现里拟合尺度固定 3mm（试过按 σ 自动放大，薄壁件 4mm 就把主面翻到侧壁，
//   见 split_planes 中的记录）；σ 只作为诊断量上报，帮助判断"点云本身是不是太差"。
//   用户想要的"过滤大点"属于**吸收/覆盖**的语义，由 absorbTol 承担。
static double cloud_noise_sigma(const std::vector<Vec3>& P, const CloudNormals& cn)
{
    const int N = (int)P.size();
    if(N < 60) return 0.5;
    const int K = CloudNormals::K;
    const int maxSamples = 3000;
    const int stride = std::max(1, N / maxSamples);
    std::vector<double> res;
    res.reserve((size_t)(N/stride) + 2);
    for(int i=0; i<N; i+=stride){
        const int* row = &cn.knn[(size_t)i*(size_t)K];
        Vec3 nf, cf;
        if(!ls_plane_p(P, row, K, nf, cf)) continue;
        double s = 0.0;
        for(int t=0;t<K;++t){
            const int j = row[t];
            if(j < 0) continue;
            s += std::fabs(dot(P[j]-cf, nf));
        }
        // 单点局部 σ：K 个邻点的平均 |r| = 0.7979·s，s = σ·sqrt(1−3/K) ≈ 0.938σ
        // ⇒ σ_local = 平均|r| / 0.7484
        res.push_back(s/(double)K/0.7484);
    }
    if(res.size() < 8) return 0.5;
    size_t q = res.size()/4;
    std::nth_element(res.begin(), res.begin()+q, res.end());
    double sigma = res[q];
    if(!(sigma > 0.0)) sigma = 0.5;          // 完美平面件：给一个数值下限
    return std::min(3.0, sigma);
}

// 沿射线 (o,dir) 找最近命中距离；量程内无命中返回 -1。
// 为什么不用"最近邻居距离"：同面邻点永远最近（沿射线投影≈0、横向 1mm），
// 会把远处的真命中挡住。所以按球收集候选，再按"沿射线投影 + 横向偏移"筛。
static double ray_free_dist(const std::vector<Vec3>& P, const KdTree& tree,
                            Vec3 o, Vec3 dir, double maxDist,
                            double hitR, double collR){
    const double minT = std::max(0.8, hitR*0.6);   // 跳过起点处的同面点
    int cand[96];
    double t = minT;
    while(t <= maxDist){
        Vec3 q = o + dir*t;
        int n = tree.radiusCollect(q, collR, cand, 96);
        double best = 1e30;
        for(int i=0;i<n;++i){
            Vec3 v = P[cand[i]] - o;
            double tp = dot(v, dir);
            if(tp < minT || tp > maxDist) continue;
            double perp2 = dot(v,v) - tp*tp;
            if(perp2 > hitR*hitR) continue;
            if(tp < best) best = tp;
        }
        if(best < 1e29) return best;
        // 近段密（判"材料侧"要 <12mm 的精度），远段疏（只需知道有没有腔体对面）。
        // 步长必须 < collR，否则球与球之间会漏。
        t += (t < 20.0) ? 1.5 : 4.5;
    }
    return -1.0;
}

// 面级朝向：射线自由空间判材料侧，结果只用于上报（不改轨迹）。
// src 说明置信度：1=射线判定，2=两侧都自由（单层壳，几何二义，按原定向兜底），
// 3=用户强制翻转。
static FaceOrient classify_face_orientation(const std::vector<Vec3>& P, const KdTree& tree,
                                            const PlaneSeg& seg, const Vec3& partC,
                                            int forceFlip){
    FaceOrient fo;
    Vec3 n = seg.pl.n;

    const double maxFree  = 800.0; // 量程要覆盖腔体对面（400×300 的盒，对面在 300~400mm）
    const double nearWall = 12.0;  // <12mm 命中 = 贴材料（薄壁另一层）；更远 = 腔体对面
    const double hitR = 1.5, collR = 5.0;

    // 面内均匀取样（上限 6 点）：中位数抗单点噪声；点数多了射线开销线性涨
    std::vector<double> dp, dm;
    const size_t N = seg.inl.size();
    if(N){
        int want = (int)std::min<size_t>(6, N);
        dp.reserve(want); dm.reserve(want);
        for(int i=0;i<want;++i){
            size_t k = (size_t)((double)i/(double)want*(double)N);
            if(k >= N) k = N-1;
            Vec3 p = P[seg.inl[k]];
            dp.push_back(ray_free_dist(P, tree, p + n*0.4,  n,          maxFree, hitR, collR));
            dm.push_back(ray_free_dist(P, tree, p - n*0.4, n*(-1.0),   maxFree, hitR, collR));
        }
        auto med = [](std::vector<double>& v){
            if(v.empty()) return -1.0;
            std::nth_element(v.begin(), v.begin()+v.size()/2, v.end());
            return v[v.size()/2];
        };
        fo.freePlus  = med(dp);
        fo.freeMinus = med(dm);
    }

    // 三态分类：0 = 材料（近挡，薄壁另一层）/ 1 = 腔体（远挡，对面被围住）/ 2 = 外部（无挡）
    // 可喷侧 = 背离材料的一侧：
    //   +n 是材料           → 翻到 -n
    //   -n 是材料           → 保持 +n
    //   +n 是腔体、-n 自由   → 保持（朝腔内）
    //   -n 是腔体、+n 自由   → 翻到 -n
    //   两侧都自由 / 都材料  → 兜底（单层壳二义：外侧没扫到，几何上判不出来）
    auto cat = [&](double d){
        if(d < 0) return 2;
        return d < nearWall ? 0 : 1;
    };
    const int cp = cat(fo.freePlus), cm = cat(fo.freeMinus);
    int src = 2;
    if(forceFlip){ src = 3; }                             // 用户强制翻转，以当前 n 为准
    else if(cp == 0 && cm != 0){ n = n*(-1.0); src = 1; } // +n 贴材料
    else if(cm == 0 && cp != 0){ src = 1; }               // -n 贴材料
    else if(cp == 1 && cm == 2){ src = 1; }               // +n 是腔体
    else if(cm == 1 && cp == 2){ n = n*(-1.0); src = 1; } // -n 是腔体
    else if(cp == 1 && cm == 1){                          // 两侧都是腔体：取更开阔的一侧
        if(fo.freeMinus > fo.freePlus + 1.0) n = n*(-1.0);
        src = 1;
    }
    else if(cp == 0 && cm == 0){                          // 两侧都贴材料
        // 只有一种几何能同时满足"两侧都在 nearWall 内命中"：薄壁双层。
        // 这个面块里同时装着壁的两层（内壁 + 外壁），法向相反，射线沿两侧
        // 各走一个壁厚就撞到对面层，所以两侧自由距离都 ≈ 壁厚、完全对称。
        // 这不是"判不出来"，而是"这一块里有两张朝向相反的面"，必须标出来，
        // 不能像以前那样按 0 默认成外壁——否则内壁会整层消失。
        src = 4;
    }

    fo.n   = n;
    fo.src = src;
    if(std::getenv("SPRAYTRAJ_FACE_DBG"))
        fprintf(stderr, "[face] c=(%.1f,%.1f,%.1f) n_in=(%.3f,%.3f,%.3f) dp=%.1f dm=%.1f n_out=(%.3f,%.3f,%.3f) src=%d\n",
                seg.pl.c.x, seg.pl.c.y, seg.pl.c.z,
                seg.pl.n.x, seg.pl.n.y, seg.pl.n.z,
                fo.freePlus, fo.freeMinus, n.x, n.y, n.z, src);
    const double c = std::max(-1.0, std::min(1.0, n.z));
    fo.angZ = std::acos(c)*180.0/3.141592653589793;      // 0~180
    fo.tilt = 90.0 - fo.angZ;                            // -90~90
    fo.side = (fo.angZ < 87.0) ? 1 : (fo.angZ > 93.0 ? -1 : 0);
    // isInner 的语义是"可喷侧朝腔体"，不能只看质心：
    // 折型板的翻边朝下、质心在它上方，按质心会把外表面误标成内壁。
    // 优先用射线结论：可喷侧是"远挡（腔体对面）"才算内壁；
    // 两侧都自由（单层壳）时才退回质心判。
    if(cp == 1 && cm == 2)      fo.isInner = 1;   // +n 侧是腔体
    else if(cp == 2 && cm == 1) fo.isInner = 1;   // -n 侧是腔体
    else if(cp == 2 && cm == 2) fo.isInner = (dot(n, seg.pl.c - partC) < 0) ? 1 : 0;
    else if(cp == 1 && cm == 1) fo.isInner = (dot(n, seg.pl.c - partC) < 0) ? 1 : 0;
    // 薄壁双层：块内同时有内壁和外壁两层，单值 isInner 表达不了。
    // 返回 -1（未分），让上层决定怎么拆；绝不能默认 0 当成外壁。
    else if(cp == 0 && cm == 0) fo.isInner = -1;
    else                        fo.isInner = 0;
    return fo;
}

// 拆层后单片定向：拆层切出来的是"沿法向有间隙的两片"，而**间隙不等于壁的两侧**
// ——台阶/两级平台/凸台同样有间隙。所以两片的符号必须各自重判一次，
// 不能再用"另一片取反"这种构造：一旦把同向的两片当成壁，大片会被翻成 -Z，
// 再被可喷窗口（>130° 不喷）静默丢掉一整片（实测煤气灶背板丢了 46472 点/16%）。
// 口径与主定向一致：先归 +Z 半球（近竖直背离质心），射线能判出可喷侧时以射线为准。
//
// ---------- 步骤 1：初始定向（与主定向同口径） ----------
// ---------- 步骤 2：射线复核，只有 src==1（一侧贴材料/一侧朝腔）才采信 ----------
// ---------- 步骤 3：重算 angZ / tilt / side / isInner ----------
static void orient_piece_after_split(const std::vector<Vec3>& P, const KdTree& tree,
                                     PlaneSeg& s, const Vec3& partC)
{
    // ---------- 步骤 1：初始定向（与主定向同口径） ----------
    Vec3 n = s.pl.n;
    if(std::fabs(n.z) > 0.05){ if(n.z < 0) n = n*(-1.0); }
    else if(dot(n, s.pl.c - partC) < 0) n = n*(-1.0);   // 近竖直：背离质心
    s.pl.n = n;
    s.pl.d = -dot(n, s.pl.c);

    // ---------- 步骤 2：射线复核可喷侧 ----------
    // src==1：一侧近挡（材料）/一侧远挡（腔体），符号可信，以它为准；
    // src==2/4：两侧都自由或两侧都贴材料，几何二义，保留步骤 1 的结果。
    const FaceOrient fo = classify_face_orientation(P, tree, s, partC, 0);
    if(fo.src == 1 && dot(fo.n, s.pl.n) < 0.0){
        s.pl.n = fo.n;
        s.pl.d = -dot(s.pl.n, s.pl.c);
    }
    s.fo = fo;
    s.fo.n  = s.pl.n;
    s.fo.src = 4;                       // 来源标记仍记"薄壁拆层"，界面按它显示内外壁

    // ---------- 步骤 3：重算角度与分类 ----------
    const double D2R = 3.141592653589793/180.0;
    s.fo.angZ = std::acos(std::max(-1.0, std::min(1.0, s.pl.n.z)))/D2R;
    s.fo.tilt = 90.0 - s.fo.angZ;
    s.fo.side = (s.fo.angZ < 87.0) ? 1 : (s.fo.angZ > 93.0 ? -1 : 0);
    // 可喷侧朝腔（指向整件质心）= 内壁；这一片已经单独定过向，
    // 不能再沿用源块的 isInner（那是两层混在一起时的结论）。
    s.fo.isInner = (dot(s.pl.n, partC - s.pl.c) > 0.0) ? 1 : 0;
}

// ==================== 枪位物理验收：最后一道方向过滤（不可绕过） ====================
// 为什么角度窗口挡不住"撞机"：
//   窗口只说"法向与 +Z 的夹角 ≤ 130°"，它不知道枪位（X + standoff·n）外有没有东西。
//   拆层/双面拆出来的假面，angZ 可能只有 91°（窗口完全合法），但那侧 40mm 就是隔壁
//   那片料——枪伸不进去，走过去就是贴件/撞机。反过来，166° 也不一定真朝下。
// 所以这里**不看 angZ、不看 sideSrc、不看 isInner**，只信点云回答一个问题：
//   "从这一面的点，沿法向抬到 standoff，这一路是不是空的？"
// 三种结论：
//   OK     —— +n 侧站得住枪，按当前 n 布行；
//   FLIP   —— +n 侧站不住、-n 侧站得住，说明符号反了，翻过来（这才是真正的方向过滤）；
//   REJECT —— 两侧都站不住枪 → 拒喷。宁可漏喷，绝不撞机。
// 阈值取 0.85·standoff 而不是 standoff：拟合残差、采样点偏离平面、去噪删掉棱边点
// 都会让实测自由距离比真值小几毫米，留 15% 余量避免误杀正常面。
// 已经超出可喷窗口（angZ > 130°）的面**不给翻转洗白的机会**，直接 REJECT：
// 它本来就不该喷，翻一下变成 0° 再喷，正是"朝下被喷"的事故来源。
enum GunSideVerdict { GUN_SIDE_OK = 0, GUN_SIDE_FLIP = 1, GUN_SIDE_REJECT = 2 };

// ---------- 步骤 1：面内均匀取样，实测 ±n 两侧自由距离 ----------
// 重新实测而不复用 fo.freePlus/freeMinus：那组值是按 classify 当时的 n 算的，
// 后面"内壁喷距法向翻转朝腔"（pl.n = fo.n）之后，两侧标签已经对不上当前 n。
//
// localNor（可选）：逐点局部法向（CloudNormals::nor，下标与 P 一致）。
// **曲面块必须传**。为什么：曲面块只有一个"整块拟合的平均平面法向" pl.n，
// 而曲型板的局部法向与它相差可达几十度 —— 沿平均法向从面上的一点打射线，
// 等于贴着曲面擦过去，1~3mm 就打回曲面自身，于是被判"枪埋在材料里 → 两侧都
// 伸不进枪 → 整块拒喷"。实测曲型板（440×546×312，弧高 300）：
//   沿平均法向 +n 中位 1.2mm（±1~3）  → 拒喷，93% 点数一行不出
//   沿逐点法向 12/12 采样点 400mm 内无命中 → 空间完全自由
// 平面块的两者一致，所以默认走 pl.n，行为与原来逐位相同。
// 冻结平面的轮廓点用于计算覆盖范围。按数组下标取样会随吸收、轮廓回筛而移动，
// 甚至集中在同一边缘。按面内位置分区，从每个有点的区域取最接近中心的实测点。
static std::vector<int> sample_fitted_face(const std::vector<Vec3>& points,
                                           const std::vector<int>& indices, const Plane& plane)
{
    const Vec3 n = plane.n;
    const Vec3 u = norm(cross(n, std::fabs(n.z) < 0.9 ? Vec3{0,0,1} : Vec3{1,0,0}));
    const Vec3 v = cross(n, u);
    double uMin = 1e300, uMax = -1e300, vMin = 1e300, vMax = -1e300;
    for(int index : indices){
        const Vec3 p = points[index] - plane.c;
        const double pu = dot(p, u), pv = dot(p, v);
        uMin = std::min(uMin, pu); uMax = std::max(uMax, pu);
        vMin = std::min(vMin, pv); vMax = std::max(vMax, pv);
    }

    // 沿长边分三格、短边分两格，每侧最多六条射线。
    const int uCells = (uMax - uMin >= vMax - vMin) ? 3 : 2;
    const int vCells = 6 / uCells;
    const double du = std::max((uMax - uMin) / uCells, 1e-9);
    const double dv = std::max((vMax - vMin) / vCells, 1e-9);
    std::array<int, 6> nearest;
    std::array<double, 6> distance2;
    nearest.fill(-1);
    distance2.fill(1e300);
    for(int index : indices){
        const Vec3 p = points[index] - plane.c;
        const double pu = dot(p, u) - uMin, pv = dot(p, v) - vMin;
        const int col = std::clamp((int)(pu / du), 0, uCells - 1);
        const int row = std::clamp((int)(pv / dv), 0, vCells - 1);
        const int cell = row * uCells + col;
        const double x = pu - (col + 0.5) * du;
        const double y = pv - (row + 0.5) * dv;
        const double h = dot(p, n);
        const double d2 = x*x + y*y + h*h;
        if(d2 < distance2[cell]){
            nearest[cell] = index;
            distance2[cell] = d2;
        }
    }

    std::vector<int> samples;
    for(int index : nearest)
        if(index >= 0) samples.push_back(index);
    return samples;
}

static void measure_gun_free(const std::vector<Vec3>& P, const KdTree& tree,
                             const PlaneSeg& s, double standoff,
                             double& freePlus, double& freeMinus,
                             const std::vector<Vec3>* localNor = nullptr)
{
    const double maxFree = std::max(400.0, standoff*4.0);
    const double hitR = 1.5, collR = 5.0;
    // 采样集：合并过的面用 core（合并前的本面点），否则用 inl。
    // 为什么不能直接用 inl（2026-11-23 实测回归）：merge 把 83~88° 的侧壁并进主面后，
    // inl 里混着侧壁点；沿**主面法向**从侧壁点打射线 1.1mm 就撞到侧壁本体，
    // 中位数被拉穿，整张主面判"枪伸不进去"而拒喷（p001 主面 8433 点、p009 15631 点）。
    // 枪位验收要回答的是"从这张面自己的表面上打出去空不空"，侧壁不回答这个问题。
    const std::vector<int>& candidates = s.core.empty() ? s.inl : s.core;
    // 合并面已有独立于轮廓的验收点集，继续使用原有采样。
    const auto spatialSamples = s.fixed && s.core.empty() && !localNor && candidates.size() > 6
                              ? sample_fitted_face(P, candidates, s.pl) : std::vector<int>{};
    const std::vector<int>& samp = spatialSamples.empty() ? candidates : spatialSamples;
    const double N = (double)samp.size();
    freePlus = freeMinus = -1.0;
    if(samp.empty()) return;

    const int want = (int)std::min<size_t>(6, samp.size());
    std::vector<double> dp, dm;
    dp.reserve((size_t)want); dm.reserve((size_t)want);
    const Vec3 n = s.pl.n;
    for(int i=0;i<want;++i){
        size_t k = (size_t)((double)i/(double)want*N);
        if(k >= samp.size()) k = samp.size()-1;
        const Vec3 p = P[(size_t)samp[k]];
        Vec3 d = n;
        if(localNor){
            const Vec3 ln = (*localNor)[(size_t)samp[k]];
            if(len(ln) > 1e-9){
                d = norm(ln);
                // PCA 法向无符号：与面级可喷侧对齐，否则 "+n/-n" 会各自乱跳，
                // 中位数就失去意义。
                if(dot(d, n) < 0.0) d = d*(-1.0);
            }
        }
        dp.push_back(ray_free_dist(P, tree, p + d*0.4,  d,         maxFree, hitR, collR));
        dm.push_back(ray_free_dist(P, tree, p - d*0.4,  d*(-1.0),  maxFree, hitR, collR));
        if(std::getenv("SPRAYTRAJ_FACE_DBG"))
            fprintf(stderr, "[gun-sample] idx=%d p=(%.2f,%.2f,%.2f) h=%.2f +n=%.2f -n=%.2f\n",
                    samp[k], p.x, p.y, p.z, dot(p, n) + s.pl.d, dp.back(), dm.back());
    }
    auto med = [](std::vector<double>& v){
        if(v.empty()) return -1.0;
        std::nth_element(v.begin(), v.begin()+v.size()/2, v.end());
        return v[v.size()/2];
    };
    freePlus  = med(dp);
    freeMinus = med(dm);
}

// ---------- 步骤 2：判定，需要时把符号翻过来 ----------
// allowFlip=false 用于曲面块：曲面行用的是逐点法向 surf.N，翻面级平均法向没有意义，
// 只能"站得住就喷、站不住就拒"，不能靠翻符号洗白。
static GunSideVerdict verify_gun_side(const std::vector<Vec3>& P, const KdTree& tree,
                                      PlaneSeg& s, const Vec3& partC, double standoff,
                                      const char* tag, bool allowFlip,
                                      const std::vector<Vec3>* localNor = nullptr)
{
    double fp = -1.0, fm = -1.0;
    measure_gun_free(P, tree, s, standoff, fp, fm, localNor);
    // -1 = 量程内无命中 = 完全自由，按大值参与比较
    const double effPlus  = (fp < 0.0) ? 1e9 : fp;
    const double effMinus = (fm < 0.0) ? 1e9 : fm;
    const double need = standoff*0.85;

    s.fo.freePlus  = fp;      // 顺手修正上报口径（与当前 pl.n 一致）
    s.fo.freeMinus = fm;

    const double D2R = 3.141592653589793/180.0;
    auto reAngle = [&](){
        s.fo.n    = s.pl.n;
        s.fo.angZ = std::acos(std::max(-1.0, std::min(1.0, s.pl.n.z)))/D2R;
        s.fo.tilt = 90.0 - s.fo.angZ;
        s.fo.side = (s.fo.angZ < 87.0) ? 1 : (s.fo.angZ > 93.0 ? -1 : 0);
        s.fo.isInner = (dot(s.pl.n, partC - s.pl.c) > 0.0) ? 1 : 0;
    };
    auto flip = [&](){
        s.pl.n = s.pl.n*(-1.0);
        s.pl.d = -dot(s.pl.n, s.pl.c);
        reAngle();
        std::swap(s.fo.freePlus, s.fo.freeMinus);
    };

    // ---------- 情形 A：+n 侧站得住枪 ----------
    if(effPlus >= need) return GUN_SIDE_OK;

    // ---------- 情形 B：+n 站不住 → 先确认它是不是"本就不该喷" ----------
    // angZ 已超窗口（朝下）→ 不给翻转洗白，直接拒。
    if(s.fo.angZ > SPRAYTRAJ_NORMAL_WIN_DEG){
        if(std::getenv("SPRAYTRAJ_FACE_DBG"))
            fprintf(stderr, "[gun-gate] %s angZ=%.1f° 超窗口且 +n 侧仅 %.1fmm，拒喷（不翻转洗白）\n",
                    tag, s.fo.angZ, fp);
        return GUN_SIDE_REJECT;
    }

    // ---------- 情形 C：-n 侧站得住 → 符号反了，翻过来 ----------
    if(allowFlip && effMinus >= need){
        flip();
        if(s.fo.angZ > SPRAYTRAJ_NORMAL_WIN_DEG){   // 翻完更朝下 → 也不许
            flip();                                  // 还原，按拒喷处理
            return GUN_SIDE_REJECT;
        }
        s.fo.src = 7;                                // 7 = 枪位验收翻转（几何自证，非标签）
        if(std::getenv("SPRAYTRAJ_FACE_DBG"))
            fprintf(stderr, "[gun-gate] %s 翻转: +n 侧仅 %.1fmm、-n 侧 %.1fmm → angZ 现 %.1f°\n",
                    tag, fp, fm, s.fo.angZ);
        return GUN_SIDE_FLIP;
    }

    // ---------- 情形 D：两侧都伸不进枪 ----------
    if(std::getenv("SPRAYTRAJ_FACE_DBG"))
        fprintf(stderr, "[gun-gate] %s 拒喷: 两侧自由 +n=%.1fmm / -n=%.1fmm 均 < %.0fmm（枪伸不进去）\n",
                tag, fp, fm, need);
    return GUN_SIDE_REJECT;
}

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
                                   double tol, double curvatureTol, int minPts, double rowSpacing,
                                   double standoff){
    std::vector<PlaneSeg> out;

    // ---- 两个尺度必须分开，原来它们共用一个 tol，是"调大就塌"的根 ----
    //
    //  fitTol    = 拟合尺度：只决定"哪些点参与解平面方程"，由点云噪声定（3σ）。
    //              它的作用只有精度，不承担任何工艺语义。见 cloud_noise_sigma 的注释：
    //              这个数一旦跟着用户的大阈值走，RANSAC 的投票就退化成随机抽样。
    //  absorbTol = 吸收尺度：用户的 flatnessTol 真正要控制的东西 —— 离这张面多远的点
    //              仍然算"这张面覆盖到"（用于轮廓/覆盖/残点归并），但**不参与解方程**。
    //
    //  absorbTol 的上限不是魔数，是净空预算：枪位 = 面点 + n·standoff。若某个被吸收的点
    //  高出平面 h，该处真实净空就是 standoff − h。取 h ≤ 0.6·standoff ⇒ 最坏仍有
    //   0.4·standoff（standoff=100 时 40mm）的余量，大于默认 gunClear(20)。
    //   所以"吸收到 0.6×standoff"这条约束本身就保证了不穿工件，不需要再去猜 mm 数。
    const double sigma = cloud_noise_sigma(P, cn);
    // 拟合尺度：默认固定 3mm。它**不是**"用户容差的上限"，而是解平面方程时的内点尺度，
    // 与 flatnessTol 的工艺语义无关（flatnessTol 现在只管吸收，见上）。
    //
    // 为什么不按噪声 3σ 自动放大（试过，已回退）：实测 g_rim（300mm 薄壁盒）在 3.5mm
    // 还正常（226 行），到 4.0mm 主面就从"朝上的开口面"翻成"侧壁"，面块 7→3、
    // 行数 226→76。RANSAC 选"哪张面是主面"对容差极敏感——容差一松，面积大的侧壁
    // 就抢过面积小的开口面。所以这一层不能跟着噪声走，只能钉住。
    // 真要改，用 SPRAYTRAJ_FITTOL 强制覆盖（排查用），默认不碰。
    const double fitTolEnv = std::getenv("SPRAYTRAJ_FITTOL")
                           ? atof(std::getenv("SPRAYTRAJ_FITTOL")) : -1.0;
    const double fitTol = (fitTolEnv > 0.0) ? fitTolEnv : 3.0;
    const double absorbTol = std::min(tol, std::max(2.0, 0.6*standoff));
    fprintf(stderr, "[split] 尺度: 噪声σ=%.2fmm 拟合fitTol=%.2fmm 吸收absorbTol=%.2fmm"
                    "（用户flatTol=%.1f 净空预算=%.1f）\n",
            sigma, fitTol, absorbTol, tol, 0.6*standoff);
    if(tol > absorbTol + 1e-9)
        fprintf(stderr, "[split] 注意: flatTol %.1f 超过净空预算 0.6×standoff=%.1f，"
                        "吸收按 %.1f 执行（再大就会顶进材料）\n", tol, 0.6*standoff, absorbTol);

    // ---- 第一步：RANSAC 拟合主平面，确定高度方向 n ----
    Plane mainPl; std::vector<int> mainInl;
    bool ransacOk;
    { ProfStage ps("2c.RANSAC主平面"); ransacOk = fit_plane(P, fitTol, mainPl, mainInl); }
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

    // ---- 吸收集：平面方程到此已定死，后面不再重拟合 ----
    // 这是"过滤大点但平面不歪"的关键。原来 mergeClusters 每吸收一块就立刻
    // 用合并后的点重解平面，于是"吸收"必然"拉平面"；而且残差是对**重拟合后的**
    // 平面算的，可以链式爬升——c.ply 在 flat=40 时把相距 105mm 的两张平行面
    // 并成一张 RMS=50mm 的假面，就是这么来的。
    // 现在吸收只扩 seg.inl（轮廓/覆盖用），不参与解方程，所以吸收多少都不会歪。
    std::vector<int> mainAbs;
    std::vector<char> mainLock(P.size(), 0);   // 1 = 已被主面吸收，后续阶段不得再分配
    {
        // 吸收两条：①朝向与主面一致（含反向平行，薄壁另一层要能被吸收）②高度在容差内。
        // 朝向这条不能省：p030/p040/p036 这类件有一张与主面夹角 ~74° 的大侧面
        // （5000+ 点），它在主面法向上的投影高度同样 < absorbTol，只按高度吸收会
        // 把它整片吞进主面 —— 结果就是那张面再也不生成轨迹，等于漏喷一整面。
        // 而 p002 那两条 73°/55° 的"假侧壁"（403/237 点、距主面 p95=6mm）是主面
        // 自己的边缘带被法向估计拉偏，它们的特征是**整块贴在主面上**，
        // 交给最终合并闸门的"贴合"分支处理，不靠吸收来兜。
        const double cosAbsN = std::cos(40.0*3.141592653589793/180.0);
        mainAbs.reserve(P.size());
        double hMax = 0.0;
        int nAngRej = 0;
        for(size_t i=0;i<P.size();++i){
            if(std::fabs(dot(cn.nor[i], n)) < cosAbsN){ ++nAngRej; continue; }
            const double h = std::fabs(dot(P[i], n) + mainPl.d);
            if(h <= absorbTol){ mainAbs.push_back((int)i); mainLock[i] = 1; hMax = std::max(hMax, h); }
        }
        fprintf(stderr, "[split] 主面: 拟合内点=%zu 吸收内点=%zu（朝向不符跳过 %d）吸收最大高差=%.1fmm"
                        " ⇒ 最坏净空=%.1fmm\n",
                mainInl.size(), mainAbs.size(), nAngRej, hMax, standoff - hMax);
    }

    // ---- 高度极差判定：极差 < curvatureTol → 单一可喷面 ----
    {
        double hmin=1e18, hmax=-1e18;
        for(const auto& p:P){ double h=dot(p,n)+mainPl.d; hmin=std::min(hmin,h); hmax=std::max(hmax,h); }
        double span = hmax-hmin;
        if(span <= curvatureTol){
            // 只有一张面，不存在"重新分配"问题：这里按旧口径只按高度吸收，
            // 让陡壁/侧裙点也进 inl，重建轮廓不缺角。
            std::vector<int> allAbs; allAbs.reserve(P.size());
            for(size_t i=0;i<P.size();++i)
                if(std::fabs(dot(P[i], n) + mainPl.d) <= absorbTol) allAbs.push_back((int)i);
            PlaneSeg seg; seg.pl = mainPl; seg.inl = std::move(allAbs); seg.fixed = true;
            seg.fitSupport = std::move(mainInl);
            fprintf(stderr, "[split] 台阶高度极差 %.2f <= curvatureTol %.2f，单一可喷面（吸收 %zu 点）\n",
                    span, curvatureTol, seg.inl.size());
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
    //
    // 只在**未被主面吸收**的点上生长。
    // 这是"吸收口径要贯穿到最后"的核心：吸收已经把"高度在容差内、朝向一致"的点
    // 判给了主面，区域生长再去这些点上切块，就是同一批点被重新分配一次——
    // 切出来的碎块下游还要靠一连串合并（8° 连续合并 / 15° 同面合并 / 60° 窄条并面 /
    // 最终合并闸门）往回捞，捞不回来的就成了界面上的碎平面。
    // 现在吸收后的点直接冻结归属，区域生长只负责"主面盖不到的那部分结构"
    // （台阶、侧壁、别的朝向的面）。
    const double cosGrow = std::cos(5.0*3.141592653589793/180.0);  // 与块种子法向夹角 < 5°
    ProfStage psGrow("2d.区域生长");
    std::vector<int> regionId(P.size(), -1);
    std::vector<Vec3> seedN;
    std::vector<std::vector<int>> regions;
    std::vector<int> stack;
    for(size_t s=0; s<P.size(); ++s){
        if(regionId[s] >= 0) continue;
        if(mainLock[s]) continue;               // 已归主面，不再参与
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
                if(mainLock[j]) continue;       // 吸收集是墙，不跨过去
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
            // 种子平面**整体**冻结（法向 + 距离），不只是法向：
            // 残差必须对种子平面判，且合并过程中绝不重拟合。
            // 否则 residCap(=absorbTol) 一变，吸收的块就变，重拟合的平面跟着变——
            // 实测 absorbTol 40→60 时主面法向漂 0.0018rad，足够把 58° 陡面的
            // 内/外壁判定翻过来，于是同一个件在 40 和 60 下空移路线完全不同。
            // 结论：吸收只许改覆盖，不许改几何。平面最后只对"真共面核心"精拟合一次。
            const Plane seed = blocks[i].pl;
            for(size_t j=i+1;j<blocks.size();++j){
                if(used[j]) continue;
                // 用"对种子法向"的夹角判定，而不是对合并后重新拟合的法向判定：
                // 否则 13° 的斜坡种子可以吸收 0.5° 的地面碎块，重新拟合后法向
                // 又"平均"回地面方向，簇的身份完全乱了。
                if(std::fabs(dot(blocks[j].pl.n, seed.n)) < cs) continue;
                double maxRes = 0.0;
                for(int k : blocks[j].inl)
                    maxRes = std::max(maxRes, std::fabs(dot(P[k]-seed.c, seed.n)));
                if(const char* d=std::getenv("SPRAYTRAJ_DUMP")){
                    if(d[0]=='1' && maxRes > residCap)
                        fprintf(stderr, "[dbg] 合并拒绝: 块%zu(%zu点, 夹角%.1f°) 最大残差%.2f > %.2f\n",
                            j, blocks[j].inl.size(),
                            std::acos(std::min(1.0,std::fabs(dot(blocks[j].pl.n,seed.n))))*180.0/3.141592653589793,
                            maxRes, residCap);
                }
                if(maxRes > residCap) continue;
                cur.inl.insert(cur.inl.end(), blocks[j].inl.begin(), blocks[j].inl.end());
                used[j] = 1;
            }
            // 精拟合只做一次，且只吃"真共面核心"（到种子平面 <= 1.5·fitTol 的点）：
            // 这样解方程用的点集与 absorbTol 无关，平面不会因为阈值调大而歪。
            if(cur.inl.size() > blocks[i].inl.size()){
                const double coreTol = std::max(1.5*fitTol, 2.0);
                std::vector<int> core;
                core.reserve(cur.inl.size());
                for(int k : cur.inl)
                    if(std::fabs(dot(P[k]-seed.c, seed.n)) <= coreTol) core.push_back(k);
                Vec3 nf, cf;
                if(core.size() >= (size_t)std::max<size_t>(8, cur.inl.size()/4)
                   && ls_plane(P, core, nf, cf)){
                    if(dot(nf, seed.n) < 0) nf = nf*(-1.0);
                    if(dot(nf, seed.n) >= std::cos(3.0*3.141592653589793/180.0))
                        cur.pl = {nf, -dot(nf,cf), cf};   // 超 3° 视为漂移，保留种子平面
                }
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
    // 合并距离（判定"两个碎块是不是同一张面"）必须是**贴合尺度**，不是 absorbTol。
    // 为什么不能跟着 absorbTol 走：absorbTol 是"这张面覆盖到哪"的工艺量，可以到几十毫米；
    // 用它判平面身份，22° 的两张面（相距几十毫米）也会被判成同一张——实测 complex 件上
    // 50k 点的块 RMS 涨到 13~28mm，再被平面性判据整块丢弃，一次丢 10 行。
    // 覆盖扩展由后面的"残点归并"和"吸收集"承担，不需要合并阶段来做。
    // 3·fitTol（默认 9mm）与噪声/拟合尺度同源，与用户填多少无关。
    const double mergeResid = std::max(3.0*fitTol, 6.0);
    std::vector<RawReg> flatMerged = mergeClusters(flatRaw, 8.0, mergeResid);
    flatMerged = mergeClusters(flatMerged, 8.0, mergeResid);
    // 侧壁的最小点数与平面度门槛：太小的弧段（圆角过渡带）不算独立可喷面，
    // 否则一堆碎面各自生成轨迹，视觉上又是"乱七八糟"。
    const int    minWallPts = std::max(60, (int)(P.size()*0.008));
    const double wallRmsCap = std::max(absorbTol*0.5, 2.5);

    // ---- 主面先入列：吸收集就是主面的最终点集，平面冻结 ----
    // 顺序上放在区域生长之后、其它面之前，是为了让下面所有"并到主面"的规则
    // 有一个明确的落点；放在区域生长之前则没有区别（区域生长看不到这些点）。
    {
        PlaneSeg seg;
        seg.pl = mainPl;
        seg.inl = mainAbs;
        seg.fitSupport = std::move(mainInl);
        seg.fixed = true;                   // 平面由吸收一次定死
        fprintf(stderr, "[split] 主面入列: %zu 点（吸收集，平面冻结）\n", seg.inl.size());
        out.push_back(std::move(seg));
    }

    for(auto& m : flatMerged){
        if((int)m.inl.size() < minPts) continue;
        // 平面性检验：一堆方向各异的斜壁段拼在一起时，最小二乘平面会"平均"回
        // 地面法向，看着像个平面其实不是（棕色槽型盘盘沿的斜坡环就这样）。
        double rss=0; for(int k : m.inl){ double e=dot(P[k],m.pl.n)+m.pl.d; rss+=e*e; }
        double rms = std::sqrt(rss/(double)std::max<size_t>(1,m.inl.size()));
        if(const char* d=std::getenv("SPRAYTRAJ_DUMP"))
            if(d[0]=='1')
                fprintf(stderr, "[dbg] 平面簇 %zu点 角度%.1f° RMS=%.2f (上限%.2f)\n",
                    m.inl.size(), m.ang, rms, std::max(absorbTol, 3.0));
        if(rms > std::max(absorbTol, 3.0)){
            fprintf(stderr, "[split] 丢弃非平面簇: 点数=%zu RMS=%.2f (>%.2f) 夹角%.1f°\n",
                m.inl.size(), rms, std::max(absorbTol,3.0), m.ang);
            continue;
        }
        PlaneSeg seg; seg.pl = m.pl; seg.inl = std::move(m.inl);
        seg.fitSupport = seg.inl;
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
            if(mainLock[i]) continue;      // 已归主面的点不再被侧壁剥离抢走
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
                    if(std::fabs(s.pl.d - pl.d) < std::max(absorbTol*0.5, 2.0)) dupAt = (int)si;
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
                seg.fitSupport = seg.inl;
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
                if(mg.empty() || (b-mg.back())*binW2 >= absorbTol) mg.push_back(b);
                else if(hist2[b]>hist2[mg.back()]) mg.back()=b;
            }
            peaks2=mg;
        }
        // 塌缩判据的 RMS 上限同样受净空预算约束：否则一个跨 100mm 高度的"单面"
        // 也能被判成缓翘曲而整体塌缩，后面所有枪位就都从这张假面抬了。
        double rmsCap = std::min(std::max(2.0*absorbTol, 6.0), std::max(6.0, 0.6*standoff));
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
            seg.collapsed = true;
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
        // ---- 并回判据：看"交界处是否连续"，不看工艺容差 ----
        // 旧判据是"块内最大高差 <= absorbTol"，这里有两个方向的错：
        //   ① absorbTol 来自用户的 flatnessTol（面吸收阈值），它是"多远的点还算这张
        //      面覆盖到"的工艺口径，不该决定"这是一张面还是两张面"。p002 实测同一块
        //      弓形板：flat=5 分出 2 面、flat=40 合成 1 面，面数随工艺参数漂。
        //   ② 它量的是块自身的高度跨度。真台阶（槽型盘凸边 5mm）在 flat=40 时
        //      5 < 40 一样会被并吞；而连续弯面（p002 弓高 21mm）反而并不过去。
        // 正确的问题只有一个：两块在交界处是"贴合"还是"错开"。
        // 做法：把主面投到面内 2mm 网格上，取块 i 落在主面覆盖区（含相邻格）的点，
        //       量它们到主面平面的距离：贴合 → 只是噪声量级；台阶 → 等于台阶高。
        const Vec3 bu = norm(cross(mp.n, std::fabs(mp.n.z) < 0.9 ? Vec3{0,0,1} : Vec3{1,0,0}));
        const Vec3 bv = cross(mp.n, bu);
        auto cellKey = [](long long a, long long b){ return a*1000003LL + b; };
        auto cellOf = [&](const Vec3& p, long long& a, long long& b){
            a = (long long)std::floor(dot(p, bu)/2.0);
            b = (long long)std::floor(dot(p, bv)/2.0);
        };
        std::vector<long long> mcell; mcell.reserve(out[mi].inl.size());
        for(int k : out[mi].inl){ long long a,b; cellOf(P[k],a,b); mcell.push_back(cellKey(a,b)); }
        std::sort(mcell.begin(), mcell.end());
        mcell.erase(std::unique(mcell.begin(), mcell.end()), mcell.end());
        auto inMainFootprint = [&](const Vec3& p){
            long long a0,b0; cellOf(p,a0,b0);
            for(int da=-1; da<=1; ++da)
                for(int db=-1; db<=1; ++db)
                    if(std::binary_search(mcell.begin(), mcell.end(), cellKey(a0+da,b0+db)))
                        return true;
            return false;
        };
        // 接缝容差：真台阶的最小工艺高度（5mm 级）的一半，同时高于现场噪声（1~2mm）
        const double seamTol = 2.5;
        std::vector<int> absorbed;
        for(size_t i=0;i<out.size();++i){
            if(i == mi) continue;
            if(std::fabs(dot(out[i].pl.n, mp.n)) < cosSame) continue;
            std::vector<double> seam; seam.reserve(out[i].inl.size()/2);
            for(int k : out[i].inl)
                if(inMainFootprint(P[k])) seam.push_back(std::fabs(dot(P[k],mp.n)+mp.d));
            // 交界带至少要有几十个点，否则谈不上"贴合"——这种块是真独立的
            if(seam.size() < 20) continue;
            std::sort(seam.begin(), seam.end());
            const double seamP95 = seam[(size_t)(0.95*(double)(seam.size()-1))];
            if(seamP95 > seamTol) continue;
            fprintf(stderr, "[split] 同向同高度并入主面: %zu点 夹角%.1f° 接缝高差p95=%.2fmm(交界点%zu)\n",
                out[i].inl.size(),
                std::acos(std::min(1.0,std::fabs(dot(out[i].pl.n,mp.n))))*180.0/3.141592653589793,
                seamP95, seam.size());
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
        PlaneSeg seg; seg.pl = mainPl; seg.inl = mainAbs;
        seg.fitSupport = std::move(mainInl);
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
        const double mergeD = std::max(absorbTol, curvatureTol);
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
        // 3× 的乘子保留（翘曲板是连续弯曲，确实要比平面容差更大的高度差），
        // 但整体必须受净空预算约束：并掉的两张面高度差一旦超过 0.6×standoff，
        // 枪位就会顶进材料。c.ply 在 flat=40 时把相距 105mm 的两张平行面并成
        // 一张 RMS=50mm 的假面，根因就是这个上限只跟 tol 走、没有物理边界。
        const double warpCap  = std::min(std::max(3.0*absorbTol, 20.0),
                                         std::max(20.0, 0.6*standoff));
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
                // 组里有冻结面（主面）时保留它的平面，只并点不重拟合
                int fixedAt = -1;
                for(int m : mem) if(out[m].fixed){ fixedAt = m; break; }
                if(fixedAt >= 0){
                    seg.pl = out[fixedAt].pl;
                    seg.fitSupport = std::move(out[fixedAt].fitSupport);
                    seg.fixed = true;
                } else {
                    Vec3 nf, cf;
                    if(ls_plane(P, seg.inl, nf, cf)){
                        if(dot(nf, seg.pl.n) < 0) nf = nf*(-1.0);
                        seg.pl = {nf, -dot(nf,cf), cf};
                        seg.fitSupport = seg.inl;
                    }
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
                    // 冻结面（主面）参与合并时只并点，平面不动
                    if(out[a].fixed || out[b].fixed){
                        const PlaneSeg& fsrc = out[a].fixed ? out[a] : out[b];
                        out[a].pl = fsrc.pl;
                        if(!out[a].fixed) out[a].fitSupport = std::move(out[b].fitSupport);
                        out[a].fixed = true;
                    } else {
                        if(dot(nf, out[a].pl.n) < 0) nf = nf*(-1.0);
                        out[a].pl = {nf, -dot(nf,cf), cf};
                        out[a].fitSupport = all;
                    }
                    out[a].inl = std::move(all);
                    out.erase(out.begin()+(long)b);
                    members[a].insert(members[a].end(), members[b].begin(), members[b].end());
                    members.erase(members.begin()+(long)b);
                    did = true;
                }
            }
        }
    }

    // ---- 窄条并面：窄平面 + 同方向 + 法向差不大 → 同一张面 ----
    // 口径（2026-11-13 用户拍板）：窄平面之间只要朝向一致、法向差不大，就是同一张面，
    // 不该为了十几度的法向差各布一套行、各来一次换面空移——尤其是"面宽小于一个行距"
    // 的细条，一条行本来就能盖住，分几张面就是同一块地方反复喷。
    // 判据（四条同时满足）：
    //   ① 窄：两块自身短向跨度都 <= narrowW（默认 = 行距）。这是本条的**前提**：
    //      "窄"保证了并完仍然是一条行能盖住的带，误差有界；
    //   ② 同向：法向夹角 <= stripAng（默认 45°，比上面的通用 15° 宽，但只对窄条放开）；
    //   ③ 相邻：KNN 图里有边横跨两面（同一次扫描的连续数据，排除同向孤岛）；
    //   ④ 并后仍可喷：重拟合后短向跨度 <= narrowW，且到重拟合平面的最大偏差
    //      <= devCap（默认 8mm，低于面级平面性门 max(3,0.12*standoff)=12mm，
    //      保证并出来的面不会被下游当"假面"拒掉）。
    // 为什么"窄"才放开角度：宽板之间的缓角是真实翘曲/折边的链式漂移，一路并会塌成
    // 一张平板（这是 8° 硬门限存在的理由）；而窄条之间的缓角几乎总是同一张面被扫描
    // 条带/噪声切开，并了以后几何误差仍然有界，代价只有表面偏离几毫米。
    if(out.size() > 1 && out.size() <= 48){
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
        const double Lrow = (rowSpacing > 1.0) ? rowSpacing : 80.0;
        double narrowW = Lrow;            // 窄条上限：一条行能盖住
        double stripAng = 60.0;           // 窄条之间允许的法向差（实测 45° 卡掉半棍子一次可并的 49.4°）
        double devCap = std::max(4.0, 0.08*Lrow);   // 并后偏离上限（默认 L=80 → 6.4mm）
        if(const char* e = std::getenv("SPRAYTRAJ_STRIP_W"))  { double v=atof(e); if(v>1.0) narrowW=v; }
        if(const char* e = std::getenv("SPRAYTRAJ_STRIP_ANG")){ double v=atof(e); if(v>0.0) stripAng=v; }
        if(const char* e = std::getenv("SPRAYTRAJ_STRIP_DEV")){ double v=atof(e); if(v>0.1) devCap=v; }
        const char* dumpEnv2 = std::getenv("SPRAYTRAJ_DUMP");
        const bool dbg2 = dumpEnv2 && dumpEnv2[0] == '1';
        const double cosS = std::cos(stripAng*3.141592653589793/180.0);
        Vec3 partC2{0,0,0};
        for(const auto& p : P) partC2 = partC2 + p;
        partC2 = partC2*(1.0/(double)std::max<size_t>(1,P.size()));
        auto orient2 = [&](Vec3 n, const Vec3& c){
            if(std::fabs(n.z) > 0.05){ if(n.z < 0) n = n*(-1.0); }
            else if(dot(n, c - partC2) < 0) n = n*(-1.0);
            return n;
        };
        // 面内两个正交方向的跨度（mm）
        auto extents = [&](const Plane& pl, const std::vector<int>& inl, double& eu, double& ev){
            Vec3 ax = (std::fabs(pl.n.x) < 0.9) ? Vec3{1.0,0.0,0.0} : Vec3{0.0,1.0,0.0};
            Vec3 u = norm(cross(ax, pl.n));
            Vec3 v = cross(pl.n, u);
            double u0=1e30,u1=-1e30,v0=1e30,v1=-1e30;
            for(int k : inl){
                const Vec3 w = P[(size_t)k] - pl.c;
                const double a = dot(w,u), b = dot(w,v);
                u0=std::min(u0,a); u1=std::max(u1,a);
                v0=std::min(v0,b); v1=std::max(v1,b);
            }
            eu = u1-u0; ev = v1-v0;
        };
        std::vector<std::vector<int>> mem(F0);
        for(size_t i=0;i<F0;++i) mem[i] = {(int)i};
        // 并后枪能不能站住：并面的法向是各成员的平均，弧面上的并面法向会偏开局部径向，
        // 100mm 枪位沿它走出去就可能擦回工件（整根管实测：并面 +n 自由只剩 46.7mm，
        // 等于 2R·cos21°，枪直接撞在对面内壁上）。所以每条合并候选先用同一套几何
        // （点云最近距离）验一遍枪位空间，过不了就不并——宁可多一张面，不能并出一张喷不了的面。
        KdTree kdStrip;
        { ProfStage psKd("2e.窄条并面-KD"); kdStrip.init(P); }
        const double soRef = (standoff > 1.0) ? standoff : 100.0;
        // 与下游 verify_gun_side 完全同一口径：6 个采样点沿 ±n 打射线，取中位自由距离，
        // 两侧取大者 >= 0.85·standoff 才算"枪站得住"。
        auto gunCanStand = [&](const std::vector<int>& all, const Vec3& n)->bool{
            const size_t N = all.size();
            if(!N) return false;
            const double maxFree = std::max(400.0, soRef*4.0);
            const int want = 6;
            std::vector<double> dp, dm;
            dp.reserve((size_t)want); dm.reserve((size_t)want);
            for(int i=0;i<want;++i){
                size_t k = (size_t)((double)i/(double)want*(double)N);
                if(k >= N) k = N-1;
                const Vec3 p = P[(size_t)all[k]];
                dp.push_back(ray_free_dist(P, kdStrip, p + n*0.4,  n,        maxFree, 1.5, 5.0));
                dm.push_back(ray_free_dist(P, kdStrip, p - n*0.4,  n*(-1.0), maxFree, 1.5, 5.0));
            }
            auto med = [](std::vector<double>& v){
                if(v.empty()) return -1.0;
                std::nth_element(v.begin(), v.begin()+v.size()/2, v.end());
                return v[v.size()/2];
            };
            const double fp = med(dp), fm = med(dm);
            const double effPlus  = (fp < 0.0) ? 1e9 : fp;
            const double effMinus = (fm < 0.0) ? 1e9 : fm;
            return std::max(effPlus, effMinus) >= soRef*0.85;
        };
        bool did2 = true;
        while(did2 && out.size() > 1){
            did2 = false;
            for(size_t a=0; a<out.size() && !did2; ++a){
                for(size_t b=a+1; b<out.size() && !did2; ++b){
                    double eua,eva,eub,evb;
                    extents(out[a].pl, out[a].inl, eua, eva);
                    extents(out[b].pl, out[b].inl, eub, evb);
                    const double wa = std::min(eua,eva), wb = std::min(eub,evb);
                    if(wa > narrowW || wb > narrowW) continue;          // ① 窄
                    const double ang = std::acos(std::max(-1.0, std::min(1.0,
                        dot(orient2(out[a].pl.n, out[a].pl.c), orient2(out[b].pl.n, out[b].pl.c)))))*180.0/3.141592653589793;
                    if(ang > stripAng) continue;                        // ② 同向
                    if(std::cos(ang*3.141592653589793/180.0) < cosS) continue;
                    bool touch = false;
                    for(int x : mem[a]){ for(int y : mem[b]) if(adj[(size_t)x*F0+y]){ touch=true; break; } if(touch) break; }
                    if(!touch) continue;                                // ③ 相邻（同一片扫描数据）
                    std::vector<int> all = out[a].inl;
                    all.insert(all.end(), out[b].inl.begin(), out[b].inl.end());
                    std::sort(all.begin(), all.end());
                    all.erase(std::unique(all.begin(), all.end()), all.end());
                    Vec3 nf, cf;
                    if(!ls_plane(P, all, nf, cf)) continue;
                    Vec3 orientN = orient2(nf, cf);
                    if(dot(orientN, out[a].pl.n) < 0) nf = nf*(-1.0);
                    Plane pf{nf, -dot(nf,cf), cf};
                    double eu,ev; extents(pf, all, eu, ev);
                    const double wm = std::min(eu,ev);
                    if(wm > narrowW) continue;                          // ④ 并后仍是一条行能盖住
                    double dev = 0.0;
                    for(int k : all) dev = std::max(dev, std::fabs(dot(P[(size_t)k]-cf, nf)));
                    if(dev > devCap){
                        if(dbg2) fprintf(stderr, "[dbg] 窄条并面放弃: 面%zu-%zu 夹角%.1f° 并后偏差%.2f > %.2f\n",
                                         a, b, ang, dev, devCap);
                        continue;
                    }
                    if(!gunCanStand(all, nf)){
                        if(dbg2) fprintf(stderr, "[dbg] 窄条并面放弃: 面%zu-%zu 夹角%.1f° 并后枪位站不住\n",
                                         a, b, ang);
                        continue;
                    }
                    fprintf(stderr, "[split] 窄条并面(法向%.1f° 宽%.1fmm 偏差%.1fmm): 面%zu(%zu点) + 面%zu(%zu点) -> %zu点\n",
                            ang, wm, dev, a, out[a].inl.size(), b, out[b].inl.size(), all.size());
                    if(!(out[a].fixed || out[b].fixed)){
                        out[a].pl = pf;
                        out[a].fitSupport = all;
                    } else {
                        const PlaneSeg& fsrc = out[a].fixed ? out[a] : out[b];
                        out[a].pl = fsrc.pl;
                        if(!out[a].fixed) out[a].fitSupport = std::move(out[b].fitSupport);
                        out[a].fixed = true;
                    }
                    out[a].inl = std::move(all);
                    out.erase(out.begin()+(long)b);
                    mem[a].insert(mem[a].end(), mem[b].begin(), mem[b].end());
                    mem.erase(mem.begin()+(long)b);
                    did2 = true;
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
    std::vector<std::vector<Vec3>> rowTurns; // 截线之间的连续喷涂连接
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

static void measure_surface_row_coverage(SurfPatch& sp, double spacing)
{
    KdTree tree;
    tree.init(sp.rowP);
    const double tol = std::max(sp.rowGap, spacing)*0.75;
    int covered = 0;
    double squaredError = 0.0;
    for(const Vec3& p : sp.P){
        KnnBuf nearest;
        tree.knn(p, 1, nearest);
        if(nearest.n == 0) continue;
        if(nearest.d2[0] <= tol*tol) ++covered;
        squaredError += nearest.d2[0];
    }
    sp.cover = (double)covered/(double)sp.P.size();
    sp.fitRms = std::sqrt(squaredError/(double)sp.P.size());
}

// 在曲面块上生成蛇形行（就地写回 sp.rowP/rowN/rowBeg）
static bool gen_surface_rows(SurfPatch& sp, double L, double pstep, double fillHole,
                             int dirMode, double speed, bool& openContours)
{
    openContours = false;
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

    // ---- 回转面：让每一行都从"上一行收尾的那个方位角"起头 ----
    // 只按高度挑种子时，相邻两行的起点会落在圈的任意方位，可能在正对面：
    // 行 k 收尾在方位 0°，行 k+1 起头在方位 180°，行间连接要横穿整个圆
    // （面点相距 200mm+），而 emitGap 判据是 "d > 3×行距就走 spray=0 的空移" ——
    // 于是整条轨迹被切成 7 段、每段之间都要抬枪跨过去。
    // 实测垃圾桶 30mm 行距：6 处空移全是这个原因（30~61mm 的枪位距离，面点 200mm+）。
    // 对齐方位后行间连接只有径向一个行距（26mm），可以连续喷涂。
    double phiPrev = 1e9;
    auto phiOf = [&](const Vec3& p)->double{
        Vec3 q = p - sp.c0;
        return std::atan2(dot(q, eIdx), dot(q, eRow));
    };

    for(size_t k=0; k<anchors.size(); ++k){
        const double tk = anchors[k];
        int seed = -1; double bestd = 1e18;
        for(int pass=0; pass<(isoMode ? 2 : 1) && seed < 0; ++pass){
            for(int i=0;i<m;++i){
                if(!isoOk[i]) continue;              // 平面点不参与等高线取种
                if(pass==0 && isoMode && phiPrev < 1e8){
                    double dd = std::fabs(std::remainder(phiOf(sp.P[i]) - phiPrev,
                                                         6.283185307179586));
                    if(dd > 0.70) continue;          // 与上一行收尾方位差 > 40° 的不用
                }
                double d = std::fabs(iv[i]-tk);
                if(d<bestd){ bestd=d; seed=i; }
            }
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
        if(isoMode && !rp.empty()) phiPrev = phiOf(rp.back());   // 下一行从这里起头
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
    // 回转面（isoMode）**跳过这一步**：等高线圈上的法向是"绕轴匀速旋转"的，
    // 相邻行在同一方位角处法向本来就一致，不需要跨行统一；而姿态场的格宽
    // 1.3×行距（39mm）在 r=100 的圈上等于 22° 方位角 —— 把连续旋转的法向
    // 压成十几段阶梯，实测法向与理想锥面的偏差中位 8~11°、p90 17°，
    // 经 100mm 放大成 14~19mm 的枪位抖动，画出来就是"折角/分叉/8 字"。
    if(sp.modelOk && !sp.rowP.empty() && !isoMode){
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

    // ---- 剔除"卷边行"：行内法向翻转 >90° ----
    // 曲型板薄壁的圆角边缘处，iso 线会顺着卷边绕到板的另一侧，同一行里前半段
    // 法向朝正面、后半段朝背面。这样的行不能按"一行一个姿态"喷：
    // 相邻两点法向几乎相反，枪位会从面的一侧瞬移到另一侧（实测 100mm 离面时
    // 跨 195mm），而这两点又都在 spray=1 的同一行里，中间没有任何空移点，
    // 折线直接穿过工件（实测折线中段最小净空 0.09mm）。
    // 判据用相邻法向点积 <0（>90°）：正常行逐点转角 <3°（实测），余量充足。
    // 整行剔除而不是拆分：拆出来的短段仍贴着卷边背面，喷它没有工艺意义，
    // 保留只会让枪绕到板背面。行距 20/60 时该行本来就不存在（iso 线离边缘更远），
    // 所以剔除它与"行距只影响布行"的预期一致。
    {
        std::vector<Vec3> rp, rn;
        std::vector<int> rb;
        rp.reserve(sp.rowP.size()); rn.reserve(sp.rowN.size());
        rb.reserve(sp.rowBeg.size()); rb.push_back(0);
        int dropped = 0;
        std::vector<double> dropMinS;
        for(int r=0;r<sp.nRow;++r){
            const int b=sp.rowBeg[r], e=sp.rowBeg[r+1];
            bool flipped = false;
            for(int k=b+1;k<e;++k){
                if(dot(sp.rowN[k-1], sp.rowN[k]) < 0.0){ flipped = true; break; }
            }
            if(flipped){
                ++dropped;
                dropMinS.push_back(dot(sp.rowP[b], eIdx));
                continue;
            }
            rp.insert(rp.end(), sp.rowP.begin()+b, sp.rowP.begin()+e);
            rn.insert(rn.end(), sp.rowN.begin()+b, sp.rowN.begin()+e);
            rb.push_back((int)rp.size());
        }
        if(dropped){
            fprintf(stderr, "[曲面] 剔除卷边行 %d 行（行内法向翻转 >90°，枪位会绕到板背面）：", dropped);
            for(size_t q=0;q<dropMinS.size();++q)
                fprintf(stderr, " s=%.1f%s", dropMinS[q], q+1<dropMinS.size()?",":"\n");
            sp.rowP.swap(rp); sp.rowN.swap(rn); sp.rowBeg.swap(rb);
            sp.nRow = (int)sp.rowBeg.size()-1;
            rowsMade = sp.nRow;
        }
    }
    if(rowsMade < 1) return false;

    // 围绕质心有点分布，并不能证明表面上存在闭合等高线。
    if(isoMode){
        const double minLength = 0.5*std::min(ext1, ext2);
        openContours = true;
        for(int r=0; r<sp.nRow; ++r){
            const int b = sp.rowBeg[r], e = sp.rowBeg[r+1];
            if(len(sp.rowP[e-1]-sp.rowP[b]) > 3.0*pstep) continue;
            double length = 0.0;
            for(int k=b+1; k<e; ++k) length += len(sp.rowP[k]-sp.rowP[k-1]);
            if(length >= minLength){ openContours = false; break; }
        }
    }

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
    measure_surface_row_coverage(sp, L);
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
static bool build_surface_mesh(const SurfPatch& sp, double fillHole,
                               std::vector<Vec3>& verts, std::vector<int>& tris,
                               std::vector<unsigned char>& filledTris, int maxVerts, int maxTris)
{
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
    std::vector<int> facingCount(NG,0);
    for(int i=0;i<m;++i){
        int iu=(int)std::lround((uv[i].x-minU)/du);
        int iv2=(int)std::lround((uv[i].y-minV)/dv);
        iu=std::min(nU-1,std::max(0,iu)); iv2=std::min(nV-1,std::max(0,iv2));
        const size_t k=(size_t)iu*nV+iv2;
        bucket[k].push_back(dot(sp.P[i]-sp.c0, sp.n0));
        if(dot(sp.N[i],sp.n0)>=0.7071067811865476) ++facingCount[k]; // 与基面法向偏差不超过 45°
    }
    std::vector<double> h(NG, 1e300);
    std::vector<unsigned char> occ(NG, 0);
    std::vector<unsigned char> rimSupport(NG,0);
    for(size_t i=0;i<NG;++i){
        std::vector<double>& bk=bucket[i];
        if(bk.empty()) continue;
        std::nth_element(bk.begin(), bk.begin()+bk.size()/2, bk.end());
        h[i]=bk[bk.size()/2];       // 中位高：单点噪声顶不出尖刺
        occ[i]=1;
        rimSupport[i] = 2*facingCount[i]>=(int)bk.size() ? 2 : 1;
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
    fill_reconstruction_holes(occ, nV, nU, dv, du, fillHole, &rimSupport);
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
    verts.clear(); tris.clear(); filledTris.clear();
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
        filledTris.push_back(occ[cc[0]]==2 || occ[cc[1]]==2 || occ[cc[2]]==2);
        tris.push_back(idx4[0]); tris.push_back(idx4[2]); tris.push_back(idx4[3]);
        filledTris.push_back(occ[cc[0]]==2 || occ[cc[2]]==2 || occ[cc[3]]==2);
    }
    int occN=0; for(size_t i=0;i<NG;++i) if(occ[i]) ++occN;
    fprintf(stderr, "[曲面] 曲面网格: %dx%d 顶点=%d 三角=%d (cell=%.1fmm 占用格=%d/%d)\n",
            nU, nV, (int)verts.size(), (int)(tris.size()/3), cell, occN, (int)NG);
    return !verts.empty() && !tris.empty();
}

// 同一网格截线遇到孔洞或断口时，拆成独立喷涂段。
template<class ClearSegment>
static bool gen_surface_mesh_rows(SurfPatch& sp, const std::vector<Vec3>& verts,
                                  const std::vector<int>& tris, const std::vector<unsigned char>& filledTris,
                                  double spacing, double step,
                                  int dirMode, double turnRadius, double fanRadius, ClearSegment clearSegment,
                                  double& scoreOut)
{
    if(verts.empty() || tris.empty() || spacing <= 0.0 || step <= 0.0) return false;
    std::vector<Vec2> footprint, rectangle;
    footprint.reserve(verts.size());
    for(const Vec3& p : verts){
        const Vec3 d = p - sp.c0;
        footprint.push_back({dot(d, sp.eL), dot(d, sp.eS)});
    }
    if(minAreaRect(footprint, rectangle) <= 0.0) return false;
    Vec2 axis = rectangle[1]-rectangle[0], other = rectangle[2]-rectangle[1];
    if((len2(axis) < len2(other)) != (dirMode == 2)) axis = other;
    const Vec3 along = norm(sp.eL*axis.x + sp.eS*axis.y);
    const Vec3 across = cross(sp.n0, along);
    std::vector<Vec2> uv;
    uv.reserve(verts.size());
    double lo = 1e18, hi = -1e18;
    for(const Vec3& p : verts){
        const Vec3 d = p - sp.c0;
        uv.push_back({dot(d, along), dot(d, across)});
        lo = std::min(lo, uv.back().y);
        hi = std::max(hi, uv.back().y);
    }
    if(hi - lo < 1e-6) return false;
    const int count = (int)std::min(2000.0, std::ceil((hi-lo)/spacing));
    const double gap = (hi-lo)/count;
    // 固定枪轴有条件：只在重建面朝向和工件实测点都支持时出喷涂点。
    constexpr double minNormalDot = 0.8660254037844386; // 入射偏角 30°
    const double supportDist = std::max(sp.mCell, 2.0*sp.rFit);
    // 沿枪轴核对材料支撑。枪位高度沿已重建的连续网格走，
    // 以免台阶附近最近点在不同高度层之间来回跳。
    std::vector<Vec3> footprintSamples;
    footprintSamples.reserve(sp.P.size());
    for(size_t i=0; i<sp.P.size(); ++i){
        if(dot(sp.N[i], sp.n0) < minNormalDot) continue;
        const Vec3 d = sp.P[i]-sp.c0;
        footprintSamples.push_back({dot(d, along), dot(d, across), 0.0});
    }
    if(footprintSamples.size() < 3) return false;
    KdTree footprintTree;
    footprintTree.init(footprintSamples);
    auto supported = [&](const Vec3& p){
        const Vec3 d = p-sp.c0;
        KnnBuf neighbors;
        footprintTree.knn({dot(d,along), dot(d,across), 0.0}, 3, neighbors);
        return neighbors.n == 3 && neighbors.d2[0] <= supportDist*supportDist;
    };
    struct Section { double x0, x1; Vec3 first, last; bool filled; };
    auto sectionRow = [&](double y, bool reverse, std::vector<Vec3>& points,
                          std::vector<int>& beginnings){
        const size_t rowBounds = beginnings.size();
        double runScore = 0.0;
        std::vector<Section> sections;
        for(size_t i=0; i<tris.size(); i+=3){
            Vec3 hits[2];
            double xs[2];
            int nHits = 0;
            for(int j=0; j<3; ++j){
                const int a = tris[i+j], b = tris[i+(j+1)%3];
                if((uv[a].y <= y) == (uv[b].y <= y)) continue;
                const double f = (y-uv[a].y)/(uv[b].y-uv[a].y);
                hits[nHits] = verts[a] + (verts[b]-verts[a])*f;
                xs[nHits++] = uv[a].x + (uv[b].x-uv[a].x)*f;
            }
            if(nHits != 2 || std::fabs(xs[1]-xs[0]) < 1e-6) continue;
            if(xs[1] < xs[0]){ std::swap(xs[0], xs[1]); std::swap(hits[0], hits[1]); }
            sections.push_back({xs[0], xs[1], hits[0], hits[1], filledTris[i/3]!=0});
        }
        std::sort(sections.begin(), sections.end(), [](const Section& a, const Section& b){
            return a.x0 < b.x0;
        });
        if(reverse){
            std::reverse(sections.begin(), sections.end());
            for(auto& s : sections){
                std::swap(s.x0, s.x1);
                std::swap(s.first, s.last);
            }
        }
        auto finishRun = [&]{
            const int start = beginnings.back();
            if((int)points.size() - start < 2){
                points.resize(start);
                return;
            }
            const double span = std::fabs(dot(points.back()-points[start], along));
            if(span < 2.0*step){
                points.resize(start);
                return;
            }
            beginnings.push_back((int)points.size());
            runScore += span - std::min(gap*0.5, 2.0*supportDist);
        };
        std::vector<Vec3> lane;
        std::vector<char> support;
        auto finishLane = [&]{
            // 同一个连续网格条带内，允许跨过短的测量缺口。孔洞与净空失败已在外层断开。
            for(size_t j=0; j<lane.size();){
                if(support[j]){ ++j; continue; }
                const size_t begin = j;
                while(j < lane.size() && !support[j]) ++j;
                if(begin > 0 && j < lane.size() &&
                   std::fabs(dot(lane[j]-lane[begin-1],along)) <= 2.0*supportDist)
                    std::fill(support.begin()+begin, support.begin()+j, 1);
            }
            for(size_t j=0; j<lane.size(); ++j){
                if(support[j]) points.push_back(lane[j]);
                else finishRun();
            }
            finishRun();
            lane.clear();
            support.clear();
        };
        double end = -1e18;
        for(const auto& s : sections){
            const bool adjacent = std::fabs(s.x0-end) <= 1e-5;
            if(!adjacent) finishLane();
            const int steps = std::max(1, (int)std::ceil(len(s.last-s.first)/step));
            for(int j=0; j<=steps; ++j){
                if(j == 0 && adjacent && !lane.empty()) continue;
                const Vec3 p = s.first + (s.last-s.first)*((double)j/steps);
                if(!clearSegment(p,p)){
                    finishLane();
                    continue;
                }
                if(!lane.empty()){
                    const Vec3 delta = p-lane.back();
                    if(std::fabs(dot(delta,sp.n0)) > std::fabs(dot(delta,along))+step ||
                       !clearSegment(lane.back(), p))
                        finishLane();
                }
                lane.push_back(p);
                support.push_back(s.filled || supported(p));
            }
            end = s.x1;
        }
        finishLane();
        // 锯齿边缘的短残段容易多出一趟绕行。仅在本行已经有完整长段、
        // 且短段仍落在相邻长段的喷幅内时将其吸收。
        const size_t firstRun = rowBounds;
        if(beginnings.size()>firstRun+1){
            size_t write=beginnings[firstRun-1];
            const size_t runCount=beginnings.size();
            const auto spans=beginnings;
            beginnings.resize(firstRun);
            for(size_t r=firstRun; r<runCount; ++r){
                const int b=spans[r-1], e=spans[r];
                const double span=len(points[e-1]-points[b]);
                const bool coveredBefore = beginnings.size()>firstRun &&
                    len(points[write-1]-points[beginnings[beginnings.size()-2]])>fanRadius &&
                    len(points[e-1]-points[write-1])<=fanRadius;
                const bool coveredAfter = r+1<runCount &&
                    len(points[spans[r+1]-1]-points[spans[r]])>fanRadius &&
                    len(points[b]-points[spans[r]])<=fanRadius;
                if(span<=fanRadius && (coveredBefore || coveredAfter)){
                    runScore -= span - std::min(gap*0.5,2.0*supportDist);
                    continue;
                }
                for(int k=b; k<e; ++k) points[write++]=points[k];
                beginnings.push_back((int)write);
            }
            points.resize(write);
        }
        return runScore;
    };
    // 整组行平移同一距离，行距保持一致；优先保留长的连续喷涂段。
    std::vector<Vec3> bestPoints;
    std::vector<int> bestBeginnings;
    std::vector<int> bestRows;
    double bestScore = -1e18;
    for(double offset : {0.0, -1.0/6, 1.0/6, -1.0/3, 1.0/3}){
        std::vector<Vec3> points;
        std::vector<int> beginnings{0};
        std::vector<int> rows;
        double score = 0.0;
        bool reverse = false;
        for(int row=0; row<count; ++row){
            const size_t before = rows.size();
            const double y=lo+(row+0.5+offset)*gap;
            double rowScore=sectionRow(y,reverse,points,beginnings);
            if(count>1 && (row==0 || row+1==count) && beginnings.size()>before+2){
                // 最外行恰好落在阶梯轮廓上时，在半喷幅内向面内移，优先用完整直行。
                const double budget=std::min(fanRadius,gap*0.25);
                const double sign=row==0 ? 1.0 : -1.0;
                size_t bestRuns=beginnings.size()-1-before;
                const double penalty=std::min(gap*0.5,2.0*supportDist);
                const double originalSpan=rowScore+penalty*bestRuns;
                for(double shift : {budget*0.5,budget}){
                    std::vector<Vec3> trialPoints;
                    std::vector<int> trialBeginnings{0};
                    const double trialScore=sectionRow(y+sign*shift,reverse,trialPoints,trialBeginnings);
                    const size_t runs=trialBeginnings.size()-1;
                    if(runs==0 || runs>=bestRuns || trialScore+penalty*runs+budget<originalSpan)
                        continue;
                    const int base=beginnings[before];
                    points.resize(base);
                    points.insert(points.end(),trialPoints.begin(),trialPoints.end());
                    beginnings.resize(before+1);
                    for(size_t k=1; k<trialBeginnings.size(); ++k)
                        beginnings.push_back(base+trialBeginnings[k]);
                    bestRuns=runs;
                    rowScore=trialScore;
                }
            }
            score += rowScore;
            rows.resize(beginnings.size()-1, row);
            if(rows.size() != before) reverse = !reverse;
        }
        for(size_t r=1; r<rows.size(); ++r){
            const int k = beginnings[r];
            score -= 0.5*std::fabs(dot(points[k]-points[k-1],along));
        }
        if(beginnings.size() > 1 && score > bestScore){
            bestPoints.swap(points);
            bestBeginnings.swap(beginnings);
            bestRows.swap(rows);
            bestScore = score;
        }
    }
    if(bestBeginnings.empty()) return false;
    scoreOut = bestScore;
    // 实体上的相邻行允许保持喷涂。板边外的同侧折返仍保留弧线，走安全净空后关枪通过。
    auto meshPoint = [&](Vec2 q, Vec3& p, bool& filled){
        for(size_t k=0; k<tris.size(); k+=3){
            const int a = tris[k], b = tris[k+1], c = tris[k+2];
            const Vec2 u = uv[b]-uv[a], v = uv[c]-uv[a], d = q-uv[a];
            const double det = cross2(u,v);
            if(std::fabs(det) < 1e-9) continue;
            const double t = cross2(d,v)/det, s = cross2(u,d)/det;
            if(t < -1e-7 || s < -1e-7 || t+s > 1.0+1e-7) continue;
            p = verts[a] + (verts[b]-verts[a])*t + (verts[c]-verts[a])*s;
            filled = filledTris[k/3]!=0;
            return (int)(k/3);
        }
        return -1;
    };
    auto toUV = [&](const Vec3& p){
        const Vec3 d = p-sp.c0;
        return Vec2{dot(d, along), dot(d, across)};
    };
    auto traceOnMesh = [&](const Vec3& a, const Vec3& b, double bulge,
                           std::vector<Vec3>& path, double side=0.0){
        const Vec2 av=toUV(a), bv=toUV(b);
        const int samples=std::max(2,(int)std::ceil((len(b-a)+2.0*(std::fabs(bulge)+std::fabs(side)))/
                                                   std::max(0.5,step*0.5)));
        path.clear();
        Vec3 prev=a;
        for(int i=1; i<samples; ++i){
            const double f=(double)i/samples;
            const double progress = bulge==0.0 ? f : f*f*(3.0-2.0*f);
            const double arch=4.0*f*(1.0-f);
            const Vec2 q=av+(bv-av)*progress+Vec2{arch*bulge,arch*arch*side};
            Vec3 p;
            bool filled=false;
            if(meshPoint(q,p,filled)<0 || len(p-prev)>2.0*step || !clearSegment(prev,p))
                return false;
            path.push_back(p);
            prev=p;
        }
        return len(b-prev)<=2.0*step && clearSegment(prev,b);
    };
    // 以三角形重心和共享边中点构建面内通路。只有真正共享边的三角片可以互通。
    const int triangleCount = (int)(tris.size()/3);
    std::vector<Vec3> nodes;
    std::vector<std::vector<int>> links(triangleCount);
    nodes.reserve(tris.size());
    for(int t=0; t<triangleCount; ++t)
        nodes.push_back((verts[tris[t*3]]+verts[tris[t*3+1]]+verts[tris[t*3+2]])*(1.0/3.0));
    std::map<std::pair<int,int>,int> edgeOwners;
    for(int t=0; t<triangleCount; ++t) for(int e=0; e<3; ++e){
        const int a=tris[t*3+e], b=tris[t*3+(e+1)%3];
        auto [owner, inserted] = edgeOwners.emplace(std::minmax(a,b), t);
        if(inserted) continue;
        const int portal = (int)nodes.size(), neighbor = owner->second;
        nodes.push_back((verts[a]+verts[b])*0.5);
        links.push_back({t,neighbor});
        links[t].push_back(portal);
        links[neighbor].push_back(portal);
    }
    auto routeOnMesh = [&](const Vec3& a, const Vec3& b, std::vector<Vec3>& path){
        Vec3 surface;
        bool filled = false;
        const int start=meshPoint(toUV(a),surface,filled), goal=meshPoint(toUV(b),surface,filled);
        if(start<0 || goal<0 || !clearSegment(a,nodes[start]) || !clearSegment(nodes[goal],b))
            return false;
        std::vector<double> distance(nodes.size(),1e30);
        std::vector<int> parent(nodes.size(),-1);
        using Entry = std::pair<double,int>;
        std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>> pending;
        distance[start]=0.0;
        pending.emplace(0.0,start);
        while(!pending.empty()){
            const auto [cost,u] = pending.top();
            pending.pop();
            if(cost>distance[u]) continue;
            if(u==goal) break;
            for(int v : links[u]){
                const double next=cost+len(nodes[v]-nodes[u]);
                if(next >= distance[v] || !clearSegment(nodes[u],nodes[v])) continue;
                distance[v]=next;
                parent[v]=u;
                pending.emplace(next,v);
            }
        }
        if(distance[goal]>=1e30) return false;
        std::vector<Vec3> waypoints{b};
        for(int u=goal; u>=0; u=parent[u]) waypoints.push_back(nodes[u]);
        waypoints.push_back(a);
        std::reverse(waypoints.begin(), waypoints.end());
        path.clear();
        // 在有效网格上尝试跳过中间节点，消除三角片边上的锯齿折返。
        std::vector<Vec3> shortcut;
        for(size_t i=0; i+1<waypoints.size();){
            size_t j=waypoints.size()-1;
            for(; j>i+1; --j){
                if(traceOnMesh(waypoints[i],waypoints[j],0.0,shortcut)){
                    path.insert(path.end(),shortcut.begin(),shortcut.end());
                    break;
                }
            }
            if(j==i+1){
                const Vec3 delta=waypoints[j]-waypoints[i];
                const int samples=std::max(1,(int)std::ceil(len(delta)/step));
                for(int k=1; k<samples; ++k)
                    path.push_back(waypoints[i]+delta*((double)k/samples));
            }
            if(j+1<waypoints.size()) path.push_back(waypoints[j]);
            i=j;
        }
        return true;
    };
    std::vector<std::vector<Vec3>> turns(bestRows.size());
    std::vector<int> starts(bestBeginnings.begin(),bestBeginnings.end()-1);
    std::vector<int> ends(bestBeginnings.begin()+1,bestBeginnings.end());
    const double trimLimit=std::min(fanRadius, gap*0.25);
    int alignedTurns=0, meshDetours=0;
    for(size_t r=1; r<bestRows.size(); ++r){
        const Vec3& a = bestPoints[ends[r-1]-1];
        const Vec3& b = bestPoints[starts[r]];
        const Vec2 av = toUV(a), bv = toUV(b);
        const double sign = (av.x >= toUV(bestPoints[starts[r-1]]).x) ? 1.0 : -1.0;
        const bool adjacent = bestRows[r]-bestRows[r-1]==1 && std::fabs(bv.x-av.x)<=1.5*gap;
        const double bulge = adjacent ? std::min(turnRadius, gap*0.5)*sign : 0.0;
        auto& turn=turns[r];
        if(adjacent && std::fabs(av.x-bv.x) <= trimLimit){
            // 先对齐相邻行的行端，最多退到喷幅仍能覆盖原端点的位置。
            bool aligned=false;
            for(double extra : {0.0,trimLimit*0.25,trimLimit*0.5}){
                const double endLine=std::min(sign*av.x,sign*bv.x)-extra;
                int ai=ends[r-1]-1, bi=starts[r];
                while(ai>starts[r-1]+1 && sign*toUV(bestPoints[ai]).x>endLine+step*0.5) --ai;
                while(bi<ends[r]-2 && sign*toUV(bestPoints[bi]).x>endLine+step*0.5) ++bi;
                if(sign*(av.x-toUV(bestPoints[ai]).x)>trimLimit ||
                   sign*(bv.x-toUV(bestPoints[bi]).x)>trimLimit ||
                   std::fabs(toUV(bestPoints[ai]).x-toUV(bestPoints[bi]).x)>step)
                    continue;
                if(traceOnMesh(bestPoints[ai],bestPoints[bi],bulge,turn) ||
                   traceOnMesh(bestPoints[ai],bestPoints[bi],bulge*0.5,turn) ||
                   traceOnMesh(bestPoints[ai],bestPoints[bi],0.0,turn)){
                    ends[r-1]=ai+1; starts[r]=bi;
                    ++alignedTurns;
                    aligned=true;
                    break;
                }
            }
            if(aligned) continue;
        }
        if(adjacent && (traceOnMesh(a,b,bulge,turn) || traceOnMesh(a,b,bulge*0.5,turn)
                     || traceOnMesh(a,b,-bulge*0.5,turn))) continue;
        if(traceOnMesh(a,b,0.0,turn)) continue;
        if(bestRows[r]==bestRows[r-1]){
            // 同一行内的小凹口，用一条与喷涂行相切的短弧绕开。
            const double reach=std::min(gap*0.5,2.0*fanRadius);
            bool bowed=false;
            for(double retreat : {0.0,trimLimit*0.5,trimLimit}){
                int ai=ends[r-1]-1, bi=starts[r];
                while(ai>starts[r-1]+1 && len(bestPoints[ai-1]-a)<=retreat) --ai;
                while(bi<ends[r]-2 && len(bestPoints[bi+1]-b)<=retreat) ++bi;
                for(double f : {0.25,0.5,0.75,1.0}){
                    const double d=reach*f;
                    if(traceOnMesh(bestPoints[ai],bestPoints[bi],0.0,turn,d) ||
                       traceOnMesh(bestPoints[ai],bestPoints[bi],0.0,turn,-d)){
                        ends[r-1]=ai+1; starts[r]=bi;
                        bowed=true;
                        scoreOut-=d+retreat;
                        break;
                    }
                }
                if(bowed) break;
            }
            if(bowed) continue;
        }
        if(routeOnMesh(a,b,turn)){
            double length=0.0;
            Vec3 previous=a;
            for(const auto& p : turn){ length+=len(p-previous); previous=p; }
            scoreOut -= gap + length+len(b-previous)-len(b-a);
            ++meshDetours;
            continue;
        }
        fprintf(stderr, "[曲面] 第 %zu 段与下一段没有满足净空的面内连接\n", r);
        return false;
    }
    double turnCost=0.0;
    for(size_t r=1; r<starts.size(); ++r){
        Vec2 prev=toUV(bestPoints[ends[r-1]-1]);
        Vec2 direction=norm2(prev-toUV(bestPoints[ends[r-1]-2]));
        auto addTurn = [&](const Vec3& p){
            const Vec2 next=toUV(p), delta=next-prev;
            if(len2(delta)<step*0.1) return;
            const Vec2 dir=norm2(delta);
            // 惩罚急拐角，自动方向优先选能接成常规折返的排布。
            turnCost+=std::max(0.0,0.5-dot2(direction,dir));
            prev=next;
            direction=dir;
        };
        for(const auto& p : turns[r]) addTurn(p);
        addTurn(bestPoints[starts[r]]);
        addTurn(bestPoints[starts[r]+1]);
    }
    scoreOut-=gap*turnCost;
    sp.rowP.clear();
    sp.rowBeg.assign(1,0);
    for(size_t r=0; r<starts.size(); ++r){
        sp.rowP.insert(sp.rowP.end(),bestPoints.begin()+starts[r],bestPoints.begin()+ends[r]);
        sp.rowBeg.push_back((int)sp.rowP.size());
    }
    sp.rowN.assign(sp.rowP.size(), sp.n0);
    sp.rowTurns.swap(turns);
    sp.nRow = (int)sp.rowBeg.size()-1;
    sp.rowGap = gap;
    if(std::getenv("SPRAYTRAJ_DUMP"))
        fprintf(stderr,"[曲面] 对齐行端 %d 处，网格绕行 %d 处，布行评分 %.1f\n",
                alignedTurns,meshDetours,scoreOut);
    return true;
}

// 用紧公差顺序 RANSAC 剥离，检查少量平面能否覆盖曲面候选块。
static bool patch_is_polyhedral(const std::vector<Vec3>& P, const std::vector<int>& blk,
                                int nmax, double tolPlane, double coverReq,
                                int& nPlaneOut, double& covOut,
                                std::vector<Vec3>* planeN = nullptr)
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
        if(planeN) planeN->push_back(bn);   // 供"法向是否成对"判据使用
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
//     max(4*curvatureTol, 2*重建尺度) 有 4 倍以上余量，不会误伤其他工件。
//   * reconGap 是重建几何尺度，**不是工艺行距**。以前这里传 scanRowSpacing，
//     行距 60 时门限被抬到 120mm，曲型板 p95=74 反而被判成"平"的，
//     整块板被拆成十几个小平面——重建结果随布行参数变化。
static bool detect_curved_patch(const std::vector<Vec3>& P, const CloudNormals& cn,
                               const Vec3& partC, double curvatureTol, double reconGap,
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
    const double thr = std::max(4.0*curvatureTol, 2.0*reconGap);
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
        // 公差取 0.2×重建尺度（12mm → 2.4mm）：既容得下现场点云的墙面噪声（~0.5mm），
        // 又远小于"浅曲面被切成薄带"的厚度，避免把 sag 30mm 级的缓曲面误判成多面体
        // （sag 30mm 的板在 2.4mm 带宽下要 12 个平面才铺满，nmax=8 只够盖 ~75%）。
        const double tolPl = std::min(3.0, std::max(1.5, 0.20*reconGap));
        int np = 0; double cov = 0.0;
        std::vector<Vec3> pn;                       // 各剥离平面的法向
        bool poly = patch_is_polyhedral(P, outInl, 8, tolPl, 0.93, np, cov, &pn);
        fprintf(stderr, "[曲面] 多平面校验: 紧公差(%.1fmm)剥离 %d 个平面覆盖 %.1f%%%s\n",
                tolPl, np, cov*100.0,
                poly ? " → 多面体，退回平面分割" : "（未达 93%）");
        if(poly) return false;

        // 闭合多面体的平面法向必然成对出现（±X/±Y/±Z），曲面薄带的法向则是单向连续铺开的。
        // 盒子/棱柱这类件面数常超过 nmax，覆盖率天然到不了 93%，光看覆盖率会把整个盒子
        // 判成"曲面"，然后用粗网格去拟合 —— 面片就会皱成一片（薄壁盒 86% 点被吃进曲面块）。
        // 这里补一条"成对轴"判据放行，而不是简单调低覆盖率门限：调低门限会把 sag 30mm 级的
        // 缓曲面（切成 12 个薄带就能盖 95%）也误判成多面体。
        if(np >= 3 && cov >= 0.75){
            int paired = 0;
            std::vector<char> used(pn.size(), 0);
            for(size_t i=0;i<pn.size();++i){
                if(used[i]) continue;
                for(size_t j=i+1;j<pn.size();++j){
                    if(used[j]) continue;
                    if(dot(pn[i], pn[j]) < -0.95){    // 反向平行 = 同一根轴的两面
                        ++paired; used[i]=used[j]=1; break;
                    }
                }
            }
            if(paired >= 2){
                fprintf(stderr, "[曲面] 多平面校验: %d 对反向平行面(≥2) → 闭合多面体，退回平面分割\n", paired);
                return false;
            }
        }
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
    //   ② 这些点在垂直于 nf 的平面上最小宽度 ≥ 3×重建尺度 —— 排除弯曲板被 2.4mm
    //      薄层切出的窄条（曲型板实测宽度只有 4~5mm）。
    // 循环最多 3 轮：一个件可能同时有顶部唇边和底面两个平盘。
    for(int pass=0; pass<3; ++pass){
        // 公差固定 3.0mm：曲型板的两个平盘（侧板）只有在这个公差下才会被完整
        // 剥出来；退回 2.4mm 时剥不干净，整块板会被当成一个曲面块，侧板丢失。
        // 该值只描述"平盘/曲面"的几何分界，与工艺行距无关。
        const double tolPl = 3.0;
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
        const bool wide = (wmin >= 3.0*reconGap);
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
    p.gunClear       = 20.0;   // 枪帽 Ø32 半径16 + 余量4；大枪/旋杯按实际改
    p.edgeExtend     = -1.0;   // 行端点外扩：-1 = 自动(1×行距)，0 = 不外扩
    p.edgeExtendAcross = -1.0; // 跨行方向外扩：-1 = 自动(0.25×行距)，0 = 首末行压边界
    p.shellThick     = 0.0;    // 薄壳双面壁厚：0 = 关闭（不推断另一面）
    // 空移让开方向：默认 0 = 沿用旧口径（沿内部 +Z 让开，等价"工件朝上、相机俯拍"）。
    // 侧装/立式悬挂工件必须由调用方填「基座 +X 在点云系里的方向」（手眼 R 第 1 行），
    // 否则会从工件顶上翻过去，一次过渡走出一两米（见 spray_traj.h 的字段注释）。
    p.liftDirX = 0.0; p.liftDirY = 0.0; p.liftDirZ = 0.0;
    // 最小可喷面宽度：短边 < 30mm 判为线状/退化面（0=关闭）。
    // 2026-11-23 从 18 提到 30：批量实测里 12~16mm 那批"面"根本是**挂具挂钩**
    // （p003 宽 12.0mm、p013 宽 16.0mm、p007 宽 12.0mm），不是工件本体。
    // 给挂钩出轨迹既浪费节拍，又会把枪伸到挂具缝隙里去。30mm 能滤掉挂钩，
    // 又保住 20mm 级窄筋（p016 那条折边 20.6mm 走的是吸收通道，不靠这道门）。
    p.minFaceWidth = 30.0;
    // 换面不关枪长度：80mm ≈ 0.8×standoff。短于它的换面直接连过去、枪不停；
    // 更长的换面仍然关枪（停枪抬枪绕行）。
    p.crossFaceGunOnLen = 80.0;
    // 换面不关枪的枪轴转角上限：60°（现场折角件常见 70~90°，按件调）。
    // 0/负 = 不判夹角（只看距离、净空、穿面）。
    p.crossFaceMaxTurnDeg = 60.0;
    // 喷幅宽度：0 = 未提供 → 退回"由行距反推"的旧口径（W = 行距/(1−搭接率)）。
    // 填了就走新口径：行距是"选出来的"，喷幅才是不随行距变的物理量。
    p.fanWidth = 0.0;
    // 有效喷距：0 = 不判（只按横向喷幅半径判覆盖，旧口径）。
    // 填了之后"远了喷不到"才真正生效：喷距 = standoff + 点相对行点沿枪轴的弧高。
    p.maxSprayDist = 0.0;
    return p;
}

// ==================== 可喷侧参考方向 up：入口旋转 / 出口逆旋转 ====================
// 为什么需要这一层：整套算法内部把 +Z 当成"工件正面朝向"，有三处硬编码依赖它：
//   1) 初始定向 "|n.z| > 0.05 → 归 +Z 半球"；
//   2) angZ/tilt 按与 +Z 的夹角算，可喷窗口 130° 也是相对 +Z；
//   3) 安全高度 zSafe = 工件最高 z + 1.3×standoff，空移沿 +Z 让开。
//
// 而输入点云**来自相机系**，那里的事实是固定的、不需要人工选：
//   · 相机系 +Z 是光轴、指向景物（远离相机）；像素/(u,v,depth) 反投影就得到这个约定；
//   · 相机只能看到朝向它的那一面，所以可见面外法向必然 n·(+Z) < 0，即 n_z < 0；
//   · 相机与枪同装在一个腕上（eye-in-hand），枪只能从相机这一侧伸过来。
//   ⇒ "可喷侧 = 朝相机侧 = 更低 Z"是「看得见才喷得到」的必然结果，不是工位参数。
// 旧代码把这个事实当成"工件朝上、相机在上方"（+Z 半球），于是可喷面被翻到背面，
// 枪位整体跑到工件后面（实测背板_X+800：轨迹 z 中位数比点云高 78mm，k 朝 −Z）。
//
// 逐处改那三处要动二十几处、每处都得重测；改成**入口做一次刚体旋转**：把 up 方向
// 转到 +Z，内部一行不动，出口再把轨迹点/枪轴/平面法向与顶点逆变换回原坐标系。
// 旋转 det=+1，不是镜像：叉积、行方向、枪口朝向都保持右手系，安全高度也自动落到
// 正确一侧（转完以后"最高 z"就是沿 up 最远的点）。
//
// 默认 up = (0,0,-1)：相机系恒定不变的那个方向。
// 仅当点云已被搬到别的坐标系（例如手眼变换后的机器人基座系，+Z 竖直向上）时，
// 才用环境变量 SPRAYTRAJ_UP="x,y,z" 覆盖 —— 那是换了个坐标系，不是换了个工位习惯。
//   相机在上方（工件朝上）       → 不设 / 0,0,1
//   相机朝工件、朝相机的方向 Z 低 → 0,0,-1
//   侧装（相机在 -Y 侧）         → 0,-1,0
// 轴角旋转（Rodrigues）：v 绕单位轴 axis 转 ang 弧度。右手系，det=+1，不镜像。
static Vec3 rot_axis(const Vec3& v, const Vec3& axis, double ang){
    const double c = std::cos(ang), s = std::sin(ang);
    return v*c + cross(axis, v)*s + axis*(dot(axis, v)*(1.0-c));
}

struct UpFrame {
    bool   active = false;
    Vec3   axis{1,0,0};
    double ang = 0.0;
    Vec3   up{0,0,1};
    Vec3 fwd(const Vec3& v) const { return active ? rot_axis(v, axis,  ang) : v; }
    Vec3 inv(const Vec3& v) const { return active ? rot_axis(v, axis, -ang) : v; }
};

static UpFrame read_up_frame(){
    UpFrame f;
    // 相机系事实：可喷侧永远是朝相机侧，也就是 −Z。不需要调用方传、也不该让人去选。
    f.up = Vec3{0,0,-1};
    if(const char* e = std::getenv("SPRAYTRAJ_UP")){
        double a=0,b=0,c=0;
        if(std::sscanf(e, "%lf,%lf,%lf", &a,&b,&c) == 3){
            Vec3 v{a,b,c};
            if(len(v) > 1e-9) f.up = norm(v);
        }
    }
    const Vec3 z{0,0,1};
    const double c0 = std::max(-1.0, std::min(1.0, dot(f.up, z)));
    if(c0 > 0.999999) return f;                       // up 已经是 +Z：点云已在"正面朝 +Z"的坐标系
    f.active = true;
    f.ang = std::acos(c0);
    const Vec3 ax = cross(f.up, z);                   // 把 up 转到 +Z 的转轴
    if(len(ax) < 1e-9){ f.axis = Vec3{1,0,0}; f.ang = 3.141592653589793; }  // up = -Z
    else f.axis = norm(ax);
    return f;
}

// 注意：结果**写进出参**，不按值返回。
// SprayTrajResult 约 289KB；C# P/Invoke 调用时托管侧已经占了一块返回缓冲，
// 若这一层再多一个按值返回的 289KB 栈帧，UI 线程的 1MB 栈会被直接顶爆
// （2026-11-12 实测：加了一层 retry 包装后，背板 10.6 万点调用 Stack overflow）。
// 2D 离群小簇过滤：布行前剔除与主面断开的极小孤立簇。
// 背景：p001 重建后 P2 里残留 19 个尾端离群点（距主面 147mm），
// minAreaRect 走凸包，19 个点就能把矩形拉成 248×112@134°，行方向随之歪掉。
// 只剔除“同时满足”两个条件的小簇：
//   1) 点数 < P2 的 smallFrac（默认 1%）——大碎片/真实小面不动；
//   2) 到主连通域的最近距离 > farGap（默认 3×重建间隙=36mm）——贴着主面的
//      缺口/小舌不动。
// 这样 p001 的 19 点（0.5%、147mm）会被清掉，而真实碎片面、窄条上的小分叉保留。
static void filterP2Outliers(std::vector<Vec2>& pts, double link, double farGap,
                             double smallFrac, double maxDropFrac, size_t faceIdx)
{
    const int n = (int)pts.size();
    if(n < 40) return;
    // 1) 按 link 连通，标记所有连通域
    std::vector<int> par(n);
    std::iota(par.begin(), par.end(), 0);
    std::function<int(int)> find = [&](int x){
        while(par[x] != x){ par[x] = par[par[x]]; x = par[x]; }
        return x;
    };
    auto uni = [&](int a, int b){
        a = find(a); b = find(b);
        if(a != b) par[a] = b;
    };
    const double cell = std::max(1e-6, link);
    std::map<std::pair<int,int>, std::vector<int>> grid;
    std::vector<int> ix(n), iy(n);
    for(int i=0;i<n;++i){
        ix[i] = (int)std::floor(pts[i].x / cell);
        iy[i] = (int)std::floor(pts[i].y / cell);
        grid[{ix[i], iy[i]}].push_back(i);
    }
    const double l2 = link * link;
    for(int i=0;i<n;++i){
        for(int dx=-1; dx<=1; ++dx)
        for(int dy=-1; dy<=1; ++dy){
            auto it = grid.find({ix[i]+dx, iy[i]+dy});
            if(it == grid.end()) continue;
            for(int j : it->second){
                if(j <= i) continue;
                const double ddx = pts[i].x - pts[j].x;
                const double ddy = pts[i].y - pts[j].y;
                if(ddx*ddx + ddy*ddy <= l2) uni(i, j);
            }
        }
    }
    // 2) 按根分组，找主连通域
    std::map<int, std::vector<int>> comps;
    for(int i=0;i<n;++i) comps[find(i)].push_back(i);
    if(comps.size() <= 1) return;
    int mainRoot = -1; size_t mainCnt = 0;
    for(auto& kv : comps)
        if(kv.second.size() > mainCnt){ mainCnt = kv.second.size(); mainRoot = kv.first; }

    // 3) 主域点建 farGap 网格，用于近似最近距离
    const double fg = std::max(1e-6, farGap);
    std::map<std::pair<int,int>, std::vector<int>> mainGrid;
    for(int i : comps[mainRoot]){
        int cx = (int)std::floor(pts[i].x / fg);
        int cy = (int)std::floor(pts[i].y / fg);
        mainGrid[{cx, cy}].push_back(i);
    }
    std::vector<char> drop(n, 0);
    size_t dropComps = 0, dropPts = 0;
    const double fg2 = fg * fg;
    for(auto& kv : comps){
        if(kv.first == mainRoot) continue;
        const std::vector<int>& c = kv.second;
        const bool small = (double)c.size() < smallFrac * (double)n;
        if(!small) continue;
        bool nearMain = false;
        for(int i : c){
            int cx = (int)std::floor(pts[i].x / fg);
            int cy = (int)std::floor(pts[i].y / fg);
            for(int dx=-1; dx<=1 && !nearMain; ++dx)
            for(int dy=-1; dy<=1 && !nearMain; ++dy){
                auto it = mainGrid.find({cx+dx, cy+dy});
                if(it == mainGrid.end()) continue;
                for(int j : it->second){
                    const double ddx = pts[i].x - pts[j].x;
                    const double ddy = pts[i].y - pts[j].y;
                    if(ddx*ddx + ddy*ddy <= fg2){ nearMain = true; break; }
                }
            }
            if(nearMain) break;
        }
        if(!nearMain){
            for(int i : c) drop[i] = 1;
            ++dropComps; dropPts += c.size();
        }
    }
    if(dropPts == 0) return;
    // 总剔除量上限：一个面里如果“小簇”加起来已经占了不少（例如无边框件
    // 边缘碎成很多小片），说明这是真实碎片/欠采样，不是孤立脏点；整面保留，
    // 交给 planeFrag / 轮廓切段。p001 只删 19/3726=0.5%，不受影响。
    if((double)dropPts > maxDropFrac * (double)n){
        fprintf(stderr, "[traj] 面%zu: P2 离群小簇合计 %zu 点（%.1f%% > %.0f%%），保留原 P2\n",
                faceIdx, dropPts, 100.0 * (double)dropPts / (double)n,
                maxDropFrac * 100.0);
        return;
    }

    std::vector<Vec2> kept; kept.reserve((size_t)n - dropPts);
    for(int i=0;i<n;++i) if(!drop[i]) kept.push_back(pts[i]);
    fprintf(stderr, "[traj] 面%zu: P2 剔除离群小簇 %zu 个 / %zu 点"
                    "（单簇<%0.1f%%、距主面 >%.0fmm）  %zu -> %zu\n",
            faceIdx, dropComps, dropPts, smallFrac * 100.0, fg,
            (size_t)n, kept.size());
    pts.swap(kept);
}

enum class TrajectoryFailure { Other, NoRows, AngleWindow, CandidatePath };

static void spray_traj_core(
    const double* xyz, int n, const SprayTrajParams* params, SprayTrajResult* out,
    TrajectoryFailure& failure, bool enableCandidateCheck)
{
    failure = TrajectoryFailure::Other;
    SprayTrajResult& res = *out;
    std::memset(out, 0, sizeof(*out));
    res.error = nullptr;
    ProfStage psTotal("0.总耗时");   // 任意返回路径都会打印

    std::string err;
    if(!xyz || n<10){ err="点云点数不足(<10)"; res.error=_strdup(err.c_str()); return; }

    // 可喷侧参考方向：非 +Z 时先把点云整体旋转到"up = +Z"，出口再转回（见 UpFrame 注释）
    const UpFrame uf = read_up_frame();
    if(uf.active)
        fprintf(stderr, "[坐标系] 点云按相机系处理：可喷侧 up=(%.2f,%.2f,%.2f) → 入口旋转 %.1f° 到 +Z，出口转回（轨迹仍在原坐标系）\n",
                uf.up.x, uf.up.y, uf.up.z, uf.ang*180.0/3.141592653589793);

    std::vector<Vec3> P(n);
    for(int i=0;i<n;++i) P[i]=uf.fwd({xyz[i*3],xyz[i*3+1],xyz[i*3+2]});
    res.input_points=n;
    fprintf(stderr, "[参数] 工艺行距=%.2fmm（只用于布行） 重建几何尺度=%.1fmm（固定，与工艺行距解耦）\n",
            params->scanRowSpacing, kReconGap);

    // ---- 体素降采样：只服务"布行 + 安全判据"，不再参与分割/重建 ----
    // 用户 2026-11-20 15:46 的口径：
    //   "重建平面应该不需要体素过滤那么严重吧，体素过滤主要还是为了轨迹规划和
    //    安全层服务"
    // 实测支持这条：p047 现场件 226114 点，全分辨率分割出 **2 个面**，
    // 2mm 体素后出 2 面 + 1 线状 + 1 碎面（碎面还得靠环级丢弃才拦住）。
    // 也就是说那个"碎面"根本不是碎面，是体素把主面撕下来的一块 —— 碎面过滤
    // 在很大程度上是在给体素擦屁股。全分辨率下它被主面正常吸收。
    //
    // 所以分工改成：
    //   分割 / 重建 / 碎面判定 / 轮廓   ← Pf（全分辨率，只做 SOR）
    //   布行（P2 的原子）/ 枪位净空 / 空移直线可达   ← 体素代表点（够用即可）
    // 轨迹阶段只用到点坐标（P2 投影、枪位抬升）与最近点查询，两者在 2mm 网格上
    // 与全分辨率逐位等价到亚毫米，但 KD 树小 3 倍、枪位过滤快 3 倍。
    //
    // 阈值 12 万点：小件（tray 1.4 万 / complex 8.4 万）完全不走这条路，行为逐位不变。
    // 可用 SPRAYTRAJ_VOXEL=<mm> 显式指定，=0 强制关闭。
    size_t nAfterVoxel = P.size();
    std::vector<int> planIdx;      // 布行/安全层用的点下标（升序，指向 Pf）
    double planLeaf = 0.0;         // 0 = 不抽稀（小件走全量，行为逐位不变）
    {
        if(const char* e = std::getenv("SPRAYTRAJ_VOXEL")) planLeaf = std::atof(e);
        else if(g_voxelLeafOverride > 1e-6) planLeaf = g_voxelLeafOverride;   // 兜底重试档
        else if(n > 120000) planLeaf = auto_voxel_leaf(params);
    }
    nAfterVoxel = P.size();

    double tSor=0, tSplit=0, tMesh=0, tTraj=0;
    // 1. 去噪
    std::vector<Vec3> Pf;
    {
        ProfStage ps("1.SOR去噪", &tSor);
        Pf = sor_filter(P, params->sorK, params->sorZ);
    }
    res.removed_points = (int)nAfterVoxel - (int)Pf.size();   // 只统计 SOR 剔除，不含体素降采样
    if(Pf.size()<10){ err="去噪后点数不足"; res.error=_strdup(err.c_str()); return; }
    {
        double zmin=Pf[0].z, zmax=Pf[0].z;
        for(auto&p:Pf){ zmin=std::min(zmin,p.z); zmax=std::max(zmax,p.z); }
        fprintf(stderr, "[去噪] 点数 %zu->%zu, Z范围[%.1f,%.1f]\n", P.size(), Pf.size(), zmin, zmax);
    }

    // ---- 布行/安全层点集：在 Pf 上取体素代表点（下标），不复制坐标 ----
    // 与 Pf 同索引空间，所以 seg.inl 无需改写；轨迹阶段按 planIdx 过滤即可。
    // 小件（planLeaf<=0）时 planIdx 留空，全部退回"全量"，行为逐位不变。
    if(planLeaf > 1e-6 && Pf.size() >= 20000){
        ProfStage ps("0b.布行层抽稀");
        planIdx = voxel_downsample_idx(Pf, planLeaf);
        // 太少就退回全量：抽稀后不足 500 点说明这件本来就稀疏，抽了没有意义。
        if(planIdx.size() < 500) planIdx.clear();
        else fprintf(stderr, "[体素] 布行层 %.1fmm 代表点 %zu/%zu (%.0f%%)；"
                             "分割/重建/碎面判定用全分辨率 %zu 点\n",
                     planLeaf, planIdx.size(), Pf.size(),
                     100.0*(double)planIdx.size()/(double)Pf.size(), Pf.size());
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
    // 法向树提到块外：面级朝向诊断要沿法向打射线，必须在"统一法向定向"之后还能用。
    // 否则诊断拿到的是定向前的随机符号，n0 和最终 nx 会对不上。
    CloudNormals cn;
    {
        ProfStage ps("2.平面分割", &tSplit);
        cn.build(Pf);

        Vec3 partC{0,0,0};
        for(const auto& p : Pf) partC = partC + p;
        partC = partC*(1.0/(double)std::max<size_t>(1,Pf.size()));

        std::vector<int> curvedInl;
        isCurved = detect_curved_patch(Pf, cn, partC, params->curvatureTol,
                                       kReconGap, surf, curvedInl);
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
                                        std::max(20, (int)(Pf.size()*0.02)), params->scanRowSpacing,
                                        params->standoff);
                for(auto& s : sub){
                    for(int& k : s.inl) k = back[k];
                    for(int& k : s.fitSupport) k = back[k];
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
            segs = split_planes(Pf, cn, params->flatnessTol, params->curvatureTol,
                                std::max(20, (int)(Pf.size()*0.02)), params->scanRowSpacing,
                                params->standoff);
        }
    }

    // 若分割结果为空（极少见），退回单平面
    if(segs.empty()){
        Plane pl; std::vector<int> inl;
        // 兜底路径同样不能用 flatnessTol 当拟合容差（理由见 split_planes 里的尺度说明）
        const double fitTolFb = std::min(8.0, std::max(3.0, 3.0*cloud_noise_sigma(Pf, cn)));
        if(!fit_plane(Pf, fitTolFb, pl, inl)){ err="无法重建有效平面"; res.error=_strdup(err.c_str()); return; }
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
    // ⚠ 而"独立朝向"豁免本身又会放过另一类东西：**两个面的交棱处的过渡带**。
    // 现场实测（p006 tray 件）：0.1° 主面 + 51°/67° 两个侧壁的交角附近，区域生长切出
    // 一条 379 点的窄带，拟合法向夹在主面与侧壁之间（51° 偏差），因此命中"独立朝向"
    // 豁免被保留。但这条带子的真实几何是 165.8 × 5.5 mm —— 宽只有 5.5mm，
    // 比重建栅格（kReconGap=12mm）还窄，重建时整条带被压成一个 10×10 的方格
    // （重建面积 100mm²，而按点密度本该 ~900mm²），随后在这个方格上布 1 行、
    // 并为此多出 2 处面间过渡（额外抬枪空移）。视觉上就是"大面上贴着一小块垃圾"。
    //
    // 判据（两条同时成立才丢，缺一不可）：
    //   ① 面内**窄向宽度** < 重建栅格 kReconGap。不是调参数字，而是"这套重建分辨率
    //      能不能表达该面"的硬下限：窄向跨度不足一个栅格时该面必然退化成单个 cell
    //      方块（本例重建面积 10×10=100mm²，而按点密度本该 ~900mm²），轮廓与布行
    //      都已失去几何意义。
    //   ② 它是**冗余**的：该面绝大多数点落在另一个更大面的平面附近（中位距离很小），
    //      且法向与那个面只差十几度 —— 即它本来就是那个面的一部分，只是被区域生长
    //      沿条带切了出来（本例：与面1法向仅差 16.2°，点到面1中位距离 1.6mm）。
    //   加②的理由：只按宽度会误杀**真实窄结构**——窄筋/窄翻边同样比栅格窄，
    //   但它们立在工件本体上（离邻近平面十几到几十毫米）、法向与邻面近乎垂直，
    //   ②不成立，于是被保住。两条 AND 起来，误删面只会来自"与邻面共面且窄于栅格"，
    //   这种面本来就喷不出独立内容。
    //   仍然只对靠 distinct 豁免才活下来的面（!big）生效：占比 ≥3% 的窄面可能是
    //   工件本体或主面，误删会整件漏喷，不放行这条。
    auto inPlaneExtents = [](const std::vector<Vec3>& P, const std::vector<int>& idx,
                             Vec3 n, double& wNarrow, double& wLong)->bool{
        if(idx.size() < 12) return false;
        n = norm(n);
        Vec3 ax = std::fabs(n.z) < 0.9 ? Vec3{0,0,1} : Vec3{1,0,0};
        Vec3 u = norm(cross(n, ax));
        if(len(u) < 0.5) return false;
        Vec3 v = cross(n, u);
        double mu=0, mv=0;
        for(int k : idx){ mu += dot(P[k],u); mv += dot(P[k],v); }
        mu /= (double)idx.size(); mv /= (double)idx.size();
        double suu=0, suv=0, svv=0;
        for(int k : idx){
            double a = dot(P[k],u)-mu, b = dot(P[k],v)-mv;
            suu += a*a; suv += a*b; svv += b*b;
        }
        const double th = 0.5*std::atan2(2.0*suv, suu-svv);
        const double c = std::cos(th), sn = std::sin(th);
        std::vector<double> pa, pb;
        pa.reserve(idx.size()); pb.reserve(idx.size());
        for(int k : idx){
            double a = dot(P[k],u)-mu, b = dot(P[k],v)-mv;
            pa.push_back( a*c + b*sn);   // 面内主轴
            pb.push_back(-a*sn + b*c);   // 面内次轴
        }
        // 2%~98% 分位跨度：抗端点毛刺（掠射边缘会甩出几个孤立点把极值带飞）
        auto spread=[](std::vector<double>& z)->double{
            const size_t lo = (size_t)(z.size()*0.02);
            const size_t hi = std::min(z.size()-1, (size_t)(z.size()*0.98));
            std::nth_element(z.begin(), z.begin()+lo, z.end());
            const double a0 = z[lo];
            std::nth_element(z.begin()+lo, z.begin()+hi, z.end());
            return z[hi]-a0;
        };
        const double a1 = spread(pa), b1 = spread(pb);
        wLong   = std::max(a1,b1);
        wNarrow = std::min(a1,b1);
        return true;
    };

    if(segs.size() > 1){
        size_t maxIdx=0; for(size_t i=1;i<segs.size();++i)
            if(segs[i].inl.size() > segs[maxIdx].inl.size()) maxIdx=i;
        double total = (double)Pf.size();
        const Plane& mainPl = segs[maxIdx].pl;
        const size_t absFloor = std::max<size_t>(80, (size_t)(total*0.005));
        // 元数据快照：下面的 keep.push_back(std::move(s)) 会把已判定的面的 inl 搬空，
        // 之后再回头查 segs[j].inl.size() 全是 0 —— "找宿主面"就永远找不到。
        struct SegMeta { size_t n; Plane pl; };
        std::vector<SegMeta> meta(segs.size());
        for(size_t i=0;i<segs.size();++i){ meta[i].n = segs[i].inl.size(); meta[i].pl = segs[i].pl; }
        // 排查用：把每张面的"占比/夹角/面内尺寸/是否 big/distinct"打出来。
        // 只看不判，判定逻辑在下面；没有这条日志时无法解释"为什么这个碎块没被丢"。
        if(const char* d=std::getenv("SPRAYTRAJ_DUMP")) if(d[0]=='1'){
            fprintf(stderr, "[dbg] 碎面判定: 参与分割总点数=%zu 面数=%zu\n", Pf.size(), segs.size());
            for(size_t i=0;i<meta.size();++i){
                double wN=0, wL=0;
                bool okE = inPlaneExtents(Pf, segs[i].inl, meta[i].pl.n, wN, wL);
                double ratio = (double)meta[i].n/(double)Pf.size();
                double angMain = std::acos(std::min(1.0,std::fabs(dot(meta[i].pl.n,mainPl.n))))*180.0/3.141592653589793;
                bool bigI = ratio >= 0.03;
                bool disI = (angMain > 40.0) && (meta[i].n >= absFloor);
                fprintf(stderr, "[dbg]   面%zu 点数=%zu 占比=%.2f%% 与主面夹角=%.1f° 面内=%s %.1f x %.1fmm big=%d distinct=%d\n",
                        i, meta[i].n, ratio*100.0, angMain, okE?"":"?", okE?wN:-1.0, okE?wL:-1.0,
                        (int)bigI, (int)disI);
            }
        }
        std::vector<PlaneSeg> keep;
        // 碎面口径（用户 2026-11-20 拍板）：命中**任一**碎面特征即丢，不再要求三条全中。
        //   F1 点数占比 < 3%
        //   F2 面内窄向宽度 < 重建栅格 kReconGap（重建时必然压成单格，轮廓与布行都失去意义）
        //   F3 贴住某个更大的同向面（夹角<25° 且 点到该面中位距离<=5mm）→ 被切出来的冗余带
        // 唯一豁免：与主面夹角>40°（独立朝向）且 窄向宽度>=栅格 且 不贴任何同向大面。
        //   依据是历史实测：cao 槽型盘 66.1° 侧壁 578 点/占比 2.75%，纯按占比会被误删 → 整面漏喷。
        //   这不是"放宽口径"，而是"独立朝向 + 能重建出几何"本身就不构成碎面。
        // 回退：SPRAYTRAJ_FRAG_AND=1 走旧的全中口径，便于 A/B。
        const bool fragAnd = (std::getenv("SPRAYTRAJ_FRAG_AND") != nullptr);
        // 碎面短边阈值：与重建栅格对齐（重建表达不了一个栅格宽以下的面）。
        // 用户 2026-11-20：独立朝向面的豁免短边阈值要可调 → SPRAYTRAJ_FRAG_SHORT 覆盖。
        double fragShort = kReconGap;
        if(const char* e = std::getenv("SPRAYTRAJ_FRAG_SHORT")){
            const double v = std::atof(e);
            if(v > 0.0) fragShort = v;
        }
        // 点数占比阈值（用户 2026-11-20 14:19 拍板：从 3% 拉到 10%）。
        // 理由：轨迹要尽量简洁，不是所有小面都值得喷——只出一两行、又和主面近共面的小面，
        // 喷了只是给机器人加一堆换行/抬枪，膜厚收益可以忽略。
        // 3% 那条在 p047 上正好被"2593/79270=3.3%"擦边放行，属于阈值定得太窄。
        double fragShare = 0.10;
        if(const char* e = std::getenv("SPRAYTRAJ_FRAG_SHARE")){
            const double v = std::atof(e);
            if(v > 0.0) fragShare = v;
        }
        // "近共面"的夹角上限：超过这个角度算独立朝向（侧壁/翻边），归豁免管。
        double fragAng = 40.0;
        if(const char* e = std::getenv("SPRAYTRAJ_FRAG_ANG")){
            const double v = std::atof(e);
            if(v > 0.0) fragAng = v;
        }
        size_t fragDropPts = 0, nFragDrop = 0;
        for(size_t i=0;i<segs.size();++i){
            PlaneSeg& s = segs[i];
            const double ratio = (double)meta[i].n/total;
            const double angMain = std::acos(std::min(1.0, std::fabs(dot(meta[i].pl.n, mainPl.n))))*180.0/3.141592653589793;
            const bool big      = ratio >= 0.03;
            const bool major    = ratio >= fragShare;                   // 10%：够一块"面"的分量
            const bool distinct = (angMain > 40.0) && (meta[i].n >= absFloor);
            double wN=0, wL=0;
            const bool okE = inPlaneExtents(Pf, s.inl, meta[i].pl.n, wN, wL);

            // F3：找"更大且近同向"的宿主面，并量点到它的中位距离
            size_t hostIdx = meta.size();
            double hostAng = 0.0, dMed = 1e18;
            for(size_t j=0;j<meta.size();++j){
                if(j==i || meta[j].n <= meta[i].n) continue;
                double a = std::acos(std::min(1.0, std::fabs(dot(meta[j].pl.n, meta[i].pl.n))))*180.0/3.141592653589793;
                if(a < 25.0){ hostIdx=j; hostAng=a; break; }
            }
            if(hostIdx < meta.size()){
                const Plane& hp = meta[hostIdx].pl;
                std::vector<double> dd; dd.reserve(s.inl.size());
                for(int k : s.inl) dd.push_back(std::fabs(dot(Pf[k], hp.n) + hp.d));
                std::nth_element(dd.begin(), dd.begin()+dd.size()/2, dd.end());
                dMed = dd[dd.size()/2];
            }
            // 三条独立判据（用户 2026-11-20 定稿口径）
            const bool d1 = okE && (wN < fragShort);                   // 短边 < 阈值（无条件丢）
            // D5：点数占比 < fragShare(默认10%) **且** 与主面近共面（夹角 < fragAng，默认40°）。
            // 加"近共面"这一维是用户 14:19 的要求：小面如果独立朝向（侧壁/翻边/窄筋），
            // 主面的喷锥打不到，必须自己走；但**和主面几乎共面**的小面，主面布行时早就
            // 连带覆盖了，再给它单独出一两行只是白加换行/抬枪。
            const bool d5 = (!major) && (angMain < fragAng);
            const bool d6 = (hostIdx < meta.size()) && (hostAng < 25.0) && (dMed <= 5.0);  // 夹角小且距离近
            const bool isMainFace = (i == maxIdx);
            // 豁免：独立朝向(>40°) 且 短边 >= 阈值 → 窄筋/翻边/侧壁，真实结构
            const bool realStruct = (angMain > 40.0) && okE && (wN >= fragShort);

            // ---- 环形 / 类环形豁免（用户 2026-11-20 14:23：要考虑环形、类环形）----
            // D5「占比小 + 近共面」会误杀环带：一圈凸筋、方框焊环、跑道环的面积确实小、
            // 法向也确实和主面平行（很多还是同向抬高的筋），但它的轨迹是**一圈闭合路径**，
            // 不是一两行直线；按碎面丢掉 = 整圈漏喷。D6「贴着更大的同向面」同样会误杀
            // 抬高的细环筋（环贴着主面、中位距离只有 1~3mm）。
            // 判据与轨迹侧的 ring/loop 同源（detect_band_uv）：中间是空的、带宽窄、
            // 方位绕成圈 → 是环，保留，交给环形布行。
            bool isRingFace = false;
            if(!isMainFace && (d5 || d6) &&
               !(std::getenv("SPRAYTRAJ_RING_KEEP") && std::atoi(std::getenv("SPRAYTRAJ_RING_KEEP")) == 0)){
                const Vec3 nb = meta[i].pl.n;
                Vec3 t1 = (std::fabs(nb.z) < 0.9) ? cross(nb, Vec3{0,0,1}) : cross(nb, Vec3{1,0,0});
                t1 = norm(t1);
                const Vec3 t2 = cross(nb, t1);
                const Vec3 c0 = meta[i].pl.c;
                std::vector<Vec2> uvq; uvq.reserve(s.inl.size());
                for(int k : s.inl){ const Vec3 q = Pf[k] - c0; uvq.push_back(Vec2{dot(q,t1), dot(q,t2)}); }
                BandUV bi;
                isRingFace = detect_band_uv(uvq, bi);
                if(isRingFace)
                    fprintf(stderr, "[split] 保留环带面: 点数=%zu (占比%.1f%%, 与主面夹角%.1f°) "
                                    "内径%.1f 外径%.1f 带宽%.1f(垂直) 方位%.0f%% → 绕环布行，不按碎面丢\n",
                            meta[i].n, ratio*100.0, angMain,
                            bi.rInMed, bi.rOutMed, bi.wPerp, bi.closure*100.0);
            }

            bool drop = false;
            if(fragAnd){                       // 旧口径（A/B 用）
                const bool f2 = okE && (wN < kReconGap);
                const bool f3 = (hostIdx < meta.size()) && (dMed <= 5.0);
                if(!big && distinct){
                    if(f2 && f3) drop = true;
                    else if(f2) fprintf(stderr, "[split] 保留窄面: 点数=%zu 面内宽=%.1fmm < 栅格%.0fmm，"
                                                "但不存在共面宿主面 → 真实窄结构（窄筋/翻边），不丢\n",
                                        meta[i].n, wN, kReconGap);
                } else drop = (!big && !distinct);
            } else {                           // 新口径：任一命中即丢
                drop = !isMainFace && !isRingFace && (d1 || d6 || (d5 && !realStruct));
            }

            if(drop){
                char why[96]; why[0]=0;
                if(d1) strcat(why, "D1短边<阈值 ");
                if(d5) strcat(why, "D5占比小且近共面 ");
                if(d6) strcat(why, "D6夹角小且距离近 ");
                fprintf(stderr, "[split] 丢弃碎面: 点数=%zu (占比%.1f%%, 与主面夹角%.1f°, "
                                "面内%.1f×%.1fmm, 短边阈值%.0fmm%s) 命中[%s]\n",
                        meta[i].n, ratio*100.0, angMain, okE?wN:-1.0, okE?wL:-1.0, fragShort,
                        (hostIdx<meta.size()) ? "" : " 无同向宿主", why);
                if(d6) fprintf(stderr, "[split]   D6 佐证: 宿主=面%zu 夹角=%.1f° 点到宿主面中位距离=%.2fmm\n",
                               hostIdx, hostAng, dMed);
                fragDropPts += meta[i].n; ++nFragDrop;
                continue;
            }
            if(!big && distinct)
                fprintf(stderr, "[split] 保留独立朝向小面: 点数=%zu (占比%.1f%%, 与主面夹角%.1f° > 40°, "
                                "面内%.1f×%.1fmm) → 真实结构\n",
                        meta[i].n, ratio*100.0, angMain, okE?wN:-1.0, okE?wL:-1.0);
            keep.push_back(std::move(s));
        }
        if(nFragDrop && fragDropPts*100.0/(double)total > 40.0)
            fprintf(stderr, "[split] ⚠ 碎面过滤丢弃 %zu 面 / 占 %.1f%% 的点（>40%%），请复核是否误删真实结构\n",
                    nFragDrop, fragDropPts*100.0/(double)total);
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

    // 面级朝向诊断（内壁/外壁）：必须在上面这轮定向之后跑，
    // 否则拿到的是 PCA 的随机符号，n0 与最终 nx 会对不上。
    // 只填 PlaneSeg::fo，不回写 pl.n —— 轨迹输出与本改动前逐字节一致。
    //
    // 曲面块例外：射线分类是按"整块平均平面法向"测两侧自由距离的，在强曲面上
    // 这条射线会贴着曲面擦过去、1~3mm 就打回自身，于是把"可喷侧"判到背面
    // （曲型板实测：+n 侧 1.1mm 命中 → 判 +n 贴材料 → 翻向背面 → angZ 135.8°
    //  > 130° 窗口 → 93% 点数的整个曲面块一行不出）。曲面块本来就有逐点法向
    // （surf.n0 已定向朝外、surf.N 逐点对齐），朝向以它为准，这里不再复核；
    // "内壁"语义对单层外蒙皮也没有意义，记 src=2（几何二义兜底）不触发内壁翻向。
    for(auto& s : segs){
        if(s.curved && isCurved){
            FaceOrient fo;
            const Vec3& n = s.pl.n;
            fo.n = n;
            measure_gun_free(Pf, cn.tree, s, params->standoff,
                             fo.freePlus, fo.freeMinus, &cn.nor);   // 逐点法向，仅供上报
            fo.angZ = std::acos(std::max(-1.0, std::min(1.0, n.z)))*180.0/3.141592653589793;
            fo.tilt = 90.0 - fo.angZ;
            fo.side = (fo.angZ < 87.0) ? 1 : (fo.angZ > 93.0 ? -1 : 0);
            fo.isInner = 0;
            fo.src = 2;
            s.fo = fo;
            if(std::getenv("SPRAYTRAJ_FACE_DBG"))
                fprintf(stderr, "[face] 曲面块 n=(%.3f,%.3f,%.3f) angZ=%.1f° 逐点法向自由 +%.1f/-%.1f（不做射线定向）\n",
                        n.x, n.y, n.z, fo.angZ, fo.freePlus, fo.freeMinus);
        } else {
            s.fo = classify_face_orientation(Pf, cn.tree, s, partC, params->normalFlip);
        }
    }

    // 内壁的喷距法向必须朝腔内。
    // 上面那轮定向（|nz|>0.05 归 +Z 半球 / 近竖直背离质心）对**单层壳内壁**会把法向
    // 翻到盒外：点在内腔壁上、质心在腔内，"背离质心"就等于朝外。于是清单里出现
    // "分类=内壁、内外反"——判成内壁，枪却停在外表面，喷的其实是外面。
    // 射线结论可靠时（src==1：一侧近挡=材料、另一侧远挡=腔体）以可喷侧 fo.n 为准，
    // 让枪真的进腔。兜底/强制/薄壁拆层的结果不动（它们的符号另有依据）。
    {
        int nFix = 0;
        for(auto& s : segs){
            if(s.fo.src != 1 || s.fo.isInner <= 0) continue;
            if(dot(s.pl.n, s.fo.n) > 0.0) continue;        // 本来就一致
            s.pl.n = s.fo.n;
            s.pl.d = -dot(s.pl.n, s.pl.c);
            if(s.curved && isCurved){
                surf.n0 = surf.n0*(-1.0);
                for(auto& v : surf.N) v = v*(-1.0);
            }
            ++nFix;
        }
        if(nFix) fprintf(stderr, "[fix] 内壁喷距法向翻转朝腔: %d 面\n", nFix);
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

    // ---- 薄壁双层拆段：一块面里同时装着内壁和外壁 ----
    // 根因：分割阶段用无向法向(|dot|)做区域生长，壁的两层法向相反却被判成同向，
    // 于是并成一块、块质心落在两层中间；射线沿两侧各走一个壁厚就撞到对面层，
    // 两侧自由距离都 ≈ 壁厚且对称，靠距离永远分不出内外。
    //
    // 必须在**这里**拆：轨迹是按 seg 生成的，只在输出层拆行的话内壁永远分不到
    // 轨迹（之前就是这样——清单里有内壁、橙色面片上一根线都没有）。拆成两个
    // 独立 seg 后，内壁和普通平面一样走 2D 布行，区别只是法向朝腔内。
    //
    // 判据一（主）：沿块法向投影的最大间隙，不用射线结论——射线是从两层混合的块里
    //   采样的，取点位置不同结论就会变（实测 5 块只命中 1 块）。间隙是几何事实。
    // 判据二（补充）：投影直方图的两个密度峰。
    //   现场扫描的开口盒几乎都有 4mm 顶沿，沿口把内外壁连成**一张连续面**，投影后
    //   两层之间被沿口点填满 → 最大间隙只剩 1.88mm（实测），判据一必然失效，内外壁
    //   被当成一块、只出外壁一条轨迹。但"两个壁面各自很密"这件事没有变，密度峰还在，
    //   所以判据一失败时退回找双峰，两峰之间的谷就是分界。
    //   双峰路径额外要求"两片各自薄而平、且互相平行"，否则会把厚实块/斜杂面误拆。
    {
        std::vector<PlaneSeg> outerSegs, innerSegs;
        outerSegs.reserve(segs.size());
        innerSegs.reserve(segs.size());
        int nSplit = 0, nSplitPeak = 0, nSameDir = 0, nDupMerge = 0;

        // 直方图双峰检测：返回分界位置（沿法向，相对块质心）。
        auto twoPeaks = [](const std::vector<double>& h, double& mid)->bool{
            if(h.size() < 120) return false;
            const double lo = h.front(), hi = h.back();
            const double span = hi - lo;
            if(span < 2.0) return false;                  // 跨度连 2mm 都不到，不可能两层
            int nb = (int)std::floor(span/0.5) + 1;       // 0.5mm 一桶
            if(nb < 5) return false;
            if(nb > 240) nb = 240;
            const double w = span/(double)nb;
            std::vector<double> d(nb, 0.0);
            for(double v : h){
                int k = (int)((v-lo)/w);
                if(k < 0) k = 0; else if(k >= nb) k = nb-1;
                d[k] += 1.0;
            }
            // [1,2,1]/4 平滑：压掉采样网格造成的碎峰（不平滑时一块面能有十几个小峰）
            std::vector<double> m(nb, 0.0);
            for(int i=0;i<nb;++i){
                const double a = d[i];
                const double l = (i>0 ? d[i-1] : d[i]);
                const double r = (i+1<nb ? d[i+1] : d[i]);
                m[i] = 0.25*l + 0.5*a + 0.25*r;
            }
            int p1=-1, p2=-1; double v1=-1.0, v2=-1.0;
            for(int i=1;i<nb-1;++i){
                if(m[i] >= m[i-1] && m[i] > m[i+1]){
                    if(m[i] > v1){ p2=p1; v2=v1; p1=i; v1=m[i]; }
                    else if(m[i] > v2){ p2=i; v2=m[i]; }
                }
            }
            if(p1<0 || p2<0) return false;
            if(p1>p2){ std::swap(p1,p2); std::swap(v1,v2); }
            if((p2-p1)*w < 2.0) return false;             // 两峰太近 = 同一层的厚尾
            if(v2 < 0.35*v1)    return false;             // 次峰太弱 = 长尾/离群
            int g=p1; double gv=m[p1];
            for(int i=p1;i<=p2;++i) if(m[i] < gv){ gv=m[i]; g=i; }
            if(gv > 0.70*std::min(v1,v2)) return false;   // 谷不够深 = 一整团，不是两片
            mid = lo + ((double)g + 0.5)*w;
            return true;
        };

        for(size_t si2=0; si2<segs.size(); ++si2){
            PlaneSeg& s = segs[si2];
            if(s.curved || s.inl.size() < 60){ outerSegs.push_back(std::move(s)); continue; }
            const Vec3 n = s.pl.n, c0 = s.pl.c;
            std::vector<double> hh; hh.reserve(s.inl.size());
            for(int k : s.inl) hh.push_back(dot(Pf[k]-c0, n));
            std::sort(hh.begin(), hh.end());
            // 只在"两侧各至少 30 点"的候选里找最大间隙：边缘离群点的间隙往往比
            // 真正的层间间隙还大（实测最大间隙会落在排序后第 2 个点上）。
            int gapAt = -1; double maxGap = 0.0;
            for(size_t t=30; t+30<=hh.size(); ++t){
                const double g = hh[t]-hh[t-1];
                if(g > maxGap){ maxGap = g; gapAt = (int)t; }
            }
            // 单层平面的 h 散布只有噪声量级（~0.5mm），1.5mm 门限不会误切
            s.layerThick = hh[(size_t)((hh.size()-1)*0.98)] - hh[(size_t)((hh.size()-1)*0.02)];

            double hMid = 0.0; bool viaPeak = false;
            if(gapAt >= 0 && maxGap >= 1.5){
                hMid = 0.5*(hh[gapAt-1]+hh[gapAt]);
            } else if(twoPeaks(hh, hMid)){
                viaPeak = true;
            } else {
                s.layerCnt = 1; outerSegs.push_back(std::move(s)); continue;
            }

            std::vector<int> lo, hi;
            lo.reserve(s.inl.size()); hi.reserve(s.inl.size());
            for(int k : s.inl) (dot(Pf[k]-c0, n) < hMid ? lo : hi).push_back(k);
            if(lo.size() < 30 || hi.size() < 30){ s.layerCnt = 1; outerSegs.push_back(std::move(s)); continue; }
            // 少数侧还要"够一张面"：一张可喷面得有面积，1~2% 的边缘点簇不是第二张面。
            // p002 实测：5417 点主面 + 92 个边缘点（1.7%）被判成"两层"，低侧单独建面
            // 后轮廓提取直接失败（"点太少/全散点"），却仍占掉一个面号——界面上就是
            // 一块平板显示成 2~3 个面。门槛取 2%：真双层件最不均也在 15% 以上
            // （norim/rim 接近 1:1），而 p037 实测真第二层占 2.8%，不能被吞掉。
            {
                const size_t nMinor = std::min(lo.size(), hi.size());
                const size_t nAll   = lo.size() + hi.size();
                if((double)nMinor < 0.02*(double)nAll){
                    if(std::getenv("SPRAYTRAJ_FACE_DBG"))
                        fprintf(stderr,"[split] seg%zu 少侧仅 %zu/%zu 点(%.1f%%)，不足一张面，按单层处理\n",
                                si2, nMinor, nAll, 100.0*(double)nMinor/(double)nAll);
                    s.layerCnt = 1; outerSegs.push_back(std::move(s)); continue;
                }
            }

            // 少数侧必须"铺得开"：面内两个方向的包围盒都要达到多数侧的 25% 以上。
            // 单层大面的边缘小团（本次用户件 36 点，仅 5.3x5.0x2.9mm）排序后也会
            // 出现 3.98mm 的"最大间隙"；若把它当成另一层壁，29 万点的主面就会被
            // 误标成内壁、翻转法向，再被"朝下 >130°"整面跳过——界面上就只剩边缘
            // 两条线。真双层两面共享同一块投影区域（norim 两面 ext 比 >90%，
            // rim 的稀疏层也铺满整个 400x300 底面）。
            // 曾提到 60% 想把 p002 那条 18mm 边带（41%）拦住，结果把 p037 的真第二层
            // 也吞掉（2.8% → 单层 → 少一排 → 换边撞枪）。改回 25%：边带由上面那条
            // 占比<2% 的闸门处理，这里只管"小团冒充壁"。
            const Vec3 bu = norm(cross(n, std::fabs(n.z) < 0.9 ? Vec3{0,0,1} : Vec3{1,0,0}));
            const Vec3 bv = cross(n, bu);
            auto ext2 = [&](const std::vector<int>& v, double& eu, double& ev){
                double u0=1e300,u1=-1e300,v0=1e300,v1=-1e300;
                for(int k : v){
                    const double a=dot(Pf[k],bu), b=dot(Pf[k],bv);
                    u0=std::min(u0,a); u1=std::max(u1,a);
                    v0=std::min(v0,b); v1=std::max(v1,b);
                }
                eu = (v.empty()? 0.0 : u1-u0);
                ev = (v.empty()? 0.0 : v1-v0);
            };
            double lu,lv,hu,hv;
            ext2(lo,lu,lv); ext2(hi,hu,hv);
            const double coverU = std::min(lu,hu) / std::max(1e-9, std::max(lu,hu));
            const double coverV = std::min(lv,hv) / std::max(1e-9, std::max(lv,hv));
            if(std::getenv("SPRAYTRAJ_FACE_DBG"))
                fprintf(stderr,"[split] seg%zu 双层候选: 少侧=%zu 多侧=%zu 面内跨度 少=%.1fx%.1f 多=%.1fx%.1f 覆盖=%.0f%%x%.0f%%\n",
                        si2, std::min(lo.size(),hi.size()), std::max(lo.size(),hi.size()),
                        std::min(lu,hu), std::min(lv,hv), std::max(lu,hu), std::max(lv,hv),
                        100.0*coverU, 100.0*coverV);
            if(std::min(coverU,coverV) < 0.25){
                if(std::getenv("SPRAYTRAJ_FACE_DBG"))
                    fprintf(stderr,"[split] seg%zu 少侧覆盖不足 %.1f%%x%.1f%%，按单层处理\n",
                            si2, 100.0*coverU, 100.0*coverV);
                s.layerCnt = 1; outerSegs.push_back(std::move(s)); continue;
            }

            // 双峰路径的防误拆校验：两片必须各自薄而平，且互相平行。
            // 判据一拆出来的不加这道校验（原来能拆的一律保持原样，避免回归漂移）。
            if(viaPeak){
                Vec3 n1,c1,n2,c2;
                bool okFit = ls_plane(Pf, lo, n1, c1) && ls_plane(Pf, hi, n2, c2);
                double thk = 0.0, r1 = 1e9, r2 = 1e9, angLR = 0.0;
                if(okFit){
                    const double R2D = 180.0/3.141592653589793;
                    angLR = std::acos(std::max(-1.0,std::min(1.0,std::fabs(dot(n1,n2)))))*R2D;
                    Vec3 nu = norm(n1);
                    thk = std::fabs(dot(c2-c1, nu));
                    auto rmsA=[&](const std::vector<int>& v, const Vec3& nn, const Vec3& cc){
                        double s2=0.0;
                        for(int k : v){ const double dd = dot(Pf[k]-cc, nn); s2 += dd*dd; }
                        return (v.empty()? 1e9 : std::sqrt(s2/(double)v.size()));
                    };
                    r1 = rmsA(lo, n1, c1); r2 = rmsA(hi, n2, c2);
                }
                // 25°：超过就不是同一个薄壁的两面了（该由区域生长分开）；
                // 每片厚度上限 max(1.5, 0.35×间距)：壁厚 4mm 时允许 1.5mm 起伏，
                // 现场点云噪声 1mm 级也过得去，但"厚实块"会被挡下来。
                const double tolR = std::max(1.5, 0.35*thk);
                // 分离度：真两层是"两片薄点、中间空"——每片沿法向只有噪声那么厚。
                // 而一张连续的面被假双峰劈开时，两片各自本身就铺开好几毫米，
                // RMS 挡不住这种情况（劈出来的两半各自非常平，r1/r2 都很小）：
                // 棕槽实测 间距 5.8mm、两片各摊 2.9mm、maxGap 只有 0.08mm，
                // 一路过了夹角/间距/RMS 三道校验，被切成 17371+21640 两半。
                auto span98 = [&](const std::vector<int>& v){
                    if(v.size() < 10) return 1e9;
                    std::vector<double> h2; h2.reserve(v.size());
                    for(int k : v) h2.push_back(dot(Pf[k]-c0, n));
                    std::sort(h2.begin(), h2.end());
                    return h2[(size_t)((h2.size()-1)*0.98)] - h2[(size_t)((h2.size()-1)*0.02)];
                };
                const double t1 = span98(lo), t2 = span98(hi);
                double kSpan = 0.35;   // 调到 1.0 等于关掉这道校验（真双层被误挡时用）
                if(const char* e = std::getenv("SPRAYTRAJ_PEAK_SPAN")){ double v=atof(e); if(v>0) kSpan=v; }
                const double tolT = std::max(1.2, kSpan*thk);
                if(!okFit || angLR > 25.0 || thk < 1.5 || r1 > tolR || r2 > tolR
                   || std::max(t1,t2) > tolT){
                    if(std::getenv("SPRAYTRAJ_FACE_DBG"))
                        fprintf(stderr,"[split] 块%zu 双峰但否决: 夹角=%.1f° 间距=%.2f r1=%.2f r2=%.2f "
                                       "片厚=%.2f/%.2f (tol=%.2f/%.2f)\n",
                                si2, angLR, thk, r1, r2, t1, t2, tolR, tolT);
                    s.layerCnt = 1; outerSegs.push_back(std::move(s)); continue;
                }
                ++nSplitPeak;
            }

            auto cen=[&](const std::vector<int>& v){
                Vec3 t{0,0,0};
                for(int k : v) t = t + Pf[k];
                return t*(1.0/(double)v.size());
            };
            const Vec3 cLo = cen(lo), cHi = cen(hi);

            // 两片先把几何立起来；法向符号一律交给 orient_piece_after_split 重判，
            // 这里给的 n 只是"源块法向"这个起点。
            PlaneSeg so;                                     // 高侧（沿源块法向更远）
            so.curved = false;
            so.collapsed = s.collapsed;
            so.pl.n = n; so.pl.c = cHi; so.pl.d = -dot(n, cHi);
            so.inl  = std::move(hi);
            so.layerCnt = 2; so.layerThick = s.layerThick;

            PlaneSeg sn;                                     // 低侧
            sn.curved = false;
            sn.collapsed = s.collapsed;
            sn.pl.n = n; sn.pl.c = cLo; sn.pl.d = -dot(n, cLo);
            sn.inl  = std::move(lo);
            sn.layerCnt = 2; sn.layerThick = s.layerThick;
            // 拆层用该侧全体点求质心，保留两片各自的几何支撑点。
            so.fitSupport = so.inl;
            sn.fitSupport = sn.inl;

            // ---- 拆层也要过方向过滤：两片各自定可喷侧，不许"一片取反" ----
            orient_piece_after_split(Pf, cn.tree, so, partC);
            orient_piece_after_split(Pf, cn.tree, sn, partC);

            // ---- 两片到底是不是"壁的两侧" ----
            // 只有可喷侧确实相反才算内外壁；同向 = 台阶/两级平台/凸台，
            // 两片各按自己的法向布行（都朝 +Z 就都喷），谁也不许被翻过来。
            const bool facing = dot(so.pl.n, sn.pl.n) < -0.5;
            if(!facing){
                // ---- 同向回并：这两片其实是同一张脸被缝隙劈开的 ----
                // 现象（2026-11-13 现场）：一张完整平面沿法向铺开 5.8mm，点排序后
                // 出现一段 1.5mm 以上的空档，判据一就把这张脸劈成高低两片；两片共用
                // 源块法向、只差一个平移，界面上就是"同一个位置两个平面、法向完全
                // 一样"，各自布行 = 同一块地方喷两遍、膜厚翻倍。
                //
                // 判据一缺的就是"两片到底分没分开"这一问：它只看空档大小，不看两片
                // 自身的厚度。真双层是"两片薄点、中间空"，而一张起伏面排序后同样有
                // 空档，但两片各自本身就铺开好几毫米。
                //
                // 回并判据（三条全中，缺一条都会把真台阶吞掉）：
                //   ① 平行：两片法向夹角 <= 5°（劈出来的两片本来就共用源块 n，必过）；
                //   ② 贴得近：沿法向间距 <= dupGap（默认 8mm，可用环境变量覆盖）；
                //   ③ 投影重合：面内 2mm 网格的 IoU >= 0.85。
                // 真台阶的两级在面内一般错开，IoU 低，保持两片；上一级完全盖住下一级
                // 时 IoU 高，但那种下级本来就喷不到，并成一片也是对的。
                bool doDup = true;
                if(const char* e = std::getenv("SPRAYTRAJ_DUP_MERGE")) doDup = (e[0] != '0');
                double dupGap = 8.0, dupIoU = 0.85;
                if(const char* e = std::getenv("SPRAYTRAJ_DUP_GAP")){ double v=atof(e); if(v>0) dupGap=v; }
                if(const char* e = std::getenv("SPRAYTRAJ_DUP_IOU")){ double v=atof(e); if(v>0) dupIoU=v; }
                const double dupThk = std::fabs(dot(sn.pl.c - so.pl.c, n));
                if(doDup && dupThk <= dupGap && dot(so.pl.n, sn.pl.n) > 0.9962){
                    auto cellsOf = [&](const std::vector<int>& v){
                        std::vector<long long> c; c.reserve(v.size());
                        for(int k : v){
                            const long long a = (long long)std::floor(dot(Pf[k],bu)/2.0);
                            const long long b = (long long)std::floor(dot(Pf[k],bv)/2.0);
                            c.push_back(a*1000003LL + b);
                        }
                        std::sort(c.begin(), c.end());
                        c.erase(std::unique(c.begin(), c.end()), c.end());
                        return c;
                    };
                    const std::vector<long long> ca = cellsOf(so.inl), cb = cellsOf(sn.inl);
                    size_t ia=0, ib=0, inter=0;
                    while(ia<ca.size() && ib<cb.size()){
                        if(ca[ia]==cb[ib]){ ++inter; ++ia; ++ib; }
                        else if(ca[ia]<cb[ib]) ++ia; else ++ib;
                    }
                    const double iou = (double)inter /
                                       (double)std::max<size_t>(1, ca.size()+cb.size()-inter);
                    // 两条口径取"或"，对应两种误拆形态：
                    //   ① 重叠型：两片占同一块地方（IoU 高）——劈在同一个面上；
                    //   ② 相邻型：一张面被从中间劈成左右两半（IoU 很低），但两半
                    //      共面，小片到大片平面的高度 p95 只有 1mm 级。
                    //      棕槽实测就是这种：质心 y 差 38mm、法向完全相同、间距 5.8mm。
                    // 真台阶的 p95 ≈ 台阶高（5mm 级），过不了 ②。
                    PlaneSeg* big = (so.inl.size() >= sn.inl.size()) ? &so : &sn;
                    PlaneSeg* sml = (big == &so) ? &sn : &so;
                    std::vector<double> hh2; hh2.reserve(sml->inl.size());
                    for(int k : sml->inl) hh2.push_back(std::fabs(dot(Pf[k]-big->pl.c, big->pl.n)));
                    double h95 = 1e9;
                    if(!hh2.empty()){
                        const size_t q = (size_t)(0.95*(double)(hh2.size()-1));
                        std::nth_element(hh2.begin(), hh2.begin()+q, hh2.end());
                        h95 = hh2[q];
                    }
                    const double dupH = 2.0;
                    if(iou >= dupIoU || h95 <= dupH){
                        // 保留点数多那片的归属，但平面要按合并后的点重拟合：
                        // 直接用某一侧的平面会让另一侧的点整体偏 dupThk，喷距跟着错。
                        big->inl.insert(big->inl.end(), sml->inl.begin(), sml->inl.end());
                        std::sort(big->inl.begin(), big->inl.end());
                        big->inl.erase(std::unique(big->inl.begin(), big->inl.end()), big->inl.end());
                        Vec3 nf2, cf2;
                        if(ls_plane(Pf, big->inl, nf2, cf2)){
                            if(dot(nf2, n) < 0) nf2 = nf2*(-1.0);
                            big->pl.n = nf2; big->pl.c = cf2; big->pl.d = -dot(nf2, cf2);
                            big->fitSupport.clear();
                        }
                        big->layerCnt = 0;                       // 合回一张面，不再是双层
                        orient_piece_after_split(Pf, cn.tree, *big, partC);
                        big->fo.isInner = 0;
                        big->fo.src = 7;                         // 7 = 同向回并
                        ++nDupMerge;
                        if(std::getenv("SPRAYTRAJ_FACE_DBG"))
                            fprintf(stderr,"[split] 同向回并: 间距=%.2fmm 重合度=%.2f 共面p95=%.2f → %zu点 "
                                           "(阈值 gap<=%.1f IoU>=%.2f h<=%.1f)\n",
                                    dupThk, iou, h95, big->inl.size(), dupGap, dupIoU, dupH);
                        outerSegs.push_back(std::move(*big));
                        continue;
                    }
                }
                // src=6 与"薄壁拆层(4)"区分开：清单里要能一眼看出这两片不是壁的两侧，
                // 而是被间隙误分成两片的一张脸，两片都该喷。
                so.fo.isInner = 0; sn.fo.isInner = 0;
                so.fo.src = 6;     sn.fo.src = 6;
            }

            if(std::getenv("SPRAYTRAJ_FACE_DBG"))
                fprintf(stderr, "[split] seg%zu 双层: 间隙=%.2fmm 高侧点=%d 低侧点=%d "
                                "高侧n=(%.2f,%.2f,%.2f) 低侧n=(%.2f,%.2f,%.2f) → %s\n",
                        si2, maxGap, (int)so.inl.size(), (int)sn.inl.size(),
                        so.pl.n.x,so.pl.n.y,so.pl.n.z, sn.pl.n.x,sn.pl.n.y,sn.pl.n.z,
                        facing ? "壁两侧(内外壁)" : "同向两面(不是壁，各自布行)");
            ++nSplit;
            if(!facing) ++nSameDir;
            // 顺序：先喷完所有外壁，再进腔喷内壁。
            // 若按"外/内/外/内"交替，枪会在工件内外来回穿，每次都要走一遍抬枪绕行。
            if(!facing){
                outerSegs.push_back(std::move(so));
                outerSegs.push_back(std::move(sn));
            } else if(so.fo.isInner > 0){
                innerSegs.push_back(std::move(so));
                outerSegs.push_back(std::move(sn));
            } else {
                outerSegs.push_back(std::move(so));
                innerSegs.push_back(std::move(sn));
            }
        }
        // 顺序：先喷完所有外壁，再进腔喷内壁。
        // 若按"外/内/外/内"交替，枪会在工件内外来回穿，每次都要走一遍抬枪绕行。
        segs.clear();
        segs.reserve(outerSegs.size()+innerSegs.size());
        for(auto& s : outerSegs) segs.push_back(std::move(s));
        for(auto& s : innerSegs) segs.push_back(std::move(s));
        if(nSplit || nDupMerge) fprintf(stderr, "[split] 薄壁拆段: %d 块 → %zu 段（内壁 %d 面，其中双峰判据 %d 块；"
                                   "同向两面 %d 块按普通面各自布行；同向回并 %d 块合回单面）\n",
                           nSplit, segs.size(), nSplit-nSameDir, nSplitPeak, nSameDir, nDupMerge);
    }

    // ---- 最终合并闸门：把"吸收口径"贯穿到最后 ----
    //
    // 现场问题（2026-11-19，面板 p001-173313-cam 38922 点，flat=40）：第一步"吸收"
    // 已经把 38211/38211 点全部归进主面，最后却仍然走出 5 个面。三条本该拦住它的
    // 防线全部落空，而且落空的原因各不相同：
    //   ① 吸收：只决定"点归谁"（轮廓/覆盖率），完全不参与"要不要拆"，所以后面的
    //      区域生长、薄壁拆段根本不受它约束——吸收率 100% 也照样拆。
    //   ② 合并：三道合并各用各的阈值（接缝 p95<=2.5 / 夹角<=8° / IoU>=0.85 或 h95<=2），
    //      而"薄壁拆段"排在它们**之后**，新劈出来的两片一条都过不去，等于把并好的
    //      面又拆开、且没人再兜底。
    //   ③ 喷幅：代码里没有喷幅参数（行距才是唯一覆盖依据），这道防线根本不存在。
    //
    // 根子是同一个：**"多远的点还算这张面"这个用户口径，全流程只用了一次**。
    // 这里用同一个 absorbTol 做最后一次裁决，不引入新魔数：
    //   ① 与主面夹角 <= 30°（与"同向同高度并入主面"共用同一条线），超过就是真折边/侧壁；
    //   ② 该面点到主面平面的高度 p95 <= absorbTol —— 用户自己说"这个距离内算一张面"；
    //   ③ 与主面在 KNN 图里直接相连 —— 是同一块料，不是工件远端恰好同向的另一块。
    // 三条全中 → 点并入主面，平面按并集重拟合（不重拟合的话喷距会整体偏一个台阶高）。
    // 回退开关：SPRAYTRAJ_MERGE_GATE=0。
    if(segs.size() > 1 && segs.size() <= 32){
        bool gateOn = true;
        if(const char* e = std::getenv("SPRAYTRAJ_MERGE_GATE")) gateOn = (e[0] != '0');
        if(gateOn){
            const double D2R = 3.141592653589793/180.0;
            const double absorbTol =
                std::min(params->flatnessTol, std::max(2.0, 0.6*params->standoff));
            size_t mi = 0;
            for(size_t i=1;i<segs.size();++i)
                if(segs[i].inl.size() > segs[mi].inl.size()) mi = i;
            // 面的归属标签（下标空间 = Pf）
            std::vector<int> lab(Pf.size(), -1);
            for(size_t i=0;i<segs.size();++i)
                for(int k : segs[i].inl)
                    if(k >= 0 && (size_t)k < Pf.size()) lab[k] = (int)i;
            // 主面的 KNN 邻居集合：有边连到主面的点即"与主面相接"
            const std::vector<int>& knn = cn.knn;
            const int KN = CloudNormals::K;
            std::vector<char> nearMain(Pf.size(), 0);
            // 反向账：主面点的 KNN 邻居落在哪个面上。只查"主面→别人"会漏掉
            // 被薄壁拆段劈出去的那一层——两层相距 2mm 时，各自的 10 近邻全是
            // 同层点，主面那侧永远看不到对方（p001 实测 573 点被这条件误杀）。
            std::vector<char> adjToMain(segs.size(), 0);
            if(knn.size() >= Pf.size()*(size_t)KN){
                for(size_t i=0;i<Pf.size();++i){
                    if(lab[i] != (int)mi) continue;
                    const int* row = &knn[i*(size_t)KN];
                    for(int t=0;t<KN;++t){
                        const int j = row[t];
                        if(j < 0 || (size_t)j >= Pf.size()) continue;
                        nearMain[j] = 1;
                        const int lj = lab[j];
                        if(lj >= 0 && (size_t)lj < segs.size()) adjToMain[lj] = 1;
                    }
                }
            }
            auto touchesMain = [&](size_t si){
                if(adjToMain[si]) return true;
                for(int k : segs[si].inl)
                    if(k >= 0 && (size_t)k < Pf.size() && nearMain[k]) return true;
                return false;
            };
            const Vec3 mn = segs[mi].pl.n;
            const double md = segs[mi].pl.d;
            // 主面的 2D 覆盖范围（沿主面平面的包围盒）。判据的物理含义是用户那句话：
            // "主平面的喷幅完全足够覆盖碎面" —— 主面的行本来就从这块碎面上方扫过去，
            // 只要碎面落在主面的投影范围内、高度又在吸收口径内，它就不是另一张面。
            // 光靠 KNN 图会漏判：薄壁拆段劈出的那一层若点稀疏（573 点铺 291x122mm，
            // 近邻间距 ~8mm），K=25 的近邻全是同层点，两层之间一条边都没有。
            const Vec3 fbu = norm(cross(mn, std::fabs(mn.z) < 0.9 ? Vec3{0,0,1} : Vec3{1,0,0}));
            const Vec3 fbv = cross(mn, fbu);
            double u0 = 1e30, u1 = -1e30, v0 = 1e30, v1 = -1e30;
            for(int k : segs[mi].inl){
                const double uu = dot(Pf[k], fbu), vv = dot(Pf[k], fbv);
                u0 = std::min(u0, uu); u1 = std::max(u1, uu);
                v0 = std::min(v0, vv); v1 = std::max(v1, vv);
            }
            const double fpMargin = 10.0;
            auto insideFootprint = [&](size_t si)->bool{
                if(segs[si].inl.empty()) return false;
                size_t cnt = 0;
                for(int k : segs[si].inl){
                    const double uu = dot(Pf[k], fbu), vv = dot(Pf[k], fbv);
                    if(uu >= u0-fpMargin && uu <= u1+fpMargin &&
                       vv >= v0-fpMargin && vv <= v1+fpMargin) ++cnt;
                }
                return (double)cnt >= 0.90*(double)segs[si].inl.size();
            };
            std::vector<int> addPts;
            std::vector<char> drop(segs.size(), 0);
            for(size_t i=0;i<segs.size();++i){
                if(i == mi) continue;
                const PlaneSeg& s = segs[i];
                const bool dbg = (std::getenv("SPRAYTRAJ_FACE_DBG") != nullptr);
                if(s.curved){
                    if(dbg) fprintf(stderr, "[merge] 面%zu(%zu点) 跳过：曲面块（法向逐点变）\n", i, s.inl.size());
                    continue;
                }
                if(s.inl.empty()) continue;
                // 先算"整块贴不贴在主面上"。这是几何尺度，不是工艺阈值：
                // 3*fitTol(3mm)=9mm 量级，和噪声/拟合同源。
                // 为什么要单独一条分支：p002 的那两条 73°/55° 窄带，法向是边缘效应
                // 拉偏的，但整块点距主面 p95 只有 6.0mm —— 按用户口径"投影下去"就该
                // 回主面。而真侧面（p030 的 -Z 面、p037 的端面）整块跨几十毫米，
                // 过不了这条，仍然独立成面。
                std::vector<double> hh; hh.reserve(s.inl.size());
                for(int k : s.inl) hh.push_back(std::fabs(dot(Pf[k], mn) + md));
                const size_t q = (size_t)(0.95*(double)(hh.size()-1));
                std::nth_element(hh.begin(), hh.begin()+q, hh.end());
                const double h95 = hh[q];
                const double flatCap = 8.0;
                const bool stuck = (h95 <= flatCap) && (touchesMain(i) || insideFootprint(i));
                if(!stuck){
                    // ---- 吸收口径（2026-11-23 用户两次定稿）----
                    // 用户原话：「不管真不真侧面，我都要吸收」→「在吸收范围内，就没有门槛」。
                    // 意思是：**夹角不参与裁决**。只要这块点落在吸收范围内（h95 <= absorbTol
                    // 且与主面相连/落在主面投影内），就归主面，不再问它是折边还是侧面。
                    // 旧口径 60° 门把 p016 那条 65.7° 的 20×203mm 折边拒在门外，结果
                    // 折边自己占了一个面序号、主面 49171 点反而整面不出轨迹。
                    // 夹角现在只进日志（>60° 特意打一行），不参与判断。
                    // 回退：SPRAYTRAJ_MERGE_SIDE=<角度> 恢复旧行为（默认不设 = 无门槛）。
                    const double cang = dot(s.pl.n, mn);
                    const double ang  = std::acos(std::max(-1.0,std::min(1.0,cang)))/D2R;
                    double sideCap = 0.0;          // 0 = 不设门槛
                    if(const char* e = std::getenv("SPRAYTRAJ_MERGE_SIDE")){
                        double v = atof(e);
                        if(v > 1.0 && v <= 179.0) sideCap = v;
                    }
                    if(sideCap > 0.0 && ang > sideCap){
                        if(dbg) fprintf(stderr, "[merge] 面%zu(%zu点) 跳过：夹角%.1f° > 手动门槛%.0f° h95=%.2f\n",
                                        i, s.inl.size(), ang, sideCap, h95);
                        continue;
                    }
                    if(ang > 90.0)
                        // 反向面被吸收必须留痕：这类面并进去后不会单独布行，
                        // 万一是真另一张面，这里就是唯一线索，不能静默。
                        fprintf(stderr, "[merge] 注意 面%zu(%zu点) 反向面被吸收（夹角%.1f°）"
                                        "→ 不再单独布行，交主面喷锥覆盖\n", i, s.inl.size(), ang);
                    else if(ang > 60.0)
                        fprintf(stderr, "[merge] 面%zu(%zu点) 侧面折边吸收：夹角%.1f° h95=%.2f → 不单独布行\n",
                                i, s.inl.size(), ang, h95);
                    if(!touchesMain(i) && !insideFootprint(i)){
                        if(dbg) fprintf(stderr, "[merge] 面%zu(%zu点) 跳过：既不与主面相接，也不在主面投影范围内\n", i, s.inl.size());
                        continue;
                    }
                    if(h95 > absorbTol){
                        if(dbg) fprintf(stderr, "[merge] 面%zu(%zu点) 跳过：高度p95=%.2fmm > 吸收 %.1fmm\n", i, s.inl.size(), h95, absorbTol);
                        continue;
                    }
                    const double ang2 = std::acos(std::min(1.0, cang))/D2R;
                    fprintf(stderr, "[merge] 面%zu(%zu点) → 主面%zu：夹角%.1f° 高度p95=%.2fmm(≤%.1f) %s\n",
                            i, s.inl.size(), mi, ang2, h95, absorbTol,
                            touchesMain(i) ? "KNN相连" : "落在主面投影内");
                } else {
                    fprintf(stderr, "[merge] 面%zu(%zu点) → 主面%zu：整块贴合 高度p95=%.2fmm(≤%.1f) 夹角%.1f° %s\n",
                            i, s.inl.size(), mi, h95, flatCap,
                            std::acos(std::max(-1.0,std::min(1.0,dot(s.pl.n,mn))))/D2R,
                            touchesMain(i) ? "KNN相连" : "落在主面投影内");
                }
                addPts.insert(addPts.end(), s.inl.begin(), s.inl.end());
                drop[i] = 1;
            }
            if(!addPts.empty()){
                // 先存一份合并前的本面点给枪位验收用（见 PlaneSeg::core 注释）。
                // 必须在 insert 之前拷，之后 inl 就混进侧面点了。
                if(segs[mi].core.empty())
                    segs[mi].core = segs[mi].inl;
                std::vector<PlaneSeg> keep;
                keep.reserve(segs.size());
                for(size_t i=0;i<segs.size();++i){
                    if(drop[i] || i == mi) continue;
                    keep.push_back(std::move(segs[i]));
                }
                segs[mi].inl.insert(segs[mi].inl.end(), addPts.begin(), addPts.end());
                std::sort(segs[mi].inl.begin(), segs[mi].inl.end());
                segs[mi].inl.erase(std::unique(segs[mi].inl.begin(), segs[mi].inl.end()),
                                   segs[mi].inl.end());
                // 主面平面由吸收一次定死（fixed）：这里只并点，不重拟合，
                // 否则"吸收完的平面"又被下游改回去，等于吸收白做。
                if(segs[mi].curved && isCurved){
                    // 重建网格与布行都读 SurfPatch，吸收后必须与面点归属同步。
                    // 曲面参考基仍沿原曲面的朝向，避免被并入的侧壁扭转。
                    surf.P.clear();
                    surf.N.clear();
                    surf.P.reserve(segs[mi].inl.size());
                    surf.N.reserve(segs[mi].inl.size());
                    for(int k : segs[mi].inl){
                        surf.P.push_back(Pf[k]);
                        const Vec3 n = cn.nor[k];
                        surf.N.push_back(dot(n, surf.n0) < 0.0 ? n*(-1.0) : n);
                    }
                    surf.tree.init(surf.P);
                    segs[mi].pl = {surf.n0, -dot(surf.n0,surf.c0), surf.c0};
                }else if(!segs[mi].fixed){
                    Vec3 nf, cf;
                    if(ls_plane(Pf, segs[mi].inl, nf, cf)){
                        PlaneSeg& main = segs[mi];
                        if(dot(nf, main.pl.n) < 0) nf = nf*(-1.0);
                        const Plane refit{nf, -dot(nf,cf), cf};
                        const double rmsCap = std::max(3.0, 0.12*params->standoff);
                        auto supportRms = [&](const PlaneSeg& s, const Plane& pl){
                            if(s.fitSupport.size() < 20) return 1e30;
                            double rss = 0.0;
                            for(int k : s.fitSupport){
                                const double e = dot(Pf[k]-pl.c, pl.n);
                                rss += e*e;
                            }
                            return std::sqrt(rss/(double)s.fitSupport.size());
                        };
                        const double refitRms = enableCandidateCheck ? supportRms(main, refit) : 0.0;
                        if(refitRms > rmsCap && supportRms(main, main.pl) <= rmsCap){
                            // 已验收的吸收点继续参与轮廓。重拟合拉偏载体时，先尝试原载体的外侧，
                            // 方位不合格时再从同方向的吸收面中寻找可喷侧载体。
                            const Vec3 originalNormal = main.pl.n;
                            size_t source = segs.size();
                            Plane carrier;
                            double window = SPRAYTRAJ_NORMAL_WIN_DEG;
                            if(const char* e = std::getenv("SPRAYTRAJ_NORMAL_WIN")){
                                const double value = atof(e);
                                if(value > 1.0 && value <= 180.0) window = value;
                            }
                            double bestAlignment = -1.0;
                            for(size_t order=0; order<segs.size(); ++order){
                                const size_t i = (mi+order)%segs.size();
                                if(i != mi && !drop[i]) continue;
                                const PlaneSeg& s = segs[i];
                                const double alignment = std::fabs(dot(originalNormal,s.pl.n));
                                if(alignment < std::cos(10.0/D2R) ||
                                   alignment <= bestAlignment || supportRms(s,s.pl) > rmsCap) continue;
                                Plane trial = s.pl;
                                if(dot(trial.n, trial.c-partC) < 0.0){
                                    trial.n = trial.n*(-1.0);
                                    trial.d = -trial.d;
                                }
                                const double angZ = std::acos(std::clamp(trial.n.z,-1.0,1.0))*D2R;
                                if(angZ > window) continue;
                                source = i;
                                carrier = trial;
                                bestAlignment = alignment;
                                if(i == mi) break;
                            }
                            if(source < segs.size()){
                                const double carrierRms = supportRms(segs[source], carrier);
                                if(source != mi){
                                    main.core = std::move(segs[source].inl);
                                    main.fitSupport = std::move(segs[source].fitSupport);
                                }
                                main.pl = carrier;
                                main.fo.n = carrier.n;
                                main.fo.angZ = std::acos(std::clamp(carrier.n.z,-1.0,1.0))*D2R;
                                main.fo.tilt = 90.0-main.fo.angZ;
                                main.fo.side = main.fo.angZ < 87.0 ? 1 : (main.fo.angZ > 93.0 ? -1 : 0);
                                main.fo.isInner = 0;
                                main.preservedCarrier = true;
                                fprintf(stderr, "[merge] 重拟合载体 RMS %.2fmm > %.2fmm，沿外侧载体面%zu（RMS %.2fmm）布行并复核枪位\n",
                                        refitRms, rmsCap, source, carrierRms);
                            }
                        }
                        if(!main.preservedCarrier) main.pl = refit;
                    }
                }
                segs[mi].layerCnt = 0;                 // 合回一张面，不再是双层
                segs[mi].fo.src = 8;                    // 8 = 最终合并闸门
                keep.push_back(std::move(segs[mi]));
                segs.swap(keep);
                int nMerged = 0;
                for(char c : drop) if(c) ++nMerged;
                fprintf(stderr, "[merge] 最终合并: %d 面并入主面（%zu 点），剩 %zu 面（吸收口径 %.1fmm）\n",
                        nMerged, addPts.size(), segs.size(), absorbTol);
            }
        }
    }

    // ---- 薄壳双面：一层点代表一张板，板有两面 ----
    // 上面所有内外壁判据都要求"点云里有两片点"（间隙或双峰）。但薄壁件（1~4mm）的
    // 激光通常只打到可及的那一面，只有一层点 —— 这时"内外"不是两片点，而是**同一张
    // 板的两个朝向**，这是几何事实，不需要第二片点来证明。
    // 判据用射线自由距离：两侧都够伸枪 → 两面都是空腔、都能喷，按 shellThick 沿法向
    // 偏移复制一个反向面；只有一侧够 → 实心件外表面，另一面埋在材料里，不复制。
    // 默认 shellThick=0 关闭，行为与之前逐点一致。
    if(params->shellThick > 0.0){
        const double gcD = (params->gunClear > 0.0) ? params->gunClear
                                                    : std::max(12.0, 0.2*params->standoff);
        const double needClear = std::max(25.0, 2.0*gcD);   // 能伸进枪的最小空腔
        const double D2R = 3.141592653589793/180.0;
        std::vector<PlaneSeg> back;
        back.reserve(segs.size());
        int nDual = 0;
        for(size_t i=0;i<segs.size();++i){
            const PlaneSeg& s = segs[i];
            if(s.curved) continue;                 // 曲面块不推断（法向逐点变，偏移没有定义）
            if(s.layerCnt == 2) continue;          // 双层已经拆出内外两片，不用再补
            if(segs.size() + back.size() >= (size_t)SPRAYTRAJ_PLANE_SLOTS) break;
            // -1 = 量程内无命中，视为无限自由
            const bool pFree = (s.fo.freePlus  < 0.0) || (s.fo.freePlus  >= needClear);
            const bool mFree = (s.fo.freeMinus < 0.0) || (s.fo.freeMinus >= needClear);
            if(!(pFree && mFree)) continue;
            PlaneSeg b = s;                        // inl 是点索引，与原面共享即可
            // 背面 = 真实可喷侧（fo.n）的**反方向**。
            // 不能用 pl.n：pl.n 是"枪从哪边来"，对侧壁它朝盒外、对顶/底面它朝腔内
            // （|nz|>0.05 强制归 +Z 半球），方向语义不统一，拿它偏移会把底面的背面
            // 算到腔里去（实测内底 z=0 的背面被算成 z=+4）。
            const Vec3 nb = s.fo.n * (-1.0);
            b.pl.n = nb;
            b.pl.c = s.pl.c + nb*params->shellThick;
            b.pl.d = -dot(b.pl.n, b.pl.c);
            b.fitSupport.clear();
            b.preservedCarrier = false;
            b.fo   = s.fo;
            b.fo.n = b.pl.n;
            b.fo.angZ = std::acos(std::max(-1.0,std::min(1.0,b.pl.n.z)))/D2R;
            b.fo.tilt = 90.0 - b.fo.angZ;
            b.fo.side = (b.fo.angZ < 87.0) ? 1 : (b.fo.angZ > 93.0 ? -1 : 0);
            b.fo.isInner = (s.fo.isInner > 0) ? 0 : 1;     // 内壁的背面就是外壁
            b.fo.src = 5;                                   // 5 = 薄壳双面（推断的另一面）
            std::swap(b.fo.freePlus, b.fo.freeMinus);
            b.layerCnt = 1;
            b.layerThick = params->shellThick;
            back.push_back(std::move(b));
            ++nDual;
        }
        for(auto& b : back) segs.push_back(std::move(b));
        if(nDual) fprintf(stderr, "[dual] 薄壳双面: %d 面按壁厚 %.1fmm 补出另一侧（总 %zu 段）\n",
                          nDual, params->shellThick, segs.size());
    }

    // ===== 碎面终检（2026-11-20）：拆段/补背面之后**新造出来**的面必须再过一遍 =====
    // 上面那轮过滤只筛"区域生长的原始面"。5166 的薄壁拆段、上方的薄壳双面都在它之后
    // 新造面（一个宽面被拆成两条窄带、或按壁厚补出背面），这些面从未被 D1/D5/D6 判过：
    // 判据看的是拆之前那整块（宽度够、占比够），拆完的小片自然全部过关。
    // 现场表现就是"短边只有几毫米、或贴着大面的小平面，照样进清单并布出一行"。
    // 这里用完全相同的口径补一次，只对**衍生面**生效（主面永远保留）。
    if(segs.size() > 1){
        size_t maxIdx2 = 0;
        for(size_t i=1;i<segs.size();++i)
            if(segs[i].inl.size() > segs[maxIdx2].inl.size()) maxIdx2 = i;
        const double total2 = (double)Pf.size();
        const Plane& mainPl2 = segs[maxIdx2].pl;
        double fragShort2 = kReconGap;
        if(const char* e = std::getenv("SPRAYTRAJ_FRAG_SHORT")){
            const double v = std::atof(e);
            if(v > 0.0) fragShort2 = v;
        }
        std::vector<double> wN2(segs.size(), -1.0), wL2(segs.size(), -1.0);
        for(size_t i=0;i<segs.size();++i){
            double a=0,b=0;
            if(inPlaneExtents(Pf, segs[i].inl, segs[i].pl.n, a, b)){ wN2[i]=std::min(a,b); wL2[i]=std::max(a,b); }
        }
        std::vector<char> drop2(segs.size(), 0);
        for(size_t i=0;i<segs.size();++i){
            if(i == maxIdx2) continue;
            const double ratio = (double)segs[i].inl.size()/total2;
            const double angMain = std::acos(std::min(1.0, std::fabs(dot(segs[i].pl.n, mainPl2.n))))*180.0/3.141592653589793;
            const bool okE = wN2[i] > 0.0;
            const bool d1 = okE && (wN2[i] < fragShort2);
            const bool d5 = ratio < 0.03;
            const bool realStruct = (angMain > 40.0) && okE && (wN2[i] >= fragShort2);
            bool d6 = false; size_t host=segs.size(); double hAng=0.0, dMed=1e18;
            for(size_t j=0;j<segs.size();++j){
                if(j==i || segs[j].inl.size() <= segs[i].inl.size()) continue;
                const double a = std::acos(std::min(1.0, std::fabs(dot(segs[j].pl.n, segs[i].pl.n))))*180.0/3.141592653589793;
                if(a < 25.0){ host=j; hAng=a; break; }
            }
            if(host < segs.size()){
                std::vector<double> dd; dd.reserve(segs[i].inl.size());
                for(int k : segs[i].inl) dd.push_back(std::fabs(dot(Pf[k], segs[host].pl.n) + segs[host].pl.d));
                std::nth_element(dd.begin(), dd.begin()+dd.size()/2, dd.end());
                dMed = dd[dd.size()/2];
                d6 = (dMed <= 5.0);
            }
            if(d1 || d6 || (d5 && !realStruct)){
                drop2[i] = 1;
                fprintf(stderr, "[split] 终检丢弃衍生面: 点数=%zu (占比%.1f%%, 与主面夹角%.1f°, "
                                "面内%.1f×%.1fmm) 命中[%s%s%s]\n",
                        segs[i].inl.size(), ratio*100.0, angMain, wN2[i], wL2[i],
                        d1?"D1短边 ":"", d5?"D5占比<3% ":"", d6?"D6贴同向大面 ":"");
                if(d6) fprintf(stderr, "[split]   终检 D6: 宿主=面%zu 夹角=%.1f° 中位距离=%.2fmm\n", host, hAng, dMed);
            }
        }
        size_t nD = 0; for(char c : drop2) if(c) ++nD;
        if(nD){
            std::vector<PlaneSeg> keep2; keep2.reserve(segs.size()-nD);
            for(size_t i=0;i<segs.size();++i) if(!drop2[i]) keep2.push_back(std::move(segs[i]));
            if(!keep2.empty()){ segs = std::move(keep2); fprintf(stderr, "[split] 终检: 衍生面 %zu → 剩 %zu 面\n", nD, segs.size()); }
            else fprintf(stderr, "[split] 终检: 衍生面 %zu 面全部命中，但会清空面集 → 放弃本次终检\n", nD);
        }
    }

    // 调试：每个面的内点单独写一份 ply（SPRAYTRAJ_FACE_PLY=<目录>）。
    // 用途：现场说"这几个小片为什么没过滤"时，先确认它到底属于哪个面 ——
    // 面清单里的一个"面"可能在重建后被切成多个环（岛），肉眼看像好几个小面。
    if(const char* fdir = std::getenv("SPRAYTRAJ_FACE_PLY")){
        for(size_t i=0;i<segs.size();++i){
            char fn[600];
            snprintf(fn,sizeof(fn),"%s/face%02zu.ply", fdir, i);
            FILE* fp = fopen(fn,"w");
            if(!fp) continue;
            fprintf(fp,"ply\nformat ascii 1.0\nelement vertex %zu\nproperty float x\nproperty float y\nproperty float z\nend_header\n",
                    segs[i].inl.size());
            for(int k : segs[i].inl) fprintf(fp,"%.3f %.3f %.3f\n", Pf[k].x, Pf[k].y, Pf[k].z);
            fclose(fp);
        }
        fprintf(stderr, "[dbg] 面级 ply 已写出: %zu 面 -> %s\n", segs.size(), fdir);
    }

    // 输出所有分割平面（供 WPF 显示重建面）
    // 为每个面计算局部基(u,v)主轴 + 内点包围盒，WPF 直接画矩形（不做 PCA 近似）
    res.plane_count = 0;
    size_t wi = 0;                     // 只统计"网格重建成功"的面，避免空面被前端画成包围盒
    std::vector<char> meshOk(segs.size(), 0);   // 没重建出面片的面不生成轨迹
    // 窄条补全记录：补全面要让轨迹也铺满整个矩形，不能只跟着点云走
    std::vector<char> planeFilled(segs.size(), 0);
    std::vector<std::array<Vec3,4>> planeFillRect(segs.size());
    // 每个面重建出来的轮廓多边形（面局部 u,v 坐标，与 planeBase 配套）。
    // 布行时用它把行端点夹在面内 —— 轮廓已经算出来了，布行却按 P0.5/P99.5 包围盒走，
    // 等于白算：包围盒会把"面在这一行没有材料"的空档也算进去，实测出框 338mm。
    std::vector<std::vector<std::vector<Vec2>>> planeLoops(segs.size());
    // 面局部基（uf,vf）与面心，用于把 (u,v) 与世界坐标互转（与 res.planes 上的 ux/uy 同源）
    std::vector<Vec3> planeBaseC(segs.size()), planeBaseU(segs.size()), planeBaseV(segs.size());
    // 重建网格是否碎片化（碎成多块且最大块 <= 50%）：这类面必须按轮廓切段，
    // 否则行会横穿两块材料之间的空档。在重建阶段算一次，轨迹阶段直接用。
    std::vector<char> planeFrag(segs.size(), 0);
    std::vector<char> sparseDisconnected(segs.size(), 0);
    std::vector<Vec3> surfaceVerts;
    std::vector<int> surfaceTris;
    std::vector<unsigned char> surfaceFilledTris;
    std::vector<int> planeSeg;          // res.planes[wi] 对应的 segs 下标（拆薄壁行时要用）
    for(size_t i=0;i<segs.size() && wi<SPRAYTRAJ_PLANE_SLOTS;++i){
        Plane& P = segs[i].pl;
        res.planes[wi].nx = P.n.x;
        res.planes[wi].ny = P.n.y;
        res.planes[wi].nz = P.n.z;
        res.planes[wi].cx = P.c.x;
        res.planes[wi].cy = P.c.y;
        res.planes[wi].cz = P.c.z;
        res.planes[wi].count = (int)segs[i].inl.size();
        {
            // 面级法向诊断（angZ 同时决定是否喷：朝下超窗口的面不生成轨迹）
            const FaceOrient& fo = segs[i].fo;
            res.planes[wi].n0x = fo.n.x; res.planes[wi].n0y = fo.n.y; res.planes[wi].n0z = fo.n.z;
            res.planes[wi].angZ      = fo.angZ;
            res.planes[wi].tiltDeg   = fo.tilt;
            res.planes[wi].freePlus  = fo.freePlus;
            res.planes[wi].freeMinus = fo.freeMinus;
            res.planes[wi].side      = fo.side;
            res.planes[wi].isInner   = fo.isInner;
            res.planes[wi].sideSrc   = fo.src;
            res.planes[wi].layerCnt  = segs[i].layerCnt;
            res.planes[wi].layerThick= segs[i].layerThick;
        }

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
                okMesh = build_surface_mesh(surf, params->fillHoleSize,
                                            verts3, tris, surfaceFilledTris,
                                            SPRAYTRAJ_PLANE_MAX_POLY_VERTS,
                                            SPRAYTRAJ_PLANE_MAX_POLY_TRIS);
            }
            res.planes[wi].polyVertCount = 0;
            res.planes[wi].polyTriCount  = 0;
            // 最小可喷宽度：曲面块也按投影到面内 (uf,vf) 的最小外接矩形短边判。
            // 线状/退化点云即使被拟合成"曲面"，在面内也只是一条细带，必须拒掉。
            double faceW = -1.0, faceL = -1.0;
            if(okMesh){
                std::vector<Vec2> wv; wv.reserve(verts3.size());
                for(const auto& p3 : verts3){
                    Vec3 d = p3 - P.c;
                    wv.push_back(Vec2{ dot(d, uf), dot(d, vf) });
                }
                faceW = minRectShortSide(wv, &faceL);
            }
            if(okMesh && params->minFaceWidth > 1e-9 &&
               (faceW < 0.0 || faceW < params->minFaceWidth)){
                fprintf(stderr,
                        "[重建] 面片%zu: 曲面宽度 %.1fmm < 最小可喷宽度 %.1fmm，判为线状/退化面，不输出该面\n",
                        i, faceW, params->minFaceWidth);
                continue;
            }
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
                fprintf(stderr, "[重建] 面片%zu: 顶点=%d 三角=%d（曲面块）宽度=%.1fmm 内点=%zu\n",
                        wi, nv, nt, faceW, segs[i].inl.size());
                meshOk[i] = 1;
                surfaceVerts = std::move(verts3);
                surfaceTris = std::move(tris);
                planeSeg.push_back((int)i);
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
            std::vector<std::vector<Vec2>> meshLoops;   // 重建轮廓（局部 u,v），供布行裁剪
            bool okMesh = false;
            {
                ProfStage ps("3.轮廓重建", &tMesh);
                // 环形窄面：重建尺度必须能分辨环宽。若用工艺行距，行距调到 30mm
                // 时尺度 30 > 环宽 32 的一半，开运算把整圈削成几段弧，内点被剔掉
                // 一半（6977→3653），轨迹只能盖住 250° 的环。
                double cellMesh = kReconGap;
                {
                    Vec2 rc; double rin=0, rout=0, rw=0; int ns=0;
                    if(detect_ring_uv(uv, kReconGap, rc, rin, rout, &rw, &ns)){
                        cellMesh = std::min(kReconGap,
                                            std::max(1.5, rw/6.0));
                        fprintf(stderr, "[重建] 面片%zu 环形(内径%.1f 外径%.1f 宽%.1f 方位%d/36): 重建格距 %.1fmm (固定尺度 %.1f)\n",
                                i, rin, rout, rw, ns, cellMesh, kReconGap);
                    } else {
                        // 类环形（方框/椭圆/跑道/开口环）：同一件窄带，格距不跟着环宽收，
                        // 开运算同样会把整圈削成几段弧 —— 环越窄、削得越狠，
                        // 后面的轨迹就只能在残存的几段弧上布行（表现为"环上只有半圈有轨迹"）。
                        BandUV bi;
                        if(detect_band_uv(uv, bi)){
                            cellMesh = std::min(kReconGap, std::max(1.5, bi.wPerp/6.0));
                            fprintf(stderr, "[重建] 面片%zu 类环形(内径%.1f 外径%.1f 带宽%.1f 方位%.0f%%): 重建格距 %.1fmm (固定尺度 %.1f)\n",
                                    i, bi.rInMed, bi.rOutMed, bi.wPerp, bi.closure*100.0, cellMesh, kReconGap);
                        }
                    }
                }
                okMesh = build_plane_mesh(uv, cellMesh, params->fillHoleSize, params->minFaceWidth,
                                          verts, tris,
                                          SPRAYTRAJ_PLANE_MAX_POLY_VERTS,
                                          SPRAYTRAJ_PLANE_MAX_POLY_TRIS,
                                          &keepMask, &fillRect, &meshLoops);
            }
            // 补全面：记下矩形四角的世界坐标，轨迹阶段据此在补全区布行
            if(!fillRect.empty() && fillRect.size() >= 4){
                planeFilled[i] = 1;
                for(int k=0;k<4;++k)
                    planeFillRect[i][k] = P.c + uf*fillRect[k].x + vf*fillRect[k].y;
            }
            // 记下该面的重建轮廓与局部基，供轨迹阶段裁剪行端点。
            // **坐标系必须统一**：build_plane_mesh 产出的轮廓是 PCA 旋转后基 (uf,vf) 下
            // 的坐标，而轨迹阶段的 P2/行坐标用的是 (uu,vv)=cross(z,n)/cross(n,uu)。
            // 两者相差一个 ang 旋转（uf=uu*ca+vv*sa, vf=vv*ca-uu*sa），这里就地转成
            // (uu,vv) 口径存下来，轨迹侧再乘 theta 旋转即可直接和 xs/xe/yc 比较。
            {
                const double ca=std::cos(ang), sa=std::sin(ang);
                std::vector<std::vector<Vec2>> conv; conv.reserve(meshLoops.size());
                for(auto& lp : meshLoops){
                    std::vector<Vec2> q; q.reserve(lp.size());
                    for(auto& t : lp)
                        q.push_back(Vec2{ t.x*ca - t.y*sa, t.x*sa + t.y*ca });  // (uf,vf) -> (uu,vv)
                    conv.push_back(std::move(q));
                }
                planeLoops[i] = std::move(conv);
            }
            planeBaseC[i] = P.c; planeBaseU[i] = uf; planeBaseV[i] = vf;
            res.planes[wi].polyVertCount = 0;
            res.planes[wi].polyTriCount  = 0;
            // 最小可喷宽度：用重建轮廓的最小外接矩形短边，而不是点云包围盒。
            // 重建轮廓已经过开运算/孔洞填补/直边规整，能代表真正要显示和布行的面；
            // 线状退化面在这里表现为一条细带（短边≈重建格距或点距）。
            double faceW = -1.0, faceL = -1.0;
            if(okMesh) faceW = minRectShortSide(verts, &faceL);
            if(okMesh && params->minFaceWidth > 1e-9 &&
               (faceW < 0.0 || faceW < params->minFaceWidth)){
                fprintf(stderr,
                        "[重建] 面片%zu: 宽度 %.1fmm < 最小可喷宽度 %.1fmm，判为线状/退化面，不输出该面\n",
                        i, faceW, params->minFaceWidth);
                continue;
            }
            // ---- 内点回筛（去掉"面外散点"）----
            // 语义：重建出来的多边形是"这个面到底覆盖哪里"的唯一权威。
            //
            // 旧写法带一条 `kept >= inl/2` 的保护，本意是怕形态学把好面削光导致丢面；
            // 但它在"这个面本身就是稀疏过渡带"时必然失效——两个平行面之间的棱带被区域
            // 生长拟合成一个虚面，开运算只留下沿边沿的碎片，回筛只保住 36% 的点，
            // 条件不成立 ⇒ 整个碎片点集原样进入布行 ⇒ 轨迹按 P0.5/P99.5 包围盒张开
            // （现场实测：675×418mm 的空框，而重建面只有 3 块共 14,467mm² 的角块）。
            //
            // 所以这里改成**无条件按掩膜回筛**：面外散点一律不参与布行。丢面风险由
            // 下面的 "kept<20 就保持原样" 兜住（点数太少的面本来也布不出行）。
            //
            // 注意：**不能按覆盖率/碎片化直接丢面**。真实面也可能是稀疏或碎片的
            // （curved 用例的斜面内点 7429 个散在 145,386mm² 上，重建碎成 3 块、
            // 最大块仅 53%，但它是实实在在的一块面），一律丢面会漏喷。
            // 真正要抑制的是"行越出面的范围"，那件事在布行阶段按重建轮廓裁剪解决
            // （build_plane_mesh 现在会回传重建轮廓 loopsOut，genRows 用它夹住行端点）。
            const size_t inputCount = segs[i].inl.size();
            if(okMesh && keepMask.size()==inputCount){
                std::vector<int> kept; kept.reserve(segs[i].inl.size());
                for(size_t t=0;t<keepMask.size();++t) if(keepMask[t]) kept.push_back(segs[i].inl[t]);
                if(kept.size() >= 20){
                    if(kept.size() != segs[i].inl.size())
                        fprintf(stderr, "[重建] 面片%zu: 剔除面外内点 %zu -> %zu（保留 %.0f%%）\n",
                                i, segs[i].inl.size(), kept.size(),
                                100.0*(double)kept.size()/(double)std::max<size_t>(1,segs[i].inl.size()));
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
                    planeFrag[i] = (nc > 1 && lf <= 0.50) ? 1 : 0;
                    sparseDisconnected[i] = (nc > 1 && segs[i].inl.size() <= inputCount/2);
                    fprintf(stderr, "[重建] 面片%zu: 顶点=%d 三角=%d（规则轮廓）连通块=%d 最大块占比=%.0f%% 面积=%.0fmm2 宽度=%.1fmm 内点=%zu%s\n",
                            wi, nv, nt, nc, lf*100.0, ar, faceW, segs[i].inl.size(),
                            planeFrag[i] ? "  ← 碎片化，轨迹按轮廓切段" : "");
                }
                // ===== 环级占比判据（2026-11-20）=====
                // D1/D5/D6 判的都是"面"，而碎片化面在重建后被切成多个互不相连的环，
                // 每环各自布一行：现场看到的就是"一个面上叠着几块小方片，各带一圈空移"。
                // 这些环单独看点数占比都很小，却从没被当成"面"审过 —— 判据在面这一层
                // 就已经放行（p047 面2：4 环、最大块 38%、整面占比 3.3%）。
                // 这里对**碎片化面**补一次逐环占比：所有环都 < ringShare 时整面不喷。
                // 只在"整面全碎"时丢，不做逐环丢（逐环丢要重建网格，且 curved 用例
                // "一块真实斜面带碎成 3 块"会被误伤，见 5722 行的教训）。
                if(planeFrag[i] && planeLoops[i].size() > 1){
                    double ringShare = 0.10;   // 与 D5 同口径（用户 2026-11-20 14:19：3% → 10%）
                    if(const char* e = std::getenv("SPRAYTRAJ_RING_SHARE")){
                        const double rv = std::atof(e);
                        if(rv >= 0.0) ringShare = rv;
                    }
                    if(ringShare > 0.0){
                        std::vector<size_t> cnt(planeLoops[i].size(), 0);
                        for(int k : segs[i].inl){
                            Vec3 d = Pf[k] - P.c;
                            Vec2 q{ dot(d,uu), dot(d,vv) };
                            for(size_t r=0;r<planeLoops[i].size();++r)
                                if(pointInPoly(q, planeLoops[i][r])){ ++cnt[r]; break; }
                        }
                        size_t mx = 0; for(size_t c2 : cnt) mx = std::max(mx, c2);
                        const double share = (double)mx / (double)std::max<size_t>(1, Pf.size());
                        if(mx > 0 && share < ringShare){
                            fprintf(stderr, "[split] 环级丢弃: 面%zu 碎片化 %zu 环，最大环仅 %zu 点 = %.2f%% < %.1f%% → 整面不喷\n",
                                    i, planeLoops[i].size(), mx, share*100.0, ringShare*100.0);
                            res.planes[wi].polyVertCount = 0;
                            res.planes[wi].polyTriCount  = 0;
                            continue;
                        }
                    }
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
                planeSeg.push_back((int)i);
                ++wi;
            } else {
                fprintf(stderr, "[重建] 面片%zu: 轮廓提取失败（点太少/全散点），不输出该面\n", i);
            }
        }
    }
    res.plane_count = (int)wi;


    // 3. 对每个分割出的平面段，分别生成 S 形轨迹并拼接
    struct Row { double xs, xe, yc; int segIdx; };  // segIdx: 该行是扫描线与面交集的第几段（环形面会产生多段）
    // ---- 行端点按"重建轮廓"夹紧 ----
    // 病灶：布行阶段用 P0.5/P99.5 分位包围盒决定"这一行铺到哪"，而重建面多边形
    // 早就算出来了却没用上。包围盒会把"面在这一行根本没有材料"的空档也算进长度，
    // 于是行直接横穿空处。辅助 gapSeg=fillHole(默认 300mm) 更放大了这个问题 —
    // 300mm 以内的空洞一律被"桥接"成连续行，等于对着空气喷 300mm。
    // 现场实测：重建面只有 3 块共 14,467mm² 的角块，轨迹却张开成 675×418mm 的框，
    // 最远超出面轮廓 338mm。
    //
    // 修法：把重建轮廓作为行范围的唯一权威。每行先按支撑点聚出初段，再与该
    // yc 水平线同轮廓的交集求交，交集之外一律不喷。轮廓已经把"孔洞该不该填"
    // 的决策做完（fillHole 以内的洞在重建时已填实、以外的是真空），所以这里
    // 不需要再按 gapSeg 猜一次，两侧口径自然一致。
    auto clipRowByLoops=[&](double yc, double xs, double xe,
                           const std::vector<std::vector<Vec2>>& loops)->std::vector<std::pair<double,double>>{
        std::vector<std::pair<double,double>> out;
        if(loops.empty()){ out.emplace_back(xs,xe); return out; }
        std::vector<std::pair<double,double>> iv;      // 该 y 高度上的轮廓内部区间（已合并）
        for(const auto& lp : loops){
            const size_t n = lp.size();
            if(n < 3) continue;
            std::vector<double> xs0;
            for(size_t k=0;k<n;++k){
                const Vec2& a = lp[k]; const Vec2& b = lp[(k+1)%n];
                if((a.y > yc) == (b.y > yc)) continue;   // 不相交 / 水平边
                const double t = (yc - a.y)/(b.y - a.y);
                xs0.push_back(a.x + (b.x - a.x)*t);
            }
            if(xs0.size() < 2) continue;
            std::sort(xs0.begin(), xs0.end());
            for(size_t k=0;k+1<xs0.size();k+=2) iv.emplace_back(xs0[k], xs0[k+1]);
        }
        if(iv.empty()) return out;                       // 该行整条不在面内
        std::sort(iv.begin(), iv.end());
        std::vector<std::pair<double,double>> mg;         // 合并重叠区间
        for(auto& r : iv){
            if(!mg.empty() && r.first <= mg.back().second + 1e-9)
                mg.back().second = std::max(mg.back().second, r.second);
            else mg.push_back(r);
        }
        const double minKeep = std::max(10.0, kReconGap);  // 短于重建尺度的残段没有喷涂价值
        for(auto& r : mg){
            const double a = std::max(xs, r.first), b = std::min(xe, r.second);
            if(b - a >= minKeep) out.emplace_back(a, b);
        }
        return out;
    };
    double L=params->scanRowSpacing, pstep=params->pointSpacing;
    // 跨行方向外扩：首末行相对点云边界各向外多少（mm）。<0 = 自动 = 0.25×行距，
    // 与行端点外扩 edgeExtend 是两个方向、两套默认值：
    //   · edgeExtend（行方向）默认 1×行距——点云在行端缺的那一截通常较大；
    //   · 跨行方向只差"网格截断 + 点云边缘薄薄一层"，0.25×行距已能覆盖。
    const double acrossExt = (params->edgeExtendAcross >= 0.0)
                           ? params->edgeExtendAcross
                           : 0.25*params->scanRowSpacing;
    double turnR = params->turnRadiusFactor * L;

    // 主方向选择 + 行端点计算：对给定 2D 点集，穷举角度(每5°)选最平滑方向
    // 行端/行范围用 P0.5/P99.5 分位数包围盒（全长覆盖，防长尾离群点，也防稀疏端部被行过滤裁掉）
    auto quant=[&](std::vector<double> v, double f)->double{
        if(v.empty()) return 0.0;
        size_t k=(size_t)((v.size()-1)*f);
        std::nth_element(v.begin(), v.begin()+k, v.end());
        return v[k];
    };
    auto genRows=[&](const std::vector<Vec2>& pts, double ang, double L, double fillHole, double acrossExt)->std::vector<Row>{
        double c=std::cos(ang), s=std::sin(ang);
        std::vector<Vec2> R(pts.size());
        std::vector<double> rx(pts.size()), ry(pts.size());
        for(size_t i=0;i<pts.size();++i){
            R[i]={c*pts[i].x+s*pts[i].y, -s*pts[i].x+c*pts[i].y};
            rx[i]=R[i].x; ry[i]=R[i].y;
        }
        double mnx=quant(rx,0.005), mxx=quant(rx,0.995);
        double mny=quant(ry,0.005), mxy=quant(ry,0.995);
        // 小面（长向本身 < 行距）允许一条“全长行”通过：41×132 的板在
        // rowSpacing=150 时，132mm 的行会被下面的 ">= L" 当碎段丢掉，
        // 整面直接输出“没有可用行”。大面仍用 L，碎片过滤不变。
        const double minSeg = std::min(L, std::max(1.0, 0.5*std::max(mxx-mnx, mxy-mny)));
        std::vector<Row> rs;
        // 支撑点阈值自适应：小平面（如窄侧边）点数少，用相对阈值；
        // 大平面用绝对阈值防噪声孤点。
        int supMin = (int)std::max(3, (int)(pts.size()*0.001));
        supMin = std::min(supMin, 20);
        // ---- 首末行位置（跨行方向的外扩策略）----
        // 从点云边界 mny/mxy 各向外跨 acrossExt，行数按行距推算，行距在两个端点
        // 之间等分——首末行严格落在外扩后的边界上。
        // 原实现是"从 mny+0.25L 起、每 L 一步、直到超过 mxy+0.25L"：末行位置由固定
        // 网格决定，最多少掉将近一整个行距（实测薄壁盒内壁顶部末行停在 z=281，真实
        // 边沿 z=300，漏 19mm；补行端点没用，因为漏的方向是跨行方向不是行方向）。
        // 等分后行距偏差 = (跨度/行数) - L，几十行时 < 1%，肉眼与工艺上都无差别。
        // 行端按"该行带内点的 x 聚类分段"：实心面 → 一整行（全长走刀）；
        // 环形/框形面（如槽型盘凸边环）→ 每行分成左右两段，中央空区不喷。
        // 点云缺陷覆盖重建：间隙 < fillHole 的小洞/缺口视为可填补——直接
        // 跨越生成连续行（重建面已把该洞补齐，行应贴重建面连续无断点）；
        // 间隙 >= fillHole 的真实孔/槽才断开为两段。
        //
        // 曾经的尝试：把阈值改成 max(L*1.2, 20) 想抑制"对着空气喷 300mm"。
        // **不能这么做**——线扫点云沿行方向的相邻点距本身就受扫描遮挡影响，
        // 稀疏面上 40~80mm 的"点间隙"随处可见（不是真空），阈值降到 1.2L 后
        // norim 行数 195→864、longstripe 47→113，轨迹碎成一地。
        // 正确的解耦方式是让"行要不要断开"跟着**重建轮廓**判断，而不是拿点间隙猜：
        // 见下方按重建轮廓修正行的处理（只对轮廓相对行严重缩水的行动手）。
        const double gapSeg = std::max(fillHole, L*0.5);   // x 间隙超过此值视为两段；以下则跨越填补
        // 跨行外扩不许越出面的真实边界：面外没有材料可喷，扩出去只是空走+撞壁风险。
        // 实测开口盒内壁（行距80、-ea 默认16）：首行被推到 z=316（越顶 16mm），行端点
        // 再沿行外扩 8mm 到距邻壁 12mm，斜向净空 sqrt(12²+16²)=20.0 恰好卡在阈值上；
        // 连接段从这一点出发判定全部失败 → 该面每个换行都抬枪绕行，轨迹碎成一截截。
        // 行中心压在真实边沿上即足够（喷幅本身有 20~30% 搭接），越界纯属多余。
        const double ryLo = *std::min_element(ry.begin(), ry.end());
        const double ryHi = *std::max_element(ry.begin(), ry.end());
        const double y0 = std::min(std::max(mny - acrossExt, ryLo), mxy);
        const double y1 = std::max(std::min(mxy + acrossExt, ryHi), mny);
        int nRow = 1; double stepA = L;
        if(y1 - y0 > 1e-9){
            // 行数按跨度/行距四舍五入（行距允许比设定值大 20% 以内——喷涂本来有
            // 20~30% 搭接余量）；超过 20% 才补一条，防止窄面上两条行被拉得太开。
            // 不能用 floor 起算的老写法：那样末行由固定网格决定，最多少掉近一个行距。
            // 跨度不足半个行距时只走一条中心行：喷幅本身有 20~30% 搭接余量，
            // 一条行足以盖住整条窄带。旧写法用 max(1,...) 保证 nIv>=1，等于强制
            // 布两行，行距被压成面宽的一半——实测 26mm 宽的窄面被布成 13mm 行距，
            // 同一块地方喷两遍，行程直接翻倍（半棍子 6 个面因此从 3.3m 涨到 6.6m）。
            // nRow==1 时下方 yc 走 (y0+y1)/2 分支，中心行压在中线上，不会偏。
            // 跨度不超过一个行距时只走一条中心行：行带按 |Δy| <= L/2 采集，一条行
            // 的喷幅正好盖住整条带，布两条只是把同一块地方喷两遍。
            // （旧写法门槛是 0.5L，0.5L~1L 这一段仍被布成两条行，行距被压成面宽。）
            int nIv = 0;
            if((y1-y0) <= L){
                nRow  = 1;
                stepA = L;
            }else{
                // ---- 首末行必须落在边界上（不要改成"内移 L/2"）----
                // 试过内移半步（喷幅边缘正好压边沿，理论上更省一条行），实测覆盖率
                // 反而掉：norim 300×300 侧壁 97%→91%（缺口 max 29→67mm）、cao 面2
                // 100%→69%。原因是这些面靠"末行压在边沿"来盖住最后一条窄带，
                // 内移后那一条带没人喷。省行程不能拿漏喷换，故恢复为边界对齐。
                nIv = (int)std::floor((y1-y0)/L + 0.5);
                if((y1-y0)/(double)nIv > 1.2*L) ++nIv;
                nRow  = nIv + 1;
                stepA = (y1-y0)/(double)nIv;
            }
        }
        std::vector<Row> raw;
        for(int ri=0; ri<nRow; ++ri){
            const double yc = (nRow>1) ? (y0 + stepA*(double)ri) : 0.5*(y0+y1);
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
                if(segCnt>=supMin && (segPrev-segStart) >= minSeg)
                    raw.push_back({segStart, segPrev, yc, segIdx});
                segIdx++;
            };
            for(size_t i=1;i<bx.size();++i){
                if(bx[i]-segPrev > gapSeg){ flushSeg(); segStart=bx[i]; segCnt=0; }
                segPrev=bx[i]; segCnt++;
            }
            flushSeg();
        }
        // 极端行补回：跨行外扩把首/末行推到点云之外时，那条行带的支撑点太少
        // （带有一半落在空处），会被上面的支撑点/段长过滤掉——恰恰是最该保留的
        // "补边沿"那条行被删了，表现为"ea 调大反而少一行"。这里按最近的已有行
        // 克隆一条到目标位置（长度沿用邻居），后面还有枪位兜底会把它退回来。
        if(nRow > 1 && acrossExt > 1e-9 && !raw.empty()){
            auto cloneAt=[&](double yTarget){
                double best = 1e18;
                for(auto& r : raw) best = std::min(best, std::fabs(r.yc - yTarget));
                if(best < 1e-6) return;                 // 该行本来就在
                std::vector<Row> add;
                for(auto& r : raw)
                    if(std::fabs(std::fabs(r.yc - yTarget) - best) < 1e-6){
                        Row q = r; q.yc = yTarget; add.push_back(q);
                    }
                for(auto& q : add) raw.push_back(q);
            };
            cloneAt(y0);
            cloneAt(y1);
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
        if(std::getenv("SPRAYTRAJ_ROW_DBG"))
            fprintf(stderr, "[rows] ang=%.1f° 跨行 %.1f..%.1f 跨度=%.1f 行数=%d 间距=%.2f 出行=%zu\n",
                    ang*180.0/3.141592653589793, mny, mxy, mxy-mny, nRow, stepA, rs.size());
        if(rs.empty()) return rs;
        return rs;
    };

    std::vector<SprayTrajPoint> allPts;
    // 每个面的"面内轨迹"在 allPts 里的区间，供最后的重排阶段使用。
    // 只记区间、不改任何生成逻辑：区间内是该面自己的点（喷涂点 + 面内换行空移），
    // 区间外是面间过渡点。重排时会丢掉旧的过渡点、按新顺序重建（见"5.面顺序重排"）。
    struct FaceSpan { size_t beg, end, si; };
    std::vector<FaceSpan> faceSpan;
    int totalRows = 0;
    int nDownSkip = 0;      // 朝下超出角度窗口、不喷的面数
    int nGunReject = 0;     // 枪位验收未通过（两侧都伸不进枪）、不布行的面数
    int nGunFlip = 0;       // 枪位验收翻符号纠正的面数
    int nPlanarSkip = 0;    // 平面性不足（误合并/强曲面），不按平面布行的面数
    ProfStage psTraj("4.轨迹生成", &tTraj);

    // ---- 面间过渡的安全高度 ----
    // 面间转移是空移，但控制器会把相邻两个点用直线插值走完。原来只补一个
    // spray=0 哨兵点（位置不变），下一点直接是另一面的首点 —— 实际运动就是一
    // 根穿过工件凸台/槽壁的斜线。工业上必须"抬枪绕行"，否则撞件。
    // 处理：沿上一面法向退出 clear → 抬到安全高度 → 安全高度上平移 → 下降接下一面首点。
    //   clear 取 standoff：喷头本来就离面 100mm，再退 100mm 即彻底离开该面。
    //
    // ---- 让开方向 liftU：这是第二十六轮的关键修正，第二十七轮（2026-11-11）改口径 ----
    // 原来一律让到 z = 常数（内部 +Z）。这条隐含假设是"可喷侧 ≈ 竖直向上"，也就是
    // 工件朝上摆、相机俯拍。侧装/立式悬挂工件上它直接崩掉：
    //   背板_X+800 实测可喷侧 ≈ 基座 +X（水平），重力向上 = 基座 +Y，
    //   而内部 +Z 经两次旋转后落到基座系的另一个**水平**轴上。
    //   于是"让开"变成在工件前面横划 400~900mm（Z 向跨度就那么大），
    //   竖直方向只让了几十毫米 —— 看着"空移绕了老远"的根因就是这个。
    // 第二十六轮改成"沿重力向上让开"，背板实测反而更糟：空移从工件顶上翻过去，
    //   一次面间过渡空移净位移 ΔY=336 / ΔZ=523、单段走 1.8 米，
    //   而工件只有 866×539×862 —— 工件上沿外面绕了一大圈。
    // 第二十七轮（用户 2026-11-11 拍板）：让开方向 = **基座 +X，机器人正面**。
    //   沿 +X 退让 = 枪先撤到工件前方（这一侧没有夹具/背架，机器人一定探得到），
    //   再从前方平移过去、降回下一面 —— 物理上最短、也最接近人的直觉。
    // 三轴全 0 时退回 +Z，行为与最早版本**逐位一致**。所有"让到安全面"的点统一由 liftTo() 产出。
    //
    // ⚠ 只给"轴"不够，还要定"往哪一侧退"：
    //   +X 这条轴有两个方向。工件厚 130mm、贴在 X∈[761,891]，正面朝机器人（X 小的一侧）。
    //   若一律退到 s 最大侧（X≈891 外面的背面），每次过渡都得先穿过工件本体才到得了
    //   那张安全面 —— 实测就是这个原因让侧壁面（918 个喷点）判成"无可行退出路径·整面不喷"。
    //   退让侧应该由**枪的位置**定：枪 = 面点 + standoff·n，哪一侧有枪就从哪一侧走。
    //   这里按「全部喷涂点 s(枪位) 的中位数」定侧，取多数面所在的那侧；
    //   退化（两侧各半，如对喷的盒内壁）时沿用最大侧，与旧行为一致。
    Vec3 liftU{0,0,1};
    bool hasLiftDir = false;
    {
        const Vec3 d{params->liftDirX, params->liftDirY, params->liftDirZ};
        if(len(d) > 1e-9){
            // 与点云同做入口旋转：点云被旋到 up=+Z，方向向量也必须跟着旋。
            liftU = norm(uf.fwd(d));
            hasLiftDir = true;
            const double cAng = std::acos(std::max(-1.0, std::min(1.0, liftU.z)));
            fprintf(stderr, "[traj] 空移让开方向 = 点云(%.2f,%.2f,%.2f) → 内部(%.2f,%.2f,%.2f)"
                            "，与可喷侧夹角 %.1f°%s\n",
                    d.x,d.y,d.z, liftU.x,liftU.y,liftU.z,
                    cAng*180.0/3.141592653589793,
                    (cAng > 60.0*3.141592653589793/180.0)
                        ? "（两者接近垂直：旧口径会横着让，本轮修正生效）" : "");
        }
    }
    // 沿 liftU 的高度：s(p) = p·liftU。liftU=+Z 时 s(p) 就是 p.z，与旧代码逐位等价
    // （p.x*0 + p.y*0 + p.z*1 的加法都是精确的，但这里仍显式走快速路径，免得
    //  编译器把 dot 展开成 mul+add 时引入 fma 导致末位不同）。
    //
    // ⚠ 不要再加"按枪位自动翻侧"这类隐式纠正：
    //   调用方给的方向就是**从工件指向枪那一侧**（用户口径：基座 +X 是正面，
    //   工件正面朝机器人，枪就从这一侧伸过来），方向本身就是对的。
    //   我加过一次"枪位中位数在负侧就翻过来"，结果把 +X 翻成了 −X ——
    //   因为 liftU 在**内部系**（点云被入口旋转 180°），而 segs[].pl.n 是**点云系**，
    //   拿 pl.n 直接点乘 liftU 等于混两个坐标系，符号整体反掉，
    //   安全面于是落到工件背面，侧壁面被判"无可行退出路径·整面不喷"。
    //   要判也必须先用 uf.fwd() 把法向转到内部系 —— 但既然方向由调用方给定，
    //   这层推断本身就不该存在。
    const bool liftIsZ = (std::fabs(liftU.x) < 1e-12 && std::fabs(liftU.y) < 1e-12
                          && liftU.z > 0);
    auto sOf = [&](const Vec3& p)->double{ return liftIsZ ? p.z : dot(p, liftU); };
    double sPartMax = -1e18;
    for(const auto& p : Pf){ sPartMax = std::max(sPartMax, sOf(p)); }
    // 把 p 沿 liftU 挪到高度 sTarget（垂直于 liftU 的平面内的分量保持不动）。
    //   liftU=+Z 时 = {p.x, p.y, sTarget}，正是旧写法的 {p.x, p.y, zSafe}。
    // 这个快速路径不是优化，是**正确性**：p + u*(s - p·u) 在浮点下有约 10%
    // 的点 z 会与 s 差 1 ulp（实测 200k 次里 19018 次），而"默认参数下与旧版本
    // 逐位一致"是这条改动的验收标准（现场要能直接对比改前/改后的 CSV）。
    auto liftTo = [&](const Vec3& p, double sTarget)->Vec3{
        if(liftIsZ) return Vec3{p.x, p.y, sTarget};
        return p + liftU*(sTarget - sOf(p));
    };
    // ---- 安全面：必须落在**枪所在的那一侧** ----
    //
    // 推导：枪位 = 面点 + standoff·n（n 指向可喷侧 = 枪那一侧），故
    //   s(枪位) = s(面点) + standoff·(n·liftU)，落在 s 的枪侧外沿之外一点点。
    // 所以"让开面"就是「工件沿 liftU 的极值 + 1.3×standoff」，取哪一侧的极值，
    // 取决于枪在哪一侧。
    //
    // 旧口径（liftU = +Z，工件朝上、相机俯拍）里枪在上方 ⇒ 枪在 s **最大**侧，
    // 于是 sSafe = sMax + 1.3×standoff，这就是旧的 zSafe，**逐位一致**（不许改）。
    //
    // 但立式悬挂工件上 liftU 变成了水平轴（基座 ±X），枪在工件**正面**，
    // 沿这根轴的 s 可能是**最小**侧：背板实测 s ∈ [474, 1300]（点云 z，越大越靠里），
    // 枪位 s ≈ 1333 落在 s 小的一端之外。这时还取 sMax 侧，等于把安全面放到工件
    // **背面**，每次面间过渡都得先穿过板体才到得了 —— 实测侧壁面（918 个喷点）
    // 就是这么被判成"无可行退出路径·整面不喷"的。
    //
    // 判据：拿每个面的枪位 s 取中位数，看它离 sMin 近还是离 sMax 近。
    // 枪侧那一端就往外让 1.3×standoff。liftU 指向哪一头都不影响结论 ——
    // 用户给的只是"轴"，物理上必须退向有枪的那一侧，否则必穿件。
    const double clear = params->standoff;
    double sSafe;
    if(!hasLiftDir){
        // ---- 默认口径（liftDir 三轴全 0）：逐位保持一致，不许改 ----
        // 这条分支必须与旧版本字节级相同：现场要靠改前/改后 CSV 直接 diff 验收。
        // 所以**不进**任何枪侧推断，哪怕推断在默认口径下也会得出同样的数 ——
        // 只要多算一个 dot/分支，浮点结果就不再保证逐位一致。
        sSafe = sPartMax + 1.3*params->standoff;
    }else{
        // ---- 指定了让开方向：安全面必须落在**枪所在的那一侧** ----
        //
        // ⚠ 关键：不能拿"全部枪位的中位数"投票。
        //   背板实测（点云系=基座系，liftU=基座 +X）：
        //     面0/1 法向≈(-0.996, 0.049, 0.069) —— 与 liftU 近乎**反平行**，
        //       这两面 9.7 万点、是工件正反两大面，枪位 X≈711（在工件外侧小端）。
        //     面2   法向≈(0.084, 0.997, 0.003) —— 与 liftU **垂直**，侧壁，
        //       它的枪位 X≈845 落在工件主体内部，对"轴的两个方向"毫无信息。
        //     面3   斜向，也无信息。
        //   若把面2/面3 的枪位混进中位数，中位数被拉到 853（工件内部），
        //   判据恒为"枪在最大侧"，安全面落到工件背面 —— 实测就是侧壁面被判
        //   "无可行退出路径·整面不喷"（丢 918 个喷点）的直接原因。
        //
        //   正确做法：只有法向与 liftU **不垂直**的面才有资格投票，权重取
        //   |n·liftU|（越平行越说明"枪沿这根轴站在哪一头"）。判据本身用
        //   "投票枪位的加权中位数离哪一端更近"。
        //
        // ⚠ 边界必须用稳健分位，不能用 min/max：
        //   点云 X∈[24,891]，但 X<760 只有 146 点（0.14%，飞点），主体在 [760,891]。
        //   用 raw min 会让 sMin=24，dist(枪位,24) 永远大于 dist(枪位,sMax)。
        std::vector<double> ss;
        ss.reserve(Pf.size());
        for(const auto& p : Pf) ss.push_back(sOf(p));
        double sPartLo = sPartMax, sPartHi = sPartMax;
        double sRawMin = sPartMax;
        if(!ss.empty()){
            std::sort(ss.begin(), ss.end());
            sRawMin = ss.front();
            const size_t lo = (size_t)(0.005 * (double)(ss.size() - 1));
            const size_t hi = (size_t)(0.995 * (double)(ss.size() - 1));
            sPartLo = ss[lo];
            sPartHi = ss[hi];
        }

        // 逐面投票：只有 |n·liftU| > 0.3 的面参与（近垂直的侧壁不投票）
        struct Vote { double s; double w; };
        std::vector<Vote> votes;
        votes.reserve(segs.size());
        for(const auto& s : segs){
            const double c = dot(s.pl.n, liftU);
            const double wgt = std::fabs(c);
            if(wgt <= 0.3) continue;                     // 侧壁：对轴向无信息
            // 注意：pl.c / pl.n 出自 split_planes(Pf)，Pf 已是**内部系**点，
            // 与 sOf/liftU 同系，**不要**再套 uf.fwd（双重旋转会让符号整体反掉）。
            votes.push_back({dot(s.pl.c + s.pl.n * params->standoff, liftU), wgt});
        }

        double gunS = 0.0;
        bool haveGun = false;
        if(!votes.empty()){
            std::sort(votes.begin(), votes.end(),
                      [](const Vote& a, const Vote& b){ return a.s < b.s; });
            double wSum = 0.0;
            for(const auto& v : votes) wSum += v.w;
            double acc = 0.0;
            gunS = votes.back().s;
            for(const auto& v : votes){
                acc += v.w;
                if(acc >= wSum * 0.5){ gunS = v.s; break; }
            }
            haveGun = true;
        }
        const bool gunOnMaxSide = !haveGun
                               || std::fabs(gunS - sPartHi) <= std::fabs(gunS - sPartLo);
        sSafe = gunOnMaxSide ? (sPartHi + 1.3*params->standoff)
                             : (sPartLo - 1.3*params->standoff);
        fprintf(stderr, "[traj] 安全面: 加权枪位 s=%.1f（%zu/%zu 面投票，侧壁不投票），"
                        "工件沿让开轴稳健 s∈[%.1f,%.1f]（原始 [%.1f,%.1f]）→ 枪在%s侧，让到 s=%.1f\n",
                gunS, votes.size(), segs.size(), sPartLo, sPartHi, sRawMin, sPartMax,
                gunOnMaxSide ? "远(最大 s)" : "近(最小 s)", sSafe);
    }
    Vec3 prevNormal{0,0,1};
    int transCount = 0;
    int nGapSafe = 0;          // 曲面行间直线会穿件→改走安全空移的次数
    int nGunOnLink = 0;        // 换面短接（枪不停）的次数
    int nApproachFix = 0;   // 面间过渡：垂直直降会穿件、改沿法向进入的次数
    int nDropFace = 0;      // 面间过渡：没有任何可行进入路径（窄缝/深腔）、整面不喷的面数

    // ---- 空移段"直线可达"判定 ----
    // 原来面间过渡一律"沿法向退出 → 抬到全局安全高度 zSafe → 平移 → 下降"，
    // 哪怕两个端点本来就在工件上方很高处。垃圾桶"锥面顶 → 顶环"就是这样：
    // 两端 z=231.9 / 289.9，直线只有 191mm，却被抬到 391 再落下，走了 507mm。
    // 更糟的是"顶环 → 内腔底面"这种**下降进入开口腔**的过渡：腔口是敞开的，
    // 根本没有障碍，还是先抬到 391 再垂直落 291mm，走了 506mm。
    //
    // 现在改成直接判几何：沿直线等距采样，逐点查最近工件点，全程都 >= moveClear
    // 就直连。这比"两端 z 都高于工件最高点"宽松得多，也严格得多——后者对
    // 上面第一种情形都会漏判（231.9 比 zPartMax+0.5·standoff=241.1 只低 9mm）。
    KdTree partTree;
    // 安全层点集：布行层代表点（planIdx 非空时）。KD 树是这里最大的一笔开销，
    // 而"枪位到工件最近距离"这种判据在 2mm 网格上与全分辨率差不到 1mm
    // （代表点本身取自实测点，几何只会变稀不会变形），所以安全层吃抽稀点。
    // 索引口径要一致：ptFace 也按同一套下标建，否则整树剪枝会错位。
    std::vector<int> ptFace;
    std::vector<Vec3>    safePts;
    const std::vector<Vec3>* safeSrc = &Pf;
    if(!planIdx.empty()){
        safePts.reserve(planIdx.size());
        for(int k : planIdx) safePts.push_back(Pf[(size_t)k]);
        safeSrc = &safePts;
    }
    ptFace.assign(safeSrc->size(), -1);
    std::vector<unsigned long long> faceSkipMask(segs.size(), 0ULL);
    if(segs.size() <= 64){
        // 原始点下标 → 抽稀后下标 的映射（同一套顺序，所以是"第几个代表点"）
        std::vector<int> pos(Pf.size(), -1);
        if(!planIdx.empty()){ for(size_t t=0;t<planIdx.size();++t) pos[(size_t)planIdx[t]] = (int)t; }
        for(size_t k=0;k<segs.size();++k){
            for(int i : segs[k].inl){
                const int t = planIdx.empty() ? i : pos[(size_t)i];
                if(t >= 0 && (size_t)t < ptFace.size()) ptFace[(size_t)t] = (int)k;
            }
        }
    }
    { ProfStage ps("4a1.点云KD树"); partTree.init(*safeSrc, segs.size() <= 64 ? &ptFace : nullptr); }
    const double moveClear = std::max(6.0, 0.25*params->standoff);
    // 枪体安全半径：内壁/内底布行 + 行端点外扩两处共用同一个口径
    const double gunClearEff = params->gunClear > 0.0
                             ? params->gunClear
                             : std::max(12.0, 0.2*params->standoff);
    KdTree candidateTree;
    bool usedCandidateCheck = false;

    // ---- 面几何"挡板"清单 ----
    // 只拿点云测距会看漏两种地方：
    //   ① 去噪把棱边/顶沿点删掉（边界点邻域只有半边支撑），壁顶一圈没有点，
    //      枪位压在壁面上仍测出"还有 20mm"，顶行外扩就会冲出板边、连接段穿壁；
    //   ② 单层点云、稀疏采样，面对面之间的空隙没有点可测。
    // 面本身是解析的（平面法向 n + 面内 2D 包围盒），点到面的距离可以直接算：
    // 枪位投影落在面范围内时，|带符号距离| 就是该面对枪位的净空。
    // 顶行外扩越过壁顶时投影落到包围盒之外，判定自动失效——正好是"从棱边上方
    // 越过去"应有的行为，不会把正常的越顶外扩拦掉。
    struct BlkFace {
        Vec3 n, c, u, v; double u0, u1, v0, v1; int si;
        bool curved = false;      // 曲面块的拟合平面只用于净空估计，不参与穿面判定
        // 面内占据栅格：矩形包围盒会覆盖面上的空洞（槽口、孔、环形中间的空区），
        // 把"枪从空洞里穿过去"误判成"压在面上"。5mm 一格 + 3×3 邻域容忍稀疏点。
        double cell = 5.0; int gu = 0, gv = 0;
        std::vector<unsigned char> occ;
        bool solidAt(double a, double b) const {
            if(occ.empty()) return true;              // 未建栅格 → 按实心处理（保守）
            const int i = (int)std::floor((a - u0)/cell);
            const int j = (int)std::floor((b - v0)/cell);
            for(int di=-1; di<=1; ++di)
                for(int dj=-1; dj<=1; ++dj){
                    const int ii = i+di, jj = j+dj;
                    if(ii<0 || jj<0 || ii>=gu || jj>=gv) continue;
                    if(occ[(size_t)jj*(size_t)gu + (size_t)ii]) return true;
                }
            return false;
        }
    };
    std::vector<BlkFace> blkFaces;
    {
    ProfStage ps("4a2.面几何栅格");
    for(size_t k=0;k<segs.size();++k){
        const auto& s = segs[k];
        if(s.inl.size() < 20) continue;
        const Vec3 n = s.pl.n;
        Vec3 ax = (std::fabs(n.x) < 0.9) ? Vec3{1.0,0.0,0.0} : Vec3{0.0,1.0,0.0};
        Vec3 u = norm(cross(ax, n));
        Vec3 v = cross(n, u);
        BlkFace bf{n, s.pl.c, u, v, 1e30, -1e30, 1e30, -1e30, (int)k};
        bf.curved = s.curved;
        // 必须用全部内点求范围：抽样会漏掉面最外侧那一圈，包围盒比真实面小
        // 1~2mm，恰好把"枪位压在棱边内侧"这种最该拦的情形放过去。
        for(size_t q=0;q<s.inl.size();++q){
            const Vec3 w = Pf[s.inl[q]] - s.pl.c;
            const double a = dot(w,u), b = dot(w,v);
            bf.u0 = std::min(bf.u0,a); bf.u1 = std::max(bf.u1,a);
            bf.v0 = std::min(bf.v0,b); bf.v1 = std::max(bf.v1,b);
        }
        if(bf.u1 < bf.u0) continue;
        bf.gu = std::max(1, (int)std::ceil((bf.u1-bf.u0)/bf.cell) + 1);
        bf.gv = std::max(1, (int)std::ceil((bf.v1-bf.v0)/bf.cell) + 1);
        bf.occ.assign((size_t)bf.gu*(size_t)bf.gv, 0);
        for(size_t q=0;q<s.inl.size();++q){
            const Vec3 w = Pf[s.inl[q]] - s.pl.c;
            const int i = (int)std::floor((dot(w,u) - bf.u0)/bf.cell);
            const int j = (int)std::floor((dot(w,v) - bf.v0)/bf.cell);
            if(i>=0 && j>=0 && i<bf.gu && j<bf.gv)
                bf.occ[(size_t)j*(size_t)bf.gu + (size_t)i] = 1;
        }
        blkFaces.push_back(std::move(bf));
    }
    }
    // 点到各面几何的最小净空（投影落在面的**实心区**才算）
    auto planeClearAt = [&](const Vec3& g, int skipSi)->double{
        double best = 1e30;
        for(const auto& f : blkFaces){
            if(f.si == skipSi) continue;
            const Vec3 w = g - f.c;
            const double sd = std::fabs(dot(w, f.n));
            if(sd >= best) continue;                 // 已经不可能更近
            const double a = dot(w, f.u), b = dot(w, f.v);
            if(a < f.u0 || a > f.u1 || b < f.v0 || b > f.v1) continue;
            if(!f.solidAt(a, b)) continue;           // 投影落在面上的空洞里 → 没有材料
            best = sd;
        }
        return best;
    };
    // 线段 A→B 是否穿过某个面（交点落在面内包围盒里）。
    // 这是"穿模"的精确判据。gunClear 只测距离，测不出穿越：内壁蛇形换行时折返弧
    // 的 cr 恰好平行于行方向，弧鼓 40mm 到面外（枪位越过 x=400 的壁到 420），
    // 墙前 20mm 与墙后 20mm 在距离上完全一样，旧判据一律放行。
    // 注意不能用"点在面背侧"来判：枪在腔内时，对**所有外壁面**都在背面（合法）；
    // 枪在件外时，对所有内壁面都在背面（也合法）。只有"跨越"才是穿模。
    auto segCrossFace=[&](const Vec3& A, const Vec3& B)->bool{
        for(const auto& f : blkFaces){
            if(f.curved) continue;      // 曲面/球面的拟合平面横穿自身，跨越它≠穿模；靠点云净空判
            const double da = dot(A - f.c, f.n), db = dot(B - f.c, f.n);
            if(da * db > 0.0) continue;                       // 同侧，不可能穿
            if(std::fabs(da - db) < 1e-12) continue;          // 与面平行且共面
            const double t = da / (da - db);
            if(t < -1e-9 || t > 1.0 + 1e-9) continue;
            const Vec3 X = A + (B - A)*t;
            const Vec3 w = X - f.c;
            const double a = dot(w, f.u), b = dot(w, f.v);
            if(a >= f.u0 - 1e-6 && a <= f.u1 + 1e-6 &&
               b >= f.v0 - 1e-6 && b <= f.v1 + 1e-6 &&
               f.solidAt(a, b)) return true;     // 穿过面上空洞不算穿模
        }
        return false;
    };
    // 点到各面**边界线**的最小净空。
    // 平面判据管"枪压在面上"，管不了"枪贴着棱边擦过去"——开口盒的壁顶就是一条
    // 棱边，去噪往往正好把这条棱边上的点删光（边界点邻域只有半边支撑），于是
    // 枪位离壁顶只有 16mm 也测出"离得很远"，顶行外扩会把枪甩出板外再横扫回来。
    // 面的边界用面内 2D 包围盒的 4 条边代表（矩形壁面是精确的；异形面是保守过
    // 近似，只会让外扩更早收敛，不会放宽）。棱边距离不依赖点云是否采到棱边。
    auto ptSegD = [](const Vec3& p, const Vec3& a, const Vec3& b)->double{
        const Vec3 ab = b - a;
        const double L2 = dot(ab, ab);
        if(L2 < 1e-12) return len(p - a);
        double t = dot(p - a, ab)/L2;
        t = std::max(0.0, std::min(1.0, t));
        return len(p - (a + ab*t));
    };
    auto edgeClearAt = [&](const Vec3& g, int skipSi)->double{
        double best = 1e30;
        for(const auto& f : blkFaces){
            if(f.si == skipSi) continue;
            const Vec3 C[4] = { f.c + f.u*f.u0 + f.v*f.v0,
                                f.c + f.u*f.u1 + f.v*f.v0,
                                f.c + f.u*f.u1 + f.v*f.v1,
                                f.c + f.u*f.u0 + f.v*f.v1 };
            for(int e=0;e<4;++e){
                const double d = ptSegD(g, C[e], C[(e+1)&3]);
                if(d < best) best = d;
            }
        }
        return best;
    };
    // 某点的枪位（沿该面法向抬 standoff）到工件的净空：点云、面、棱边三者取小
    auto gunClearAt = [&](const Vec3& surfPt, const Vec3& n, int skipSi = -1)->double{
        Vec3 g = surfPt + n*params->standoff;
        // 只用于和 gunClearEff(+1) 比较，半径上界取 gunClearEff+2 足够；
        // 超出上界的点直接返回上界，不改变任何 "< 阈值" 判定。
        const unsigned long long sk = (skipSi >= 0 && skipSi < (int)faceSkipMask.size())
                                    ? faceSkipMask[(size_t)skipSi] : 0ULL;
        double d = partTree.nearestWithinExcept(g, gunClearEff + 2.0, sk);
        d = std::min(d, planeClearAt(g, skipSi));
        d = std::min(d, edgeClearAt(g, skipSi));
        return d;
    };
    auto segClearOK = [&](const Vec3& a, const Vec3& b)->bool{
        const double d = len(b-a);
        if(d < 1e-6) return true;
        const int nsep = std::max(2, (int)std::ceil(d/std::max(1.0, pstep)));   // 步距 = 点距（原 2×点距，细长凸台会被跨过）
        for(int i=0;i<=nsep;++i){
            Vec3 q = a + (b-a)*((double)i/(double)nsep);
            if(partTree.nearestWithin(q, moveClear) < moveClear) return false;
        }
        return true;
    };
    // 段到工件的“最小净空”打分：只用于出口/落点改道搜索（正常路段用 segClearOK 快速判）。
    // 返回 0 表示该段穿过某个面（几何穿模），搜索时优先级最低。
    auto segMinClear = [&](const Vec3& a, const Vec3& b)->double{
        const double d = len(b-a);
        const int ns = std::max(2, (int)std::ceil(d/std::max(1.0, pstep)));
        double mn = 1e30;
        for(int k=0;k<=ns;++k){
            Vec3 q = a + (b-a)*((double)k/(double)ns);
            KnnBuf kb; partTree.knn(q, 1, kb);
            if(kb.n < 1) return 0.0;
            mn = std::min(mn, std::sqrt(kb.d2[0]));
            if(segCrossFace(a, q)) return 0.0;
        }
        return mn;
    };
    // 按**指定**净空阈值的擦件判定。与 segClearOK 的区别只在阈值：
    // segClearOK 用的是 moveClear = 0.25×standoff（空移的安全余量，standoff 150 时
    // 是 37.5mm），而换面短接要问的是"枪体真的会不会碰到"——那个阈值是 gunClearEff
    // （枪帽半径+余量，默认 20mm）。用 moveClear 判短换面等于凭空多要 17.5mm 余量，
    // p037 实测 22.36mm 明明不碰枪，却因 <37.5 被判"擦件"，用户看到的就是
    // "长度阈值填 800 还是关枪"（2026-11-20 报障）。
    auto segClearGe = [&](const Vec3& a, const Vec3& b, double thr)->bool{
        const double d = len(b-a);
        if(d < 1e-6) return true;
        const int nsep = std::max(2, (int)std::ceil(d/std::max(1.0, pstep)));
        for(int i=0;i<=nsep;++i){
            Vec3 q = a + (b-a)*((double)i/(double)nsep);
            if(partTree.nearestWithin(q, thr) < thr) return false;
        }
        return true;
    };

    // 原始点云、解析面和边界共同复查实际枪位。半步余量覆盖采样点之间的净空。
    auto candidateSegmentClear = [&](const Vec3& a, const Vec3& b, double need)->bool{
        if(segCrossFace(a, b)) return false;
        const double length = len(b-a);
        const int steps = std::max(1, (int)std::ceil(length/std::min(pstep, 2.0)));
        const double margin = need + 0.5*length/steps;
        for(int j=0; j<=steps; ++j){
            const Vec3 q = a + (b-a)*((double)j/steps);
            if(candidateTree.nearestWithin(q, margin) < margin ||
               planeClearAt(q, -1) < margin || edgeClearAt(q, -1) < margin) return false;
        }
        return true;
    };

    // ================= 换面短接：不够长就不停枪 =================
    // 现场口径（2026-11-19）：换面时只要枪**不用位移出来**（不退枪、不抬到安全面），
    // 就该直接连过去继续喷 —— 停一次枪要减速、断喷、再起喷，既费节拍，接缝还容易
    // 补枪过量。只有换面空移**太长**（或必须绕行）才关枪。
    //
    // 这一段只负责"短接"这一种情形；不满足时调用方仍走原来的关枪 + 抬枪绕行。
    // 输出点的 kind=1：查看器据此画橙色（= 换面过渡，但枪没关），和红色布行分开。
    // 诊断：为什么这段换面没有按"不关枪"处理。阈值填了却不生效时，光看"没短接"
    // 是查不出来的 —— 四条判据（长度/夹角/净空/穿面）任意一条拒掉都是同一个结果。
    // SPRAYTRAJ_TRANS_DBG=1 时把被拒理由和当时的实测值打出来。
    const bool xfDbg = std::getenv("SPRAYTRAJ_TRANS_DBG") != nullptr;
    auto tryShortFaceLink=[&](const Vec3& pEnd, const Vec3& nEnd,
                              const Vec3& pStart, const Vec3& nStart,
                              std::vector<SprayTrajPoint>& out)->bool{
        const double d = len(pStart - pEnd);
        const double nd = dot(nEnd, nStart);
        const double angDeg = std::acos(std::max(-1.0,std::min(1.0,nd)))*180.0/3.141592653589793;
        auto rej=[&](const char* why, double meas = -1.0, const char* unit = "")->bool{
            if(xfDbg){
                fprintf(stderr, "[xface] 拒绝: %-14s 距离 %.1f / 阈值 %.1f，两端夹角 %.1f°",
                        why, d, params->crossFaceGunOnLen, angDeg);
                if(meas >= 0.0) fprintf(stderr, "，实测 %s %.2f%s", why, meas, unit);
                fprintf(stderr, "\n");
            }
            return false;
        };
        if(params->crossFaceGunOnLen <= 0.0) return rej("阈值=0（关闭）");
        if(d < 1e-6) return true;                       // 两点重合：本来就不用连
        if(d > params->crossFaceGunOnLen) return rej("超过阈值"); // 太长 → 交给调用方关枪
        // 姿态差太大时插值法向会过零（90° 时模长只剩 0.707×standoff），枪位会陷进
        // 工件；真机上也是"一边喷一边猛拧腕"。上限由 crossFaceMaxTurnDeg 给（默认 60°），
        // 想放行折角件就调大 —— 不给用户一个能看的旋钮，他只会看到"阈值填了没用"。
        if(params->crossFaceMaxTurnDeg > 0.0){
            const double lim = std::cos(params->crossFaceMaxTurnDeg*3.141592653589793/180.0);
            if(nd < lim) return rej("两端夹角超上限", angDeg, "°");
        }
        // 擦件判据要报**实测最小净空**，不能只说"擦件"：现场拿这个数决定
        // 是"调小 gunClear 让它过"还是"这段必须关枪"（凹折角件枪位本来就要
        // 从棱角外侧绕，直线过棱一定贴件，那种情况怎么调都没用）。
        // 阈值用 gunClearEff 而不是 moveClear：这里判的是"枪体会不会碰到"，
        // 不是"空移够不够安全"。虚增的余量会让明明是近距换面的段被判成擦件。
        if(!segClearGe(pEnd, pStart, gunClearEff)){
            const double mc = segMinClear(pEnd, pStart);
            if(mc < 2.0) return rej("直线穿件", mc, "mm");
            return rej("直线擦枪", mc, "mm");
        }
        if(segCrossFace(pEnd, pStart)) return rej("直线穿面");
        const int np = std::max(2, (int)std::ceil(d/std::max(1.0, pstep)));
        // 球面插值而不是"线性插值再归一化"：后者在夹角大时角速度前慢后快，
        // 82° 那一段中途会突然抢转，腕部加速度比匀速大出一截。夹角越小两者越接近，
        // 但这里既然专门为"大夹角也要能过"开了口子，就必须按球面走。
        // 硬上限 150°：再往上就接近"翻到反面去喷"，slerp 两端权重同时趋零，
        // 归一化会放大浮点噪声（数值上等于除零），必须拒掉。
        const double th = std::acos(std::max(-1.0, std::min(1.0, nd)));
        if(th > 150.0*3.141592653589793/180.0) return rej("两端夹角>150°（形同翻面）");
        const double sth = std::sin(th);
        for(int t=1;t<np;++t){
            const double f = (double)t/(double)np;
            Vec3 nrm;
            if(sth < 1e-6) nrm = norm(nEnd*(1.0-f) + nStart*f);   // 夹角≈0 退化
            else{
                const double wa = std::sin((1.0-f)*th)/sth, wb = std::sin(f*th)/sth;
                nrm = norm(nEnd*wa + nStart*wb);
            }
            const Vec3 q   = pEnd + (pStart - pEnd)*f;
            out.push_back({q.x,q.y,q.z, -nrm.x,-nrm.y,-nrm.z, params->speed, 1, 1});
        }
        return true;   // 计数由调用方做：buildTrans 在"代价预估"里也会被调用，不能在这儿加
    };

    // ================= 单带件：整件横向跨度 <= 一个行距 → 全件只布一条行 =================
    // 现象（2026-11-19 现场 p004，531×110 的折带）：分割阶段按法向台阶把这块料切成
    // 6 片窄面（宽 49.1 / 36.0 / 30.0 / 22.0 / 19.8 / 36.2mm），布行阶段再**逐片独立**
    // 决策——每片都因为"自己只有 20~50mm 宽、比行距窄"而只布 1 行，于是同一块材料上
    // 落下 3 条互相平行、只差半个喷幅的行。用户看到的就是"一条料被来回喷了三遍"。
    //
    // 病根不在布行公式（每片的 nRow=1 单独看都是对的），而在**没有任何环节问过
    // "这些碎片加起来有多宽"**。分割是逐块生长的，它只保证每块内部平面性，不保证
    // 块与块之间是不同零件——折带被法向台阶切碎，正是这种情形。
    //
    // 判据：把"有望可喷"的面合并，在合并面上取最小外接矩形；短边 <= 行距 ⇒
    // 一条行的喷幅（按 20~30% 搭接反推约 1.3~1.4 倍行距）必然盖住整件，
    // 无论上游切出几片，都只该有一条行。
    //
    // 落点：选法向最贴近合并面平均法向的那片当"载体"。行仍落在真实面上、按真实法向
    // 抬 standoff（枪位与朝向不被这个规则改动），只是把同带其余面的点并入它的 2D
    // 支撑集，让这条行落在整件的横向中心、长度覆盖全长；其余面不再布行。
    //
    // 安全阀：载体若被后面的枪位/朝向/平面性门拒掉，carrierEmitted 保持 false，
    // 各面按原逻辑各自布行——任何情况下都不会比改动前少喷。
    // 关掉：SPRAYTRAJ_SINGLEBAND=0；放宽/收紧倍数：SPRAYTRAJ_SINGLEBAND_K。
    bool singleBand = false;
    size_t singleBandCarrier = (size_t)-1;
    std::vector<size_t> bandFaces;
    {
        bool enable = true;
        if(const char* e = std::getenv("SPRAYTRAJ_SINGLEBAND")) enable = (e[0] != '0');
        double kBand = 1.0;
        if(const char* e = std::getenv("SPRAYTRAJ_SINGLEBAND_K")){ double v=atof(e); if(v>0.05) kBand=v; }
        // 带厚上限（相对带宽）。折带是"一张料弯出来的"，厚度 = 弯曲矢高，实测
        // 6.6/85.7=0.08、13.8/76.4=0.18；而实心件（complex 162×596×110）的"厚度"
        // 就是它自己的截面尺寸，实测 43.7/134.5=0.33。0.25 两边都有余量。
        // 这条是"折带"与"实心件"的分水岭——单看短边分不开（complex 134.5 也 < 150），
        // 单看法向差更分不开（实测所有件都在 70~88°）。
        double kThick = 0.25;
        if(const char* e = std::getenv("SPRAYTRAJ_SINGLEBAND_T")){ double v=atof(e); if(v>0.02) kThick=v; }
        if(enable){
            for(size_t i=0;i<segs.size();++i){
                if(!meshOk[i]) continue;
                if(segs[i].curved) continue;                        // 曲面块走自己的行追踪
                if(segs[i].fo.isInner == 1) continue;               // 内壁不参与合并
                if(segs[i].fo.angZ > SPRAYTRAJ_NORMAL_WIN_DEG) continue;
                if(segs[i].inl.size() < 60) continue;
                bandFaces.push_back(i);
            }
            if(bandFaces.size() >= 2){
                std::vector<int> all;
                for(size_t i : bandFaces)
                    all.insert(all.end(), segs[i].inl.begin(), segs[i].inl.end());
                Vec3 nU{0,0,1}, cU{0,0,0};
                double wShort = 1e30, maxAng = 0.0, wLong = 0.0;
                double thick95 = 0.0;
                if(all.size() >= 120 && ls_plane(Pf, all, nU, cU)){
                    nU = norm(nU);
                    const Vec3 ax = (std::fabs(nU.x) < 0.9) ? Vec3{1,0,0} : Vec3{0,1,0};
                    const Vec3 ae1 = norm(cross(ax, nU));
                    const Vec3 ae2 = cross(nU, ae1);
                    std::vector<Vec2> uv; uv.reserve(all.size());
                    for(int k : all) uv.push_back({dot(Pf[k]-cU, ae1), dot(Pf[k]-cU, ae2)});
                    std::vector<Vec2> rect;
                    if(minAreaRect(uv, rect) > 1e-9 && rect.size() == 4){
                        const double a0 = len2(rect[1]-rect[0]);
                        const double b0 = len2(rect[2]-rect[1]);
                        wShort = std::min(a0, b0); wLong = std::max(a0, b0);
                    }
                    // 最大法向差：只做诊断输出，**不进判据**。
                    // 注意它算的是"各面 与 混合平均面 nU"的夹角，**不是两面互夹**。
                    // p048 实测这个数 89.2°，而两面真实互夹只有 44.6° —— 差一倍，
                    // 曾据它得出"两面近正交所以适用单带件"的错误结论。
                    const double R2D = 180.0/3.141592653589793;
                    for(size_t i : bandFaces){
                        const double d = std::fabs(dot(norm(segs[i].pl.n), nU));
                        maxAng = std::max(maxAng, std::acos(std::min(1.0, d))*R2D);
                    }
                }

                // ---- 载体选择：以**点数最多（≈面积最大）的面**为准 ----
                // 历史实现选"法向最贴近合并面平均法向"的那片，理由是让枪轴居中。
                // 但折带的两片常常是近正交的，此时两片的 |dot(n, nU)| 都接近 0，
                // "最贴近"退化成比谁的噪声小：
                //   p048 → 选中的是面0（面积 5856mm²、1427 点、带宽 48mm 的边角小片），
                //          主平面面2（面积 27771mm²、13744 点）反而没被选。
                // 行于是按**小片**的几何铺出来、按小片的轮廓夹端点，
                // 现场看到的就是"只有边边那一小块有轨迹，主平面整片光着"
                // （用户 2026-11-23 贴图提问的原话）。
                // 口径与"侧面折边一律吸收"一致：小片是折边，主面才是要喷的那张面。
                // 所以按内点数取最大面；点数接近（<=5% 差）时才用原法向贴近度决胜。
                size_t bestFace = (size_t)-1; double bestDot = 0.0;
                {
                    size_t cand = (size_t)-1;
                    for(size_t i : bandFaces)
                        if(cand == (size_t)-1 || segs[i].inl.size() > segs[cand].inl.size())
                            cand = i;
                    if(cand != (size_t)-1){
                        const size_t nMax = segs[cand].inl.size();
                        bestFace = cand;
                        bestDot = std::fabs(dot(norm(segs[cand].pl.n), nU));   // 诊断用真值
                        for(size_t i : bandFaces){
                            if(segs[i].inl.size()*20 < nMax*19) continue;   // 明显更小的片不参与
                            const double d = std::fabs(dot(norm(segs[i].pl.n), nU));
                            if(d > bestDot + 1e-9 ||
                               (std::fabs(d-bestDot) <= 1e-9 && segs[i].inl.size() > segs[bestFace].inl.size())){
                                bestDot = d; bestFace = i;
                            }
                        }
                    }
                }

                // ---- 带厚：必须按**载体面自己的法向**量，不能按混合平均面 nU ----
                // 这是 p048 误判成单带件的根因（用户 2026-11-23 追问"两个面不是都可以喷吗"）。
                //
                // 折带的"厚度"= 弯曲矢高，前提是**整条带共用一张面**，所以沿那张面的
                // 法向投影才有意义。可一旦两个面是**前后串联**（p048：沿行方向差 186mm、
                // 两面互夹 44.6°），nU 就是一个谁都不代表的中庸方向 —— 两个面都在 86°
                // 附近、几乎垂直于它，沿 nU 投影当然"很薄"（18.6mm，刚好卡进 19.1mm 门限）。
                // 但行实际是沿**载体面法向**抬 standoff 的，那个方向上两面的高差是
                // p048 实测 **60.0mm** —— 远超 19.1mm，本来就不该走单带件。
                //
                // 一句话：带厚要和"行实际用的方向"同口径，否则量的是投影长度，不是厚度。
                // 回退旧口径：SPRAYTRAJ_THICK_BY_NU=1
                double thick95_nU = 0.0;   // 旧口径，留作诊断对照
                if(bestFace != (size_t)-1 && !all.empty()){
                    const Plane& cpl = segs[bestFace].pl;
                    const Vec3 cn = norm(cpl.n);
                    std::vector<double> hh; hh.reserve(all.size());
                    for(int k : all) hh.push_back(std::fabs(dot(Pf[k]-cpl.c, cn)));
                    const size_t q95 = (size_t)(0.95*(double)(hh.size()-1));
                    std::nth_element(hh.begin(), hh.begin()+(ptrdiff_t)q95, hh.end());
                    thick95 = hh[q95];
                    // 旧口径（按混合平均面 nU）始终算一份，只作诊断对照 ——
                    // 判据未过时要能一眼看出"是哪个口径把门卡住的"。
                    {
                        std::vector<double> h2; h2.reserve(all.size());
                        for(int k : all) h2.push_back(std::fabs(dot(Pf[k]-cU, nU)));
                        const size_t q2 = (size_t)(0.95*(double)(h2.size()-1));
                        std::nth_element(h2.begin(), h2.begin()+(ptrdiff_t)q2, h2.end());
                        thick95_nU = h2[q2];
                    }
                    if(std::getenv("SPRAYTRAJ_THICK_BY_NU")) thick95 = thick95_nU;
                }
                if(bestFace != (size_t)-1 && wShort <= kBand*L &&
                   thick95 <= kThick*std::max(1.0, wShort)){
                    singleBand = true;
                    singleBandCarrier = bestFace;
                    fprintf(stderr, "[traj] 单带件: 可喷面 %zu 片，合并最小外接矩形 %.1f x %.1fmm "
                                    "（短边<=行距 %.1fmm）、带厚p95=%.1fmm（<=%.2f×短边，按载体面%zu法向）、"
                                    "最大法向差%.1f°(各自对nU) → 全件只布一条行（载体=面%zu，法向贴合度 %.2f）\n",
                            bandFaces.size(), wLong, wShort, L, thick95, kThick, bestFace,
                            maxAng, bestFace, bestDot);
                } else if(std::getenv("SPRAYTRAJ_FACE_DBG") || std::getenv("SPRAYTRAJ_BAND_DBG")) {
                    fprintf(stderr, "[traj] 单带件判据未过: %zu 片，合并矩形 %.1f x %.1fmm"
                                    "（短边门限 %.1fmm）、带厚p95=%.1fmm（门限 %.1fmm，按载体面%zu法向）、"
                                    "旧口径(按nU)=%.1fmm → 各面照常布行\n",
                            bandFaces.size(), wLong, wShort, kBand*L,
                            thick95, kThick*std::max(1.0, wShort), bestFace, thick95_nU);
                }
            }
        }
    }

    // ---- 面顺序：生成阶段不动 ----
    // 逐面生成仍按 split_planes 的发现顺序（内点数 / 面积降序），面间过渡也照旧就地接上。
    // 真正的顺序优化放在**最后的重排阶段**（见"5.面顺序重排"）：到那时每个面的真实
    // 首末枪位已经躺在轨迹里了，用真实端点去算真实过渡代价才可靠 —— 而过渡的固定开销
    // （退枪 → 抬到安全面 → 平移 → 落枪）远大于两端点的直线距离，用面心/包围盒角点
    // 这类几何代理估，实测在 box_traj 上空移反而 +34%。
    // 例外：单带件的载体提到最前，先让它出结果——它过了门才谈得上"其余面不必布行"。
    std::vector<size_t> faceOrder(segs.size());
    for(size_t i=0;i<faceOrder.size();++i) faceOrder[i] = i;
    if(singleBand){
        for(size_t i=0;i<faceOrder.size();++i)
            if(faceOrder[i] == singleBandCarrier){ std::swap(faceOrder[0], faceOrder[i]); break; }
    }
    bool carrierEmitted = false;

    // ---- 喷幅覆盖：枪在某个行点开火时，喷锥照到的面上不再单独布行 ----
    //
    // 口径纠正（用户 2026-11-20 16:19）：
    //   「喷幅应该跟行距不是一回事吧，行距设多少，喷幅只跟出粉量有关，
    //     还有关的是里面距离，毕竟远了喷不到」
    //
    // 历史实现把喷幅从行距**反推**：R = L/(2(1−o))，行距 120/搭接 30% → R=85.7mm。
    // 物理上反了：喷幅是"枪 + 出粉量 + 雾化气压 + 喷距"定出来的，行距则是**为覆盖**
    // 而选出来的（L = W(1−o)）。行距调大不代表枪能喷更宽。
    //
    // 现在喷幅 W 由参数 fanWidth 独立给出；没给才退回旧的由行距反推口径
    // （保住既有回归逐位不变）。
    //
    // 覆盖判据必须**两个方向分开**——这是"远了喷不到"的落点：
    //   · 横向（垂直枪轴）：点到最近行点的垂距 <= R = W/2  → 喷锥截面上够得着
    //   · 轴向（沿枪轴）：  该行的喷距 standoff + Δ 落在 [0.5·standoff, maxSprayDist]
    //       Δ = 点沿**该行枪轴**方向的偏移（凸起为正、凹陷为负）。
    //     以前只看 3D 欧氏距离，一个高出主面 150mm 的凸台到行的距离算出来 85mm，
    //     被判"喷得到"——而它的真实喷距 250mm 早就超出有效喷距了。
    //   · 面级夹角 <= 60°：喷锥边缘 cos 衰减还能接受。
    bool fanCoverOn = true;
    if(const char* e=std::getenv("SPRAYTRAJ_FAN_COVER")) fanCoverOn = (e[0] != '0');
    double fanOverlap = 0.30;
    if(const char* e=std::getenv("SPRAYTRAJ_FAN_OVERLAP")){
        double v = atof(e); if(v >= 0.0 && v < 0.9) fanOverlap = v;
    }
    // 喷幅半径：优先用独立参数；未给则退回"由行距反推"（旧口径）。
    double fanR;
    const double fanWSrc = (params->fanWidth > 1.0) ? params->fanWidth
                                                    : params->scanRowSpacing/std::max(0.2, 1.0-fanOverlap);
    fanR = 0.5*fanWSrc;
    if(const char* e=std::getenv("SPRAYTRAJ_FAN_R")){
        double v = atof(e); if(v > 1.0) fanR = v;
    }
    // 有效喷距：<=0 表示不判轴向（旧口径）。上限同时给出"太近"的下限 0.5·standoff，
    // 因为枪贴到面上会把漆吹散、冲稀，膜厚同样不达标。
    const double maxSprayD = (params->maxSprayDist > 1.0) ? params->maxSprayDist : 0.0;
    const double minSprayD = 0.5*params->standoff;
    std::vector<Vec3> fanP, fanN;                       // 已布行面的"面点"与所属面法向
    // 网格边长：要能容纳"横向很近、轴向很远"的候选点。启用有效喷距后，一个点可能
    // 横向只差 5mm 而轴向差 150mm，用 fanR 当格距时 3×3 邻域根本搜不到它。
    // 取 max(fanR, maxSprayD) 后，3×3 邻域一定覆盖判定所需的全部候选。
    const double fanCellSz = std::max(fanR, maxSprayD > 0.0 ? maxSprayD : fanR);
    std::map<std::array<int,3>, std::vector<int>> fanGrid;
    auto fanCell = [&](const Vec3& p)->std::array<int,3>{
        return { (int)std::floor(p.x/fanCellSz), (int)std::floor(p.y/fanCellSz), (int)std::floor(p.z/fanCellSz) };
    };
    auto fanAdd = [&](const Vec3& p, const Vec3& n){
        const int idx = (int)fanP.size();
        fanP.push_back(p); fanN.push_back(n);
        fanGrid[fanCell(p)].push_back(idx);
    };
    // 找最近行点，并把距离**拆成横向/轴向两个分量**（沿最近行点的枪轴 fanN[k]）：
    //   横向 h = |Δ − (Δ·n)n| ：能不能被喷锥截面扫到，与 fanR 比
    //   轴向 s = Δ·n           ：点相对"面点"沿枪轴抬高多少（凸起为正、凹陷为负）
    //     枪位在 fanP + n·standoff，所以真实喷距 = standoff − s：
    //       凸起 s>0 → 点更靠近枪口 → 喷距变短（太近会吹散漆）
    //       凹陷 s<0 → 点更远离枪口 → 喷距变长（太远喷不到）
    //     注意 s 的正负号别搞反：写成 standoff + s 会把凸起当"更远"，判据整个反过来。
    // 只在"先按 3D 距离筛出候选"时用欧氏距离做粗筛，精判一律用分量。
    auto fanNear = [&](const Vec3& p)->std::pair<double,int>{
        double bd = 1e18; int bi = -1;
        const auto c = fanCell(p);
        for(int dx=-1;dx<=1;++dx) for(int dy=-1;dy<=1;++dy) for(int dz=-1;dz<=1;++dz){
            auto it = fanGrid.find({c[0]+dx, c[1]+dy, c[2]+dz});
            if(it == fanGrid.end()) continue;
            for(int k : it->second){
                const double d = len(p - fanP[k]);
                if(d < bd){ bd = d; bi = k; }
            }
        }
        return {bd, bi};
    };
    // 把"到最近行点的偏移"拆成横向 / 轴向。bi<0 时返回 (1e18, 0)。
    auto fanSplit = [&](const Vec3& p, int bi, double& h, double& s){
        if(bi < 0){ h = 1e18; s = 0.0; return; }
        const Vec3 dv = p - fanP[(size_t)bi];
        const Vec3& n = fanN[(size_t)bi];
        s = dot(dv, n);
        h = len(dv - n*s);
    };
    // 把本面已写入的**喷涂点**反算回面点收进覆盖集（枪位 = 面点 + n·standoff）。
    auto fanHarvest = [&](size_t beg, size_t end, const Vec3& n){
        for(size_t q=beg; q<end; ++q){
            if(allPts[q].spray == 0) continue;
            fanAdd(Vec3{allPts[q].x, allPts[q].y, allPts[q].z} - n*params->standoff, n);
        }
    };
    size_t maxFacePts = 0;
    for(const auto& s : segs) maxFacePts = std::max(maxFacePts, s.inl.size());
    if(fanCoverOn){
        // 口径自检：喷幅是被行距反推出来的（旧），还是独立给的（新）
        const bool fanIndep = (params->fanWidth > 1.0);
        fprintf(stderr, "[fan] 喷幅覆盖判据: R=%.1fmm（喷幅%.1fmm %s）", fanR, 2.0*fanR,
                fanIndep ? "独立给定" : "※由行距反推，建议用 -fanw 直接给喷幅");
        if(maxSprayD > 0.0)
            fprintf(stderr, " 有效喷距[%.1f, %.1f]mm（standoff=%.1f）", minSprayD, maxSprayD, params->standoff);
        else
            fprintf(stderr, " 有效喷距未启用（只判横向，远了也可能判成喷得到——建议填 maxSprayDist）");
        fprintf(stderr, "\n");
        if(!fanIndep && params->scanRowSpacing > 40.0)
            fprintf(stderr, "[fan] ⚠ 行距 %.1fmm 反推出喷幅 %.1fmm，等同于假设一把枪能喷这么宽。"
                            "喷幅应由出粉量/雾化/喷距决定，请用 -fanw 显式给出。\n",
                    params->scanRowSpacing, 2.0*fanR);
    }

    // 布行层代表点掩膜：P2 只吃这些点（planIdx 为空时全 1，行为逐位不变）。
    // 注意只影响"布行的原子"，不改 seg.inl —— 面点数/覆盖率/统计仍按全分辨率口径，
    // 否则日志里的"内点=NNNN"会随体素档位漂移，反而不好追。
    std::vector<unsigned char> isPlan(Pf.size(), 1);
    if(!planIdx.empty()){
        std::fill(isPlan.begin(), isPlan.end(), (unsigned char)0);
        for(int k : planIdx) isPlan[(size_t)k] = 1;
    }

    for(size_t oi=0; oi<faceOrder.size(); ++oi){
        const size_t si = faceOrder[oi];
        if(!meshOk[si]) continue;          // 轮廓都重建不出来的面，不喷
        // 单带件：载体已经布出一条行，其余碎片面的材料都落在这一行的喷幅带内，
        // 再各布一条就是对同一处反复喷（膜厚翻倍、行程翻倍）。
        if(singleBand && carrierEmitted && si != singleBandCarrier){
            fprintf(stderr, "[traj] 面%zu 属单带件（短边<=行距）且载体面%zu 已布行，跳过不布行\n",
                    si, singleBandCarrier);
            continue;
        }
        auto& seg = segs[si];
        Plane& pl = seg.pl;
        const auto& inl = seg.inl;
        if(inl.size() < 20) continue;

        // 喷幅覆盖检查：前面（更大）的面已经把行走到它上方/旁边时，这里不必再布一遍。
        if(fanCoverOn && !fanP.empty() && si != faceOrder[0] &&
           inl.size()*2 < maxFacePts && inl.size() >= 8){
            std::vector<double> hh, dd;          // hh=横向垂距  dd=3D 欧氏距离
            hh.reserve(inl.size()); dd.reserve(inl.size());
            std::vector<double> angs; angs.reserve(inl.size());
            std::vector<unsigned char> reach(inl.size(), 0);   // 该点是否被喷锥真正覆盖
            double sprayDmin = 1e18, sprayDmax = -1e18;        // 真实喷距范围（诊断用）
            for(size_t ii=0; ii<inl.size(); ++ii){
                const int k = inl[ii];
                const auto pr = fanNear(Pf[k]);          // 最近行点（3D 距离 + 下标）
                double h = 1e18, sn = 0.0;
                fanSplit(Pf[k], pr.second, h, sn);
                hh.push_back(pr.second >= 0 ? h : fanR*4.0);   // 没搜到 = 肯定够不着，记个哨兵值便于统计
                dd.push_back(pr.second >= 0 ? pr.first : fanR*4.0);
                // 真实喷距：枪位在面点 + n·standoff，点在面点 + n·s ⇒ 喷距 = standoff − s。
                // 凸起 s>0 → 变近（太近会吹散漆）；凹陷 s<0 → 变远（太远喷不到）。
                const double sprayD = params->standoff - sn;
                if(pr.second >= 0){
                    sprayDmin = std::min(sprayDmin, sprayD);
                    sprayDmax = std::max(sprayDmax, sprayD);
                }
                bool ok;
                if(maxSprayD > 0.0){
                    // 新口径：横向、轴向分开判（用户 2026-11-20 16:19）
                    ok = (h <= fanR) && (sprayD >= minSprayD) && (sprayD <= maxSprayD);
                } else {
                    // 旧口径：3D 欧氏距离 <= R。**必须保留**，否则"横向近、轴向远"的
                    // 凸台会被判成喷得到，等于把漏喷当成省行。
                    ok = (pr.first <= fanR);
                }
                reach[ii] = ok ? 1 : 0;
                if(pr.second >= 0){
                    const double c = dot(seg.pl.n, fanN[(size_t)pr.second]);
                    angs.push_back(std::acos(std::max(-1.0, std::min(1.0, c)))*180.0/3.141592653589793);
                }
            }
            std::sort(hh.begin(), hh.end());
            std::sort(dd.begin(), dd.end());
            std::sort(angs.begin(), angs.end());
            const double h50  = hh[hh.size()/2];
            const double h95  = hh[(size_t)(0.95*(double)(hh.size()-1))];
            const double d50  = dd[dd.size()/2];
            const double d95  = dd[(size_t)(0.95*(double)(dd.size()-1))];
            const double a50 = angs.empty() ? 180.0 : angs[angs.size()/2];
            size_t cov = 0; for(unsigned char rch : reach) if(rch) ++cov;
            const double covFrac = (double)cov/(double)reach.size();
            // 口径决定用哪个"够得着"的量：
            //   新口径 → 横向垂距 h95 <= R（再叠轴向喷距，已在 reach 里判过）
            //   旧口径 → 3D 欧氏距离 d95 <= R（与本参数引入前逐位一致）
            const bool near95OK = (maxSprayD > 0.0) ? (h95 <= fanR) : (d95 <= fanR);
            if(covFrac >= 0.90 && near95OK && a50 <= 60.0){
                fprintf(stderr, "[fan] 面%zu(%zu点) 已被前面布行的喷锥覆盖：R=%.1fmm "
                                "横向p50=%.1f p95=%.1f 3D p95=%.1f 喷距%.1f~%.1f 覆盖=%.0f%% 夹角中位=%.1f° → 跳过布行\n",
                        si, inl.size(), fanR, h50, h95, d95, sprayDmin, sprayDmax, covFrac*100.0, a50);
                continue;
            }
            if(std::getenv("SPRAYTRAJ_FACE_DBG"))
                fprintf(stderr, "[fan] 面%zu(%zu点) 未通过喷锥覆盖：横向p50=%.1f p95=%.1f "
                                "3D p50=%.1f p95=%.1f(R=%.1f) 喷距%.1f~%.1f%s 覆盖=%.0f%% 夹角中位=%.1f° → 照常布行\n",
                        si, inl.size(), h50, h95, d50, d95, fanR, sprayDmin, sprayDmax,
                        maxSprayD > 0.0 ? "" : "(未判轴向)", covFrac*100.0, a50);
        }
        if(std::getenv("SPRAYTRAJ_PROF"))
            fprintf(stderr, "[prof] ==== 面%zu 进入: 内点=%zu curved=%d isCurved=%d ====\n",
                    si, inl.size(), (int)seg.curved, (int)isCurved);

        // 面级方向和枪位快检。可信平面可用后续实际行复核，角度窗口仍必过。
        // 位置很关键：必须在**所有**定向/劈面/拆层/内壁修正之后，在布行之前。
        // 上游（区域生长 → 初始定向 → 射线定向 → 拆层 → 内壁翻朝腔）每一步都可能把
        // 符号弄错，而它们的错误只有两种下场：漏喷（可接受）或从背面伸枪（撞机，不可接受）。
        // 这里不再相信任何标签，只问"枪能不能站到这一侧"：
        //   站得住 → 按当前法向布行；
        //   站不住但反侧站得住（且没超窗口）→ 翻过来；
        //   快检未过 → 可信平面进入候选行复核，其余面拒喷。
        bool needsCandidateCheck = seg.preservedCarrier;
        {
            char tag[48];
            std::snprintf(tag, sizeof(tag), "面%zu", si);
            const bool isPlanarFace = !(seg.curved && isCurved);
            const bool allowFlip = isPlanarFace && !seg.preservedCarrier;
            // 曲面块用逐点法向测自由距离：整块平均法向在强曲面上会把射线擦回自身，
            // 见 measure_gun_free 注释（曲型板曾因此整块拒喷）。平面块仍传 nullptr。
            const std::vector<Vec3>* ln = isPlanarFace ? nullptr : &cn.nor;
            ProfStage psGate("4b0.枪位验收");
            const GunSideVerdict gv = verify_gun_side(Pf, cn.tree, seg, partC,
                                                      params->standoff, tag, allowFlip, ln);
            if(gv == GUN_SIDE_REJECT && enableCandidateCheck && isPlanarFace &&
               (!seg.fitSupport.empty() || seg.collapsed) && seg.fo.angZ <= SPRAYTRAJ_NORMAL_WIN_DEG){
                needsCandidateCheck = true;
                ++nGunReject;
                fprintf(stderr, "[gun] 面%zu 面级枪位验收未通过（+n=%.1f / -n=%.1fmm）→ 复核候选行\n",
                        si, seg.fo.freePlus, seg.fo.freeMinus);
            } else if(gv == GUN_SIDE_REJECT){
                ++nGunReject;
                fprintf(stderr, "[skip] 面%zu 枪位验收未通过（按%s测得 +n=%.1fmm / -n=%.1fmm，负值表示量程内无命中），"
                                "不布行（内点数 %zu）\n", si, ln ? "逐点法向" : "面法向",
                        seg.fo.freePlus, seg.fo.freeMinus, inl.size());
                continue;
            }
            if(gv == GUN_SIDE_FLIP){
                ++nGunFlip;
                fprintf(stderr, "[fix] 面%zu 枪位验收翻转法向：angZ 现 %.1f°（原向一侧枪伸不进去）\n",
                        si, seg.fo.angZ);
            }
        }

        // ---- 可喷角度窗口：法向朝下的面不生成轨迹 ----
        // 与面级清单用同一条判据（SPRAYTRAJ_NORMAL_WIN_DEG = 130°）。以前窗口只在
        // driver 打印清单时执行，轨迹照喷，于是开口盒的"外底面"（法向 180°、朝下、
        // 贴地面）也有轨迹——清单里已经滤掉，这里必须同口径滤掉。
        // 注意用 fo.angZ（真实可喷侧法向），不是 pl.n：内底朝上（0°）在这条之外，照喷。
        // SPRAYTRAJ_NORMAL_WIN 是**外壁**的朝上窗口，默认 130° = 外壁也喷（用户 2026-10-30
        // 确认：整个工件都要涂层，侧面漏不得；加倾斜角后侧壁枪轴与竖直是 63° 左右，不是 90°，
        // 因为安全层会绕表面点把枪轴倾斜上去，喷距不变）。
        // 想只喷接近朝上的面再调小（如 60°）；**内壁不受这个窗口管**——箱体内壁法向本来就是
        // 水平，必须水平喷，一刀切到 60° 会把开口盒从 118 行砍到 22 行。
        double normalWin = SPRAYTRAJ_NORMAL_WIN_DEG;
        if(const char* e = std::getenv("SPRAYTRAJ_NORMAL_WIN")){
            double v = atof(e);
            if(v > 1.0 && v <= 180.0) normalWin = v;
        }
        const bool isInnerWall = (seg.fo.isInner == 1);
        const double effWin = isInnerWall ? SPRAYTRAJ_NORMAL_WIN_DEG : normalWin;
        if(seg.fo.angZ > effWin){
            ++nDownSkip;
            fprintf(stderr, "[skip] 面%zu %s朝向 %.1f° > %.0f°，不喷（内点数 %zu）\n",
                    si, isInnerWall ? "内壁" : "外壁", seg.fo.angZ, effWin, inl.size());
            continue;
        }

        // ---- 面级平面性门：拟合出来的"平面"必须真的是平面 ----
        // 区域生长在薄壁件上会把整圈壁吞成一个假面：实测 4mm 薄壁开口盒
        // （g_norim）生成了一个 193956 点、法向 (0,0,1)、质心 z=115.8 的"水平面"，
        // 而点云 z 实际从 0 铺到 300 —— 沿法向 RMS 约 85mm。在假面上布行，
        // 轨迹是一条横穿所有壁的水平线，实测枪位距件 0.06mm（等于贴件走过）。
        // 平面性门按 RMS 卡：正常扫描平面 1~3mm，超限的要么是真曲面（该走曲面块）、
        // 要么是误合并（根本不该喷），两种情况都不能按平面布行。
        {
            double rss=0.0, emax=0.0;
            for(int k : inl){
                const double e = dot(Pf[k]-pl.c, pl.n);
                rss += e*e;
                emax = std::max(emax, std::fabs(e));
            }
            const double fullRms = std::sqrt(rss/(double)inl.size());
            double frms = fullRms;
            const bool useFitSupport = enableCandidateCheck && !seg.fitSupport.empty();
            if(useFitSupport){
                double fitRss = 0.0;
                for(int k : seg.fitSupport){
                    const double e = dot(Pf[k]-pl.c, pl.n);
                    fitRss += e*e;
                }
                frms = std::sqrt(fitRss/(double)seg.fitSupport.size());
            }
            const double rmsCap = std::max(3.0, 0.12*params->standoff);
            // 曲面块本来就不是平面（pl 只是平均平面），不能按平面性拒掉——
            // 它后面走 gen_surface_rows 在真实曲面上追行，尺度完全不同。
            const bool isSurfBlock = seg.curved && isCurved;
            const bool skipPlanarity = isSurfBlock || seg.collapsed;
            if(!skipPlanarity && frms > rmsCap){
                ++nPlanarSkip;
                fprintf(stderr, "[skip] 面%zu 平面性不足 RMS=%.1fmm > %.1fmm（%zu 点），"
                                "疑为误合并/强曲面，不按平面布行\n", si, frms, rmsCap, useFitSupport ? seg.fitSupport.size() : inl.size());
                continue;
            }
            // 本面点都在枪位下方 standoff 处。查询点本身也可能偏离拟合平面
            // （真实内点最多 emax），故最近距离下界是 standoff - 2*emax。
            // 取 emax < 0.2*standoff，保证下界 >= 60mm，远大于 20mm 安全半径；
            // 只有强起伏/离群点才退回全局查询。典型平面 emax≈5mm，仍走剪枝。
            if(enableCandidateCheck && !skipPlanarity && fullRms > rmsCap){
                needsCandidateCheck = true;
                fprintf(stderr, "[fit] 面%zu 载体 RMS=%.2fmm，覆盖点 RMS=%.2fmm → 复核候选行\n",
                        si, frms, fullRms);
            }
            if(!needsCandidateCheck && !isSurfBlock && segs.size() <= 64 && si < 64 &&
               emax < 0.2*params->standoff)
                faceSkipMask[si] = 1ULL << si;
        }

        enum class RayTarget { NoHit, OnTarget, Blocked };
        std::vector<unsigned char> target;
        if(needsCandidateCheck){
            if(candidateTree.nodes.empty()) candidateTree.init(P);
            target.resize(Pf.size(), 0);
            for(int k : inl) target[(size_t)k] = 1;
            for(int k : seg.core) target[(size_t)k] = 1;
            for(int k : seg.fitSupport) target[(size_t)k] = 1;
        }
        auto candidateHit = [&](const Vec3& gunPos, bool continuousFace = false)->RayTarget{
            constexpr double rayRadius = 1.5; // 与面级射线的命中半径一致。
            const double standoff = params->standoff;
            const double extra = std::max(0.0, std::min(params->flatnessTol, std::max(2.0, 0.6*standoff)));
            const double rayLength = maxSprayD > 0.0 ? maxSprayD : standoff + extra;
            const Vec3 dir = pl.n*(-1.0);
            const int hit = candidateTree.firstRayPoint(gunPos, dir, rayLength, rayRadius);
            if(hit < 0) return RayTarget::NoHit;
            KnnBuf near;
            cn.tree.knn(P[(size_t)hit], 1, near);
            const double depth = dot(P[(size_t)hit]-gunPos, dir);
            const bool belongs = near.n > 0 && target[(size_t)near.id[0]];
            if(!belongs || (maxSprayD > 0.0 && depth < minSprayD)) return RayTarget::Blocked;
            if(near.d2[0] <= rayRadius*rayRadius) return RayTarget::OnTarget;
            if(!continuousFace) return RayTarget::Blocked;

            // 射线半径限定横向命中，不能同时用来判定去噪删掉的棱边点属于哪张面。
            // 规则连续面用同一重建格内、深度相符的多个本面点恢复归属。
            const double neighborhood = 0.5*kReconGap;
            KnnBuf support;
            cn.tree.knn(P[(size_t)hit], 6, support);
            int matching = 0;
            for(int j=0; j<support.n; ++j){
                if(support.d2[j] > neighborhood*neighborhood) continue;
                const int id = support.id[j];
                if(!target[(size_t)id]) return RayTarget::Blocked;
                if(std::fabs(dot(Pf[(size_t)id]-P[(size_t)hit], dir)) <= 2.0*rayRadius)
                    ++matching;
            }
            return matching >= 3 ? RayTarget::OnTarget : RayTarget::Blocked;
        };

        // ================= 曲面块：曲面上的行追踪 =================
        // 不走平面的 2D 布行——曲面投影成 2D 会在陡峭段挤压，行距失控。
        // 行点/法向已由 gen_surface_rows 在曲面上直接算好，这里只做抬升与串接。
        if(seg.curved && isCurved){
            bool meshRows = false;
            if(surf.nRow <= 0){
                ProfStage ps("4b.曲面行追踪", &tTraj);
                bool openContours = false;
                const bool traced = gen_surface_rows(surf, L, pstep, params->fillHoleSize,
                                                     params->scanDirMode, params->speed, openContours);
                const bool hasFilledSurface = std::find(surfaceFilledTris.begin(),surfaceFilledTris.end(),1)
                                              != surfaceFilledTris.end();
                if(!traced || openContours || hasFilledSurface){
                    const Vec3 lift = surf.n0*params->standoff;
                    auto clearSpraySegment = [&](const Vec3& a, const Vec3& b){
                        const Vec3 ga = a + lift, gb = b + lift;
                        return partTree.nearestWithin(gb, gunClearEff) >= gunClearEff &&
                               segClearGe(ga, gb, gunClearEff) && !segCrossFace(ga,gb);
                    };
                    double score = 0.0;
                    meshRows = gen_surface_mesh_rows(surf, surfaceVerts, surfaceTris, surfaceFilledTris, L, pstep,
                                                     params->scanDirMode, turnR, fanR, clearSpraySegment, score);
                    if(params->scanDirMode == 0 && meshRows){
                        auto points = std::move(surf.rowP), normals = std::move(surf.rowN);
                        auto starts = std::move(surf.rowBeg);
                        auto turns = std::move(surf.rowTurns);
                        const double gap = surf.rowGap;
                        double crossScore = 0.0;
                        if(!gen_surface_mesh_rows(surf, surfaceVerts, surfaceTris, surfaceFilledTris, L, pstep,
                                                  2, turnR, fanR, clearSpraySegment, crossScore) ||
                           crossScore <= score){
                            surf.rowP = std::move(points); surf.rowN = std::move(normals);
                            surf.rowBeg = std::move(starts); surf.rowTurns = std::move(turns);
                            surf.nRow = (int)surf.rowBeg.size()-1;
                            surf.rowGap = gap;
                        }
                    }else if(params->scanDirMode == 0 && !meshRows){
                        meshRows = gen_surface_mesh_rows(surf, surfaceVerts, surfaceTris, surfaceFilledTris, L, pstep,
                                                         2, turnR, fanR, clearSpraySegment, score);
                    }
                    if(meshRows){
                        measure_surface_row_coverage(surf, L);
                        fprintf(stderr, "[曲面] 蛇形截线 %d 段，按网格连续喷涂\n", surf.nRow);
                    }
                }
                if((!traced || openContours || hasFilledSurface) && !meshRows){
                    fprintf(stderr, "[曲面] 未能生成满足净空的连续轨迹，跳过该块\n");
                    continue;
                }
            }
            const int nr = surf.nRow;
            if(nr < 1 || surf.rowBeg.size() < (size_t)nr+1) continue;

            auto emitPt=[&](const Vec3& X, const Vec3& nrm, int spray, int kind=0){
                Vec3 pos = X + nrm*params->standoff;
                allPts.push_back({pos.x,pos.y,pos.z, -nrm.x,-nrm.y,-nrm.z, params->speed, spray, kind});
            };
            // 行间连接：
            //   · 正常相邻行（枪位只挪一个行距、法向一致、直线净空足够）
            //     → 按点距插值，短距保持喷涂不关枪；
            //   · 其余情况（跨到远端的卷边行 / 缺口两侧）→ 关枪，沿法向退出、
            //     抬到安全高度绕行。原来一律线性插值，曲型板实测有一段 648mm
            //     的跨行直线从板中穿过（折线最小净空 1.7mm）；两端法向相反时
            //     插值法向还会在中间过零，枪位直接落到工件表面（0.09mm）。
            auto emitGap=[&](const Vec3& a, const Vec3& na, const Vec3& b, const Vec3& nb){
                double d = len(b-a);
                if(!meshRows && d <= pstep*0.5) return;
                const Vec3 ga = a + na*params->standoff;
                const Vec3 gb = b + nb*params->standoff;
                // 关枪判据必须看**枪位行程**，不是面点跳距。
                // 枪位 = 面点 + standoff·n，两行端点法向不同，法向旋转会把面点的大跳距
                // 在枪位上大幅抵消：垃圾桶 L=12 时面点跳 93mm、枪位只走 13mm。
                // 用面点跳距判会白关 8 次枪（每次都要停-关-开-加速，很费时间），
                // 而实际喷嘴只挪了十几毫米，根本不需要断。
                const double dg = len(gb - ga);
                const double nd = dot(na, nb);
                if(meshRows && dg > 1e-6 && nd > 0.5 && segClearOK(ga,gb) && !segCrossFace(ga,gb)){
                    const Vec3 midN = norm(na+nb);
                    allPts.push_back({(ga.x+gb.x)*0.5, (ga.y+gb.y)*0.5, (ga.z+gb.z)*0.5,
                                      -midN.x,-midN.y,-midN.z,params->speed,0});
                } else if(!meshRows && dg <= 3.0*L && nd > 0.5 && segClearOK(ga,gb) && !segCrossFace(ga,gb)){
                    int np = std::max(2, (int)std::ceil(d/pstep));
                    for(int t=1;t<np;++t){
                        double f = (double)t/(double)np;
                        // 法向必须归一化：两端法向夹角大时线性插值的模长 <1，
                        // 枪位会向工件下陷（90° 夹角时只剩 0.707×standoff）。
                        Vec3 n = norm(na*(1.0-f) + nb*f);
                        // 换行过渡：kind=0。换行是同一张面内行与行之间的正常连续喷涂，
                        // 枪没关、点也在覆盖范围内，**不该标色**（用户 2026-11-20 明确：
                        // 「换面才需要换颜色，换行不需要」）。只有换面才有"关枪/不关枪"
                        // 这个需要区分的语义，见 kind=1。
                        emitPt(a + (b-a)*f, n, 1, 0);
                    }
                } else {
                    // 安全空移：沿法向退出 → 抬到工件最高点之上 → 目标上方 →
                    // （垂直下降不干净时沿目标法向退 clear 再进）→ 目标枪位。
                    ++nGapSafe;
                    auto move=[&](const Vec3& p, const Vec3& nrm){
                        allPts.push_back({p.x,p.y,p.z, -nrm.x,-nrm.y,-nrm.z, params->speed, 0});
                    };
                    // 出口两腿同样要查：腔宽≈clear 时，沿法向退到对面内壁再竖
                    // 直抬枪会贴着壁顶擦过（bin 实测 0.06mm）。改道顺序：
                    // 原地竖直 → 缩短法向退出 → 枪体安全距离，取两腿净空最大者。
                    Vec3 L1 = ga + na*clear;
                    const Vec3 L1top = liftTo(L1, sSafe);
                    if(!segClearOK(ga, L1) || segCrossFace(ga, L1) ||
                       !segClearOK(L1, L1top) || segCrossFace(L1, L1top)){
                        double best = std::min(segMinClear(ga, L1), segMinClear(L1, L1top));
                        Vec3 bestW = L1;
                        auto considerExit=[&](const Vec3& W){
                            const Vec3 Wtop = liftTo(W, sSafe);
                            double s = std::min(segMinClear(ga, W), segMinClear(W, Wtop));
                            if(s > best){ best = s; bestW = W; }
                        };
                        considerExit(ga);
                        for(double f : {0.75, 0.5, 0.25}) considerExit(ga + na*(clear*f));
                        considerExit(ga + na*gunClearEff);
                        L1 = bestW;
                        if(std::getenv("SPRAYTRAJ_TRANS_DBG"))
                            fprintf(stderr, "[trans] 曲面空移出口改道 L1=(%.1f,%.1f,%.1f) 最小净空 %.1f\n",
                                    L1.x,L1.y,L1.z, best);
                    }
                    move(L1, na);
                    move(liftTo(L1, sSafe), na);
                    const Vec3 Top = liftTo(gb, sSafe);
                    move(Top, nb);
                    if(!segClearOK(Top,gb) || segCrossFace(Top,gb)){
                        const Vec3 L2 = gb + nb*clear;
                        move(liftTo(L2, sSafe), nb);
                        move(L2, nb);
                    }
                    move(gb, nb);
                }
            };

            // ---- 端点选择（同平面口径）：首行首点 / 末行末点，取离上一面终点更近、且出得来的那端 ----
            // 反走 = 行序倒过来 + 每行点序倒过来，行间连接段（emitGap）的起止端点随之取反；
            // 走的是同一批线段，净空/穿面结论不变。出口那一端还要过一遍"退枪→抬到安全面"
            // 的净空判据（与平面路径同口径）：反走会把出口换到另一端，换成腔里就出不去了。
            bool revSurf = false;
            {
                const int lastK = surf.rowBeg[nr] - 1;
                if(!allPts.empty() && lastK > 0){
                    const Vec3 pe{allPts.back().x, allPts.back().y, allPts.back().z};
                    const int k0 = surf.rowBeg[0];
                    const Vec3 pF = surf.rowP[k0] + surf.rowN[k0]*params->standoff;
                    const Vec3 pL = surf.rowP[lastK] + surf.rowN[lastK]*params->standoff;
                    auto exitOK=[&](const Vec3& e, const Vec3& nn)->bool{
                        auto leg=[&](const Vec3& W)->double{
                            const Vec3 Wtop = liftTo(W, sSafe);
                            return std::min(segMinClear(e, W), segMinClear(W, Wtop));
                        };
                        double best = leg(e);
                        for(double f : {0.25, 0.5, 0.75, 1.0}) best = std::max(best, leg(e + nn*(clear*f)));
                        return best >= 0.9*gunClearEff;
                    };
                    // 正走以 pL 为出口；反走以 pF 为出口
                    const bool fOK = exitOK(pL, surf.rowN[lastK]);
                    const bool bOK = exitOK(pF, surf.rowN[k0]);
                    if(fOK && bOK)   revSurf = (len(pL - pe) < len(pF - pe) - 1.0);
                    else if(bOK)     revSurf = true;
                }
            }
            const int kEntry = revSurf ? (surf.rowBeg[nr]-1) : surf.rowBeg[0];
            Vec3 firstX = surf.rowP[kEntry];
            Vec3 firstN = surf.rowN[kEntry];
            Vec3 firstPos = firstX + firstN*params->standoff;

            // 面间过渡：从上一面终点抬枪绕行到曲面块首点（与平面路径同一套）
            {
                if(!allPts.empty()){
                    Vec3 pEnd{allPts.back().x, allPts.back().y, allPts.back().z};
                    const Vec3 nEnd{-allPts.back().i, -allPts.back().j, -allPts.back().k};
                    if(!usedCandidateCheck && tryShortFaceLink(pEnd, nEnd, firstPos, firstN, allPts)){
                        // 短换面：枪不停，直接连过去（点已按点距铺好，kind=1 橙色）
                        // 计数统一放在最后的过渡重建里，这里不重复加（默认会重建）
                    } else if(segClearOK(pEnd, firstPos) && !segCrossFace(pEnd, firstPos)){
                        // 直线全程离开件 >= moveClear，直连（不必抬到安全高度）。
                        // spray 是"段属性"（两端都为 1 才是喷涂段），所以必须在
                        // 中间插一个 spray=0 才能把线断开；插中点而不插端点，
                        // 三点共线、物理路径仍是这条直线，也不会产生零长段。
                        Vec3 mid = pEnd + (firstPos - pEnd)*0.5;
                        allPts.push_back({mid.x,mid.y,mid.z,
                                          -firstN.x,-firstN.y,-firstN.z, params->speed, 0});
                    } else {
                        Vec3 L1 = pEnd + prevNormal*clear;
                        const Vec3 L1top = liftTo(L1, sSafe);
                        if(!segClearOK(pEnd, L1) || segCrossFace(pEnd, L1) ||
                           !segClearOK(L1, L1top) || segCrossFace(L1, L1top)){
                            // 与曲面行间空移同样的出口改道：原地抬升 / 缩短法向退出
                            double best = std::min(segMinClear(pEnd, L1), segMinClear(L1, L1top));
                            Vec3 bestW = L1;
                            auto considerExit=[&](const Vec3& W){
                                const Vec3 Wtop = liftTo(W, sSafe);
                                double s = std::min(segMinClear(pEnd, W), segMinClear(W, Wtop));
                                if(s > best){ best = s; bestW = W; }
                            };
                            considerExit(pEnd);
                            for(double f : {0.75, 0.5, 0.25}) considerExit(pEnd + prevNormal*(clear*f));
                            considerExit(pEnd + prevNormal*gunClearEff);
                            L1 = bestW;
                        }
                        allPts.push_back({L1.x,L1.y,L1.z, prevNormal.x,prevNormal.y,prevNormal.z, params->speed, 0});
                        const Vec3 Lt = liftTo(L1, sSafe);
                        allPts.push_back({Lt.x,Lt.y,Lt.z, prevNormal.x,prevNormal.y,prevNormal.z, params->speed, 0});
                        const Vec3 Top = liftTo(firstPos, sSafe);
                        allPts.push_back({Top.x,Top.y,Top.z, -firstN.x,-firstN.y,-firstN.z, params->speed, 0});
                        // 直降穿件时沿首点法向退 clear，再降到侧上方进入
                        if(!segClearOK(Top, firstPos) || segCrossFace(Top, firstPos)){
                            const Vec3 L2 = firstPos + firstN*clear;
                            const Vec3 L2t = liftTo(L2, sSafe);
                            allPts.push_back({L2t.x,L2t.y,L2t.z, -firstN.x,-firstN.y,-firstN.z, params->speed, 0});
                            allPts.push_back({L2.x,L2.y,L2.z, -firstN.x,-firstN.y,-firstN.z, params->speed, 0});
                        }
                    }
                    ++transCount;
                }
                const size_t faceBeg = allPts.size();     // 本块面内轨迹的起点（供最后重排）
                Vec3 lastP{0,0,0}, lastN{0,0,1}; bool hasLast = false;
                for(int rr=0; rr<nr; ++rr){
                    const int r = revSurf ? (nr-1-rr) : rr;
                    const int b = surf.rowBeg[r], e = surf.rowBeg[r+1];
                    if(e-b < 2) continue;
                    const int k0 = revSurf ? (e-1) : b;
                    const int kstep = revSurf ? -1 : 1;
                    if(hasLast){
                        const auto* turn = meshRows ? &surf.rowTurns[revSurf ? r+1 : r] : nullptr;
                        if(turn){
                            const auto& points = *turn;
                            for(size_t j=0; j<points.size(); ++j)
                                emitPt(points[revSurf ? points.size()-1-j : j], surf.n0, 1);
                        }else{
                            emitGap(lastP, lastN, surf.rowP[k0], surf.rowN[k0]);
                        }
                    }
                    for(int k=k0;; k+=kstep){
                        emitPt(surf.rowP[k], surf.rowN[k], 1);
                        if(revSurf ? (k <= b) : (k+1 >= e)) break;
                    }
                    const int kEnd = revSurf ? b : (e-1);
                    lastP = surf.rowP[kEnd]; lastN = surf.rowN[kEnd];
                    hasLast = true;
                    totalRows++;
                }
                if(hasLast) prevNormal = lastN;
                if(allPts.size() > faceBeg){
                    if(fanCoverOn) fanHarvest(faceBeg, allPts.size(), surf.n0);
                    faceSpan.push_back({faceBeg, allPts.size(), si});
                }
            }
            {
                double angN = std::acos(std::min(1.0,std::fabs(dot(surf.n0,Vec3{0,0,1}))))*180.0/3.141592653589793;
                fprintf(stderr, "[traj] 面%zu(曲面): 内点=%zu 行数=%d 与Z夹角=%.1f° 行距=%.2fmm 覆盖率=%.1f%%\n",
                        si, inl.size(), nr, angN, surf.rowGap, surf.cover*100.0);
            }
            continue;
        }

        // ---- 枪位可达性：内壁/内底在角落处，枪体会顶到相邻壁 ----
        // 内壁的可喷侧朝腔内，枪得伸进腔里 standoff 才喷得到面。走到靠近相邻内壁
        // 的边缘时，枪位已经嵌进那块壁里——实测开口盒内壁轨迹两端各约 4mm 宽、
        // 最深 2mm 侵入；更糟的是内壁首点贴着壁时，面间过渡的"从安全高度垂直
        // 下降"会整整穿过 300mm 高的壁（实测一段 423mm 的直线穿墙）。
        // 判据：面点沿法向抬到枪位，查枪位到工件的净空（点云 + 面 + 棱边），太近就不布行。
        // 原来只对内壁/内底生效，理由是"外壁枪位在件外 standoff 处，最近点必然
        // >= standoff"。这个理由在**拐角/棱边处不成立**：外壁那段行一直铺到壁的
        // 边缘，枪位虽然离本面 100mm，离相邻的垂直壁只有几毫米（实测薄壁盒
        // g_norim 枪位距邻壁 0.1mm，轨迹整条贴着壁走）。现在几何判据不依赖
        // 点云是否采到棱边，对所有面统一执行——单独一块平板没有邻居面，
        // 判据自然放行，不会误删。
        if(!inl.empty()){
            ProfStage ps("4b.内点枪位过滤");
            // 留 1mm 余量：布行端点最终就落在过滤阈值上，若阈值正好等于 gunClear，
            // 行端点距邻壁刚好 20.0mm，点云拟合残差（box 实测 0.3~0.6mm）一叠加就
            // 变成 19.4 < 20，蛇形换行的连接段从端点出发必然判失败 → 整面退化成
            // 抬枪绕行。多留 1mm 让连接段有可用余量，覆盖损失 1mm 可忽略。
            const double gunClear = gunClearEff + 1.0;
            // 逐点独立，且只读点云/面几何（KD 查询是 const），可直接分块并行。
            // 大点云（68 万点）上这一段原本是整条流程 87% 的耗时。
            std::vector<unsigned char> keepFlag(inl.size(), 1);
            std::atomic<int> dropN{0};
            parallel_for((int)inl.size(), [&](int ii){
                const int i = inl[(size_t)ii];
                // 非代表点不参与判据：它们进不了 P2（isPlan 已在投影处挡掉），
                // 判了也是对结果零贡献，白白多花 3 倍查询。
                if(!isPlan[(size_t)i]) return;
                Vec3 q = Pf[i] + pl.n*params->standoff;
                // 只判 "d < gunClear"，用半径上界查询；远点 O(log n) 剪枝。
                double d = partTree.nearestWithinExcept(q, gunClear, faceSkipMask[si]);
                // 绝对兜底：**不排除任何面**，测枪位到整件的最近距离。
                // 为什么必须补这一条：上面那次查询会剪掉本面（faceSkipMask[si]），
                // 法向一旦反了，枪位就正好藏在本面材料里，而本面点又被剪掉 → 查不到 → 放行。
                // 正常面枪位离本面就是 standoff(100mm)，这条恒过、不误杀；法向反了会立刻掉到个位数。
                //
                // 这里**必须**用 nearestWithin 而不是精确 knn(q,1)：
                // 枪位天然在工件外 standoff(100mm) 处，精确 KNN 要把半径 100mm 球内的
                // 分支全部展开，薄壳点云上等于每次遍历 3 万+ 点（实测单次 ~200μs），
                // 整个过滤退化成 O(N × 3万)。而本处只判 "d < gunClear"，
                // nearestWithin 返回 min(真实最近距离, gunClear)，对这条判据逐位等价。
                d = std::min(d, partTree.nearestWithin(q, gunClear));
                d = std::min(d, planeClearAt(q, (int)si));   // 面几何（棱边点可能被去噪删掉）
                d = std::min(d, edgeClearAt(q, (int)si));    // 棱边/边界线
                if(d < gunClear){ keepFlag[(size_t)ii] = 0; ++dropN; }
            });
            std::vector<int> keep; keep.reserve(inl.size());
            for(size_t ii=0; ii<inl.size(); ++ii)
                if(keepFlag[ii]) keep.push_back(inl[ii]);
            const int drop = dropN.load();
            if(drop > 0 && keep.size() >= 20){
                seg.inl.swap(keep);        // inl 是 seg.inl 的引用，内容随之更新
                if(std::getenv("SPRAYTRAJ_FACE_DBG"))
                    fprintf(stderr, "[gun] 面%zu%s: 枪位距件 <%.0fmm，剔除 %d 点，剩 %zu\n",
                            si, seg.fo.isInner > 0 ? " 内壁" : "", gunClear, drop, seg.inl.size());
            }
        }

        // 该平面的局部基
        Vec3 zaxis{0,0,1};
        Vec3 u = norm(cross(zaxis, pl.n));
        if(len(u)<1e-6) u = Vec3{1,0,0};
        Vec3 v = cross(pl.n, u);

        // 单带件载体的"整条带轮廓"，行系 2D。空 = 沿用本面轮廓。
        // 只在这个面是单带件载体、且凸包有效时才被填（见下方夹端点那段的说明）。
        std::vector<std::vector<Vec2>> bandLoops;

        // 投影到 2D
        // 用布行层代表点建 P2：全分辨率点在这里只是"更密的同一批坐标"，
        // 对 minAreaRect / genRows / 分位数包围盒没有增益（体素 2mm 远小于行距），
        // 却让后面每个点的枪位过滤都要多算 3 倍。
        std::vector<Vec2> P2; P2.reserve(inl.size());
        for(int i:inl){ if(!isPlan[(size_t)i]) continue; Vec3 q=Pf[i]-pl.c; P2.push_back({dot(q,u), dot(q,v)}); }
        // P2 是布行的唯一输入。重建后的 keepMask 只按“栅格是否被占用”筛点，
        // 尾端孤立小簇仍可能留在 P2 里（p001 残留 19 点，距主面 147mm）：
        // minAreaRect 走凸包，19 个点就能把矩形拉成 248×112@134°。
        // 只剔“极小且远离主面”的簇；真实碎片/小分叉（点数>=1% 或距主面<=36mm）保留，
        // 继续交给 planeFrag / 轮廓切段处理。
        filterP2Outliers(P2, kReconGap, 3.0 * kReconGap, 0.01, 0.02, si);

        // 单带件载体：把同带其余面的点并入支撑集。
        // 只当"支撑点"用——面的平面、法向、枪位全部不变，变的只是行落点与行长度：
        // 不并入的话这条行会压在本片碎面自己的中线上（实测 49mm 宽的一条侧棱），
        // 整件 110mm 的带子只盖住一半；并入后行落在整件横向中心，两端也够长。
        if(singleBand && si == singleBandCarrier){
            size_t nAdded = 0;
            for(size_t j : bandFaces){
                if(j == si) continue;
                for(int k : segs[j].inl){
                    if(!isPlan[(size_t)k]) continue;      // 与主循环同口径：只吃代表点
                    const Vec3 q = Pf[k]-pl.c;
                    P2.push_back({dot(q,u), dot(q,v)});
                    ++nAdded;
                }
            }
            if(nAdded)
                fprintf(stderr, "[traj] 面%zu 单带件支撑合并: +%zu 点（来自同带其余 %zu 片）\n",
                        si, nAdded, bandFaces.size()-1);
        }

        // 补全面：把补全矩形作为"虚拟支撑点"并入点集，让 genRows 在补全区也布行。
        // 否则会出现"面补齐了、轨迹还是只走有数据的那一段"。
        //
        // 关键：虚拟点必须**沿矩形的两条边**铺，不能铺成轴对齐包围盒（AABB）。
        // 矩形在轨迹系里与行方向存在几度偏差时，990×28.5 的长条取 AABB 会变成
        // 990×107（4.5° 偏差就能撑出 78mm），于是窄边被排出 9 行、轨迹张成扇形。
        // 同时把该面的行方向**强制取矩形长边**，让矩形在旋转后严格轴对齐。
        double forcedTheta = -1.0;
        const size_t nRealP2 = P2.size();     // inl 部分到此为止，后面是补全虚拟点
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
        // ---- 补全虚拟点的枪位过滤 ----
        // 虚拟支撑点是沿矩形边"外推"出来的，**不经过前面那段 inl 枪位过滤**
        // （那段在 P2 构建之前、只作用于真实点云点）。补全区恰恰是"没有数据、
        // 靠外推猜"的地方，最可能猜错：实测 complex 面4 补全 42mm 宽，行落在
        // 补全区里，枪位距工件只剩 1.2mm —— gunClear 硬违约 823 个点。
        // 这里对补全新增的那部分 P2 再过滤一次；删掉的点就是不该喷的外推区。
        if(P2.size() > nRealP2){
            std::vector<Vec2> keepP2; keepP2.reserve(P2.size());
            keepP2.insert(keepP2.end(), P2.begin(), P2.begin()+(ptrdiff_t)nRealP2);
            int dropF = 0;
            for(size_t k=nRealP2; k<P2.size(); ++k){
                Vec3 X = pl.c + u*P2[k].x + v*P2[k].y;
                if(gunClearAt(X, pl.n, (int)si) < gunClearEff + 1.0){ ++dropF; continue; }
                keepP2.push_back(P2[k]);
            }
            if(dropF){
                P2.swap(keepP2);
                fprintf(stderr, "[gun] 面%zu: 补全区 %d/%zu 点枪位不可达（<%.0fmm），已剔除\n",
                        si, dropF, P2.size()+(size_t)dropF-1, gunClearEff+1.0);
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
        // 近水平面（法向与内部 +Z 夹角 < 2.9°）的"行方向被离群点带歪"纠偏。
        // minAreaRect 取凸包边，起停扫的稀疏斜向线尾 / 重建栅格量化出的斜边，
        // 会以极少数点把这条边带歪（双长板实测 minAreaRect 与数据主轴差 21.5°，
        // 行方向映回 3D 便歪 21.25°，还被 clipRowByLoops 夹成 183mm 的斜短线）。
        // PCA 由全部面内点投票，对这类离群点稳健。
        // ★ 只做"异常纠偏"而非"替换"：偏差 <= 10° 时保持 minAreaRect 原值，
        //   正常件（含斜面、侧壁、环形面）逐位不变；只有明显被污染时才切到 PCA，
        //   避免把全量回归里的行/点数扰动带进现场。
        if(haveRectAng && P2.size() >= 3 && len(cross(zaxis, pl.n)) < 0.05){
            double mx=0.0, my=0.0;
            for(const auto& q : P2){ mx += q.x; my += q.y; }
            mx /= (double)P2.size(); my /= (double)P2.size();
            double sxx=0.0, sxy=0.0, syy=0.0;
            for(const auto& q : P2){
                double dx=q.x-mx, dy=q.y-my;
                sxx += dx*dx; sxy += dx*dy; syy += dy*dy;
            }
            double thPCA = 0.5*std::atan2(2.0*sxy, sxx-syy);
            double d = std::atan2(std::sin(thLong-thPCA), std::cos(thLong-thPCA));
            if(d >  1.5707963267948966) d -= 3.141592653589793;
            if(d < -1.5707963267948966) d += 3.141592653589793;
            // 只纠"细长面"：长短轴方差比 > 16（长宽比 > 4）。这类面的长轴是唯一的、
            // 工艺上也必须沿长边走；近方形 / L 形面 minAreaRect 与 PCA 差 70°+ 属正常
            // （两条边等长，选哪条都对），一旦纠偏反而把既有布局整片换掉。
            double disc = std::sqrt((sxx-syy)*(sxx-syy) + 4.0*sxy*sxy);
            double lam1 = 0.5*(sxx+syy+disc), lam2 = 0.5*(sxx+syy-disc);
            bool slender = lam2 > 1e-9 && lam1 > 16.0*lam2;
            bool fix = slender && std::fabs(d) > 10.0*3.141592653589793/180.0;
            if(std::getenv("SPRAYTRAJ_FACE_DBG"))
                fprintf(stderr, "[rect] 面%zu 近水平面: minAreaRect %.1f° / PCA %.1f° 偏差 %.1f° 细长=%d -> %s\n",
                        si, thLong*180.0/3.141592653589793, thPCA*180.0/3.141592653589793,
                        std::fabs(d)*180.0/3.141592653589793, (int)slender,
                        fix ? "采用 PCA" : "保持不动");
            if(fix) thLong = thPCA;
        }
        if(std::getenv("SPRAYTRAJ_FACE_DBG"))
            fprintf(stderr, "[rect] 面%zu P2=%zu haveRect=%d thLong=%.1f forced=%.1f  c=(%.1f,%.1f,%.1f) n=(%.2f,%.2f,%.2f)\n",
                    si, P2.size(), (int)haveRectAng, thLong*180.0/3.141592653589793,
                    forcedTheta*180.0/3.141592653589793,
                    pl.c.x, pl.c.y, pl.c.z, pl.n.x, pl.n.y, pl.n.z);
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
            bestRows = genRows(P2, theta, L, params->fillHoleSize, acrossExt);
        }
        else if(haveRectAng){
            // mode1/0 = 沿长边；mode2 = 沿短边（长边方向转 90°）
            // ---- 候选 A：矩形长边（几何基准，先算出来）----
            theta = normAng((params->scanDirMode == 2) ? thLong + 1.5707963267948966 : thLong);
            double e1=0,e2=0; extentAlong(theta,e1,e2);
            bool wantShort = (params->scanDirMode == 2);
            if((!wantShort && e1 < e2) || (wantShort && e1 > e2)){
                theta = normAng(theta + 1.5707963267948966);
                std::swap(e1,e2);
            }
            const std::vector<Row> rowsLong = genRows(P2, theta, L, params->fillHoleSize, acrossExt);

            // ---- 候选 B：竖直面取面内水平 ----
            // 这条规则治的是"同一工件上方形侧壁(300×300)一横一竖交叉成乱麻"：
            // minAreaRect 的两条边等长时，边序由凸包起点决定，x=0/x=400 壁排成沿 z、
            // y=0/y=300 壁排成沿 x，用户反馈"就不能是规则的轨迹嘛"。工艺上侧壁也
            // 确实是水平往复（枪沿水平走、靠升降换行，重力方向重叠更均匀）。
            //
            // **但水平只是打破平局的手段，不能无条件覆盖几何最优方向。** 对非等边的
            // 面，水平方向并不等价：242×40 的立板长边在面内 88.5°、水平在 20°，锁水平
            // 等于沿 40mm 短边布行，行与重建轮廓斜交 68.5°，clipRowByLoops 把这类斜行
            // 整条判为面外丢掉（行距 60 时 3 行丢 2 行，只剩一条 91mm 斜线，两侧各出界
            // 22mm）。同理 400×150 若被锁到偏 30° 的水平向，行数也会翻两三倍。
            //
            // 因此判据改为"**两个候选排出来的行布局比较**"，取水平的前提是它与长边
            // 方向基本等价：行数可差 1（行数是取整量）、行均长不低于 85%。超出这个
            // 等价带宽就说明水平方向明显更碎（换行更多、行更短），坚持它等于拿
            // "方向一致"换"越喷越碎"。
            // 判据只读两个候选自身的行布局，与面的绝对尺寸、长宽比、行距、坐标系
            // 都无关 —— 不是"长宽比 >= N"这类按尺寸拍的阈值。
            //
            // **只对真·竖直面生效**（|n.z| < 0.3，即与竖直夹角 >72°）：斜面（曲型板
            // 那个 52° 面，n.z=0.61）长边本来沿坡向，取水平会把 4 行 300mm 切成 13 行
            // ~90mm（实测总行数 18→27），本就不在候选范围内。
            bool forcedHoriz = false;
            if(std::fabs(pl.n.z) < 0.3){
                Vec3 hz = cross(Vec3{0.0,0.0,1.0}, pl.n);
                const double hl = len(hz);
                if(hl > 1e-9){
                    hz = hz * (1.0/hl);
                    double thH = std::atan2(dot(hz,v), dot(hz,u));
                    thH = normAng((params->scanDirMode == 2)
                                  ? thH + 1.5707963267948966 : thH);
                    const std::vector<Row> rowsH = genRows(P2, thH, L, params->fillHoleSize, acrossExt);
                    if(!rowsH.empty() && !rowsLong.empty()){
                        auto meanLen = [](const std::vector<Row>& rs){
                            double s = 0.0;
                            for(const Row& r : rs) s += (r.xe - r.xs);
                            return rs.empty() ? 0.0 : s/(double)rs.size();
                        };
                        const double mL = meanLen(rowsLong), mH = meanLen(rowsH);
                        // 行数容差取"1 行"与"10%"的较大者：行数是取整量，等长的两个
                        // 方向因行端裁剪/边界吸附会差 1~2 行（box.ply 的 300×300 侧壁
                        // 实测 22 行 vs 24 行，差 9%），用固定 ±1 会把等边面误判成"更碎"，
                        // 又退回一横一竖。
                        const double rl = (double)rowsLong.size();
                        const double tolRows = std::max(1.0, std::ceil(0.10*rl));
                        const bool sameRows = ((double)rowsH.size() <= rl + tolRows);
                        const bool sameLen  = (mL <= 1e-9) || (mH >= 0.85*mL);
                        if(sameRows && sameLen){
                            theta   = thH;
                            bestRows = rowsH;
                            forcedHoriz = true;
                        }
                    }
                }
            }
            if(!forcedHoriz) bestRows = rowsLong;
            extentAlong(theta, e1, e2);
            // 只有"面足够宽、却只排出不到 3 行"才说明方向选错了、需要重选；
            // 窄面（垂直行向的宽度撑不满 3 行）本来就只能铺 1~2 行，必须保留，
            // 否则会掉进下面的穷举搜索，选到沿窄边的方向 → 上百条 20mm 的碎行。
            if(bestRows.size() < 3 && e2 >= 3.0*L) bestRows.clear();
        }
        if(bestRows.empty()){
            double bestCost=1e18;
            for(int deg=0; deg<180; deg+=5){
                double ang=deg*3.141592653589793/180.0;
                auto rs=genRows(P2,ang,L,params->fillHoleSize, acrossExt);
                if(rs.empty()) continue;
                // 单行候选的合法条件：这个方向**跨行跨度本来就放不下第二行**
                // （<= 1.2 倍行距）。跨度明明够两行以上、却只排出 1 行，说明多行被
                // 支撑点/碎段过滤吃掉了，属退化，不能当最优方向。
                // 旧写法是 size()<2 一刀切淘汰——恰好把"面确实窄、只该有一行"的
                // 正确方向（行距 >= 面宽时沿长边那条）直接扔掉，于是被迫选斜向多行，
                // 再被轮廓裁剪砍成一条出界斜线。判定用跨度比 L，与具体尺寸无关。
                double s1=0.0, s2=0.0; extentAlong(ang, s1, s2);
                const bool oneRow = (rs.size() < 2);
                if(oneRow && s2 > 1.2*L) continue;
                double jag=0, rowLen=0.0;
                for(size_t i=0;i<rs.size();++i){
                    rowLen += std::max(0.0, rs[i].xe - rs[i].xs);
                    if(i) jag += std::fabs(rs[i].xs-rs[i-1].xs) + std::fabs(rs[i].xe-rs[i-1].xe);
                }
                if(rs.size() > 1) jag /= (double)(rs.size()-1);
                double cost = jag + rs.size()*L*0.02;
                // 单行候选之间比行长：只能走一行时，工艺上就该沿面最长方向走，
                // 行越长覆盖面越大、空走越少。系数取 0.02（行数代价同一量级），
                // 不会把多行候选的整体排序带偏。
                if(oneRow) cost -= rowLen*0.02;
                if(cost<bestCost){ bestCost=cost; theta=ang; bestRows=rs; }
            }
        }
        // 注意：这里**不能**直接 continue。窄环带（环宽 < 行距）本来就可能一行都排不出来
        // （行被短段过滤吃光），但它是环，应该由下面的同心圆/沿轮廓/射线环带走，
        // 提前返回等于把整圈漏喷。真正的"无行"判定放到环带判据之后。
        const bool noRows0 = bestRows.empty();
        if(noRows0)
            fprintf(stderr, "[traj] 段(点数=%zu) 没有可用平行行，先试环形/类环形判据\n", inl.size());
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
        // ---- 框形/环形面：沿轮廓走同心偏置环（用户 2026-11-20 14:19 拍板）----
        // 现场是长条形框料（p047 面1 = 712×616 的框、带宽约 31mm、占空 19%）。
        // 平行扫描会把一条连续环带切成几十行，而且每行还要被中央空洞断开成左右两段，
        // 换行/抬枪比有效行程还多；环带本身只值几圈。
        // 这里改成"沿外轮廓与内孔之间等分出 k 圈闭合环"：
        //   外轮廓每个顶点找内孔上最近的一点，按 t=j/(k+1) 线性插值出第 j 圈。
        // 两条边界近似平行（框料就是等宽的），插值出的环自然也是等距内缩的。
        // 判据（缺一条都会误伤实心盘 / 小垫片 / 带一个小孔的板）：
        //   ① 有内孔（loops>=2）；② 占空 < 55%（实心盘 97%）；
        //   ③ 带宽 <= 3×行距（宽了平行扫描本来就正常）；
        //   ④ 跨度 >= 3×带宽（排除小垫片，一圈周长还没行距长）。
        bool loopMode = false;
        std::vector<std::vector<Vec2>> loopPaths;
        // 每条环是否闭合：轮廓环（loopMode）恒闭合；射线环带的开口弧（C/U 形）不闭合。
        // 下游按"闭合折线"采样时会把首尾自动接上，开口弧必须显式声明，否则
        // 弧尾会直接跨回到弧头，在缺口上画出一道穿腔直线（实测 ring_cut 多跑 175mm）。
        std::vector<char> loopOpen;
        // noRows0：连一行平行行都排不出来的段，**不走轮廓环**。轮廓环依赖重建轮廓
        // 分出内外两圈，而"排不出行"多半意味着点太少/轮廓退化（实测 p028 这类
        // 59 点、208mm 长、带宽 1.1mm 的稀疏条会被判成"沿轮廓 1 环"，只出 14 个点，
        // 把整件从 57 点压到 14 点）。这类段只保留环形/类环形的机会——
        // 那两条看的是点分布的极坐标结构，不看轮廓。
        if(!ringMode && !noRows0 && !std::getenv("SPRAYTRAJ_LOOP_OFF")){
            const auto& loopsUU = planeLoops[si];
            if(loopsUU.size() >= 2){
                const double ct=std::cos(theta), stt=std::sin(theta);
                std::vector<std::vector<Vec2>> lps; lps.reserve(loopsUU.size());
                for(auto& lp0 : loopsUU){
                    std::vector<Vec2> q; q.reserve(lp0.size());
                    for(auto& t : lp0) q.push_back(Vec2{ct*t.x + stt*t.y, -stt*t.x + ct*t.y});
                    lps.push_back(std::move(q));
                }
                size_t iOut=0; double aMax=0;
                for(size_t k=0;k<lps.size();++k){
                    const double a=std::fabs(polyArea2(lps[k]));
                    if(a>aMax){ aMax=a; iOut=k; }
                }
                double perOut=0;
                for(size_t k=0;k<lps[iOut].size();++k)
                    perOut += len2(lps[iOut][(k+1)%lps[iOut].size()] - lps[iOut][k]);
                double fill=aMax, ahMax=0, perHole=0; size_t iHole=lps.size();
                for(size_t k=0;k<lps.size();++k){
                    if(k==iOut) continue;
                    const double a=std::fabs(polyArea2(lps[k]));
                    fill -= a;
                    double pf=0;
                    for(size_t t2=0;t2<lps[k].size();++t2)
                        pf += len2(lps[k][(t2+1)%lps[k].size()] - lps[k][t2]);
                    if(a>ahMax){ ahMax=a; iHole=k; perHole=pf; }
                }
                double ux0=1e18,ux1=-1e18,vy0=1e18,vy1=-1e18;
                for(auto& lp0:lps) for(auto& q:lp0){
                    ux0=std::min(ux0,q.x); ux1=std::max(ux1,q.x);
                    vy0=std::min(vy0,q.y); vy1=std::max(vy1,q.y);
                }
                double bboxA = std::max(1e-9,(ux1-ux0)*(vy1-vy0));
                const double occ = fill/bboxA;
                const double bandW = fill/std::max(1.0, 0.5*(perOut+perHole));
                const double span = std::max(ux1-ux0, vy1-vy0);
                // "孔"必须真的在"外圈"里面。轮廓重建给出的多个 loop 里，既可能是
                // 外圈+内孔（真框形），也可能是**同一张面被切成几块**（实测 p043：
                // 208 点、20mm 宽的长条被切成 2 块，第 2 块被当成"内孔"，于是判成
                // 沿轮廓 1 环，69 个点的正常单行退化成 20 个点）。
                // 用"内孔顶点是否落在外圈多边形内"区分，这两件事在几何上完全不同。
                const bool holeInside = (iHole < lps.size()) && !lps[iHole].empty()
                                      && pointInPoly(lps[iHole][0], lps[iOut]);
                if(iHole<lps.size() && holeInside && occ<0.55 && bandW>=3.0 && bandW<=3.0*L && span>=3.0*bandW){
                    const int kk = std::max(1, (int)std::floor(bandW/L + 0.5));
                    const auto& O = lps[iOut]; const auto& H = lps[iHole];
                    loopPaths.reserve((size_t)kk);
                    loopOpen.assign((size_t)kk, 0);
                    for(int j=1;j<=kk;++j){
                        const double t = (double)j/(double)(kk+1);
                        std::vector<Vec2> lp; lp.reserve(O.size());
                        for(auto& o : O){
                            double db=1e18; Vec2 hh{0,0};
                            for(auto& h0 : H){ const double d2=len2(h0-o); if(d2<db){ db=d2; hh=h0; } }
                            lp.push_back(Vec2{o.x + (hh.x-o.x)*t, o.y + (hh.y-o.y)*t});
                        }
                        loopPaths.push_back(std::move(lp));
                    }
                    loopMode = !loopPaths.empty();
                    fprintf(stderr, "[traj] 面%zu 框形判据: 占空%.0f%% 带宽%.1fmm 跨度%.0fmm → 沿轮廓 %d 环（平行扫描要 %zu 行）\n",
                            si, occ*100.0, bandW, span, kk, bestRows.size());
                }
            }
        }
        // ---- 环形 / 类环形：沿环带走轨迹（用户 2026-11-20 14:23「要考虑环形、类环形」）----
        // 与上两条互补：
        //   ringMode  正圆同心圆 —— 只在"正圆 + 同心 + 方位铺满 34/36"时成立；
        //   loopMode  有内孔的框形轮廓 —— 要求重建轮廓分出 >=2 个环；
        // 这里补上两者都盖不住的环带：椭圆环、跑道环、圆角矩形环、方框（内孔被
        // fillHoleSize 糊住或轮廓没分出内孔）、以及开口环（C/U 形）。
        // 几何来自 detect_band_uv 的射线模型：rIn/rOut 逐方位给出环带的内外边界，
        // 第 j 圈就取 rIn + (rOut-rIn)·j/(k+1)，方位角连起来即是贴合的闭合/开口曲线。
        bool bandMode = false;
        if(!ringMode && !loopMode && forcedTheta < 0 && !std::getenv("SPRAYTRAJ_BAND_OFF")){
            BandUV bi;
            const bool bandOK = detect_band_uv(P2, bi);
            if(std::getenv("SPRAYTRAJ_BAND_DBG"))
                fprintf(stderr, "[band] 面%zu N=%d 种子(%.1f,%.1f) 方位%.0f%% 多段%.0f%% "
                                "带宽%.1f⊥(内%.1f 外%.1f) 行距%.1f → %s\n",
                        si, bi.N, bi.c.x, bi.c.y, bi.closure*100.0, bi.multi*100.0,
                        bi.wPerp, bi.rInMed, bi.rOutMed, L, bandOK ? "是环带" : "否");
            if(bandOK){
                const int kk = std::max(1, (int)std::floor(bi.wPerp/L + 0.5));
                const double ct=std::cos(theta), st=std::sin(theta);
                auto toRow=[&](Vec2 q)->Vec2{ return Vec2{ ct*q.x + st*q.y, -st*q.x + ct*q.y }; };
                // 小缺口（<=2 个方位，约 7~15mm）按邻近方位线性插值补上：这些缺口来自
                // 稀疏/遮挡，不是结构断口。不补的话一圈会被切成 2~3 段，每段都要单独
                // 起停 + 抬枪（p047 那块闭合框带 73% 覆盖，被 1 个方位的缺口切成 2 段）。
                {
                    const int MAXGAP = 2;
                    int a = 0;
                    while(a < bi.N){
                        if(bi.rIn[(size_t)a] >= 0.0f){ ++a; continue; }
                        int g=0; while(a+g < bi.N && bi.rIn[(size_t)(a+g)] < 0.0f) ++g;
                        if(g <= MAXGAP && a > 0 && a+g < bi.N &&
                           bi.rIn[(size_t)(a-1)] >= 0.0f && bi.rIn[(size_t)(a+g)] >= 0.0f){
                            for(int t=0;t<g;++t){
                                const double f = (double)(t+1)/(double)(g+1);
                                bi.rIn[(size_t)(a+t)]  = (float)(bi.rIn[(size_t)(a-1)]*(1.0-f) + bi.rIn[(size_t)(a+g)]*f);
                                bi.rOut[(size_t)(a+t)] = (float)(bi.rOut[(size_t)(a-1)]*(1.0-f) + bi.rOut[(size_t)(a+g)]*f);
                            }
                        }
                        a += std::max(g,1);
                    }
                }
                // 连续有材料的方位段；首尾都覆盖时合并成一条跨 0 的闭合环
                std::vector<std::vector<int>> idxRuns;
                {
                    std::vector<int> cur;
                    for(int a=0;a<bi.N;++a){
                        if(bi.rIn[(size_t)a] >= 0.0f) cur.push_back(a);
                        else if(!cur.empty()){ idxRuns.push_back(cur); cur.clear(); }
                    }
                    if(!cur.empty()){
                        if(!idxRuns.empty() && idxRuns.front().front()==0 && cur.back()==bi.N-1){
                            std::vector<int> merged = cur;
                            merged.insert(merged.end(), idxRuns.front().begin(), idxRuns.front().end());
                            idxRuns.front() = std::move(merged);
                        } else idxRuns.push_back(cur);
                    }
                }
                if(std::getenv("SPRAYTRAJ_BAND_DBG")){
                    fprintf(stderr, "[band] 面%zu 环带分段: %zu 段", si, idxRuns.size());
                    for(const auto& r0 : idxRuns)
                        fprintf(stderr, " [%d..%d]/%zu", r0.front(), r0.back(), r0.size());
                    fprintf(stderr, "\n");
                }
                for(int j=1;j<=kk;++j){
                    const double t = (double)j/(double)(kk+1);
                    for(const auto& run : idxRuns){
                        if(run.size() < 3) continue;
                        std::vector<Vec2> lp; lp.reserve(run.size()+1);
                        for(int a : run){
                            const double r = (double)bi.rIn[(size_t)a]
                                           + ((double)bi.rOut[(size_t)a]-(double)bi.rIn[(size_t)a])*t;
                            const double aa = 2.0*3.14159265358979*(double)a/(double)bi.N;
                            lp.push_back(toRow(Vec2{ bi.c.x + r*std::cos(aa), bi.c.y + r*std::sin(aa) }));
                        }
                        // 开口弧：不补首点（补了也只会多出一条 0 长边，首尾照样被当成
                        // 闭合边接上），而是登记 loopOpen，由下游采样时跳过收尾连接段。
                        const bool closed = (run.size() == (size_t)bi.N);
                        if(!closed && run.size() >= 2 && len2(lp.back() - lp.front()) < 1e-6)
                            lp.pop_back();
                        loopPaths.push_back(std::move(lp));
                        loopOpen.push_back(closed ? 0 : 1);
                    }
                }
                bandMode = !loopPaths.empty();
                if(bandMode)
                    fprintf(stderr, "[traj] 面%zu 类环形判据: 内径%.1f 外径%.1f 带宽%.1f⊥ 方位%.0f%% 缺口%.0f%% "
                                    "→ 沿环带 %d 圈/%zu 段（平行扫描要 %zu 行）\n",
                            si, bi.rInMed, bi.rOutMed, bi.wPerp, bi.closure*100.0, bi.multi*100.0,
                            kk, loopPaths.size(), bestRows.size());
                loopMode = bandMode;
            }
        }
        // 三条环带判据都没成立、平行行也没有 → 这个段确实布不出轨迹
        if(!ringMode && !loopMode && bestRows.empty()){
            fprintf(stderr, "[traj] 段(点数=%zu) 没有可用行，跳过\n", inl.size());
            continue;
        }
        // ---- 按重建轮廓裁剪行端点（抑制"轨迹出框"的关键一环）----
        // bestRows 的 xs/xe 来自 P0.5/P99.5 分位包围盒，会把"面在这一行没有材料"
        // 的空档也算成行长。这里用该面的重建轮廓多边形把每行夹回面内。
        // 轮廓已把"fillHole 以内该填、以外是真空"的决策做完了，所以这里不需要
        // 再按 gapSeg 猜一次——两侧口径一致，不会出现"面没补、轨迹却喷过去"。
        // 环形面不走这里（同心圆本来就贴着环心/环宽算，没有包围盒问题）。
        // ---- 单带件载体：轮廓换成"整条带"的轮廓 ----
        // 这一层是**保险**，不是 p048 的根因（A/B 实测：把载体选择修好后，
        // 开关这一层都得 143 点 / 1 行，贡献为 0）。p048 真正的根因是载体选错，
        // 见前面"载体选择"那段 —— 选了边角小片，行就按小片几何铺、按小片轮廓夹。
        //
        // 它防的是另一种情形：带子由 3 片以上拼成、载体只覆盖其中一段。
        // 那种情况下同带其余面的点已并进 P2（行按整条带铺），但夹端点用的
        // loops 仍只有载体那一片的范围，行会被削掉另外几段。
        // 口径与侧面吸收一致（用户 2026-11-23：「吸收了就成了主平面啊，
        // 按照重建平面来布行，重建平面怎么样就怎么样布行」）：合成了整条带，
        // 就按整条带的实际范围布行。
        // 做法：用 P2 全部支撑点在行系下的凸包作唯一轮廓。
        // 凸包会把带子上的凹陷填平——对单带件这是对的：带子横向跨度 <= 一个行距，
        // 一条行的喷幅本来就把整条盖住，凹陷处也是喷幅覆盖范围内的料
        // （与"波浪瓦楞板按平面处理"同一口径）。
        // 关掉：SPRAYTRAJ_BANDCLIP_OFF=1
        if(singleBand && si == singleBandCarrier && !P2.empty() && !std::getenv("SPRAYTRAJ_BANDCLIP_OFF")){
            const double cr=std::cos(theta), sr=std::sin(theta);
            std::vector<Vec2> hullIn; hullIn.reserve(P2.size());
            for(auto& q : P2) hullIn.push_back(Vec2{ cr*q.x + sr*q.y, -sr*q.x + cr*q.y });
            std::vector<Vec2> hull = convexHull(hullIn);
            if(hull.size() >= 3 && std::fabs(polyArea2(hull)) > 4.0*kReconGap*kReconGap){
                // **不外扩**：行端点本来就是从同一份 P2 的分位包围盒出来的，
                // 已经落在凸包内。这里只负责"别再用小片轮廓去夹"，
                // 不负责把行抻长——抻长是后面"行端点外扩"那一环的事（带枪位校验）。
                bandLoops.push_back(std::move(hull));
                fprintf(stderr, "[traj] 面%zu 单带件载体: 行端点改用整条带轮廓（%zu 顶点）"
                                "代替本面 %zu 环轮廓\n",
                        si, bandLoops[0].size(), planeLoops[si].size());
            }
        }
        std::vector<std::vector<Vec2>> rowFootprint;
        bool routeAroundVoid = false;
        if(!ringMode && !loopMode && !bestRows.empty() && !std::getenv("SPRAYTRAJ_NO_ROWCLIP")){
            const std::vector<std::vector<Vec2>>& loopsUU = planeLoops[si];
            if(!loopsUU.empty()){
                // 轮廓存的是 (uu,vv) 口径，行坐标是再转 theta 之后的口径，这里统一到行系。
                const double ct=std::cos(theta), st=std::sin(theta);
                std::vector<std::vector<Vec2>> loops; loops.reserve(loopsUU.size());
                for(auto& lp : loopsUU){
                    std::vector<Vec2> q; q.reserve(lp.size());
                    for(auto& t : lp) q.push_back(Vec2{ ct*t.x + st*t.y, -st*t.x + ct*t.y });
                    loops.push_back(std::move(q));
                }
                // 单带件载体：换成整条带的轮廓（见上方 bandLoops 的注释）
                if(!bandLoops.empty()) loops = bandLoops;
                rowFootprint = loops;
                // ---- 只做"端点夹紧"，不做切段 ----
                // 反复试过三种策略，只有这一种在所有用例上单调改善：
                //   A) 按轮廓切段 → 框形面对了，但稀疏/锯齿轮廓被切成碎屑
                //      （norim 行数 195→268、zhe 面1 行数 5→2），得不偿失；
                //   B) 面级判"轮廓可信" → 判据在真框形面和锯齿面之间没有清晰分界，
                //      同一批用例两边都误判；
                //   C) 只夹端点（本实现）→ 出框的根本原因是"行端点被稀疏离群点
                //      或分位包围盒撑到面外"，把端点收回到轮廓交集里就消除 90% 的
                //      越界；行内部的真空档交给上游 gapSeg（它本来就管这件事）。
                // 关键：只取"与当前行重叠最长的那一段"并把 xs/xe 收进去，
                // 段的取舍尽量少改动行数——宁可行短一点，也不要碎成几十段。
                // ---- 只处理"轮廓相对行严重缩水"的行，其余原样保留 ----
                // 反复试过四种策略，只有这一种在所有用例上单调改善：
                //   A) 按轮廓切段 → 框形面对了，但稀疏/锯齿轮廓被切成碎屑
                //      （norim 行数 195→268、zhe 面1 行数 5→2），得不偿失；
                //   B) 面级判"轮廓可信" → 判据在真框形面和锯齿面之间没有清晰分界；
                //   C) 无条件夹端点 → 会漏喷。实测诊断：complex 面1 夹掉 1091mm 里
                //      58% 有材料支撑、面4 是 67%、cao 面0 是 94% —— 点云在那些位置
                //      确实有料，只是重建轮廓因形态学开运算"缩"了进去。
                //   D) 本实现：只对"轮廓相比行缩水严重"的行动手。
                // 判据：轮廓给出的可喷长度 / 行原长 < 40% 才夹/丢。这个门槛对应
                // "行有六成以上落在重建面外"，那已经不是边界误差而是整段悬空
                // （bigplane 面0 是 5%、norim 面2 是 35%），而正常面的边界误差
                // 都在 90% 以上（complex 面1 = 91%、cao 面0 = 97%）。
                std::vector<Row> kept; kept.reserve(bestRows.size());
                int nClip = 0, nDrop = 0;
                double cutLen = 0.0;
                // 碎片化面（重建网格碎成多块、最大块 <= 50%）：这类面本身就是
                // "两块材料被一个空档连着"——行必须按轮廓切段，否则一定横穿空档。
                // bigplane 面0 就是这种：164 个点被重建成 3 块共 14,467mm²，
                // 行却按 675×418 的包围盒铺，实测对着空气喷 211mm。
                // 判据放在面级（不按行判），所以不会像逐行切段那样在锯齿面上碎掉：
                // 真·锯齿面（norim 面9）连通块 4 但最大块 52%、面积 46,376mm²，
                // 已经 >= 50% 门槛，走"只夹端点"分支。
                const bool fragmented = planeFrag[si] || sparseDisconnected[si];
                if(fragmented && std::getenv("SPRAYTRAJ_ROWCLIP_DBG"))
                    fprintf(stderr, "[traj] 面%zu: 重建网格碎片化 → 行按轮廓切段\n", si);
                // ---- 首末行豁免"缩水"判据 ----
                // 首末行按设计就压在点云边沿上，行带 |Δy|<=L/2 天生只有一半有料：
                // 它的行长是"边沿那排点"给出的，本来就比中间行短，被 40% 缩水判据
                // 一夹就断成半截（残白 242×40 立板 row=12：末行 179→58mm，接着蛇形
                // 过渡斜跨 117mm，3D 里就是一条孤零零的斜线——用户读成"轨迹歪了"）。
                // 首末行不做压缩水夹取：轮廓在这里缩进是常态（形态学开运算），
                // 夹它只会制造空洞。整行落在轮廓外的丢弃判断也不适用于它们
                // （行端最远只超出 extRow，且 extRow 已被面尺度夹过）。
                double ycLo=1e18, ycHi=-1e18;
                for(const Row& q : bestRows){ ycLo=std::min(ycLo,q.yc); ycHi=std::max(ycHi,q.yc); }
                const bool multiRow = bestRows.size() > 1;
                for(const Row& r : bestRows){
                    const bool edgeRow = multiRow &&
                        (std::fabs(r.yc - ycLo) < 1e-6 || std::fabs(r.yc - ycHi) < 1e-6);
                    auto iv = clipRowByLoops(r.yc, r.xs, r.xe, loops);
                    bool materialGap = false;
                    for(size_t k=1; k<iv.size(); ++k)
                        if(iv[k].first-iv[k-1].second >= kReconGap){ materialGap = true; break; }
                    if(edgeRow && !fragmented && !materialGap){ kept.push_back(r); continue; }
                    const double rowLen = std::max(1e-9, r.xe - r.xs);
                    if(std::getenv("SPRAYTRAJ_ROWCLIP_DBG")){
                        fprintf(stderr, "[clip] yc=%8.1f row=[%8.1f,%8.1f] len=%7.1f iv=%zu",
                                r.yc, r.xs, r.xe, rowLen, iv.size());
                        for(auto& s : iv) fprintf(stderr, " [%.1f,%.1f]", s.first, s.second);
                        fprintf(stderr, "\n");
                    }
                    if(iv.empty()){
                        // 整行落在轮廓外：只有"行本身很短"才认为它是噪声（长行整条在面外
                        // 往往说明轮廓被开运算削掉了，删了就是漏喷）
                        if(rowLen < 2.0*L || fragmented){ ++nDrop; continue; }
                        kept.push_back(r); continue;
                    }
                    double cov = 0.0;
                    for(auto& s : iv) cov += std::max(0.0, std::min(s.second, r.xe) - std::max(s.first, r.xs));
                    if(fragmented || materialGap){
                        // 碎片化面：每个区间都是独立的一块材料，逐段保留
                        if(iv.size() > 1 || cov / rowLen < 0.99){
                            int keptSeg = 0;
                            for(size_t k=0;k<iv.size();++k){
                                const double minLength = (sparseDisconnected[si] || materialGap)
                                    ? std::max(10.0, 2.0*pstep) : std::max(10.0, L*0.5);
                                if(iv[k].second - iv[k].first < minLength) continue;
                                Row q = r;
                                q.xs = iv[k].first; q.xe = iv[k].second;
                                q.segIdx = r.segIdx*8 + (int)k;
                                kept.push_back(q);
                                ++keptSeg;
                                cutLen += rowLen - (q.xe - q.xs);
                            }
                            if(keptSeg) ++nClip;
                            routeAroundVoid = routeAroundVoid || (materialGap && keptSeg>1);
                            continue;
                        }
                    }else if(cov / rowLen >= 0.40){
                        kept.push_back(r); continue;   // 基本在面内，不动
                    }
                    // 缩水严重：取重叠最长一段夹端点
                    size_t bi = 0; double bl = -1;
                    for(size_t k=0;k<iv.size();++k){
                        const double ov = std::min(iv[k].second, r.xe) - std::max(iv[k].first, r.xs);
                        if(ov > bl){ bl = ov; bi = k; }
                    }
                    Row q = r;
                    q.xs = std::max(r.xs, iv[bi].first);
                    q.xe = std::min(r.xe, iv[bi].second);
                    if(q.xe - q.xs < std::max(10.0, L*0.5)){ ++nDrop; continue; }
                    ++nClip; cutLen += rowLen - (q.xe - q.xs);
                    kept.push_back(q);
                }
                if(!kept.empty() && (nClip || nDrop))
                    fprintf(stderr, "[traj] 面%zu: 按重建轮廓夹端点 %d 行、丢弃 %d 行（行数 %zu -> %zu，%zu 环，夹掉 %.0fmm）\n",
                            si, nClip, nDrop, bestRows.size(), kept.size(), loops.size(), cutLen);
                if(!kept.empty() || sparseDisconnected[si]) bestRows.swap(kept);
            }
        }
        if(!ringMode && !loopMode && bestRows.empty()) continue;
        // 载体真的出了行，才允许"其余同带面不布行"；否则（被枪位/朝向/平面性门拒掉）
        // 保持原行为，各面各自布行 —— 宁可多喷，不能漏喷。
        if(!needsCandidateCheck && singleBand && si == singleBandCarrier &&
           (ringMode ? (nRing > 0) : (loopMode ? !loopPaths.empty() : !bestRows.empty())))
            carrierEmitted = true;

        const int faceRows = ringMode ? nRing : (loopMode ? (int)loopPaths.size() : (int)bestRows.size());
        if(!needsCandidateCheck) totalRows += faceRows;
        {
            double angN = std::acos(std::min(1.0,std::fabs(dot(pl.n,Vec3{0,0,1}))))*180.0/3.141592653589793;
            if(ringMode)
                fprintf(stderr, "[traj] 面%zu: 内点=%zu 行数=%d(同心圆) 与Z夹角=%.1f°\n",
                    si, inl.size(), nRing, angN);
            else if(loopMode)
                fprintf(stderr, "[traj] 面%zu: 内点=%zu 行数=%zu(%s) 与Z夹角=%.1f°\n",
                    si, inl.size(), loopPaths.size(),
                    bandMode ? "沿环带" : "沿轮廓环", angN);
            else
                fprintf(stderr, "[traj] 面%zu: 内点=%zu 行数=%zu 与Z夹角=%.1f° 行方向=%.1f°\n",
                    si, inl.size(), bestRows.size(), angN, theta*180.0/3.141592653589793);
        }

        double ca=std::cos(theta), sa=std::sin(theta);
        auto rotInv=[&](Vec2 p)->Vec2{ return {ca*p.x - sa*p.y, sa*p.x + ca*p.y}; };

        // ---- 跨行外扩的枪位兜底 ----
        // 首/末行按 acrossExt 外扩出去后，还要过一遍 gunClear：内壁底部那条首行的
        // 枪位本来就被 gunClear 卡在临界（内点已按 20mm 裁过），再往外 3mm 就侵入
        // 安全余量。这里只把"外扩出来的那一条"按 1mm 往回退，退到安全为止，
        // 最多退光 acrossExt（回到点云边界）；行距因此会有 ≤3mm 的不均匀，可接受。
        if(!ringMode && !loopMode && acrossExt > 1e-9 && !bestRows.empty()){
            const double rmn = [&]{ double m=1e18; for(auto& q:P2) m=std::min(m, ca*q.x+sa*q.y); return m; }();
            const double rmx = [&]{ double m=-1e18; for(auto& q:P2) m=std::max(m, ca*q.x+sa*q.y); return m; }();
            auto rowGunOK=[&](const Row& r)->bool{
                const double xs3[3] = { r.xs, 0.5*(r.xs+r.xe), r.xe };
                for(int t=0;t<3;++t){
                    Vec2 w = rotInv(Vec2{xs3[t], r.yc});
                    Vec3 X = pl.c + u*w.x + v*w.y;
                    if(gunClearAt(X, pl.n, (int)si) < gunClearEff) return false;
                }
                return true;
            };
            const int maxStep = (int)std::ceil(acrossExt);
            int nAdj = 0;
            size_t i = 0;
            while(i < bestRows.size()){
                size_t j = i;
                while(j+1 < bestRows.size() && bestRows[j+1].segIdx == bestRows[i].segIdx) ++j;
                for(int side=0; side<2; ++side){
                    if(side==1 && j==i) break;
                    Row& r = (side==0) ? bestRows[i] : bestRows[j];
                    const double dir = (side==0) ? +1.0 : -1.0;   // 0=低边向内退，1=高边向内退
                    int k=0;
                    while(k < maxStep && !rowGunOK(r)){ r.yc += dir; ++k; }
                    if(k) ++nAdj;
                }
                i = j+1;
            }
            if(nAdj)
                fprintf(stderr, "[traj] 面%zu 跨行外扩: %d 条首/末行因枪位不足退回（点云边界 %.1f..%.1f）\n",
                        si, nAdj, rmn, rmx);
        }


        // ---- 行端点外扩（布行边界外扩策略，可配置）----
        // 病灶：点云在面边沿常常"缺一截"——扫描遮挡，或者边缘点被判给了相邻面——
        // 于是行端点落在比真实边沿更靠里的位置，边沿那一条带就漏喷。实测薄壁开口盒
        // 内壁顶部漏 9mm，就是行端点比 z=300 的顶沿短了一截。这和枪体安全半径
        // （gunClear）是两个独立机制：gunClear 管"这里该不该布行"，这里管"行该走多长"。
        // 做法：每条行的两端沿行方向各外扩 extRow mm，1mm 一步试探，每步都查枪位——
        // 外扩点抬 standoff 后若距件 < gunClear（再往外枪头就顶进相邻壁了）就地停住。
        // 于是"往空处外扩"照扩（开口盒内壁顶部朝空，扩满），"往实体里外扩"自动拦住
        // （内壁底部朝内底，扩一两毫米就停），不需要事先知道面的真实边界在哪。
        // 配置：params->edgeExtend，<0 = 自动 = 0.25×行距（行距 80 时 20mm），0 = 不外扩。
        // 自动值原来取 1×行距（80mm）：行端点整整多喷一个行距到空处，既浪费涂料，
        // 3D 视图里又是一圈"出框"的线，看着不像长在工件上的轨迹。半条喷幅以内
        // 就足够把边沿那一条带盖住（喷幅本身有 20~30% 搭接）。
        double extRow = params->edgeExtend >= 0.0 ? params->edgeExtend : 0.25*params->scanRowSpacing;
        if(sparseDisconnected[si] || routeAroundVoid) extRow = 0.0;
        // ---- 盲扩量还要受"面自身尺度"约束 ----
        // 行距是工艺参数（喷幅/搭接），与这块面有多大无关；拿它决定"往面外伸多长"
        // 在大行距下会离谱：242×40 的板用 rowSpacing=150，extRow=37.5mm，两端各伸
        // 37mm，行总长 242→304 —— 3D 里就是两条悬在空处的线，用户直接看成"轨迹偏了"。
        // 而遮挡缺口本身只有一个扫线间距的量级（现场薄壁盒顶部实测 9mm）。所以再加一道
        // 与面尺度挂钩的上限：每端盲扩 <= 面内行向跨度的 4%（下限 3mm 保证小面也有补边能力）。
        // 判据只依赖 P2 在行向的跨度，与零件尺寸、坐标系、行距都无关。
        {
            double gxLo=1e18, gxHi=-1e18;
            for(const auto& q : P2){ const double t=ca*q.x+sa*q.y; gxLo=std::min(gxLo,t); gxHi=std::max(gxHi,t); }
            const double spanRow = gxHi - gxLo;
            if(spanRow > 1e-6) extRow = std::min(extRow, std::max(3.0, 0.04*spanRow));
        }
        const bool extDbg = std::getenv("SPRAYTRAJ_CONN_DBG") != nullptr;
        if(!ringMode && !loopMode && extRow > 1e-9 && !bestRows.empty()){
            int nExt=0, nBlk=0; double extSum=0.0, extMax=0.0;
            for(auto& r : bestRows){
                for(int side=0; side<2; ++side){
                    const double sgn = (side==0) ? +1.0 : -1.0;   // +1 = xe 端，-1 = xs 端
                    double add = 0.0;
                    for(double d=1.0; d<=extRow+1e-9; d+=1.0){
                        Vec2 q2{ (sgn>0 ? r.xe : r.xs) + sgn*d, r.yc };
                        Vec2 w = rotInv(q2);
                        Vec3 X = pl.c + u*w.x + v*w.y;
                        const double clr = gunClearAt(X, pl.n, (int)si);
                        if(extDbg && (d<=3 || clr < gunClearEff+5))
                            fprintf(stderr, "[ext] 面%zu row yc=%.1f side=%d d=%.0f X=(%.1f,%.1f,%.1f) clr=%.1f\n",
                                    si, r.yc, side, d, X.x, X.y, X.z, clr);
                        if(clr < gunClearEff){ ++nBlk; break; }
                        add = d;
                    }
                    if(add > 0.0){
                        if(sgn>0) r.xe += add; else r.xs -= add;
                        ++nExt; extSum += add; extMax = std::max(extMax, add);
                    }
                }
            }
            if(nExt || nBlk)
                fprintf(stderr, "[traj] 面%zu 行端点外扩: %d 端成功(均值 %.1fmm / 最大 %.1fmm)，"
                                "%d 端起扩即被枪位拦停（上限 %.1fmm）\n",
                        si, nExt, nExt? extSum/nExt : 0.0, extMax, nBlk, extRow);
        }

        // 生成该平面的 2D 轨迹。元素 = (旋转系 2D 坐标, spray, 额外抬升 lift)：
        // 3D 枪位 = 面点 X(p2) + pl.n·(standoff + lift)。lift=0 即普通贴面轨迹点；
        // lift>0 用于"抬枪绕行"的过渡点（净空不足的连接段降级，见下方 connSafe）。
        // 任意 3D 点都能表达：h = dot(P-pl.c, pl.n)，p2 = rot(投影到面内 2D)，
        // lift = h - standoff —— (2D, lift) → 3D 是仿射映射，直线段仍是直线段。
        struct T2 { Vec2 p; int spray; double lift; int kind; };
        std::vector<T2> traj2;
        // kind：0=布行本体（含换行掉头，换行不标色）；1=换面不关枪（橙色）
        auto sampleLine=[&](Vec2 a, Vec2 b, int spray, int kind=0){
            double LL=len2(b-a);
            int np=std::max(2,(int)std::ceil(LL/pstep)+1);
            for(int t=0;t<np;++t){ double f=(double)t/(np-1); traj2.push_back({a+(b-a)*f, spray, 0.0, kind}); }
        };
        auto sampleArc=[&](Vec2 a, Vec2 b, Vec2 apex, int kind=0){
            const int M=60; std::vector<Vec2> crv(M+1);
            for(int i=0;i<=M;++i){ double t=(double)i/M; crv[i]=a*((1-t)*(1-t)) + apex*(2*(1-t)*t) + b*(t*t);}
            std::vector<double> s(M+1,0);
            for(int i=1;i<=M;++i) s[i]=s[i-1]+len2(crv[i]-crv[i-1]);
            double arc=s[M];
            int np=std::max(4,(int)std::ceil(arc/pstep));
            for(int i=0;i<=np;++i){
                double target=(double)i/np*arc;
                int lo=0,hi=M; while(lo<hi){int mid=(lo+hi)/2; if(s[mid]>=target)hi=mid;else lo=mid+1;}
                traj2.push_back({crv[lo], 1, 0.0, kind});
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
                        // 圈间过渡：kind=0（正常连续喷涂，不标色）。只有换面才标橙。
                        for(int t=1;t<nI;++t)
                            traj2.push_back({prevEnd + (first-prevEnd)*((double)t/(double)nI), 1, 0.0, 0});
                    }
                }
                for(int t=0;t<=np;++t){
                    int ti = rev ? (np - t) : t;
                    double aa = 2.0*3.14159265358979*(double)ti/(double)np;
                    Vec2 q{ ringC.x + rad*std::cos(aa), ringC.y + rad*std::sin(aa) };
                    traj2.push_back({toRow(q), 1, 0.0});
                }
                prevEnd = traj2.back().p; hasPrev = true;
            }
        }
        if(loopMode){
            // 沿轮廓环：每圈首尾相接，相邻圈反向走（呼应蛇形），圈间在对应顶点处短接——
            // 两圈只隔 bandW/(k+1) ≈ 8mm，这一段是正常连续喷涂（kind=0，不标色 —— 换行不需要换色）。
            Vec2 prevEnd{0,0}; bool hasPrev=false;
            for(size_t k=0;k<loopPaths.size();++k){
                const auto& lp = loopPaths[k];
                if(lp.size() < 2) continue;
                const bool openPath = (k < loopOpen.size() && loopOpen[k] != 0);
                // 闭合环：最后一条边是 lp[n-1]→lp[0]；开口弧：只累计到 lp[n-1] 为止。
                // 这里必须分开算，否则开口弧会被自动接回起点，在缺口上画一道直线。
                std::vector<double> cum(lp.size(), 0.0);
                for(size_t t=0;t+1<lp.size();++t)
                    cum[t+1] = cum[t] + len2(lp[t+1]-lp[t]);
                const double per = openPath
                    ? cum.back()
                    : cum.back() + len2(lp[0]-lp.back());
                cum.push_back(per);
                if(per < 4.0*pstep) continue;
                auto at=[&](double s)->Vec2{
                    if(s<0) s=0; if(s>per) s=per;
                    size_t k2=0; while(k2+1<cum.size() && cum[k2+1]<s) ++k2;
                    if(k2>=lp.size()) k2=lp.size()-1;
                    const double seg = cum[k2+1]-cum[k2];
                    const double f = seg>1e-9 ? (s-cum[k2])/seg : 0.0;
                    const Vec2 a=lp[k2], b=lp[(k2+1)%lp.size()];
                    return Vec2{a.x+(b.x-a.x)*f, a.y+(b.y-a.y)*f};
                };
                const int np = std::max(8, (int)std::ceil(per/pstep));
                const bool rev = ((k&1)!=0);
                const Vec2 first = rev ? at(per) : at(0.0);
                if(hasPrev){
                    const int nI = std::max(1, (int)std::ceil(len2(first-prevEnd)/pstep));
                    // 环间过渡：kind=0（正常连续喷涂，不标色）。
                    for(int t=1;t<nI;++t)
                        traj2.push_back({prevEnd + (first-prevEnd)*((double)t/(double)nI), 1, 0.0, 0});
                }
                for(int t=0;t<=np;++t){
                    const double s = per*(double)t/(double)np;
                    traj2.push_back({at(rev ? (per-s) : s), 1, 0.0});
                }
                prevEnd = rev ? at(0.0) : at(per);
                hasPrev = true;
            }
        }
        if(!ringMode && !loopMode)
        for(size_t i=0;i<bestRows.size();++i){
            Vec2 a{bestRows[i].xs, bestRows[i].yc}, b{bestRows[i].xe, bestRows[i].yc};
            if(i%2==1) std::swap(a,b);
            segs2.push_back({a,b});
        }
        const bool connDbg = std::getenv("SPRAYTRAJ_CONN_DBG") != nullptr;
        if(connDbg){
            fprintf(stderr, "[conn] 面%zu 行数=%zu turnR=%.1f L=%.1f\n", si, bestRows.size(), turnR, L);
            for(size_t i=0;i<bestRows.size();++i)
                fprintf(stderr, "[conn]   row%zu xs=%.1f xe=%.1f yc=%.1f seg=%d\n",
                        i, bestRows[i].xs, bestRows[i].xe, bestRows[i].yc, bestRows[i].segIdx);
        }
        // ---- 连接段净空检查与降级（第三十二轮）----
        // 病灶（实测开口盒，行距 80 / 折返系数 0.5）：末行被跨行外扩顶到板边外
        // 16mm、再被行端点外扩推出板边最远 80mm，蛇形连接段从"板内 20mm 处"连到
        // "板外 80mm 处"，在相邻壁的半高处直接穿墙——枪位距壁最近 0.8mm。此前
        // 连接段（直线/圆弧）没有任何净空检查，行端点检查和布行前过滤都管不到它。
        // 处理：连接段逐点（≤2mm 步距）查枪位到工件点云的净空，不足时按代价从低
        // 到高降级：弧(原侧) → 弧(另一侧鼓出) → 直连 → 抬到全局安全高度绕行(spray=0)。
        // 前三档保持喷涂连续；只有确实绕不过去才抬枪断喷。
        auto clearMin3=[&](const Vec3& g)->double{
            // 只用于和 gunClearEff(-0.5) 比较；上界取 gunClearEff 保持判定一致。
            double d = partTree.nearestWithinExcept(g, gunClearEff, faceSkipMask[si]);
            // 与 gunClearAt 同口径：点云 + 面（投影在面内）+ 面边界线。只用点云时，
            // 棱边上的点一旦稀疏/被去噪删掉就会把"贴着壁面飞"判成安全。
            d = std::min(d, planeClearAt(g, (int)si));
            d = std::min(d, edgeClearAt(g, (int)si));
            return d;
        };
        auto clearOK3=[&](const Vec3& g)->bool{
            // 0.5mm 容差：枪位恰好落在安全带边界上（例如行端点正好距邻壁 20.0）时，
            // 浮点累积/降采样误差会让它变成 19.999 而误判——否则整面退化成抬枪绕行。
            return clearMin3(g) >= gunClearEff - 0.5;
        };
        auto gunOf2=[&](const Vec2& p2)->Vec3{
            Vec2 w = rotInv(p2);
            Vec3 X = pl.c + u*w.x + v*w.y;
            return X + pl.n*params->standoff;
        };
        auto polyOK=[&](const std::vector<Vec2>& ps)->bool{
            if(ps.empty()) return true;
            auto bad=[&](const char* why, const Vec3& g, double d){
                if(connDbg) fprintf(stderr,"[why] %s d=%.2f at (%.1f,%.1f,%.1f)\n",
                                     why, d, g.x,g.y,g.z);
                return false;
            };
            const double lim = gunClearEff - 0.5;
            for(size_t k=0;k<ps.size();++k){
                if(k>0){
                    Vec2 dd=ps[k]-ps[k-1]; double dl=len2(dd);
                    int ns=(int)std::ceil(dl/(pstep*2.0));
                    Vec3 prevG = gunOf2(ps[k-1]);
                    for(int t=1;t<ns;++t){
                        Vec2 m=ps[k-1]+dd*((double)t/(double)ns);
                        Vec3 g = gunOf2(m);
                        if(segCrossFace(prevG, g)) return bad("穿面", g, 0.0);
                        double d = clearMin3(g);
                        if(d < lim) return bad("净空", g, d);
                        prevG = g;
                    }
                    Vec3 gEnd = gunOf2(ps[k]);
                    if(segCrossFace(prevG, gEnd)) return bad("穿面", gEnd, 0.0);
                    double d2 = clearMin3(gEnd);
                    if(d2 < lim) return bad("净空", gEnd, d2);
                    continue;
                }
                double d0 = clearMin3(gunOf2(ps[k]));
                if(d0 < lim) return bad("净空", gunOf2(ps[k]), d0);
            }
            return true;
        };
        auto arcPoly=[&](Vec2 a, Vec2 b, Vec2 apex)->std::vector<Vec2>{
            const int M=24; std::vector<Vec2> c(M+1);
            for(int k=0;k<=M;++k){ double t=(double)k/M; c[k]=a*((1-t)*(1-t))+apex*(2*(1-t)*t)+b*(t*t); }
            return c;
        };

        bool alignedZigzag = false;
        if(needsCandidateCheck && !sparseDisconnected[si] && !routeAroundVoid && segs2.size() > 1){
            // 矩形单片的行端收齐到共有范围。保留原行距和蛇形顺序，勿跨孔或延长扫描行。
            double left = -1e30, right = 1e30, outerLeft = 1e30, outerRight = -1e30;
            double minSpacing = L;
            bool alternating = true;
            for(size_t i=0; i<segs2.size(); ++i){
                const auto& row = segs2[i];
                const double lo = std::min(row.a.x,row.b.x), hi = std::max(row.a.x,row.b.x);
                left = std::max(left,lo); right = std::min(right,hi);
                outerLeft = std::min(outerLeft,lo); outerRight = std::max(outerRight,hi);
                if(i > 0){
                    const auto& prev = segs2[i-1];
                    minSpacing = std::min(minSpacing, row.a.y - prev.a.y);
                    alternating = alternating && (prev.b.x-prev.a.x)*(row.b.x-row.a.x) < 0.0;
                }
            }
            const double endTolerance = 0.25*minSpacing;
            alignedZigzag = alternating && minSpacing > 0.0 && right-left >= std::max(10.0,L*0.5) &&
                            std::max(left-outerLeft,outerRight-right) <= endTolerance;
            if(alignedZigzag){
                for(const auto& row : segs2){
                    if(!candidateSegmentClear(gunOf2({left,row.a.y}),gunOf2({right,row.a.y}),gunClearEff)){
                        alignedZigzag = false;
                        break;
                    }
                }
            }
            if(alignedZigzag){
                for(size_t i=1; i<segs2.size(); ++i){
                    const double x = segs2[i-1].b.x > segs2[i-1].a.x ? right : left;
                    const Vec2 a{x,segs2[i-1].a.y}, b{x,segs2[i].a.y};
                    if(!polyOK({a,b}) || !candidateSegmentClear(gunOf2(a),gunOf2(b),gunClearEff)){
                        alignedZigzag = false;
                        break;
                    }
                }
            }
            if(alignedZigzag){
                for(auto& row : segs2){
                    const bool forward = row.b.x > row.a.x;
                    row.a.x = forward ? left : right;
                    row.b.x = forward ? right : left;
                }
                fprintf(stderr, "[traj] 面%zu: %zu 条蛇形扫描行端点对齐，横向直线换行\n", si, segs2.size());
            }
        }

        // ---- 行端点硬约束裁剪（gunClear 是硬约束，不是事后补救）----
        // 布行过滤发生在"面点沿法向抬 standoff"阶段，但有三种情况会让**实际发射点**
        // 偏离被检查过的位置：①补全面是矩形外推、不是真实点；②平面拟合有残差
        // （complex 面4 RMS 9.2mm，实际是斜曲面）；③行端点外扩。于是出现过
        // "过滤说合格、发射时离工件 1.2mm"的现场（complex 面4 实测 823 点）。
        // 这里在发射前逐点复核，只保留**连续合格的最长区间**：
        //   · 该行整段合格 → 原样保留（绝大多数情况，不改任何已有结果）
        //   · 部分合格 → 裁掉两端的不可达段（不足半行距就整行丢掉）
        //   · 全不合格 → 整行丢弃，宁可少喷一行，也不生成贴着工件走的轨迹
        if(!alignedZigzag && !ringMode && !segs2.empty()){
            const double lim = gunClearEff - 0.5;
            std::vector<Seg2> kept; kept.reserve(segs2.size());
            int nCut=0, nDrop=0;
            for(const auto& sg : segs2){
                const double Ls = len2(sg.b - sg.a);
                if(Ls < 1e-6){ kept.push_back(sg); continue; }
                const double checkStep = needsCandidateCheck ? std::min(pstep, 2.0) : pstep;
                const int ns = std::max(2, (int)std::ceil(Ls/checkStep));
                std::vector<char> okv((size_t)ns+1, 1);
                bool allOK = true;
                int nFail=0;
                for(int k=0;k<=ns;++k){
                    const Vec3 g = gunOf2(sg.a + (sg.b-sg.a)*((double)k/(double)ns));
                    if(clearMin3(g) < lim ||
                       (needsCandidateCheck && candidateHit(g, seg.preservedCarrier) == RayTarget::Blocked)){
                        okv[k]=0; allOK=false;
                        if(connDbg && nFail<2){
                            fprintf(stderr,"[gun]   失败点 #%d 枪位 (%.1f,%.1f,%.1f) d=%.2f 面法向 (%.2f,%.2f,%.2f) 面心 (%.1f,%.1f,%.1f)\n",
                                    k, g.x,g.y,g.z, clearMin3(g), pl.n.x,pl.n.y,pl.n.z,
                                    pl.c.x,pl.c.y,pl.c.z);
                        }
                        ++nFail;
                    }
                }
                if(allOK){ kept.push_back(sg); continue; }
                int bestA=-1, bestB=-1, cur=-1;
                for(int k=0;k<=ns;++k){
                    if(okv[k]){ if(cur<0) cur=k; if(bestA<0 || (k-cur) > (bestB-bestA)){ bestA=cur; bestB=k; } }
                    else cur=-1;
                }
                int nOK=0; for(int k=0;k<=ns;++k) nOK += okv[k]?1:0;
                if(connDbg)
                    fprintf(stderr, "[gun]   行 L=%.1f 合格 %d/%d 最长段 %d\n", Ls, nOK, ns+1,
                            (bestA<0?0:(bestB-bestA+1)));
                if(bestA<0){ ++nDrop; continue; }
                Seg2 q{ sg.a + (sg.b-sg.a)*((double)bestA/(double)ns),
                        sg.a + (sg.b-sg.a)*((double)bestB/(double)ns) };
                if(len2(q.b-q.a) < std::max(L*0.5, 10.0)){ ++nDrop; continue; }
                kept.push_back(q); ++nCut;
            }
            if(nCut || nDrop){
                segs2.swap(kept);
                fprintf(stderr, "[gun] 面%zu: 行端点裁剪 %d 行 / 整行丢弃 %d 行（枪位硬约束 %.0fmm）\n",
                        si, nCut, nDrop, gunClearEff);
            }
        }
        // ---- 转角半径自适应 ----
        // gunClear 是**硬约束**：轨迹必须始终离工件 ≥ 该距离，不能因为达不到就抬枪断喷
        // ——那是把"安全距离"当成了"放弃喷涂"的借口。正确做法是让转角收进安全域：
        // 行端点已经按 gunClear(+1mm) 收缩过，落在安全域边界附近；折返弧向行外侧鼓
        // turnR（=0.5×行距，行距 80 时 40mm），在腔内根本鼓不出去（端点距邻壁只剩
        // 20.8mm，鼓 40mm 就穿墙 19mm），于是旧代码判定失败 → 抬枪绕行 → 轨迹碎片化。
        // 这里改为在 [0, turnR] 上二分求**最大的可行半径**：能鼓多少鼓多少（半径越大
        // 越平滑），鼓不动就一路缩到 0（R=0 即两个端点的中心直线，仍是连续喷涂）。
        auto arcR=[&](Vec2 a, Vec2 b, Vec2 cr)->double{
            const Vec2 mid=(a+b)*0.5;
            if(polyOK(arcPoly(a,b,mid+cr*turnR))) return turnR;   // 满半径可过，最省
            double lo=0.0, hi=turnR;
            for(int it=0; it<6; ++it){                            // 6 次 → 分辨率 turnR/64
                const double md=0.5*(lo+hi);
                if(polyOK(arcPoly(a,b,mid+cr*md))) lo=md; else hi=md;
            }
            return lo;                                            // 0 = 只剩直连可用
        };
        // 抬到安全高度绕行：A → A 的抬升位 → B 的抬升位 → B，全部 spray=0。
        // 只用于**空移段**（跨段大跳跃/面间过渡）：那本来就不喷，抬起来不打断喷涂。
        // 喷涂段内部（同一面的换行连接）绝不用它——gunClear 不足要靠收缩转角解决。
        // 安全高度在工件沿 liftU 的最远点之上，横向平移段必然安全；两段起降与面间过渡同款。
        // 抬升位置由 liftTo 给出：liftU=+Z 时就是旧的 {x, y, zSafe}。
        auto hopZ=[&](const Vec2& a2, const Vec2& b2){
            Vec3 A = gunOf2(a2), B = gunOf2(b2);
            Vec3 P1 = liftTo(A, sSafe), P2 = liftTo(B, sSafe);
            auto emit3=[&](const Vec3& P, int sp){
                const double h = dot(P - pl.c, pl.n);
                Vec3 X = P - pl.n*h;
                Vec3 qq = X - pl.c;
                Vec2 q2{dot(qq,u), dot(qq,v)};
                Vec2 p2{ca*q2.x + sa*q2.y, -sa*q2.x + ca*q2.y};
                traj2.push_back({p2, sp, h - params->standoff});
            };
            emit3(P1,0); emit3(P2,0); emit3(B,0);
        };
        if(routeAroundVoid){
            // 包含真孔的面统一按同一轮廓收尾，避免未被切段的首末行停在轮廓外。
            std::vector<Seg2> inside;
            for(const auto& seg : segs2){
                const bool reverse = seg.a.x > seg.b.x;
                const auto iv=clipRowByLoops(seg.a.y,std::min(seg.a.x,seg.b.x),
                                             std::max(seg.a.x,seg.b.x),rowFootprint);
                for(size_t k=0; k<iv.size(); ++k){
                    const auto& span=iv[reverse ? iv.size()-1-k : k];
                    const Vec2 a{span.first,seg.a.y}, b{span.second,seg.a.y};
                    inside.push_back(reverse ? Seg2{b,a} : Seg2{a,b});
                }
            }
            segs2.swap(inside);
        }
        auto inFootprint = [](Vec2 p, const std::vector<Vec2>& loop){
            if(pointInPoly(p,loop)) return true;
            for(size_t i=0; i<loop.size(); ++i){
                const Vec2 a=loop[i], d=loop[(i+1)%loop.size()]-a;
                const double dd=dot2(d,d);
                if(dd<1e-12) continue;
                const double t=std::clamp(dot2(p-a,d)/dd,0.0,1.0);
                if(len2(p-(a+d*t))<1e-5) return true;
            }
            return false;
        };
        auto componentAt = [&](Vec2 p){
            for(size_t i=0; i<rowFootprint.size(); ++i)
                if(inFootprint(p,rowFootprint[i])) return (int)i;
            return -1;
        };
        auto connectOnPlane = [&](Vec2 a, Vec2 b, const std::vector<Vec2>& loop){
            auto valid = [&](Vec2 p, Vec2 q){
                const int n=std::max(1,(int)std::ceil(len2(q-p)/std::max(0.5,pstep*0.5)));
                for(int k=1; k<n; ++k)
                    if(!inFootprint(p+(q-p)*((double)k/n),loop)) return false;
                return polyOK({p,q});
            };
            if(valid(a,b)){ sampleLine(a,b,1); return true; }
            std::vector<Vec2> nodes{a,b};
            nodes.insert(nodes.end(),loop.begin(),loop.end());
            std::vector<double> distance(nodes.size(),1e30);
            std::vector<int> parent(nodes.size(),-1);
            std::vector<char> visited(nodes.size(),0);
            distance[0]=0.0;
            for(size_t k=0; k<nodes.size(); ++k){
                int i=-1;
                for(size_t j=0; j<nodes.size(); ++j)
                    if(!visited[j] && distance[j]<1e30 && (i<0 || distance[j]<distance[i])) i=(int)j;
                if(i<0 || i==1) break;
                visited[i]=1;
                for(size_t j=1; j<nodes.size(); ++j){
                    const double d=distance[i]+len2(nodes[j]-nodes[i]);
                    if(visited[j] || d>=distance[j] || !valid(nodes[i],nodes[j])) continue;
                    distance[j]=d; parent[j]=i;
                }
            }
            if(parent[1]<0) return false;
            std::vector<int> order;
            for(int i=1; i>=0; i=parent[i]) order.push_back(i);
            for(size_t k=order.size()-1; k>0; --k)
                sampleLine(nodes[order[k]],nodes[order[k-1]],1);
            return true;
        };
        int nHop=0, nShrink=0;
        bool missingConnection = false;
        for(size_t i=0;i<segs2.size();++i){
            sampleLine(segs2[i].a, segs2[i].b, 1);
            if(i+1<segs2.size()){
                double gap = len2(segs2[i+1].a - segs2[i].b);
                // 关枪判据和曲面块统一：一个行距级别的跨接（<= 3L）都还在工艺允许
                // 的重叠范围内，保持喷涂更省时间；只有真正的"大跳跃"（左右段切换、
                // 中间有孔洞）才关枪。原来只看 turnR*3（与行距无关），L=30 时
                // 一个 20mm 的普通换行也被关枪，白停一次。
                const double gapHard = std::max(turnR*3.0, 3.0*L);
                if(connDbg)
                    fprintf(stderr, "[conn]   conn%zu gap=%.1f hard=%.1f A=(%.1f,%.1f) B=(%.1f,%.1f)\n",
                            i, gap, gapHard, segs2[i].b.x, segs2[i].b.y, segs2[i+1].a.x, segs2[i+1].a.y);
                const Vec2 A = segs2[i].b, B = segs2[i+1].a;
                const bool strict = sparseDisconnected[si] || routeAroundVoid;
                if(alignedZigzag){
                    sampleLine(A, B, 1);
                } else if(strict || gap > gapHard){
                    const int from=componentAt(A), to=componentAt(B);
                    if(from>=0 && from==to){
                        if(connectOnPlane(A,B,rowFootprint[from])) continue;
                        missingConnection = true;
                        break;
                    }
                    if(strict && (from<0 || to<0)){
                        missingConnection = true;
                        break;
                    }
                    // 跨段大跳跃：本来就是空移（spray=0），先试直连，穿不过才抬高。
                    // 空移抬枪不打断喷涂，所以这里保留 hopZ 是合理的。
                    if(polyOK({A, B})) sampleLine(A, B, 0);
                    else { hopZ(A, B); ++nHop; }
                } else {
                    // ---- 喷涂段内部连接：必须连续，绝不抬枪 ----
                    // 转角半径能鼓多少鼓多少（受 gunClear 安全域约束），
                    // 鼓不动就缩到 0 → 退化为端点间直线，仍是 spray=1 的连续轨迹。
                    // kind=0：换行就是正常连续喷涂，不标色（用户 2026-11-20 明确：
                    // 「换面才需要换颜色，换行不需要」）。标成过渡色会让现场把
                    // 一行行的正常覆盖当成"没喷到位"，反而误导。
                    const Vec2 dir=norm2(A - segs2[i].a);
                    const Vec2 cr{-dir.y, dir.x};
                    double R = arcR(A, B, cr);
                    if(R >= 0.0){
                        if(R > 0.5) sampleArc(A, B, (A+B)*0.5 + cr*R, 0);
                        else        sampleLine(A, B, 1, 0);
                        if(R < turnR - 1e-6) ++nShrink;
                        if(connDbg)
                            fprintf(stderr, "[conn]     conn%zu 转角半径 %.1f/%0.1f\n", i, R, turnR);
                    } else {
                        // 另一侧鼓出（蛇形换行的另一条臂），同样 kind=0
                        double R2 = arcR(A, B, cr*(-1.0));
                        if(R2 > 0.5){ sampleArc(A, B, (A+B)*0.5 - cr*R2, 0); ++nShrink; }
                        else { sampleLine(A, B, 1, 0); ++nShrink; }
                    }
                }
            }
        }
        if(missingConnection){
            fprintf(stderr, "[skip] 面%zu 找不到满足净空的连续面内连接\n", si);
            continue;
        }
        if(nShrink)
            fprintf(stderr, "[traj] 面%zu: %d 处换行转角按枪体安全半径收缩（保持连续喷涂）\n",
                    si, nShrink);
        if(nHop)
            fprintf(stderr, "[traj] 面%zu: %d 处跨段空移抬高绕行（spray=0，非喷涂断点）\n", si, nHop);

        // 转 3D（用该平面自己的法向抬升，避免撞击凸起/侧壁）
        Vec3 gun = pl.n*(-1.0);

        if(needsCandidateCheck){
            const bool continuousRow = alignedZigzag || (seg.preservedCarrier && segs2.size() == 1);
            auto position = [&](const T2& t)->Vec3{
                const Vec2 q = rotInv(t.p);
                return pl.c + u*q.x + v*q.y + pl.n*(params->standoff + t.lift);
            };
            // 局部遮挡只关掉该段的喷枪，位置保持原样，关枪段同样复查净空。
            std::vector<unsigned char> blocked(traj2.size(), 0), hits(traj2.size(), 0);
            for(size_t k=0; k<traj2.size(); ++k){
                if(!traj2[k].spray) continue;
                const Vec3 g = position(traj2[k]);
                const bool fromSpray = k > 0 && traj2[k-1].spray;
                const Vec3 prev = fromSpray ? position(traj2[k-1]) : g;
                const int samples = std::max(1, (int)std::ceil(len(g-prev)/2.0));
                for(int j=1; j<=samples; ++j){
                    const Vec3 q = prev + (g-prev)*((double)j/samples);
                    const RayTarget hit = candidateHit(q, continuousRow);
                    if(hit == RayTarget::Blocked){
                        blocked[k] = 1;
                        if(fromSpray) blocked[k-1] = 1;
                        break;
                    }
                    if(hit == RayTarget::OnTarget) hits[k] = 1;
                }
            }
            size_t muted = 0;
            for(size_t k=0; k<traj2.size(); ++k){
                if(traj2[k].spray && blocked[k]){
                    traj2[k].spray = 0;
                    ++muted;
                }
            }
            bool targetHit = false;
            size_t k = 0;
            for(; k<traj2.size(); ++k){
                const size_t previous = k == 0 ? 0 : k-1;
                const bool spraySegment = traj2[k].spray && traj2[previous].spray;
                const double need = spraySegment ? gunClearEff : std::max(moveClear, gunClearEff);
                if(!candidateSegmentClear(position(traj2[previous]), position(traj2[k]), need)) break;
                targetHit = targetHit || (spraySegment && hits[k]);
            }
            if(k < traj2.size() || !targetHit){
                fprintf(stderr, "[skip] 面%zu 候选行复核未通过：%s（轨迹点 %zu）\n",
                        si, k < traj2.size() ? "枪位净空不足或穿面" : "未命中目标覆盖区", k);
                continue;
            }
            if(muted && continuousRow){
                fprintf(stderr, "[skip] 面%zu 连续蛇形候选未通过：%zu 点的首个命中或喷距不符\n", si, muted);
                continue;
            }
            if(muted)
                fprintf(stderr, "[gun] 面%zu 候选轨迹 %zu 点因首个命中或喷距不符而关枪\n", si, muted);
        }

        // ---- 端点选择（空移最短）：首/末两端，取"进得去、也出得来"的那端作起点 ----
        // 蛇形路线的两头都是合法起点：反走等于"行序与行向一起反过来"，走的是同一批
        // 线段（只是方向相反），所以净空/穿面检查的结论完全不变，不需要重算。
        // 原来固定从"第 0 行的 xs 端"进枪，上一面若正好停在另一端，就得先横穿整个
        // 面再开始喷 —— 实测 nc_complex 面3 因此凭空多走 273mm（空移 536mm）。
        //
        // ⚠ 但"入口近"不能单独决定方向：反走会把**出口**换到另一端，而出口未必走得掉。
        // 喷完本面后要"沿本面法向退 → 抬到安全面"，这一步只取决于出口那一端的几何，
        // 与下一面是谁无关，所以可以在这里先算。腔里的面尤其危险：box_traj 的两个
        // 内半环，1→0 这个方向出得去、反向出不去，只按"离上一面近"反走会让后面的
        // 面被判"无可行退出路径"整面丢掉（实测丢 2305 点）。
        // 规则：两端都出得去 → 选入口近的；只有一端出得去 → 必须让那一端当出口；
        //       两端都出不去 → 保持原方向（无论怎么反都会丢下一个面，不动更省事）。
        if(!allPts.empty() && traj2.size() > 1){
            const Vec3 pe{allPts.back().x, allPts.back().y, allPts.back().z};
            auto gunOfT=[&](const T2& t)->Vec3{
                Vec2 w = rotInv(t.p);
                Vec3 X = pl.c + u*w.x + v*w.y;
                return X + pl.n*(params->standoff + t.lift);
            };
            // 出口能力评分：与面间过渡"退出腿"同口径（法向退 + 抬到安全面的最小净空）
            auto exitScore=[&](const Vec3& e)->double{
                auto leg=[&](const Vec3& W)->double{
                    const Vec3 Wtop = liftTo(W, sSafe);
                    return std::min(segMinClear(e, W), segMinClear(W, Wtop));
                };
                double best = leg(e);
                for(double f : {0.25, 0.5, 0.75, 1.0}) best = std::max(best, leg(e + pl.n*(clear*f)));
                return best;
            };
            const Vec3 gF = gunOfT(traj2.front()), gB = gunOfT(traj2.back());
            const double dF = len(gF - pe), dB = len(gB - pe);
            const double need = 0.9*gunClearEff;
            const bool fOK = exitScore(gB) >= need;     // 正走：从末端出
            const bool bOK = exitScore(gF) >= need;     // 反走：从首端出
            bool rev = false;
            if(fOK && bOK) rev = (dB < dF - 1.0);
            else if(bOK)   rev = true;
            if(rev){
                std::reverse(traj2.begin(), traj2.end());
                if(std::getenv("SPRAYTRAJ_TRANS_DBG"))
                    fprintf(stderr, "[trans] 面%zu 反走：首端 %.0fmm/末端 %.0fmm（出口 fOK=%d bOK=%d）\n",
                            si, dF, dB, (int)fOK, (int)bOK);
            }
        }

        // 面间过渡：从上一面终点出发，抬枪绕行到本面首点。
        // 全部为 spray=0 空移点，WPF 与上位机据此断开喷涂段。
        const size_t nPtsBefore = allPts.size();   // 无可行进入路径时用来越界回滚整面
        bool approachHopeless = false;
        if(!allPts.empty() && !traj2.empty()){
            Vec3 pEnd{ allPts.back().x, allPts.back().y, allPts.back().z };
            Vec2 q0 = rotInv(traj2[0].p);
            Vec3 X0 = pl.c + u*q0.x + v*q0.y;
            Vec3 pStart = X0 + pl.n*params->standoff;

            // 姿态必须跟着"腿"走：退出腿仍是上一面姿态，进入腿才是本面姿态。
            // 原来所有过渡点一律用 gun（本面），退枪点就提前换了姿态 —— 开口盒子实测
            // 退枪段 0~100mm 内枪轴翻 90°~180°（#13654→#13655 位置不动、姿态翻 180°），
            // 真机上等于原地猛拧手腕。空移是关枪过渡，姿态应在安全高度上空慢慢转。
            auto pushMove=[&](const Vec3& p, const Vec3& nrm){
                allPts.push_back({p.x,p.y,p.z, -nrm.x,-nrm.y,-nrm.z, params->speed, 0});
            };
            // ① 直线可达就直接连（见 segClearOK 的说明）。
            //    spray 是"段属性"（两端都为 1 才是喷涂段），所以必须在中间插一个
            //    spray=0 把线断开；插中点而非端点：三点共线，物理路径不变，
            //    也不会在落点产生一个重复点（零长段，部分控制器会报错）。
            //    除点云净空外还要查**穿面**：点云稀疏/被去噪删掉的壁上，
            //    segClearOK 只看最近点会一路放行（complex 实测一段 204mm 直线穿过
            //    凸台，最差 0.38mm）。
            if(!needsCandidateCheck && !usedCandidateCheck &&
               tryShortFaceLink(pEnd, {-allPts.back().i,-allPts.back().j,-allPts.back().k},
                                pStart, pl.n, allPts)){
                // 短换面：枪不停，直接连过去（点已按点距铺好，kind=1 橙色过渡）
                // 计数统一放在最后的过渡重建里，这里不重复加（默认会重建）
                if(std::getenv("SPRAYTRAJ_TRANS_DBG"))
                    fprintf(stderr, "[trans] 面%zu 短换面不关枪 %.0fmm\n",
                            si, len(pStart - pEnd));
            } else if(segClearOK(pEnd, pStart) && !segCrossFace(pEnd, pStart)){
                if(std::getenv("SPRAYTRAJ_TRANS_DBG"))
                    fprintf(stderr, "[trans] 面%zu 直连 pEnd=(%.1f,%.1f,%.1f) pStart=(%.1f,%.1f,%.1f)\n",
                            si, pEnd.x,pEnd.y,pEnd.z, pStart.x,pStart.y,pStart.z);
                pushMove(pEnd + (pStart - pEnd)*0.5, prevNormal);
            } else {
                const Vec3 Top = liftTo(pStart, sSafe);
                auto legOK=[&](const Vec3& w0, const Vec3& w1){
                    return !segCrossFace(w0, w1) && segClearOK(w0, w1);
                };
                // ② 沿上一面法向退出（离开该面，不擦面）
                Vec3 L1 = pEnd + prevNormal*clear;
                const Vec3 L1top = liftTo(L1, sSafe);
                if(!legOK(pEnd, L1) || !legOK(L1, L1top)){
                    // 出口两条腿也必须查：腔宽与退出距离相当时，沿法向退 100mm
                    // 恰好顶到对面内壁，再抬枪就贴着壁顶擦上去（bin 实测 0.06mm）。
                    // 先试"原地抬枪"，再试缩短法向退出；取两腿最小净空最大者。
                    double exitScore = std::min(segMinClear(pEnd, L1), segMinClear(L1, L1top));
                    Vec3 exitW = L1;
                    auto considerExit=[&](const Vec3& W){
                        const Vec3 Wtop = liftTo(W, sSafe);
                        double s = std::min(segMinClear(pEnd, W), segMinClear(W, Wtop));
                        if(s > exitScore){ exitScore = s; exitW = W; }
                    };
                    considerExit(pEnd);
                    for(double f : {0.75, 0.5, 0.25}) considerExit(pEnd + prevNormal*(clear*f));
                    considerExit(pEnd + prevNormal*gunClearEff);
                    L1 = exitW;
                    if(std::getenv("SPRAYTRAJ_TRANS_DBG"))
                        fprintf(stderr, "[trans] 面%zu 出口改道 L1=(%.1f,%.1f,%.1f) 最小净空 %.1f\n",
                                si, L1.x,L1.y,L1.z, exitScore);
                    if(exitScore < 0.9*gunClearEff){
                        approachHopeless = true;
                        fprintf(stderr, "[skip] 面%zu 无可行退出路径（最佳 %.1fmm < 枪体安全 %.0fmm），整面不喷\n",
                                si, exitScore, gunClearEff);
                    }
                }
                pushMove(L1, prevNormal);              // ② 退出腿：保持上一面姿态（不再提前转向）
                // ③ 让到安全面（沿 liftU；liftU=+Z 时与旧行为一致）
                pushMove(liftTo(L1, sSafe), prevNormal);
                // ④ 在安全高度上平移到本面首点正上方，再降落。
                //    垂直降落会穿件的情形：首点正上方有凸台、或首点在腔内/凸台下
                //    （complex 实测 204mm 直降穿过工件，最差 0.38mm，且这一段
                //    spray=0→1，正是"落枪第一下就撞"）。改沿本面法向 n（指向可喷侧）
                //    偏移后降落，再沿法向进到首点。
                //    注意**两条腿都要查**：下降腿 Top→Side 和进入腿 Side→pStart。
                //    原来只查了下降腿，complex 上进入腿沿法向走 100mm，中途擦过
                //    一个凸台（实测最差 0.35mm，12 个点 <20mm），因为没人查。
                if(legOK(Top, pStart)){
                    pushMove(Top, pl.n);               // ④ 到顶空才转到本面姿态
                } else {
                    // 候选落点：沿法向由近到远 3 档 + 面内 8 个横向方位，逐条两腿验。
                    // 都不干净时取"最小净空最大"的那个 —— 宁可绕远，也不留一段 0.2mm
                    // 的擦壁直降（那是"落枪第一下就撞"）。
                    const double base = std::max(clear, gunClearEff);
                    Vec3 bestW = Top; double bestScore = -1.0;
                    auto consider=[&](const Vec3& W){
                        const Vec3 Wtop = liftTo(W, sSafe);
                        double s = std::min(segMinClear(Wtop, W), segMinClear(W, pStart));
                        if(s > bestScore){ bestScore = s; bestW = W; }
                    };
                    consider(Top);
                    for(double f : {0.25, 0.5, 1.0}) consider(pStart + pl.n*(base*f));
                    // 横向方位基：要和抬升轴正交，否则"面内 8 方位"里有两个方向
                    // 其实是沿 liftU 的 —— 抬上去再原地下降，白走一个来回。
                    // 退化保护：liftU ∥ pl.n 时 cross 为零，此时任取面内方向都行，
                    // 但必须**真的在面内**（与 pl.n 垂直）。原来这里回退到固定
                    // {1,0,0}，若 pl.n 恰好也是 ±X，inPlane 与法向平行、inPlane2
                    // 归零，8 个方位退化成同一个点，落点搜索直接失效。
                    Vec3 inPlane = cross(liftU, pl.n);
                    if(len(inPlane) < 1e-6) inPlane = cross(Vec3{0,0,1}, pl.n);
                    if(len(inPlane) < 1e-6) inPlane = cross(Vec3{1,0,0}, pl.n);
                    if(len(inPlane) < 1e-6) inPlane = Vec3{1,0,0};   // pl.n 不会是零向量，仅作兜底
                    inPlane = norm(inPlane);
                    const Vec3 inPlane2 = cross(pl.n, inPlane);
                    for(int a=0;a<8;++a){
                        const double th = a*3.141592653589793/4.0;
                        consider(pStart + (inPlane*std::cos(th) + inPlane2*std::sin(th))*base);
                    }
                    const Vec3 Wtop = liftTo(bestW, sSafe);
                    pushMove(Wtop, pl.n);
                    if(len(bestW - pStart) > 1e-6) pushMove(bestW, pl.n);
                    if(std::getenv("SPRAYTRAJ_TRANS_DBG"))
                        fprintf(stderr, "[trans] 面%zu 降落改道 W=(%.1f,%.1f,%.1f) 最小净空 %.1f\n",
                                si, bestW.x, bestW.y, bestW.z, bestScore);
                    ++nApproachFix;
                    // 连"最好的"落法都撑不到安全半径 → 这个面根本进不去（窄缝/深腔）。
                    // 与其生成一段 1.3mm 的擦壁进入（落枪就撞），不如整面不喷，
                    // 与"gunClear 是生成期硬约束"一致；日志里明说，便于工艺端补喷。
                    if(bestScore < 0.9*gunClearEff){
                        approachHopeless = true;
                        fprintf(stderr, "[skip] 面%zu 无可行进入路径（最佳 %.1fmm < 枪体安全 %.0fmm），整面不喷\n",
                                si, bestScore, gunClearEff);
                    }
                }
            }
            ++transCount;
        }
        if(approachHopeless){
            allPts.resize(nPtsBefore);      // 回滚该面已写入的过渡点
            ++nDropFace;
            continue;
        }
        if(needsCandidateCheck){
            usedCandidateCheck = true;
            const int emitted = (ringMode || loopMode) ? faceRows : (int)segs2.size();
            totalRows += emitted;
            if(singleBand && si == singleBandCarrier && emitted > 0) carrierEmitted = true;
            fprintf(stderr, "[gun] 面%zu 候选行复核通过，保留 %d 行\n", si, emitted);
        }
        prevNormal = pl.n;

        const size_t faceBeg = allPts.size();          // 本面面内轨迹的起点（供最后重排）
        for(auto& t2 : traj2){
            Vec2 q = rotInv(t2.p);
            Vec3 X = pl.c + u*q.x + v*q.y;
            Vec3 pos = X + pl.n*(params->standoff + t2.lift);
            allPts.push_back({pos.x,pos.y,pos.z, gun.x,gun.y,gun.z, params->speed, t2.spray, t2.kind});
        }
        // 同一拟合平面上可能有多个不相连的实体区域，分别保留它们的轨迹区间。
        for(size_t k=faceBeg; k<allPts.size();){
            if(!allPts[k].spray){ ++k; continue; }
            const size_t begin = k++;
            while(k<allPts.size() && allPts[k].spray) ++k;
            if(fanCoverOn) fanHarvest(begin,k,pl.n);
            if(!needsCandidateCheck) faceSpan.push_back({begin,k,si});
        }
        // 保留已复核的面内关枪区间，避免重排把关枪段重连成喷涂段。
        if(needsCandidateCheck) faceSpan.push_back({faceBeg,allPts.size(),si});
    }
    if(nGapSafe)
        fprintf(stderr, "[traj] 曲面行间安全空移: %d 处（直线会穿件，spray=0）\n", nGapSafe);
    if(transCount)
        fprintf(stderr, "[traj] 面间过渡: %d 处，安全高度 s=%.1f（工件沿抬升方向最远 %.1f）",
                transCount, sSafe, sPartMax);
    if(nApproachFix)
        fprintf(stderr, "，其中 %d 处直降会穿件→改沿面法向进入", nApproachFix);
    if(nDropFace)
        fprintf(stderr, "，%d 面因无可行进入路径整面不喷", nDropFace);
    if(transCount || nDropFace) fprintf(stderr, "\n");

    // 面间过渡：从 (pEnd,nEnd) 退枪、抬到安全面、绕到 (pStart,nStart) 上方落下。
    // 姿态必须跟着"腿"走：退出腿仍是上一面姿态，进入腿才是本面姿态 —— 原来所有过渡点
    // 一律用本面姿态，开口盒子实测退枪段 0~100mm 内枪轴翻 90°~180°（位置不动、姿态翻
    // 180°），真机上等于原地猛拧手腕。空移是关枪过渡，姿态应在安全高度上空慢慢转。
    // 返回值：0=可行，1=退出腿不可行，2=进入腿不可行（后两者按"整面不喷"处理）。
    // nFixOut 只在真正落盘时累加，预估代价时传 nullptr。
    auto buildTrans = [&](const Vec3& pEnd, const Vec3& nEnd,
                          const Vec3& pStart, const Vec3& nStart,
                          std::vector<SprayTrajPoint>& out, int* nFixOut)->int{
        out.clear();
        int code = 0;
        auto pushMove=[&](const Vec3& p, const Vec3& nrm){
            out.push_back({p.x,p.y,p.z, -nrm.x,-nrm.y,-nrm.z, params->speed, 0});
        };
        // ⓪ 短换面：枪位一跳就到、直线不擦件 → 不停枪，直接连过去。
        //    kind=1 让查看器画橙色（"这是换面过渡，不是布行，但没停枪"）。
        //    代价预估那一趟 nFixOut=nullptr，不计数。
        if(!usedCandidateCheck && tryShortFaceLink(pEnd, nEnd, pStart, nStart, out)){
            if(nFixOut) ++nGunOnLink;
            return 0;
        }
        // ① 直线可达就直接连（见 segClearOK 的说明）。spray 是"段属性"（两端都为 1
        //    才是喷涂段），所以必须在中间插一个 spray=0 把线断开；插中点而非端点：
        //    三点共线、物理路径仍是这条直线，也不会在落点产生一个零长段。
        //    除点云净空外还要查**穿面**：点云稀疏/被去噪删掉的壁上 segClearOK 只看
        //    最近点会一路放行（complex 实测一段 204mm 直线穿过凸台，最差 0.38mm）。
        if(segClearOK(pEnd, pStart) && !segCrossFace(pEnd, pStart)){
            pushMove(pEnd + (pStart - pEnd)*0.5, nEnd);
            return 0;
        }
        const Vec3 Top = liftTo(pStart, sSafe);
        auto legOK=[&](const Vec3& w0, const Vec3& w1){
            return !segCrossFace(w0, w1) && segClearOK(w0, w1);
        };
        // ② 沿上一面法向退出（离开该面，不擦面）
        Vec3 L1 = pEnd + nEnd*clear;
        const Vec3 L1top = liftTo(L1, sSafe);
        if(!legOK(pEnd, L1) || !legOK(L1, L1top)){
            // 出口两条腿也必须查：腔宽与退出距离相当时，沿法向退 100mm 恰好顶到
            // 对面内壁，再抬枪就贴着壁顶擦上去（bin 实测 0.06mm）。
            // 先试"原地抬枪"，再试缩短法向退出；取两腿最小净空最大者。
            double exitScore = std::min(segMinClear(pEnd, L1), segMinClear(L1, L1top));
            Vec3 exitW = L1;
            auto considerExit=[&](const Vec3& W){
                const Vec3 Wtop = liftTo(W, sSafe);
                double s = std::min(segMinClear(pEnd, W), segMinClear(W, Wtop));
                if(s > exitScore){ exitScore = s; exitW = W; }
            };
            considerExit(pEnd);
            for(double f : {0.75, 0.5, 0.25}) considerExit(pEnd + nEnd*(clear*f));
            considerExit(pEnd + nEnd*gunClearEff);
            L1 = exitW;
            if(exitScore < 0.9*gunClearEff) code = 1;
        }
        pushMove(L1, nEnd);                        // ② 退出腿：保持上一面姿态
        pushMove(liftTo(L1, sSafe), nEnd);         // ③ 让到安全面
        // ④ 在安全高度上平移到本面首点正上方，再降落。
        //    垂直降落会穿件的情形：首点正上方有凸台、或首点在腔内/凸台下
        //    （complex 实测 204mm 直降穿过工件，最差 0.38mm，且这一段
        //    spray=0→1，正是"落枪第一下就撞"）。改沿本面法向 n 偏移后降落，
        //    再沿法向进到首点。**两条腿都要查**：下降腿 Top→Side 和进入腿
        //    Side→pStart —— 原来只查了下降腿，complex 上进入腿沿法向走 100mm
        //    中途擦过一个凸台（最差 0.35mm，12 个点 <20mm），因为没人查。
        if(legOK(Top, pStart)){
            pushMove(Top, nStart);                 // ④ 到顶空才转到本面姿态
        } else {
            // 候选落点：沿法向由近到远 3 档 + 面内 8 个横向方位，逐条两腿验。
            // 都不干净时取"最小净空最大"的那个 —— 宁可绕远，也不留一段 0.2mm
            // 的擦壁直降（那是"落枪第一下就撞"）。
            const double base = std::max(clear, gunClearEff);
            Vec3 bestW = Top; double bestScore = -1.0;
            auto consider=[&](const Vec3& W){
                const Vec3 Wtop = liftTo(W, sSafe);
                double s = std::min(segMinClear(Wtop, W), segMinClear(W, pStart));
                if(s > bestScore){ bestScore = s; bestW = W; }
            };
            consider(Top);
            for(double f : {0.25, 0.5, 1.0}) consider(pStart + nStart*(base*f));
            // 横向方位基：要和抬升轴正交，否则"面内 8 方位"里有两个方向
            // 其实是沿 liftU 的 —— 抬上去再原地下降，白走一个来回。
            // 退化保护：liftU ∥ nStart 时 cross 为零，此时任取面内方向都行，
            // 但必须**真的在面内**（与 nStart 垂直）。原来这里回退到固定
            // {1,0,0}，若 nStart 恰好也是 ±X，inPlane 与法向平行、inPlane2
            // 归零，8 个方位退化成同一个点，落点搜索直接失效。
            Vec3 inPlane = cross(liftU, nStart);
            if(len(inPlane) < 1e-6) inPlane = cross(Vec3{0,0,1}, nStart);
            if(len(inPlane) < 1e-6) inPlane = cross(Vec3{1,0,0}, nStart);
            if(len(inPlane) < 1e-6) inPlane = Vec3{1,0,0};   // nStart 不会是零向量，仅作兜底
            inPlane = norm(inPlane);
            const Vec3 inPlane2 = cross(nStart, inPlane);
            for(int a=0;a<8;++a){
                const double th = a*3.141592653589793/4.0;
                consider(pStart + (inPlane*std::cos(th) + inPlane2*std::sin(th))*base);
            }
            const Vec3 Wtop = liftTo(bestW, sSafe);
            pushMove(Wtop, nStart);
            if(len(bestW - pStart) > 1e-6) pushMove(bestW, nStart);
            if(nFixOut) ++(*nFixOut);
            // 连"最好的"落法都撑不到安全半径 → 这个面根本进不去（窄缝/深腔）。
            // 与其生成一段 1.3mm 的擦壁进入（落枪就撞），不如整面不喷，
            // 与"gunClear 是生成期硬约束"一致；日志里明说，便于工艺端补喷。
            if(bestScore < 0.9*gunClearEff) code = 2;
        }
        return code;
    };

    // 过渡折线的实际长度（含最后一段"最后一个航点 → 面首点"）。
    auto transLen=[&](const Vec3& pEnd, const Vec3& pStart,
                      const std::vector<SprayTrajPoint>& out)->double{
        double s = 0.0; Vec3 c = pEnd;
        for(const auto& p : out){
            Vec3 q{p.x,p.y,p.z};
            s += len(q - c); c = q;
        }
        return s + len(pStart - c);
    };

    // ================= 5.面顺序重排 + 面间空移重建（最后处理） =================
    // 到这里所有面都已经按自己的几何生成完了，allPts 里躺着的是"面内轨迹 + 面间过渡"。
    // 生成阶段的面间过渡是按 split_planes 的发现顺序接的：面按点数/面积降序排下来，
    // 与空间邻接无关，于是空移在工件两头来回穿梭（实测 nc_complex 面0(-63)→面2(+101)
    // →面3(-212)→面4(-121)，X 方向折三趟）。
    //
    // 这一阶段只做一件事：把这批面按"真实过渡代价"重新排序、按新顺序重建过渡。
    // 为什么必须等到这里：过渡的代价不是一个距离，而是"退枪 → 抬到安全面 → 平移 → 落枪"，
    // 里面还含"这一头出不出得去、那一头进不进得来"。用面心 / 包围盒角点这类几何代理去估，
    // 实测在 box_traj 上空移反而 +34%（5171→6921）—— 代理估错一次就是横穿工件的一趟。
    // 现在每个面的真实首末枪位就在 allPts 里（faceSpan 记的区间），用 buildTrans 生成
    // **真实过渡折线**量长度，代价就是同一量纲、可以直接比。
    //
    // 顺序的求法：开放式路径（起点、终点都自由），每个面还有两个朝向（正走 / 反走）。
    // 反走 = 整块点序取反：每个点的姿态是它所在位置的局部法向，取反后仍然正确；行间
    // 安全空移的航点一起取反，走的仍是同一批线段（方向相反），净空/穿面结论不变。
    // 面数 ≤ 13 时按 (mask, 面, 朝向) 做精确 DP；更多时退化为多起点最近邻。
    {
        const int N = (int)faceSpan.size();
        if(N >= 2 && !std::getenv("SPRAYTRAJ_NO_REORDER")){
            auto PT=[&](size_t k){ const auto& p = allPts[k]; return Vec3{p.x,p.y,p.z}; };
            auto NR=[&](size_t k){ const auto& p = allPts[k]; return Vec3{-p.i,-p.j,-p.k}; };
            auto P0=[&](int i){ return PT(faceSpan[(size_t)i].beg); };
            auto N0=[&](int i){ return NR(faceSpan[(size_t)i].beg); };
            auto P1=[&](int i){ return PT(faceSpan[(size_t)i].end-1); };
            auto N1=[&](int i){ return NR(faceSpan[(size_t)i].end-1); };

            // ---- 代价矩阵：状态 = 面 × 朝向（2N 个），边权 = 真实过渡折线长度 ----
            // 不可行的过渡（退出/进入腿撑不到枪体安全半径）记 BIG，于是"腔里的面出不来"
            // 这类硬约束会被排到末尾（放在末尾就不需要退出），不必再手写"腔面收尾"。
            const int S2 = 2*N;
            const float BIGC = 1e5f;
            std::vector<float> C((size_t)S2*(size_t)S2, BIGC);
            std::vector<SprayTrajPoint> tmp;
            for(int i=0;i<N;++i) for(int oi=0;oi<2;++oi){
                const Vec3 pe = oi ? P0(i) : P1(i);
                const Vec3 ne = oi ? N0(i) : N1(i);
                for(int j=0;j<N;++j){
                    if(j==i) continue;
                    for(int oj=0;oj<2;++oj){
                        const Vec3 ps = oj ? P1(j) : P0(j);
                        const Vec3 ns = oj ? N1(j) : N0(j);
                        const int code = buildTrans(pe, ne, ps, ns, tmp, nullptr);
                        C[(size_t)(2*i+oi)*(size_t)S2 + (size_t)(2*j+oj)] =
                            (code==0) ? (float)std::min(transLen(pe, ps, tmp), 9e4) : BIGC;
                    }
                }
            }

            std::vector<int> ord; std::vector<char> ordRev;
            if(N <= 13){
                const int FULL = 1<<N;
                std::vector<float> f((size_t)FULL*(size_t)S2, 1e9f);
                std::vector<int>   pre((size_t)FULL*(size_t)S2, -1);
                for(int i=0;i<N;++i){
                    f[(size_t)(1<<i)*(size_t)S2 + 2*i]     = 0.0f;   // 起点自由
                    f[(size_t)(1<<i)*(size_t)S2 + 2*i + 1] = 0.0f;
                }
                for(int m=1;m<FULL;++m){
                    for(int st=0;st<S2;++st){
                        const float cur = f[(size_t)m*(size_t)S2 + (size_t)st];
                        if(cur >= 1e9f) continue;
                        for(int j=0;j<N;++j){
                            if(m & (1<<j)) continue;
                            const float* row = &C[(size_t)st*(size_t)S2];
                            for(int oj=0;oj<2;++oj){
                                const int t = 2*j+oj, nm = m | (1<<j);
                                const float nc = cur + row[t];
                                float& d = f[(size_t)nm*(size_t)S2 + (size_t)t];
                                if(nc < d){ d = nc; pre[(size_t)nm*(size_t)S2 + (size_t)t] = st; }
                            }
                        }
                    }
                }
                int bs = -1; float bv = 1e9f;
                for(int st=0;st<S2;++st){
                    const float v = f[(size_t)(FULL-1)*(size_t)S2 + (size_t)st];
                    if(v < bv){ bv = v; bs = st; }
                }
                int m = FULL-1, st = bs;
                while(st >= 0){
                    ord.push_back(st>>1); ordRev.push_back((char)(st & 1));
                    const int pst = pre[(size_t)m*(size_t)S2 + (size_t)st];
                    m ^= (1<<(st>>1));
                    st = pst;
                }
                std::reverse(ord.begin(), ord.end());
                std::reverse(ordRev.begin(), ordRev.end());
            } else {
                // 极端面数：多起点最近邻（边权已是真实代价），取最好的那条
                double bv = 1e30;
                for(int s0=0; s0<S2; ++s0){
                    std::vector<char> used((size_t)N, 0);
                    std::vector<int> t; std::vector<char> r;
                    int st = s0; double tot = 0.0;
                    t.push_back(st>>1); r.push_back((char)(st&1)); used[(size_t)(st>>1)] = 1;
                    for(int k=1;k<N;++k){
                        int bj = -1, bo = 0; float bc = 1e9f;
                        for(int j=0;j<N;++j){
                            if(used[(size_t)j]) continue;
                            for(int oj=0;oj<2;++oj){
                                const float c = C[(size_t)st*(size_t)S2 + (size_t)(2*j+oj)];
                                if(c < bc){ bc = c; bj = j; bo = oj; }
                            }
                        }
                        if(bj < 0) break;
                        tot += bc; st = 2*bj+bo;
                        t.push_back(bj); r.push_back((char)bo); used[(size_t)bj] = 1;
                    }
                    if(tot < bv){ bv = tot; ord = t; ordRev = r; }
                }
            }
            // 兜底：求解没覆盖全部面（理论上不会发生）就退回原顺序，绝不能让面凭空消失
            if((int)ord.size() != N){
                ord.clear(); ordRev.clear();
                for(int i=0;i<N;++i){ ord.push_back(i); ordRev.push_back(0); }
            }

            // ---- 旧顺序的面间空移：相邻两块之间那些过渡点的折线长 ----
            double airOld = 0.0;
            for(int i=1;i<N;++i){
                Vec3 c = P1(i-1);
                for(size_t k=faceSpan[(size_t)i-1].end; k<faceSpan[(size_t)i].beg; ++k){
                    const Vec3 q = PT(k);
                    airOld += len(q - c); c = q;
                }
                airOld += len(P0(i) - c);
            }

            // ---- 重建：按新顺序 [过渡] + [面内轨迹]，反走时连区间一起倒序 ----
            std::vector<SprayTrajPoint> rebuilt;
            rebuilt.reserve(allPts.size());
            std::vector<SprayTrajPoint> buf;
            Vec3 pe{0,0,0}, ne{0,0,1}; bool hasPrev = false;
            double airNew = 0.0; int nFix = 0, nBad = 0;
            for(size_t k=0;k<ord.size();++k){
                const int i = ord[k], o = ordRev[k];
                const size_t b = faceSpan[(size_t)i].beg, e = faceSpan[(size_t)i].end;
                const Vec3 eP = o ? P1(i) : P0(i);
                const Vec3 eN = o ? N1(i) : N0(i);
                if(hasPrev){
                    const int code = buildTrans(pe, ne, eP, eN, buf, &nFix);
                    airNew += transLen(pe, eP, buf);
                    if(code != 0){
                        // 代价表里已经尽量避开了；真避不开时按"最优落法"照走，
                        // 但要报出来（工艺端据此加过渡/改夹持），不静默丢面。
                        ++nBad;
                        fprintf(stderr, "[warn] 面%zu 过渡净空不足（%s），已按最佳落法生成\n",
                                faceSpan[(size_t)i].si, code==1 ? "上一面退不出来" : "本面进不去");
                    }
                    rebuilt.insert(rebuilt.end(), buf.begin(), buf.end());
                }
                if(o) rebuilt.insert(rebuilt.end(), allPts.rbegin() + (long)(allPts.size()-e),
                                     allPts.rbegin() + (long)(allPts.size()-b));
                else  rebuilt.insert(rebuilt.end(), allPts.begin() + (long)b, allPts.begin() + (long)e);
                pe = o ? P0(i) : P1(i);
                ne = o ? N0(i) : N1(i);
                hasPrev = true;
            }
            allPts.swap(rebuilt);

            std::string s;
            char b2[48];
            for(size_t k=0;k<ord.size();++k){
                std::snprintf(b2, sizeof(b2), k ? "→%zu%s" : "%zu%s",
                              faceSpan[(size_t)ord[k]].si, ordRev[k] ? "(反)" : "");
                s += b2;
            }
            fprintf(stderr, "[traj] 面顺序重排（真实过渡代价）: %s   面间空移 %.0f → %.0fmm",
                    s.c_str(), airOld, airNew);
            if(nFix) fprintf(stderr, "，%d 处降落改道", nFix);
            if(nBad) fprintf(stderr, "，%d 处净空不足", nBad);
            fprintf(stderr, "\n");
        }
    }

    // 换面短接统计：只有**最终重建**过的过渡才算数（上面 7900+ 那一趟会重建全部
    // 面间过渡，生成阶段写的那批被丢弃），所以这条日志必须放在重排之后。
    if(nGunOnLink)
        fprintf(stderr, "[traj] 换面短接: %d 处 <= %.0fmm，枪未关（spray=1，查看器画橙色）\n",
                nGunOnLink, params->crossFaceGunOnLen);

    // ================= 转角平滑 + 零长段清理（第十八轮） =================
    // 病灶（实测 bin.ply）：喷涂段内 >60° 的转角 56 处、>120° 的 23 处、最大 165°。
    // 两类来源：
    //   ① 往返式扫描的端点掉头：一条线走到底后**原路退回**，只横向错开 3.5mm
    //      （不是行距 12mm —— 平环上扫的是"弦"，两弦在圆周端点自然收敛）。
    //   ② 等高线圈的圈间过渡：本圈末点直接跳到下圈起点，方位差+径向偏移+z 跳变。
    // 两者的等效曲率半径都 < 1mm：200mm/s 时需要 > 80 m/s² 向心加速度，而喷涂
    // 机器人末端通常只有 10~15 m/s²。控制器只能二选一——
    //   · 降到 70~80mm/s 通过：膜厚近似正比 1/v，转角处厚 2~3 倍（流挂/色差）；
    //   · 用逼近过渡削角（KUKA $APO.CDIS / ABB zonedata）：轨迹整体偏 3~8mm，
    //     窄面上直接漏喷或撞相邻面，而算法这边毫不知情。
    // 处理：把尖角换成二次贝塞尔，控制点沿"更空"的一侧偏移 e，曲线从 A 鼓到 B。
    //   · 对端点掉头，这等价于工业上的 overspray 掉头（枪出界再绕回）；
    //   · 对圈间过渡，是把"跳过去"变成"滑过去"。
    // 安全性：偏移方向由 partTree 选"离工件更远"的一侧，即往空处绕，不会撞件；
    //         两侧都太贴件时放弃该点（保留原尖角，宁可尖也不撞）。
    {
        // ① 删零长段：实测 box 240 个、complex 62、cao 44、tr 42、bin 18。
        //    控制器的速度规划在这里会除零，KUKA/ABB 直接报 "motion not possible"。
        //    只删"与上一点 spray 状态相同"的，保留状态切换点，否则会破坏空移语义。
        std::vector<SprayTrajPoint> q;
        q.reserve(allPts.size());
        int nDup = 0;
        for(const auto& p : allPts){
            if(!q.empty()){
                const auto& r = q.back();
                double dx=p.x-r.x, dy=p.y-r.y, dz=p.z-r.z;
                if(dx*dx+dy*dy+dz*dz < 1e-12 && p.spray == r.spray){ ++nDup; continue; }
            }
            q.push_back(p);
        }

        // ② 密集区等弧长重采样
        // 现象：内孔/窄边处点数没随周长减少，点距坍缩到 0.03mm（实测 bin.ply
        //   idx673~685 连续 13 个点只走了 1.3mm，坐标在 ±0.05mm 上抖），
        //   切向完全由点云噪声决定 → 等效曲率半径 0.03~0.1mm。
        //   机器人其实"看不出"0.05mm 的位置抖动，但控制器的速度规划会把它当成
        //   "原地掉头"，直接降速或报位置超差。
        // ★ 这不是圆角能救的：横向错开 δ 只有 0.03mm 时，任何平滑曲线的掉头
        //   半径上限就是 δ/2。只能重采样 —— 把密集链按点距重新等分，采样点
        //   虽然还落在抖动曲线上，但密度降下来，等效曲率就回到正常量级。
        //   只动"连续点距 < 0.5·pstep 且链长 > 1.2·pstep"的密集链，首尾点保持
        //   原样，所以真实棱边（正常点距处的大转角）不会被削掉。
        int nRe = 0;
        {
            auto d3 = [](const SprayTrajPoint& a, const SprayTrajPoint& b){
                const double dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
                return std::sqrt(dx*dx+dy*dy+dz*dz);
            };
            std::vector<SprayTrajPoint> out;
            out.reserve(q.size());
            const double denseThr = 0.5*pstep;
            size_t i = 0;
            while(i < q.size()){
                size_t b = i;
                while(b+1 < q.size() && d3(q[b], q[b+1]) < denseThr) ++b;
                if(b > i){
                    double tot = 0.0;
                    for(size_t k=i;k<b;++k) tot += d3(q[k], q[k+1]);
                    if(tot > 1.2*pstep && (int)(b - i + 1) > 3){
                        out.push_back(q[i]);
                        double sAcc = 0.0, target = pstep;
                        for(size_t k=i;k<b;++k){
                            const double d = d3(q[k], q[k+1]);
                            if(d < 1e-12) continue;
                            while(sAcc + d >= target){
                                const double t = (target - sAcc)/d;
                                double ii=q[k].i + (q[k+1].i-q[k].i)*t;
                                double jj=q[k].j + (q[k+1].j-q[k].j)*t;
                                double kk=q[k].k + (q[k+1].k-q[k].k)*t;
                                const double nn=std::sqrt(ii*ii+jj*jj+kk*kk);
                                if(nn>1e-9){ ii/=nn; jj/=nn; kk/=nn; }
                                else       { ii=q[k].i; jj=q[k].j; kk=q[k].k; }
                                out.push_back({ q[k].x + (q[k+1].x-q[k].x)*t,
                                                q[k].y + (q[k+1].y-q[k].y)*t,
                                                q[k].z + (q[k+1].z-q[k].z)*t,
                                                ii, jj, kk, q[k].speed, q[k].spray, q[k].kind });
                                target += pstep;
                                ++nRe;
                            }
                            sAcc += d;
                        }
                        if(d3(out.back(), q[b]) > 1e-9) out.push_back(q[b]);
                        i = b + 1;
                        continue;
                    }
                }
                out.push_back(q[i]);
                ++i;
            }
            q.swap(out);
        }

        // ③ 探出点消除：行末点常被算到轮廓外一点点、下一行又从旁边 0.3mm 处起步，
        //    于是形成"出去 0.3mm 再回来"的无害往返 —— 转角 165°，但 A、B 几乎重合，
        //    等效掉头半径 ≤ 0.15mm，怎么平滑都是原地折返。直接删掉这个探出点、
        //    让前后直连。判据用 len(B-A)（往返幅度）而不是段长，因为重采样后
        //    段长已经是均匀的 pstep，靠段长看不出来。
        int nSpur = 0, nChamf = 0;
        // 迭代删除直到收敛：删掉一个拐点后，前后两点会接成新的短段，可能又构成
        // 新的发卡弯（实测单遍只抓到 1 处，实际有十几处）。每轮至少删一个就继续，
        // 6 轮足够收敛，也避免病态数据上死循环。
        for(int iter=0; iter<6 && q.size()>3; ++iter){
            std::vector<SprayTrajPoint> cur;
            cur.reserve(q.size());
            cur.push_back(q[0]);
            int removed = 0;
            for(size_t i=1;i+1<q.size();++i){
                const SprayTrajPoint& A=q[i-1];
                const SprayTrajPoint& P=q[i];
                const SprayTrajPoint& B=q[i+1];
                if(!(A.spray==1 && P.spray==1 && B.spray==1)){ cur.push_back(P); continue; }
                Vec3 a{P.x-A.x, P.y-A.y, P.z-A.z};
                Vec3 b{B.x-P.x, B.y-P.y, B.z-P.z};
                const double la=len(a), lb=len(b);
                if(la<1e-9 || lb<1e-9){ cur.push_back(P); continue; }
                const double abx=B.x-A.x, aby=B.y-A.y, abz=B.z-A.z;
                const double abLen=std::sqrt(abx*abx+aby*aby+abz*abz);
                double cg = std::max(-1.0, std::min(1.0, dot(a,b)/(la*lb)));
                const double angDeg = std::acos(cg)*57.29577951308232;
                if(angDeg > 150.0) fprintf(stderr,"[DBGJ] it=%d i=%d ang=%.1f ab=%.2f la=%.2f lb=%.2f\n",iter,(int)i,angDeg,abLen,la,lb);
                if(angDeg > 150.0 && abLen < 0.6*pstep){ ++nSpur; ++removed; continue; }

                // ④ 发卡弯倒角：行末的往返掉头，A、B 相距 1~3 个点距、转角 150~175°。
                //    这种弯的横向错开量 δ 就是 abLen，任何平滑曲线的掉头半径上限
                //    都是 δ/2，圆角救不了（试过三次贝塞尔，实测把 <3mm 从 500 处
                //    推到 606 处）。只能**删掉拐点让 A→B 直连**：转角立刻从 165°
                //    降到 ~40°，等效半径从 1.4mm 提到 5mm 以上。
                //    代价是拐点那一个点不喷（约一个点距 4mm），对 12mm 行距的覆盖
                //    重叠没有影响（实测覆盖率不变）。
                //    撞件校验：直连段中点离件不足 moveClear 就不倒角。
                if(angDeg > 150.0 && abLen < 3.0*pstep){
                    Vec3 mid{(A.x+B.x)*0.5, (A.y+B.y)*0.5, (A.z+B.z)*0.5};
                    KnnBuf kb; partTree.knn(mid, 1, kb);
                    const double dd = (kb.n>0)? std::sqrt(kb.d2[0]) : 1e9;
                    if(dd >= moveClear){ ++nChamf; ++removed; continue; }
                }
                cur.push_back(P);
            }
            cur.push_back(q.back());
            // 删除后可能出现重复点，顺手压掉（同 spray 才压，保状态切换点）
            std::vector<SprayTrajPoint> ded;
            ded.reserve(cur.size());
            for(const auto& p : cur){
                if(!ded.empty()){
                    const auto& r = ded.back();
                    const double dx=p.x-r.x, dy=p.y-r.y, dz=p.z-r.z;
                    if(dx*dx+dy*dy+dz*dz < 1e-12 && p.spray==r.spray) continue;
                }
                ded.push_back(p);
            }
            q.swap(ded);
            if(removed == 0) break;
        }
        if(nDup || nRe || nSpur || nChamf)
            fprintf(stderr, "[转角] 删零长段 %d 个，密集区重采样 %d 点，删探出点 %d 个，发卡弯倒角 %d 处\n",
                    nDup, nRe, nSpur, nChamf);
        allPts.swap(q);
    }

    // ================= 空移段整段重规划（第二十一轮，2026-10-29） =================
    // 用户口径：空移段不是"两点直线、交给控制器插补"，而是由本算法把整条路径插补出来，
    //   姿态沿弧长从上一面末点平滑转到下一面首点（每步转角有上限）；
    //   插补出来的点要能验安全（点云净空 + 不穿面），不安全就换一条途经重新验，
    //   全都不行才退回原路径并告警（留给 safety 层兜底）。
    // 候选途径（按代价从低到高）：直连 / 反向退出 / 抬升三档 / 左右各 3 档侧绕。
    // 绕行长度上限 = 原路径长度 × 2（用户 2026-10-29 拍板），超限视为不可达，不做无限绕远。
    {
        double mvStep = pstep;
        if(const char* e = std::getenv("SPRAYTRAJ_MOVE_STEP")){
            double v = atof(e);
            if(v > 0.5) mvStep = v;
        }
        const double kStepDeg   = 5.0;   // 空移相邻点姿态转角上限（度/点）
        double kMaxDetour = 2.0;         // 换途经的长度上限系数（用户 2026-10-29 拍板 2 倍）
        if(const char* e = std::getenv("SPRAYTRAJ_MOVE_DETOUR")){
            double v = atof(e);
            if(v >= 1.0) kMaxDetour = v;   // 现场可放宽（短段被 2 倍卡住时可调到 2.5~3）
        }
        // 几何安全门槛用「枪体安全半径」而不是 moveClear：moveClear=0.25×standoff=25mm 是
        // "直线可达"的宽松判据，而开口盒子这类内腔空移的几何极限就是 20mm 出头，
        // 用 25mm 当门槛会把本来合格的路径判成"无可行途经"（实测 4 段全误报）。
        // 现场可用 SPRAYTRAJ_MOVE_NEED 抬高门槛（越大越保守、越容易换途经）。
        double kSafeNeed = std::max(6.0, gunClearEff);
        if(const char* e = std::getenv("SPRAYTRAJ_MOVE_NEED")){
            double v = atof(e);
            if(v > 0.5) kSafeNeed = v;
        }

        auto pathLen = [&](const std::vector<Vec3>& v)->double{
            double s = 0.0;
            for(size_t k = 1; k < v.size(); ++k) s += len(v[k] - v[k-1]);
            return s;
        };
        // 整条折线的几何安全分：返回最小净空；任一段穿面直接 0。
        auto pathScore = [&](const std::vector<Vec3>& v)->double{
            double mn = 1e30;
            for(size_t k = 1; k < v.size(); ++k){
                if(segCrossFace(v[k-1], v[k])) return 0.0;
                mn = std::min(mn, segMinClear(v[k-1], v[k]));
            }
            return mn;
        };
        // 单位向量球面插值（姿态）：两端近似反向时绕固定轴转，避免过零退化。
        auto slerpU = [&](const Vec3& a, const Vec3& b, double t)->Vec3{
            double c = std::max(-1.0, std::min(1.0, dot(a, b)));
            const double th = std::acos(c);
            Vec3 r;
            if(th < 1e-6){
                r = a;
            } else if(th > 3.141592653589793 - 1e-3){
                Vec3 ax = norm(cross(a, (std::fabs(a.z) < 0.9) ? Vec3{0,0,1} : Vec3{1,0,0}));
                const double an = 3.141592653589793 * t;
                r = a*std::cos(an) + ax*std::sin(an);
            } else {
                const double s = std::sin(th);
                r = a*(std::sin((1.0 - t)*th)/s) + b*(std::sin(t*th)/s);
            }
            return norm(r);
        };
        // 沿折线按弧长插值；**逐段插值、强制保留每个折线顶点**。
        // 不能做"整段等距重采样"：等距点会跳过直角拐点，把"退出→抬升→平移→下降"
        // 折成一条斜线，等于把避障路径拉直（可能直接穿件）。
        auto resample = [&](const std::vector<Vec3>& v, double step,
                            std::vector<Vec3>& out, std::vector<double>& ts)->double{
            const size_t nv = v.size();
            std::vector<double> cum(nv, 0.0);
            for(size_t k = 1; k < nv; ++k) cum[k] = cum[k-1] + len(v[k] - v[k-1]);
            const double L = cum.back();
            if(L < 1e-9){ out.assign(1, v.front()); ts.assign(1, 0.0); return 0.0; }
            out.clear(); ts.clear();
            out.push_back(v.front()); ts.push_back(0.0);
            for(size_t k = 1; k < nv; ++k){
                const double sl = cum[k] - cum[k-1];
                if(sl < 1e-9){ out.push_back(v[k]); ts.push_back(cum[k] / L); continue; }
                const int m = std::max(1, (int)std::ceil(sl / std::max(0.5, step)));
                for(int t = 1; t <= m; ++t){
                    const double f = (double)t / (double)m;
                    out.push_back(v[k-1] + (v[k] - v[k-1]) * f);
                    ts.push_back((cum[k-1] + sl*f) / L);
                }
            }
            return L;
        };

        std::vector<SprayTrajPoint> merged;
        merged.reserve(allPts.size()*2);
        int nSeg = 0, nReplan = 0, nFail = 0;
        double minClearAll = 1e30, maxStepDeg = 0.0, maxDetourRatio = 1.0;
        long long nIns = 0;

        size_t i = 0;
        while(i < allPts.size()){
            if(allPts[i].spray != 0){ merged.push_back(allPts[i]); ++i; continue; }
            size_t j = i;
            while(j < allPts.size() && allPts[j].spray == 0) ++j;
            const bool hasPrev = (i > 0);
            const bool hasNext = (j < allPts.size());
            if(!hasPrev || !hasNext){                     // 首/末段没有两端喷涂点：原样保留
                for(size_t k = i; k < j; ++k) merged.push_back(allPts[k]);
                i = j; continue;
            }
            ++nSeg;
            const SprayTrajPoint& P0 = allPts[i-1];       // 起点 = 上一面喷涂末点
            const SprayTrajPoint& P1 = allPts[j];         // 终点 = 下一面喷涂首点
            const Vec3 A{P0.x,P0.y,P0.z}, B{P1.x,P1.y,P1.z};
            const Vec3 ijkA{P0.i,P0.j,P0.k}, ijkB{P1.i,P1.j,P1.k};
            const Vec3 nA{-ijkA.x,-ijkA.y,-ijkA.z}, nB{-ijkB.x,-ijkB.y,-ijkB.z};

            std::vector<Vec3> orig;
            orig.push_back(A);
            for(size_t k = i; k < j; ++k) orig.push_back({allPts[k].x, allPts[k].y, allPts[k].z});
            orig.push_back(B);
            const double L0 = std::max(1.0, pathLen(orig));
            const double c0 = pathScore(orig);

            std::vector<Vec3> best = orig;
            double bestClear = c0, bestLen = pathLen(orig);
            const char* tag = "原路径";

            // 原路径几何不合格 → 换途经：短→长依次试，取"净空最大者，同净空取最短"。
            if(c0 < kSafeNeed){
                std::vector<std::vector<Vec3>> C;
                C.push_back({A, B});                                   // 0 直连
                {                                                     // 1 反向退出（沿工件外方向退）
                    const Vec3 E = A - nA*clear;
                    const Vec3 Et = liftTo(E, std::max(sSafe, sOf(E)));
                    C.push_back({A, E, Et, liftTo(B, sOf(Et)), B});
                }
                for(double h : {0.0, 0.6*params->standoff, 1.5*params->standoff}){  // 2~4 抬升三档
                    const double sH = sSafe + h;
                    const Vec3 L1 = A + nA*clear;
                    const Vec3 L1t = liftTo(L1, std::max(sH, sOf(L1)));
                    C.push_back({A, L1, L1t, liftTo(B, sOf(L1t)), B});
                }
                // 侧绕方向 = 横向（去掉 liftU 分量）× liftU。liftU=+Z 时与旧的
                // cross(hd,{0,0,1}) 一致；侧装时不会退化成一个"沿让开方向的偏移"。
                const Vec3 dv = B - A;
                const Vec3 hz = dv - liftU*sOf(dv);
                const double span = std::max(1.0, len(hz));
                const Vec3 hd = (len(hz) < 1e-6) ? Vec3{1,0,0} : norm(hz);
                Vec3 side = cross(hd, liftU);
                side = (len(side) < 1e-6) ? Vec3{1,0,0} : norm(side);
                for(double sgn : {-1.0, 1.0})
                    for(double off : {0.6, 1.2, 2.0}){
                        const Vec3 Mid = (A + B)*0.5 + side*(sgn*span*off);
                        const Vec3 M = liftTo(Mid, sSafe);
                        C.push_back({A, M, liftTo(B, sOf(M)), B});
                    }
                const bool dbg = std::getenv("SPRAYTRAJ_MOVE_DBG") != nullptr;
                if(dbg) fprintf(stderr, "[mvdbg] A=(%.1f,%.1f,%.1f) B=(%.1f,%.1f,%.1f) 原路径 L=%.1f 净空=%.1f 门槛=%.1f\n",
                                A.x,A.y,A.z,B.x,B.y,B.z,L0,c0,kSafeNeed);
                for(const auto& cand : C){
                    const double Lc = pathLen(cand);
                    const double sc = pathScore(cand);
                    if(dbg) fprintf(stderr, "  [mvdbg] cand L=%.1f 净空=%.1f%s\n", Lc, sc,
                                    (Lc > kMaxDetour * L0) ? " (超长丢弃)" : "");
                    if(Lc > kMaxDetour * L0) continue;                // 用户口径：不做无限绕远
                    if(sc < kSafeNeed) continue;
                    const bool better = (sc > bestClear + 1e-6) ||
                                        (std::fabs(sc - bestClear) <= 1e-6 && Lc < bestLen);
                    if(better){ best = cand; bestClear = sc; bestLen = Lc; tag = "换途经"; }
                }
                if(std::string(tag) == "换途经") ++nReplan; else ++nFail;
            }
            minClearAll = std::min(minClearAll, bestClear);
            maxDetourRatio = std::max(maxDetourRatio, bestLen / L0);

            // 整段插补：位置按折线等距，姿态按弧长比例插值，单步转角超限则加密点。
            // 姿态插值策略（用户 2026-10-28 口径 + 2026-10-30 复查）：
            //   两端夹角 ≤90°（常规行间/面间）→ 直接球面插值，最短路径；
            //   >90°（侧壁↔侧壁这种反向）→ **经 -Z 极点转身**：先把枪抬到竖直朝下，
            //   再从竖直落到目标。纯 slerp 在这种情况下走的是水平面里的半圈（实测
            //   box 方位扫 179°、cao 180°、g_rim 179°），真机上就是枪身横着扫一大圈。
            //   两端本来就接近极点（仰角>80°）时不绕，直接插值，避免过极点抖动。
            std::vector<Vec3> rs; std::vector<double> ts;
            const double L = resample(best, mvStep, rs, ts);
            const Vec3 kPole{0.0, 0.0, -1.0};
            auto deg = [](double rad){ return rad * 57.29577951308232; };
            const double thAB = deg(std::acos(std::max(-1.0, std::min(1.0, dot(ijkA, ijkB)))));
            const double elA = deg(std::asin(std::max(-1.0, std::min(1.0, -ijkA.z))));
            const double elB = deg(std::asin(std::max(-1.0, std::min(1.0, -ijkB.z))));
            const bool usePole = (thAB > 90.0) && (elA <= 80.0) && (elB <= 80.0);
            const double thPole = deg(std::acos(std::max(-1.0, std::min(1.0, dot(ijkA, kPole)))))
                                + deg(std::acos(std::max(-1.0, std::min(1.0, dot(kPole, ijkB)))));
            const double thTot = usePole ? thPole : thAB;

            // 转身分配（用户 2026-10-30 口径）：平移段**保持起点姿态**，只在「接近目标的
            // 最后一段」把姿态转完，避免枪横着平移。
            // 转身段 = 路径最后一腿；太短（<15% 总长）就并入上一腿；再不够就按 5°/点
            // 所需长度补齐（否则短距离内塞 90° 会顶到 J5 限位）。
            size_t kTurn = best.size() >= 2 ? best.size() - 2 : 0;
            auto lenFromK = [&](size_t k0)->double{
                double s = 0.0;
                for(size_t q = k0 + 1; q < best.size(); ++q) s += len(best[q] - best[q-1]);
                return s;
            };
            double turnLen = lenFromK(kTurn);
            while(kTurn > 0 && turnLen < 0.15 * L){ --kTurn; turnLen = lenFromK(kTurn); }
            const double turnMin = (thTot / kStepDeg) * mvStep;
            if(turnLen < turnMin) turnLen = turnMin;
            turnLen = std::min(turnLen, std::max(1e-6, L));

            auto poseAt = [&](double tt)->Vec3{
                const double distFromEnd = (1.0 - tt) * L;
                double s = (distFromEnd >= turnLen) ? 0.0
                                                   : (turnLen - distFromEnd) / turnLen;
                s = std::max(0.0, std::min(1.0, s));
                if(!usePole) return slerpU(ijkA, ijkB, s);
                const double a = deg(std::acos(std::max(-1.0, std::min(1.0, dot(ijkA, kPole)))));
                const double thT = std::max(1e-9, thPole);
                const double sf = s * thT;
                if(sf <= a) return slerpU(ijkA, kPole, a > 1e-6 ? sf/a : 0.0);
                return slerpU(kPole, ijkB, (thT - a) > 1e-6 ? (sf - a)/(thT - a) : 1.0);
            };
            const int need = (int)std::ceil(thTot / kStepDeg) + 1;
            if(need > (int)rs.size() && L > 1e-6)
                resample(best, std::max(0.3, L / (double)(need - 1)), rs, ts);

            if(rs.size() <= 2){                                   // 太短：保留原折线，不抹掉拐点
                for(size_t k = i; k < j; ++k) merged.push_back(allPts[k]);
                i = j; continue;
            }
            Vec3 prevG = ijkA;
            for(size_t k = 1; k + 1 < rs.size(); ++k){
                const Vec3 g = poseAt(ts[k]);
                const double dg = std::acos(std::max(-1.0, std::min(1.0, dot(prevG, g)))) * 57.29577951308232;
                maxStepDeg = std::max(maxStepDeg, dg);
                prevG = g;
                merged.push_back({rs[k].x, rs[k].y, rs[k].z, g.x, g.y, g.z, P0.speed, 0});
                ++nIns;
            }
            i = j;
        }

        if(nSeg > 0){
            fprintf(stderr, "[空移规划] %d 段：换途经 %d 段，无可行途经 %d 段；"
                            "最小净空 %.1fmm，最大单步转角 %.2f°，最长绕行 %.2f×原路径，插补 %lld 点\n",
                    nSeg, nReplan, nFail, minClearAll, maxStepDeg, maxDetourRatio, nIns);
            allPts.swap(merged);
        }
    }

    // ================= 空移段等距细分（第二十轮） =================
    // 病灶：空移段只保留**拐点**，两点之间靠控制器直线插补走完。实测（行距 80）：
    //   开口盒子 8 个空移段，A→B 跨 583mm 的那段只有 3 个点，相邻间距 600mm；
    //   复杂点云有一整段只有 **1 个哨兵点**、跨 634mm。
    //   ⇒ 上位机/模拟器看到的空移是"瞬移"，中间到底撞不撞件完全无从判断，
    //     robot_safety 的逐点检查（IK / 连杆 / 自碰撞 / 间隙）也覆盖不到这些区间。
    // 处理：凡"两端至少有一端是 spray=0"的相邻点对，沿 A→B 直线按点距等距插点，
    //   一律 spray=0（落枪/抬枪的进出位移也算空移，所以同样铺点）。
    //   **只加标记、不改几何**：插点严格落在这条直线上，物理路径与原来逐点一致。
    //   枪方向用 nlerp 插值（两点反向时退化为各自端值，避免过零向量）；
    //   速度沿用前一点，不引入新的速度语义。
    // 顺序很关键：必须在转角平滑/去重/重采样**之后**做 —— 平滑要按原始拐点算曲率，
    //   先铺密点会让平滑误判"这里有一串点"，把整条空移线拉弯。
    // 可用 SPRAYTRAJ_MOVE_STEP 覆盖空移点距（mm），默认 = 喷涂点距。
    {
        double mvStep = pstep;
        if(const char* e = std::getenv("SPRAYTRAJ_MOVE_STEP")){
            double v = atof(e);
            if(v > 0.5) mvStep = v;
        }
        const int kMaxInsPerGap = 4000;      // 单段插点上限，防病态数据炸内存
        std::vector<SprayTrajPoint> dense;
        dense.reserve(allPts.size() * 2);
        long long nIns = 0;
        int nGap = 0;
        for(size_t k = 0; k < allPts.size(); ++k){
            const SprayTrajPoint& A = allPts[k];
            dense.push_back(A);
            if(k + 1 >= allPts.size()) break;
            const SprayTrajPoint& B = allPts[k + 1];
            if(A.spray != 0 && B.spray != 0) continue;      // 纯喷涂段：本来就密，不动
            const double dx = B.x - A.x, dy = B.y - A.y, dz = B.z - A.z;
            const double d = std::sqrt(dx*dx + dy*dy + dz*dz);
            if(!(d > 1.5 * mvStep)) continue;               // 已经够密（或零长）
            int m = (int)std::floor(d / mvStep);
            if(m > kMaxInsPerGap) m = kMaxInsPerGap;
            if(m < 1) continue;
            ++nGap;
            // 枪轴插值必须用 slerp：nlerp 在两端反向（cosθ<=0）时退化成"全部沿用 A 端"，
            // 于是 180° 的翻转被压到最后一个 1.5mm 的点上（开口盒子 A1 段实测 180°/1.46mm），
            // 与修复前的"原地翻腕"等价。反向(antipodal)时用固定轴（A 的任意垂线）转过去，
            // 保证每一步转角 = 180°/插点数，而不是一步到位。
            double cosAB = A.i*B.i + A.j*B.j + A.k*B.k;
            if(cosAB > 1.0) cosAB = 1.0;
            if(cosAB < -1.0) cosAB = -1.0;
            const double thAB = std::acos(cosAB);
            const bool antipodal = (thAB > 3.141592653589793 - 1e-3);
            Vec3 ax{0,0,0};
            if(antipodal){
                const Vec3 a{A.i, A.j, A.k};
                const Vec3 ref = (std::fabs(a.z) < 0.9) ? Vec3{0,0,1} : Vec3{1,0,0};
                ax = norm(cross(a, ref));
            }
            for(int t = 1; t <= m; ++t){
                const double f = (double)t * mvStep / d;
                if(f >= 1.0 - 1e-9) break;                  // 不要盖住 B 本身
                double ii, jj, kk;
                if(antipodal){
                    const double a0 = 3.141592653589793 * f;
                    ii = A.i*std::cos(a0) + ax.x*std::sin(a0);
                    jj = A.j*std::cos(a0) + ax.y*std::sin(a0);
                    kk = A.k*std::cos(a0) + ax.z*std::sin(a0);
                } else if(thAB > 1e-6){
                    const double s = std::sin(thAB);
                    const double w0 = std::sin((1.0 - f) * thAB) / s;
                    const double w1 = std::sin(f * thAB) / s;
                    ii = A.i*w0 + B.i*w1;
                    jj = A.j*w0 + B.j*w1;
                    kk = A.k*w0 + B.k*w1;
                } else {
                    ii = A.i; jj = A.j; kk = A.k;
                }
                const double nn = std::sqrt(ii*ii + jj*jj + kk*kk);
                if(nn > 1e-9){ ii /= nn; jj /= nn; kk /= nn; }
                else         { ii = A.i; jj = A.j; kk = A.k; }
                dense.push_back({ A.x + dx * f, A.y + dy * f, A.z + dz * f,
                                  ii, jj, kk, A.speed, 0 });
                ++nIns;
            }
        }
        if(nIns > 0){
            fprintf(stderr, "[空移] 铺满 %d 段空移路径，补 %lld 个 spray=0 点（点距 %.1fmm），"
                            "总点数 %zu → %zu\n",
                    nGap, nIns, mvStep, allPts.size(), dense.size());
            allPts.swap(dense);
        }
    }

    if(allPts.empty()){
        failure = nDownSkip > 0 ? TrajectoryFailure::AngleWindow : TrajectoryFailure::NoRows;
        err = "无有效扫描行";
        if(nGunReject > 0 || nDownSkip > 0 || nPlanarSkip > 0){
            err += "（";
            if(nGunReject > 0)
                err += "枪位验收未通过 " + std::to_string(nGunReject) + " 面";
            if(nDownSkip > 0){
                if(nGunReject > 0) err += "；";
                err += "超出可喷角度窗口 " + std::to_string(nDownSkip) + " 面";
            }
            if(nPlanarSkip > 0){
                if(nGunReject > 0 || nDownSkip > 0) err += "；";
                err += "平面性不足 " + std::to_string(nPlanarSkip) + " 面";
            }
            err += "）";
        }
        if(nGunReject > 0)
            err += "。请检查点云中喷涂侧的遮挡及面吸收阈值。";
        res.error=_strdup(err.c_str()); return;
    }
    if(usedCandidateCheck){
        for(size_t k=0; k<allPts.size(); ++k){
            const auto& a = allPts[k == 0 ? 0 : k-1];
            const auto& b = allPts[k];
            const double need = a.spray && b.spray ? gunClearEff : std::max(moveClear, gunClearEff);
            if(candidateSegmentClear({a.x,a.y,a.z}, {b.x,b.y,b.z}, need)) continue;
            fprintf(stderr, "[skip] 最终轨迹原始点云净空复核失败：第 %zu 段，要求 %.1fmm\n", k, need);
            failure = TrajectoryFailure::CandidatePath;
            res.error = _strdup("最终轨迹原始点云净空复核失败");
            return;
        }
    }
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

    // ---- 出口：把内部（up=+Z）坐标系逆变换回世界系 ----
    // 轨迹点/枪轴、面法向/质心/局部基/重建面顶点、总法向与质心，全部要转。
    // angZ/tilt 保留内部值（相对可喷侧 up 的夹角，0 = 正对可喷方向），不重算：
    // 世界系的 n 与 +Z 夹角对这类工位恒接近 180°，拿来显示和判窗口都是错的。
    if(uf.active){
        for(int i=0;i<cnt;++i){
            SprayTrajPoint& q = res.points[i];
            const Vec3 p = uf.inv({q.x,q.y,q.z});
            const Vec3 d = uf.inv({q.i,q.j,q.k});
            q.x=p.x; q.y=p.y; q.z=p.z;
            q.i=d.x; q.j=d.y; q.k=d.z;
        }
        for(int pi=0; pi<res.plane_count && pi<SPRAYTRAJ_PLANE_SLOTS; ++pi){
            SprayTrajPlane& q = res.planes[pi];
            const Vec3 nn = uf.inv({q.nx,q.ny,q.nz});
            const Vec3 n0 = uf.inv({q.n0x,q.n0y,q.n0z});
            const Vec3 cc = uf.inv({q.cx,q.cy,q.cz});
            q.nx=nn.x; q.ny=nn.y; q.nz=nn.z;
            q.n0x=n0.x; q.n0y=n0.y; q.n0z=n0.z;
            q.cx=cc.x; q.cy=cc.y; q.cz=cc.z;
            const Vec3 uu = uf.inv({q.ux,q.uy,q.uz});
            const Vec3 vv = uf.inv({q.vx,q.vy,q.vz});
            q.ux=uu.x; q.uy=uu.y; q.uz=uu.z;
            q.vx=vv.x; q.vy=vv.y; q.vz=vv.z;
            for(int k=0;k<q.polyVertCount*3;k+=3){
                const Vec3 v = uf.inv({q.polyVerts[k],q.polyVerts[k+1],q.polyVerts[k+2]});
                q.polyVerts[k]=v.x; q.polyVerts[k+1]=v.y; q.polyVerts[k+2]=v.z;
            }
        }
        const Vec3 rn = uf.inv({res.normal_x,res.normal_y,res.normal_z});
        const Vec3 rc = uf.inv({res.plane_cx,res.plane_cy,res.plane_cz});
        res.normal_x=rn.x; res.normal_y=rn.y; res.normal_z=rn.z;
        res.plane_cx=rc.x; res.plane_cy=rc.y; res.plane_cz=rc.z;
    }

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
}

// ---- 对外入口：整件无轨迹时按失败原因重试 ----
//
// 全程只写 out，不产生按值返回的大结构体：调用方（P/Invoke 或栈上局部）
// 只需要一块 289KB 的结果区，栈占用与改造前一致。
static void run_with_retry(const double* xyz, int n, const SprayTrajParams* params,
                           SprayTrajResult* out)
{
    TrajectoryFailure failure;
    auto runAttempt = [&](){
        spray_traj_core(xyz, n, params, out, failure, true);
        if(failure == TrajectoryFailure::CandidatePath){
            // 新候选加入后整条路径未通过，恢复原有验收流程，保住原来能喷的面。
            fprintf(stderr, "[gun] 新候选路径复核失败，按原有验收流程重新布行\n");
            spray_traj_free(out);
            spray_traj_core(xyz, n, params, out, failure, false);
        }
    };
    runAttempt();
    // 自动体素档位可能将薄壳切成两层，退一档重试。保留用户显式设置的体素档位。
    if(out->error && params && n > 120000 && !std::getenv("SPRAYTRAJ_VOXEL")
       && failure == TrajectoryFailure::NoRows){
        double leaf = auto_voxel_leaf(params);
        for(int k=0;k<2;++k){
            leaf *= 1.6;                                  // 7.5 → 12 → 19
            if(leaf > 25.0) break;
            fprintf(stderr, "[体素] 整件无有效行（薄壳/档位巧合），退档重试 %.1fmm\n", leaf);
            spray_traj_free(out);
            g_voxelLeafOverride = leaf;
            runAttempt();
            if(!out->error) break;
        }
        g_voxelLeafOverride = -1.0;
    }
}

SPRAYTRAJ_API SprayTrajResult spray_traj_generate(
    const double* xyz, int n, const SprayTrajParams* params)
{
    // 返回值直接用 sret 槽（NRVO），不额外复制：这里不能再有第二个大结构体局部变量。
    SprayTrajResult r;
    run_with_retry(xyz, n, params, &r);
    return r;
}

SPRAYTRAJ_API void spray_traj_free(SprayTrajResult* res){
    if(!res) return;
    if(res->points) free(res->points);
    if(res->error)  free((void*)res->error);
    res->points=nullptr; res->error=nullptr; res->count=0;
}

// ---- C# 友好入口（2026-11-05 接入 SprayVisionSimple 时新增）----
// 为什么不直接用 spray_traj_generate 的"返回结构体"：SprayTrajResult 约 145KB，
// 托管侧 P/Invoke 返回大结构体会在调用线程栈上分配，70 万点/多线程场景有栈溢出风险；
// 改为调用方（C#）用 Marshal.AllocHGlobal 提供缓冲区，内部浅拷贝所有权。
SPRAYTRAJ_API int spray_traj_generate_ex(const double* xyz, int n,
                                         const SprayTrajParams* params, SprayTrajResult* out){
    if(!out) return -1;
    std::memset(out, 0, sizeof(*out));
    run_with_retry(xyz, n, params, out);      // 直接写进调用方缓冲区，不产生栈上大副本
    return out->error ? 1 : 0;
}

SPRAYTRAJ_API void spray_traj_abi_sizes(int* paramsSize, int* planeSize, int* resultSize){
    if(paramsSize) *paramsSize = (int)sizeof(SprayTrajParams);
    if(planeSize)  *planeSize  = (int)sizeof(SprayTrajPlane);
    if(resultSize) *resultSize = (int)sizeof(SprayTrajResult);
}

SPRAYTRAJ_API void spray_traj_result_info(const SprayTrajResult* res, SprayTrajInfo* out){
    if(!out) return;
    std::memset(out, 0, sizeof(*out));
    if(!res) return;
    out->count = res->count;
    out->rows = res->rows;
    out->input_points = res->input_points;
    out->removed_points = res->removed_points;
    out->inlier_points = res->inlier_points;
    out->plane_count = res->plane_count;
    out->plane_rms = res->plane_rms;
    out->normal_x = res->normal_x; out->normal_y = res->normal_y; out->normal_z = res->normal_z;
    out->plane_cx = res->plane_cx; out->plane_cy = res->plane_cy; out->plane_cz = res->plane_cz;
    out->seg_len_mean = res->seg_len_mean;
    out->seg_len_max = res->seg_len_max;
}

SPRAYTRAJ_API int spray_traj_get_plane(const SprayTrajResult* res, int idx, SprayTrajPlane* out){
    if(!res || !out || idx < 0 || idx >= res->plane_count || idx >= SPRAYTRAJ_PLANE_SLOTS) return -1;
    *out = res->planes[idx];
    return 0;
}

SPRAYTRAJ_API const char* spray_traj_result_error(const SprayTrajResult* res){
    return res ? res->error : nullptr;
}
} // extern "C"
