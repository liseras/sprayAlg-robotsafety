// main.cpp - 扁平类工件喷粉轨迹生成（分区工艺版）
//
// 用法:
//   spray_traj --input cloud.ply --outdir out \
//     [--pass-spacing 100] [--standoff 200] [--margin 12] [--step 5] \
//     [--min-run 25] [--tilt 0] [--speed 400] [--flat 0] [--cell 3] \
//     [--max-hole-area 5000] [--ring-offset 30] [--ring-count 1] \
//     [--slope-angle 25] [--slope-spacing 0]
//
// 输出（写入 --outdir）:
//   trajectory.csv    机器人轨迹（位置 + 四元数 + 速度 + 工艺分区），坐标系与点云一致
//   trajectory.ply    轨迹点云
//   preview.bmp       俯视预览（分区着色 + 轨迹）
//
// 分区工艺:
//   平面区 -> 横向之字（沿长边走, v 方向按行距进给）
//   孔周区 -> 沿孔整圆走一圈（半径 = 孔半径 + ring-offset）
//   斜面区 -> 工具垂直斜面, 沿斜面长度方向走直线
// 工艺铁律: 进入工作距离后不关枪, 全部路点均为 spray, 仅 zone 区分
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include "vecmath.h"
#include "ply_io.h"
#include "plane_fit.h"
#include "coverage.h"
#include "heightfield.h"
#include "regions.h"
#include "pathgen.h"
#include "ordering.h"
#include "preview.h"

static std::string getArg(int argc, char** argv, const std::string& key,
                          const std::string& def) {
    for (int i = 1; i + 1 < argc; ++i)
        if (key == argv[i]) return argv[i + 1];
    return def;
}

int main(int argc, char** argv) {
    std::string input = getArg(argc, argv, "--input", "");
    std::string outdir = getArg(argc, argv, "--outdir", "out");
    SprayParams P;
    P.passSpacing   = atof(getArg(argc, argv, "--pass-spacing", "100").c_str());
    P.standoff      = atof(getArg(argc, argv, "--standoff", "200").c_str());
    P.margin        = atof(getArg(argc, argv, "--margin", "12").c_str());
    P.step          = atof(getArg(argc, argv, "--step", "5").c_str());
    P.minRunLen     = atof(getArg(argc, argv, "--min-run", "25").c_str());
    P.tiltDeg       = atof(getArg(argc, argv, "--tilt", "0").c_str());
    P.speed         = atof(getArg(argc, argv, "--speed", "400").c_str());
    P.flat          = atoi(getArg(argc, argv, "--flat", "0").c_str()) != 0;
    P.slopeAngle    = atof(getArg(argc, argv, "--slope-angle", "25").c_str());
    P.slopeSpacing  = atof(getArg(argc, argv, "--slope-spacing", "0").c_str());
    P.slopeMargin   = atof(getArg(argc, argv, "--slope-margin", "6").c_str());
    P.ringOffset    = atof(getArg(argc, argv, "--ring-offset", "30").c_str());
    P.ringCount     = atoi(getArg(argc, argv, "--ring-count", "1").c_str());
    P.ringSpacing   = atof(getArg(argc, argv, "--ring-spacing", "0").c_str());
    double cell     = atof(getArg(argc, argv, "--cell", "3").c_str());
    double maxHoleArea = atof(getArg(argc, argv, "--max-hole-area", "5000").c_str());
    int minCellPts  = atoi(getArg(argc, argv, "--min-cell-pts", "3").c_str());
    double planeTh  = atof(getArg(argc, argv, "--plane-thresh", "6").c_str());
    // 段顺序优化
    int orderMode   = atoi(getArg(argc, argv, "--order", "2").c_str());
    int orderMetric = atoi(getArg(argc, argv, "--order-metric", "1").c_str());
    double orderZonePenalty = atof(getArg(argc, argv, "--order-zone-penalty", "120").c_str());
    double orderAngleW      = atof(getArg(argc, argv, "--order-angle-w", "0.25").c_str());
    double orderStartU = atof(getArg(argc, argv, "--start-u", "nan").c_str());
    double orderStartV = atof(getArg(argc, argv, "--start-v", "nan").c_str());

    if (input.empty()) {
        fprintf(stderr, "用法: spray_traj --input cloud.ply --outdir out [参数]\n");
        fprintf(stderr, "  --pass-spacing 行距mm  --standoff 喷距mm  --margin 边界内缩mm\n");
        fprintf(stderr, "  --step 点间距mm  --min-run 最短行程mm  --tilt 前倾角deg\n");
        fprintf(stderr, "  --speed 速度mm/s  --flat 1=按纯平面  --cell 栅格mm\n");
        fprintf(stderr, "  --max-hole-area 小孔填充上限mm2  --ring-offset 环线离孔边mm\n");
        fprintf(stderr, "  --ring-count 每个孔环线数  --slope-angle 斜面判定角deg\n");
        fprintf(stderr, "  --slope-spacing 斜面线间距mm  --slope-margin 斜面端点内缩mm\n");
        fprintf(stderr, "  --order 段顺序优化 0=关闭 1=最近邻+2opt 2=自动(默认)\n");
        fprintf(stderr, "  --order-metric 0=直线距离 1=掩膜绕行最短路(默认)\n");
        fprintf(stderr, "  --order-zone-penalty 换区惩罚mm(默认120)  --order-angle-w 姿态转角权重(默认0.25)\n");
        fprintf(stderr, "  --start-u/--start-v 机器人起始接近点(主平面坐标, 默认自由起点)\n");
        return 1;
    }

    // ---------- 1. 读点云 ----------
    PointCloud cloud;
    std::string err;
    if (!loadPLY(input, cloud, err)) { fprintf(stderr, "错误: %s\n", err.c_str()); return 1; }
    printf("[1/8] 点云读取: %zu 点\n", cloud.pts.size());

    // ---------- 2. 粗略平面 + 覆盖掩膜 + 工件点归属 ----------
    std::vector<int> inliers;
    Plane pl = fitPlaneRobust(cloud.pts, 4, 2.0, &inliers);
    printf("[2/8] 稳健平面: 内点 %zu (%.1f%%)  法向 (%.4f %.4f %.4f)\n",
           inliers.size(), 100.0 * inliers.size() / cloud.pts.size(),
           pl.n.x, pl.n.y, pl.n.z);

    double C0[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i : inliers) {
        Vec3 d = cloud.pts[i] - pl.c;
        C0[0][0]+=d.x*d.x; C0[0][1]+=d.x*d.y; C0[0][2]+=d.x*d.z;
        C0[1][1]+=d.y*d.y; C0[1][2]+=d.y*d.z; C0[2][2]+=d.z*d.z;
    }
    C0[1][0]=C0[0][1]; C0[2][0]=C0[0][2]; C0[2][1]=C0[1][2];
    double eig0[3]; Vec3 vec0[3];
    jacobiEigen3(C0, eig0, vec0);
    Vec3 vAxis0 = normalize(cross(pl.n, vec0[2]));
    std::vector<int> allIdx(cloud.pts.size());
    for (size_t i = 0; i < cloud.pts.size(); ++i) allIdx[i] = (int)i;
    Mask mk = buildMask(cloud.pts, allIdx, pl, vec0[2], vAxis0, cell, minCellPts);
    keepLargestComponent(mk);
    std::vector<int> panelIdx;
    for (int i : allIdx) {
        Vec3 d = cloud.pts[i] - pl.c;
        int a = (int)((dot(d, vec0[2]) - mk.umin) / mk.cell);
        int b = (int)((dot(d, vAxis0) - mk.vmin) / mk.cell);
        if (mk.inside(a, b) && mk.at(a, b)) panelIdx.push_back(i);
    }
    if (panelIdx.size() < 100) { fprintf(stderr, "错误: 连通域过滤后点太少\n"); return 1; }

    // ---------- 3. RANSAC 主面 + OBB 面内轴向 ----------
    pl = fitPlaneRANSAC(cloud.pts, panelIdx, planeTh, 400);
    {
        size_t cnt = 0;
        for (int i : panelIdx)
            if (std::fabs(pl.d(cloud.pts[i])) <= planeTh) cnt++;
        printf("      主面一致集 %zu (%.1f%%), 法向 (%.4f %.4f %.4f)\n",
               cnt, 100.0 * cnt / panelIdx.size(), pl.n.x, pl.n.y, pl.n.z);
    }
    double C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i : panelIdx) {
        Vec3 d = cloud.pts[i] - pl.c;
        C[0][0]+=d.x*d.x; C[0][1]+=d.x*d.y; C[0][2]+=d.x*d.z;
        C[1][1]+=d.y*d.y; C[1][2]+=d.y*d.z; C[2][2]+=d.z*d.z;
    }
    C[1][0]=C[0][1]; C[2][0]=C[0][2]; C[2][1]=C[1][2];
    double eig[3]; Vec3 vec2_[3];
    jacobiEigen3(C, eig, vec2_);
    Vec3 u0 = vec2_[2], v0 = normalize(cross(pl.n, u0));
    double bestTheta = 0, bestArea = 1e18;
    const double PI = 3.14159265358979323846;
    for (int deg = 0; deg < 90; ++deg) {
        double th = deg * PI / 180.0;
        Vec3 ua = normalize(u0 * std::cos(th) + v0 * std::sin(th));
        Vec3 va = normalize(cross(pl.n, ua));
        double umn = 1e18, umx = -1e18, vmn = 1e18, vmx = -1e18;
        for (size_t k = 0; k < panelIdx.size(); k += 3) {
            Vec3 d = cloud.pts[panelIdx[k]] - pl.c;
            double uu = dot(d, ua), vv = dot(d, va);
            umn = std::min(umn, uu); umx = std::max(umx, uu);
            vmn = std::min(vmn, vv); vmx = std::max(vmx, vv);
        }
        double area = (umx - umn) * (vmx - vmn);
        if (area < bestArea) { bestArea = area; bestTheta = th; }
    }
    Vec3 uAxis = normalize(u0 * std::cos(bestTheta) + v0 * std::sin(bestTheta));
    Vec3 vAxis = normalize(cross(pl.n, uAxis));
    double eu = 0, ev = 0;
    for (size_t k = 0; k < panelIdx.size(); k += 3) {
        Vec3 d = cloud.pts[panelIdx[k]] - pl.c;
        eu = std::max(eu, std::fabs(dot(d, uAxis)));
        ev = std::max(ev, std::fabs(dot(d, vAxis)));
    }
    if (ev > eu) uAxis = normalize(vAxis);   // u 取长边, 之字沿长边往复
    vAxis = normalize(cross(pl.n, uAxis));
    printf("[3/8] 面内轴向 u=(%.4f %.4f %.4f)  v=(%.4f %.4f %.4f)  OBB旋转 %.1fdeg\n",
           uAxis.x, uAxis.y, uAxis.z, vAxis.x, vAxis.y, vAxis.z, bestTheta * 180.0 / PI);

    // ---------- 4. 重投影栅格：识别大孔 / 填小孔 / 腐蚀 ----------
    mk = buildMask(cloud.pts, panelIdx, pl, uAxis, vAxis, cell, minCellPts);
    keepLargestComponent(mk);
    std::vector<HoleInfo> holes = findOpenHoles(mk, maxHoleArea);
    int holeCount = 0;
    std::vector<uint8_t> smallHoles = fillEnclosedHoles(mk, maxHoleArea, &holeCount);
    Mask mkFull = mk;                 // 未腐蚀（含已填小孔）, 用于预览
    Mask mkTrav = mk;
    erodeMask(mkTrav, P.margin);      // 可喷涂/可通行区域（距边界 >= margin）
    printf("[4/8] 覆盖栅格 %dx%d (cell %.1fmm)  填小孔 %d 处  大孔 %zu 处  内缩 %.0fmm\n",
           mk.W, mk.H, cell, holeCount, holes.size(), P.margin);
    for (size_t i = 0; i < holes.size(); ++i)
        printf("      大孔#%zu: 中心(%.1f, %.1f)  等效半径 %.1f  面积 %.0fmm2\n",
               i, holes[i].cu, holes[i].cv, holes[i].rEq, holes[i].area);

    // ---------- 5. 斜面区识别（粗尺度高度场梯度） ----------
    SlopeZone slope;
    HeightField hc;
    hc.build(cloud.pts, panelIdx, pl, uAxis, vAxis, 10.0);
    hc.boxBlur();                      // 有效窗口约 20mm, 压低冲压波纹
    // 诊断: 倾角分布
    {
        int c15 = 0, c25 = 0, c35 = 0, tot = 0;
        double d = hc.cell;
        for (int j = 0; j < mk.H; ++j)
            for (int i = 0; i < mk.W; ++i) {
                if (!mk.at(i, j)) continue;
                double u = mk.u(i), v = mk.v(j);
                double gu = (hc.sample(u + d, v) - hc.sample(u - d, v)) / (2 * d);
                double gv = (hc.sample(u, v + d) - hc.sample(u, v - d)) / (2 * d);
                double ang = std::atan(std::hypot(gu, gv)) * 180.0 / PI;
                tot++;
                if (ang > 15) c15++;
                if (ang > 25) c25++;
                if (ang > 35) c35++;
            }
        printf("      倾角分布: >15deg %.1f%%  >25deg %.1f%%  >35deg %.1f%% (共 %d 格)\n",
               100.0 * c15 / std::max(1, tot), 100.0 * c25 / std::max(1, tot),
               100.0 * c35 / std::max(1, tot), tot);
    }
    // 斜面检测用更深腐蚀的掩膜: 边界卷边(尤其左右侧边)梯度会被误判为陡峭,
    // 与真实斜面带连通后污染拟合。多内缩 10mm 隔断边缘卷边。
    Mask mkSteep = mk;
    erodeMask(mkSteep, P.margin + 10.0);
    std::vector<uint8_t> steep = steepFromHeight(hc, mkSteep, P.slopeAngle);
    maskOpen(steep, mk.W, mk.H);
    // 只保留最大连通域（去掉零星陡峭噪点）
    {
        std::vector<uint8_t> best;
        std::vector<uint8_t> vis(steep.size(), 0);
        int bestCnt = 0;
        std::vector<int> comp;
        for (int j = 0; j < mk.H; ++j)
            for (int i = 0; i < mk.W; ++i) {
                size_t k = (size_t)j * mk.W + i;
                if (!steep[k] || vis[k]) continue;
                comp.clear();
                std::queue<int> q; q.push((int)k); vis[k] = 1;
                while (!q.empty()) {
                    int c = q.front(); q.pop(); comp.push_back(c);
                    int ci = c % mk.W, cj = c / mk.W;
                    const int di[4] = {1,-1,0,0}, dj[4] = {0,0,1,-1};
                    for (int t = 0; t < 4; ++t) {
                        int ni = ci + di[t], nj = cj + dj[t];
                        if (ni < 0 || ni >= mk.W || nj < 0 || nj >= mk.H) continue;
                        size_t nk = (size_t)nj * mk.W + ni;
                        if (steep[nk] && !vis[nk]) { vis[nk] = 1; q.push((int)nk); }
                    }
                }
                if ((int)comp.size() > bestCnt) {
                    bestCnt = (int)comp.size();
                    best.assign(steep.size(), 0);
                    for (int c : comp) best[(size_t)c] = 1;
                }
            }
        if (bestCnt > 0) steep.swap(best);
    }
    int steepCnt = countMask(steep);
    double steepArea = steepCnt * cell * cell;
    double maskArea = countMask(mk.m) * cell * cell;
    // 诊断: 陡峭格的 v 分布
    {
        int bins[20] = {0};
        double vLo = 1e18, vHi = -1e18;
        for (int j = 0; j < mk.H; ++j)
            for (int i = 0; i < mk.W; ++i)
                if (mk.at(i, j)) { vLo = std::min(vLo, mk.v(j)); vHi = std::max(vHi, mk.v(j)); }
        for (int j = 0; j < mk.H; ++j)
            for (int i = 0; i < mk.W; ++i)
                if (steep[(size_t)j * mk.W + i]) {
                    double v = mk.v(j);
                    int b = (int)((v - vLo) / std::max(1.0, (vHi - vLo) / 20));
                    if (b < 0) b = 0; if (b > 19) b = 19;
                    bins[b]++;
                }
        printf("      [诊断] steep v分布(20桶, %.0f..%.0f):", vLo, vHi);
        for (int b = 0; b < 20; ++b) printf(" %d", bins[b]);
        printf("\n");
    }
    printf("[5/8] 斜面区: %d 格 (%.0fmm2, 占工件 %.1f%%)\n",
           steepCnt, steepArea, 100.0 * steepArea / std::max(1.0, maskArea));
    if (steepArea > 0.45 * maskArea) {
        printf("      警告: 陡峭区占比过大, 判定为整体曲面, 本次按平面处理\n");
        std::fill(steep.begin(), steep.end(), (uint8_t)0);
        steepArea = 0;
    }
    if (steepArea > 1500) {
        if (buildSlopeZone(cloud.pts, steep, mk, pl, uAxis, vAxis, slope)) {
            printf("      斜面拟合: 法向 (%.4f %.4f %.4f)  与主面夹角 %.1fdeg\n",
                   slope.plane.n.x, slope.plane.n.y, slope.plane.n.z, slope.angleToMain);
            printf("      斜面尺寸: 长 %.0fmm  宽 %.0fmm (沿表面)\n",
                   slope.lenMax - slope.lenMin, slope.widMax - slope.widMin);
            {   // 诊断: 斜面区点相对主面的残差(判断是凸起带还是内凹卷边)
                std::vector<double> rs;
                for (int i : slope.ptIdx) rs.push_back(pl.d(cloud.pts[i]));
                std::sort(rs.begin(), rs.end());
                printf("      斜面残差: p10 %.1f  p50 %.1f  p90 %.1f  (相对主面mm)\n",
                       rs[rs.size() / 10], rs[rs.size() / 2], rs[rs.size() * 9 / 10]);
            }
        } else {
            printf("      斜面点太少, 按平面处理\n");
        }
    }

    // ---------- 6. 生成各分区段 ----------
    HeightField hf;
    hf.build(cloud.pts, panelIdx, pl, uAxis, vAxis, 5.0);
    std::vector<std::string> ringLog;
    std::vector<Run> ringRuns = buildRingRuns(holes, mkTrav, mk, pl, uAxis, vAxis, hf, P, &ringLog);
    for (auto& s : ringLog) printf("      环线: %s\n", s.c_str());

    std::vector<Run> slopeRuns = buildSlopeRuns(slope, mkFull, pl, uAxis, vAxis, P);

    // 平面之字区域 = 可喷涂区 - 斜面区(膨胀) - 孔周环带
    Mask mkMain = mkTrav;
    {
        // 斜面区从平面之字中剔除, 膨胀约 1/4 行距（过渡带让给斜面喷涂覆盖）
        std::vector<uint8_t> st = steep;
        dilateMask(st, mk.W, mk.H, std::max(2, (int)(P.passSpacing * 0.25 / cell)));
        for (size_t k = 0; k < st.size(); ++k) if (st[k]) mkMain.m[k] = 0;
    }
    for (size_t hi = 0; hi < holes.size(); ++hi) {
        double rMax = -1;
        for (auto& r : ringRuns) {
            if (r.group != (int)hi) continue;
            for (size_t k = 0; k < r.p.size(); ++k) {
                double du = r.uu[k] - holes[hi].cu, dv = r.vv[k] - holes[hi].cv;
                rMax = std::max(rMax, std::hypot(du, dv));
            }
        }
        if (rMax < 0) continue;   // 没放环线就不挖空, 让之字照常覆盖
        double rExcl = rMax + P.passSpacing * 0.5;
        for (int j = 0; j < mkMain.H; ++j)
            for (int i = 0; i < mkMain.W; ++i) {
                if (!mkMain.at(i, j)) continue;
                if (std::hypot(mkMain.u(i) - holes[hi].cu, mkMain.v(j) - holes[hi].cv) < rExcl)
                    mkMain.at(i, j) = 0;
            }
    }
    std::vector<Run> planeRuns = buildPlaneRuns(mkMain, pl, uAxis, vAxis, hf, P);

    // 组装：平面之字 + 孔周整圆 + 斜面直线，全部交给段顺序优化统一排序。
    // 原来的写法是"平面全在前 → 孔周按单端贪心插入 → 斜面完全不排序"，
    // 段一多就会出现跨工件的往返空移（p029 实测过渡占全长 41%）。
    std::vector<Run> allRuns;
    allRuns.reserve(planeRuns.size() + ringRuns.size() + slopeRuns.size());
    for (auto& r : planeRuns) allRuns.push_back(std::move(r));
    for (auto& r : ringRuns)  allRuns.push_back(std::move(r));
    for (auto& r : slopeRuns) allRuns.push_back(std::move(r));

    printf("[6/8] 分区段: 平面 %zu 段, 孔周 %zu 段, 斜面 %zu 段\n",
           planeRuns.size(), ringRuns.size(), slopeRuns.size());

    // ---------- 6.5 段顺序与方向优化（开放路径：段可反向、起点自由） ----------
    if (orderMode > 0 && allRuns.size() >= 2) {
        OrderOpts OO;
        OO.mode = orderMode;
        OO.metric = orderMetric;
        OO.zonePenalty = orderZonePenalty;
        OO.angleWeight = orderAngleW;
        OO.standoff = P.standoff;
        OO.startU = orderStartU;
        OO.startV = orderStartV;
        RunOrder ro = optimizeRunOrder(allRuns, mkTrav, pl, uAxis, vAxis, OO);
        applyRunOrder(allRuns, ro);
        printf("      段顺序: %zu 段，空移 %.0f → %.0f mm（%+.0f%%）%s\n",
               ro.order.size(), ro.costBefore, ro.costAfter,
               100.0 * (ro.costAfter - ro.costBefore) / std::max(1.0, ro.costBefore),
               ro.changed ? "" : "  [无改善，保持原序]");
        if (getenv("SPRAY_ORDER_DBG")) {
            printf("      [dbg] 执行顺序:");
            for (size_t k = 0; k < ro.order.size(); ++k) {
                const Run& r = allRuns[k];
                printf(" %s%d%s", r.zone.c_str(), r.group, ro.rev[k] ? "(反)" : "");
            }
            printf("\n");
        }
    }

    // ---------- 7. 串接 + 姿态 ----------
    std::vector<Waypoint> wp = sequenceRuns(allRuns, mkTrav, pl, uAxis, vAxis, hf, P);
    finalizeTrajectory(wp, P);
    if (wp.empty()) { fprintf(stderr, "错误: 未生成轨迹\n"); return 1; }

    double lenPlane = 0, lenRing = 0, lenSlope = 0, lenTravel = 0;
    for (size_t k = 1; k < wp.size(); ++k) {
        double L = norm(wp[k].pos - wp[k - 1].pos);
        const std::string& z = wp[k].zone;
        if (z == "plane") lenPlane += L;
        else if (z == "ring") lenRing += L;
        else if (z == "slope") lenSlope += L;
        else lenTravel += L;
    }
    double hmin = 1e18, hmax = -1e18;
    for (auto& w : wp) {
        if (w.zone == "travel") continue;
        hmin = std::min(hmin, w.h); hmax = std::max(hmax, w.h);
    }
    double total = lenPlane + lenRing + lenSlope + lenTravel;
    printf("[7/8] 轨迹: %zu 点  总长 %.0fmm  估时 %.1fs\n", wp.size(), total, total / P.speed);
    printf("      平面 %.0fmm / 孔周 %.0fmm / 斜面 %.0fmm / 过渡 %.0fmm\n",
           lenPlane, lenRing, lenSlope, lenTravel);
    printf("      表面贴合: %.1f .. %.1f mm(相对主面)\n", hmin, hmax);

    // ---------- 8. 输出 ----------
    std::filesystem::create_directories(outdir);
    std::string csvPath = outdir + "/trajectory.csv";
    FILE* f = fopen(csvPath.c_str(), "w");
    if (!f) { fprintf(stderr, "错误: 无法写入 %s\n", csvPath.c_str()); return 1; }
    fprintf(f, "seq,zone,type,x_mm,y_mm,z_mm,qw,qx,qy,qz,toolz_x,toolz_y,toolz_z,speed_mm_s,pass,group,u_mm,v_mm,h_mm,su_mm,sv_mm\n");
    for (size_t k = 0; k < wp.size(); ++k) {
        const auto& w = wp[k];
        Vec3 s = w.pos - w.nref * P.standoff;   // 表面投影（斜面段的枪投影会偏, 用表面点对照掩膜）
        Vec3 rel = s - pl.c;
        fprintf(f, "%zu,%s,%s,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.1f,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f\n",
                k, w.zone.c_str(), w.type.c_str(), w.pos.x, w.pos.y, w.pos.z,
                w.qw, w.qx, w.qy, w.qz, w.toolZ.x, w.toolZ.y, w.toolZ.z,
                w.speed, w.pass, w.group, w.u, w.v, w.h,
                dot(rel, uAxis), dot(rel, vAxis));
    }
    fclose(f);

    std::vector<Vec3> tpts;
    for (auto& w : wp) tpts.push_back(w.pos);
    saveTrajPLY(outdir + "/trajectory.ply", tpts);

    // 预览（2 px/mm）: 分区底色 + 分区着色的轨迹
    double scale = 2.0;
    int W = (int)((mk.W * cell) * scale) + 20, H = (int)((mk.H * cell) * scale) + 20;
    RGBImage img(W, H, 255, 255, 255);
    auto toPix = [&](double u, double v) -> std::pair<double,double> {
        // v 向上为正（面板顶部画在图上方）
        return {10 + (u - mk.umin) * scale,
                10 + (mk.vmin + mk.H * cell - v) * scale};
    };
    for (int j = 0; j < mkFull.H; ++j)
        for (int i = 0; i < mkFull.W; ++i) {
            if (!mkFull.at(i, j)) continue;
            size_t k = (size_t)j * mkFull.W + i;
            uint8_t cr, cg, cb;
            if (smallHoles[k] && mkTrav.at(i, j)) { cr = 252; cg = 206; cb = 120; }   // 忽略的小孔
            else if (steep[k])                    { cr = 216; cg = 190; cb = 236; }   // 斜面区
            else if (!mkTrav.at(i, j))            { cr = 205; cg = 208; cb = 216; }   // 边界内缩带
            else if (!mkMain.at(i, j))            { cr = 200; cg = 224; cb = 240; }   // 孔周环带
            else                                  { cr = 255; cg = 224; cb = 220; }
            auto [px, py] = toPix(mkFull.u(i), mkFull.v(j));
            for (int dy = 0; dy < (int)scale; ++dy)
                for (int dx = 0; dx < (int)scale; ++dx)
                    img.set((int)px + dx, (int)py + dy, cr, cg, cb);
        }
    auto zoneColor = [](const std::string& z, uint8_t& r, uint8_t& g, uint8_t& b) {
        if (z == "ring")       { r = 30;  g = 90;  b = 220; }
        else if (z == "slope") { r = 150; g = 40;  b = 190; }
        else if (z == "travel"){ r = 40;  g = 165; b = 95;  }
        else                   { r = 220; g = 40;  b = 40;  }
    };
    // 画"表面投影"(枪位置往回退 standoff), 斜面段的枪 2D 投影会偏到面板中部, 画枪投影会误导
    auto surfPix = [&](const Waypoint& w) -> std::pair<double,double> {
        Vec3 s = w.pos - w.nref * P.standoff;
        Vec3 rel = s - pl.c;
        return toPix(dot(rel, uAxis), dot(rel, vAxis));
    };
    std::pair<double,double> last(-1, -1);
    for (size_t k = 0; k < wp.size(); ++k) {
        auto p = surfPix(wp[k]);
        if (last.first >= 0) {
            uint8_t r, g, b;
            zoneColor(wp[k].zone, r, g, b);
            img.line(last.first, last.second, p.first, p.second, r, g, b);
        }
        last = p;
    }
    if (!wp.empty()) {
        auto p = surfPix(wp.front());
        img.disc(p.first, p.second, 5, 20, 160, 60);
    }
    saveBMP(outdir + "/preview.bmp", img);

    printf("[8/8] 输出:\n  %s\n  %s/trajectory.ply\n  %s/preview.bmp\n",
           csvPath.c_str(), outdir.c_str(), outdir.c_str());
    return 0;
}
