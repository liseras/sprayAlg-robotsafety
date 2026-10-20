// safety_driver.cpp —— 机器人安全层命令行驱动
//
// 用法：
//   spray_safety <traj.csv> <cloud.ply> <out.csv> [选项]
//
// 选项：
//   -standoff <mm>   喷距（必须与生成轨迹时一致，默认 100）
//   -safety <mm>     安全间隙（默认 20）
//   -search <mm>     判定带（默认 120；超过只报下界）
//   -tool <mm>       法兰→喷嘴（默认 150）
//   -tilt <deg>      非朝上面基准倾角（默认 20）
//   -tiltmax <deg>   倾角上限（默认 60）
//   -rate <deg>      每点倾角最大变化（默认 3）
//   -sprayang <deg>  工艺喷角上限（默认 30）
//   -base <x> <y> <z> <rx> <ry> <rz>   基座位姿（mm / 度，ZYX 欧拉）；不给则自动假设
//   -quiet           只打印结论
//
// 退出码：0 = 无 P0（可执行）；1 = 有 P0（撞了 / 不可达）；2 = 参数或读文件错误
//
// 为什么单独一个 exe 而不是改 driver.cpp：driver.cpp 是轨迹生成的测试入口，
// 它的输出决定 9 个非曲面样本的 CSV 基线。安全层不能有一丝机会碰到那条链。

#include "robot_safety.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------- PLY

struct PlyCloud {
    std::vector<double> xyz;
    int n = 0;
    std::string error;
};

bool readPly(const char* path, PlyCloud& out) {
    FILE* fp = fopen(path, "rb");
    if (!fp) {
        out.error = "打不开文件";
        return false;
    }
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    if (sz <= 0) {
        fclose(fp);
        out.error = "文件为空";
        return false;
    }
    std::string buf;
    buf.resize((size_t)sz);
    size_t got = fread(&buf[0], 1, (size_t)sz, fp);
    fclose(fp);
    buf.resize(got);

    bool ascii = false, binLE = false;
    long nVert = 0;
    size_t dataStart = 0;
    std::vector<std::pair<std::string, int> > props;   // (name, bytes)

    size_t pos = 0;
    while (pos < buf.size()) {
        size_t nl = buf.find('\n', pos);
        if (nl == std::string::npos) nl = buf.size();
        std::string line = buf.substr(pos, nl - pos);
        if (!line.empty() && line[line.size() - 1] == '\r') line.erase(line.size() - 1);
        if (line.rfind("format", 0) == 0) {
            ascii = line.find("ascii") != std::string::npos;
            binLE = line.find("binary_little_endian") != std::string::npos;
        } else if (line.rfind("element vertex", 0) == 0) {
            sscanf(line.c_str(), "element vertex %ld", &nVert);
        } else if (line.rfind("property", 0) == 0) {
            char type[64] = {0}, name[64] = {0};
            if (sscanf(line.c_str(), "property %63s %63s", type, name) == 2) {
                std::string ts(type);
                int bytes = 4;
                if (ts == "char" || ts == "uchar" || ts == "int8" || ts == "uint8") bytes = 1;
                else if (ts == "short" || ts == "ushort" || ts == "int16" || ts == "uint16") bytes = 2;
                else if (ts == "int" || ts == "uint" || ts == "int32" || ts == "uint32" || ts == "float" ||
                         ts == "float32") bytes = 4;
                else if (ts == "double" || ts == "float64") bytes = 8;
                props.push_back(std::make_pair(std::string(name), bytes));
            }
        } else if (line.rfind("end_header", 0) == 0) {
            dataStart = (nl < buf.size()) ? nl + 1 : buf.size();
            break;
        }
        pos = nl + 1;
    }
    if (nVert <= 0 || props.empty()) {
        out.error = "PLY 头解析失败（缺 vertex 或 property）";
        return false;
    }

    int ix = -1, iy = -1, iz = -1;
    for (size_t i = 0; i < props.size(); ++i) {
        if (props[i].first == "x") ix = (int)i;
        else if (props[i].first == "y") iy = (int)i;
        else if (props[i].first == "z") iz = (int)i;
    }
    if (ix < 0 || iy < 0 || iz < 0) {
        out.error = "PLY 缺 x/y/z 属性";
        return false;
    }

    out.xyz.reserve((size_t)nVert * 3);
    if (ascii) {
        const char* q = buf.data() + dataStart;
        const char* end = buf.data() + buf.size();
        int nProp = (int)props.size();
        while (q < end && out.n < nVert) {
            double v[3] = {0, 0, 0};
            bool ok = true;
            for (int k = 0; k < nProp; ++k) {
                while (q < end && (*q == ' ' || *q == '\t' || *q == '\r' || *q == '\n')) ++q;
                if (q >= end) {
                    ok = false;
                    break;
                }
                char* e2 = nullptr;
                double d = std::strtod(q, &e2);
                if (e2 == q) {
                    ok = false;
                    break;
                }
                if (k == ix) v[0] = d;
                else if (k == iy) v[1] = d;
                else if (k == iz) v[2] = d;
                q = e2;
            }
            if (!ok) break;
            out.xyz.push_back(v[0]);
            out.xyz.push_back(v[1]);
            out.xyz.push_back(v[2]);
            ++out.n;
            while (q < end && *q != '\n') ++q;
            if (q < end) ++q;
        }
    } else if (binLE) {
        int stride = 0;
        int off[3] = {0, 0, 0};
        for (int i = 0; i < (int)props.size(); ++i) {
            if (props[i].first == "x") off[0] = stride;
            else if (props[i].first == "y") off[1] = stride;
            else if (props[i].first == "z") off[2] = stride;
            stride += props[i].second;
        }
        // 记录坐标类型宽度（浮点也可能是 double）
        int w[3] = {4, 4, 4};
        for (int i = 0; i < (int)props.size(); ++i) {
            if (props[i].first == "x") w[0] = props[i].second;
            else if (props[i].first == "y") w[1] = props[i].second;
            else if (props[i].first == "z") w[2] = props[i].second;
        }
        const unsigned char* base = (const unsigned char*)buf.data() + dataStart;
        size_t avail = buf.size() - dataStart;
        for (long i = 0; i < nVert; ++i) {
            if ((size_t)(i * stride + stride) > avail) break;
            double v[3];
            for (int k = 0; k < 3; ++k) {
                const unsigned char* p = base + (size_t)i * stride + off[k];
                if (w[k] == 4) {
                    float f;
                    memcpy(&f, p, 4);
                    v[k] = f;
                } else if (w[k] == 8) {
                    double d;
                    memcpy(&d, p, 8);
                    v[k] = d;
                } else {
                    v[k] = 0.0;
                }
            }
            out.xyz.push_back(v[0]);
            out.xyz.push_back(v[1]);
            out.xyz.push_back(v[2]);
            ++out.n;
        }
    } else {
        out.error = "只支持 ascii / binary_little_endian PLY";
        return false;
    }

    if (out.n == 0) {
        out.error = "没读到任何点";
        return false;
    }
    return true;
}

// ---------------------------------------------------------------- 轨迹 CSV

struct TrajCsv {
    std::vector<SprayTrajPoint> pts;
    std::vector<std::string> header;
    std::vector<std::vector<std::string> > rawRows;   // 原始 token，用于原样带出额外列
    std::string error;
};

bool readTraj(const char* path, TrajCsv& out) {
    FILE* fp = fopen(path, "rb");
    if (!fp) {
        out.error = "打不开轨迹文件";
        return false;
    }
    std::string all;
    char chunk[65536];
    size_t got;
    while ((got = fread(chunk, 1, sizeof(chunk), fp)) > 0) all.append(chunk, got);
    fclose(fp);
    // 去 BOM
    if (all.size() >= 3 && (unsigned char)all[0] == 0xEF && (unsigned char)all[1] == 0xBB &&
        (unsigned char)all[2] == 0xBF)
        all.erase(0, 3);

    std::vector<std::string> lines;
    {
        size_t pos = 0;
        while (pos <= all.size()) {
            size_t nl = all.find('\n', pos);
            std::string ln = (nl == std::string::npos) ? all.substr(pos) : all.substr(pos, nl - pos);
            while (!ln.empty() && (ln[ln.size() - 1] == '\r' || ln[ln.size() - 1] == ' ')) ln.erase(ln.size() - 1);
            if (!ln.empty()) lines.push_back(ln);
            if (nl == std::string::npos) break;
            pos = nl + 1;
        }
    }
    if (lines.size() < 2) {
        out.error = "轨迹文件为空或只有表头";
        return false;
    }

    auto split = [](const std::string& s, char sep) {
        std::vector<std::string> v;
        size_t pos = 0;
        while (true) {
            size_t p = s.find(sep, pos);
            if (p == std::string::npos) {
                v.push_back(s.substr(pos));
                break;
            }
            v.push_back(s.substr(pos, p - pos));
            pos = p + 1;
        }
        return v;
    };

    out.header = split(lines[0], ',');
    for (auto& h : out.header) {
        size_t a = h.find_first_not_of(" \t");
        size_t b = h.find_last_not_of(" \t");
        h = (a == std::string::npos) ? std::string() : h.substr(a, b - a + 1);
    }
    int idx[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
    const char* want[8] = {"x", "y", "z", "i", "j", "k", "speed", "spray"};
    for (int w = 0; w < 8; ++w) {
        for (size_t c = 0; c < out.header.size(); ++c) {
            std::string h = out.header[c];
            for (auto& ch : h) ch = (char)tolower((unsigned char)ch);
            if (h == want[w]) {
                idx[w] = (int)c;
                break;
            }
        }
    }
    for (int w = 0; w < 6; ++w) {
        if (idx[w] < 0) {
            out.error = std::string("轨迹缺列 ") + want[w];
            return false;
        }
    }

    for (size_t r = 1; r < lines.size(); ++r) {
        std::vector<std::string> tok = split(lines[r], ',');
        if (tok.size() < 6) continue;
        SprayTrajPoint p;
        memset(&p, 0, sizeof(p));
        p.x = atof(tok[idx[0]].c_str());
        p.y = atof(tok[idx[1]].c_str());
        p.z = atof(tok[idx[2]].c_str());
        p.i = atof(tok[idx[3]].c_str());
        p.j = atof(tok[idx[4]].c_str());
        p.k = atof(tok[idx[5]].c_str());
        p.speed = idx[6] >= 0 && idx[6] < (int)tok.size() ? atof(tok[idx[6]].c_str()) : 200.0;
        p.spray = idx[7] >= 0 && idx[7] < (int)tok.size() ? atoi(tok[idx[7]].c_str()) : 1;
        out.pts.push_back(p);
        out.rawRows.push_back(tok);
    }
    if (out.pts.empty()) {
        out.error = "轨迹没有数据行";
        return false;
    }
    return true;
}

// 最近碰撞部位名（与控制台/CSV 共用；robot_safety.cpp 内部同名函数不对外）
std::string segName(int kind, int j0, int j1) {
    char b[48];
    switch (kind) {
        case RS_SEG_GRIP: return "工具握把";
        case RS_SEG_PIPE: return "工具细管";
        case RS_SEG_NOZZLE: return "喷嘴";
        case RS_SEG_BASE_COL: return "底座立柱";
        case RS_SEG_BASE_NECK: return "底座颈";
        case RS_SEG_JOINT: snprintf(b, sizeof(b), "关节J%d", j0); return b;
        case RS_SEG_LINK: snprintf(b, sizeof(b), "连杆J%d-J%d", j0, j1); return b;
        default: return "—";
    }
}

const char* violText(unsigned int v, std::string& buf) {
    buf.clear();
    for (int i = 0; i < RS_VIOL_COUNT; ++i) {
        if (v & (1u << i)) {
            const char* nm = robot_safety_rule_name(1u << i);
            if (!nm) continue;
            if (!buf.empty()) buf += "|";
            buf += nm;
        }
    }
    return buf.c_str();
}

}  // namespace

int main(int argc, char** argv) {
    setvbuf(stdout, NULL, _IONBF, 0);
    // --abi：打印结构体尺寸/偏移，供 C# 侧 Marshal.SizeOf 对照。
    // 跨语言 ABI 出错的表现是"字段全错位、数值离谱但不崩"，很难查，所以留一个自检入口。
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--abi")) {
            printf("RobotSafetyParams  size=%d\n", (int)sizeof(RobotSafetyParams));
            printf("RobotSafetyPoint   size=%d\n", (int)sizeof(RobotSafetyPoint));
            printf("RobotSafetyResult  size=%d\n", (int)sizeof(RobotSafetyResult));
            printf("  Point.x=%d tiltDeg=%d theta=%d spray=%d violation=%d\n",
                   (int)offsetof(RobotSafetyPoint, x), (int)offsetof(RobotSafetyPoint, tiltDeg),
                   (int)offsetof(RobotSafetyPoint, theta), (int)offsetof(RobotSafetyPoint, spray),
                   (int)offsetof(RobotSafetyPoint, violation));
            printf("  Res.byRule=%d basePose=%d error=%d\n", (int)offsetof(RobotSafetyResult, byRule),
                   (int)offsetof(RobotSafetyResult, basePose), (int)offsetof(RobotSafetyResult, error));
            return 0;
        }
    }
    if (argc < 4) {
        printf("用法: spray_safety <traj.csv> <cloud.ply> <out.csv> [选项]\n");
        printf("      -standoff -safety -search -tool -tilt -tiltmax -rate -sprayang\n");
        printf("      -threads <n>  -stride <n>  -voxel <mm>  -base <x y z rx ry rz>  -quiet\n");
        return 2;
    }
    const char* trajPath = argv[1];
    const char* cloudPath = argv[2];
    const char* outPath = argv[3];

    RobotSafetyParams p = robot_safety_default_params();
    double basePose[16];
    bool hasBase = false;
    bool quiet = false;

    for (int i = 4; i < argc; ++i) {
        if (!strcmp(argv[i], "-standoff") && i + 1 < argc) p.standoff = atof(argv[++i]);
        else if (!strcmp(argv[i], "-safety") && i + 1 < argc) p.safety = atof(argv[++i]);
        else if (!strcmp(argv[i], "-search") && i + 1 < argc) p.search = atof(argv[++i]);
        else if (!strcmp(argv[i], "-tool") && i + 1 < argc) p.toolLen = atof(argv[++i]);
        else if (!strcmp(argv[i], "-tilt") && i + 1 < argc) p.tiltNonUpDeg = atof(argv[++i]);
        else if (!strcmp(argv[i], "-tiltmax") && i + 1 < argc) p.tiltMaxDeg = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rate") && i + 1 < argc) p.tiltRateDeg = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sprayang") && i + 1 < argc) p.maxSprayAngleDeg = atof(argv[++i]);
        else if (!strcmp(argv[i], "-threads") && i + 1 < argc) p.threads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-stride") && i + 1 < argc) p.linkStride = std::max(1, atoi(argv[++i]));
        else if (!strcmp(argv[i], "-voxel") && i + 1 < argc) p.voxelMm = atof(argv[++i]);
        else if (!strcmp(argv[i], "-quiet")) quiet = true;
        else if (!strcmp(argv[i], "-base") && i + 6 < argc) {
            double x = atof(argv[i + 1]), y = atof(argv[i + 2]), z = atof(argv[i + 3]);
            double rx = atof(argv[i + 4]) * 3.14159265358979323846 / 180.0;
            double ry = atof(argv[i + 5]) * 3.14159265358979323846 / 180.0;
            double rz = atof(argv[i + 6]) * 3.14159265358979323846 / 180.0;
            i += 6;
            // R = Rz·Ry·Rx（与 DENSO P 型数据姿态约定一致）
            double cz = cos(rz), sz_ = sin(rz), cy = cos(ry), sy_ = sin(ry), cx = cos(rx), sx = sin(rx);
            double R[3][3];
            R[0][0] = cz * cy;
            R[0][1] = cz * sy_ * sx - sz_ * cx;
            R[0][2] = cz * sy_ * cx + sz_ * sx;
            R[1][0] = sz_ * cy;
            R[1][1] = sz_ * sy_ * sx + cz * cx;
            R[1][2] = sz_ * sy_ * cx - cz * sx;
            R[2][0] = -sy_;
            R[2][1] = cy * sx;
            R[2][2] = cy * cx;
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) basePose[r * 4 + c] = R[r][c];
            }
            basePose[3] = basePose[7] = basePose[11] = 0.0;
            basePose[12] = basePose[13] = basePose[14] = 0.0;
            basePose[15] = 1.0;
            basePose[3] = x;
            basePose[7] = y;
            basePose[11] = z;
            hasBase = true;
        }
    }

    TrajCsv traj;
    if (!readTraj(trajPath, traj)) {
        printf("轨迹读取失败: %s (%s)\n", traj.error.c_str(), trajPath);
        return 2;
    }
    PlyCloud cloud;
    if (!readPly(cloudPath, cloud)) {
        printf("点云读取失败: %s (%s)\n", cloud.error.c_str(), cloudPath);
        return 2;
    }
    if (!quiet) {
        printf("轨迹 %d 点 · 点云 %d 点\n", (int)traj.pts.size(), cloud.n);
        printf("参数: standoff=%.1f safety=%.1f search=%.1f tool=%.1f tilt=%.1f/%.1f rate=%.1f sprayAngle<=%.1f\n",
               p.standoff, p.safety, p.search, p.toolLen, p.tiltNonUpDeg, p.tiltMaxDeg, p.tiltRateDeg,
               p.maxSprayAngleDeg);
    }

    RobotSafetyResult r = robot_safety_check(traj.pts.data(), (int)traj.pts.size(), cloud.xyz.data(), cloud.n,
                                             hasBase ? basePose : NULL, &p);
    if (r.error) {
        printf("失败: %s\n", r.error);
        robot_safety_free(&r);
        return 2;
    }

    // ---- 输出 CSV：前 8 列按 master 格式重写，中间额外列原样带出，新列追加
    FILE* fc = fopen(outPath, "w");
    if (!fc) {
        printf("写不了 %s\n", outPath);
        robot_safety_free(&r);
        return 2;
    }
    int nCols = (int)traj.header.size();
    fprintf(fc, "X,Y,Z,I,J,K,speed,spray");
    for (int c = 8; c < nCols; ++c) fprintf(fc, ",%s", traj.header[c].c_str());
    fprintf(fc, ",tilt_deg,clearance_mm,clearance_tool_mm,clearance_link_mm,worst_part,link_part,"
                "self_gap_mm,margin_deg,standoff_err_mm,lifted_mm,cond,reverted,"
                "J1_deg,J2_deg,J3_deg,J4_deg,J5_deg,J6_deg,violation\n");
    std::string vbuf;
    for (int i = 0; i < r.count; ++i) {
        const RobotSafetyPoint& q = r.points[i];
        fprintf(fc, "%.3f,%.3f,%.3f,%.4f,%.4f,%.4f,%.1f,%d", q.x, q.y, q.z, q.i, q.j, q.k, q.speed, q.spray);
        for (int c = 8; c < nCols; ++c) {
            const std::vector<std::string>& tok = traj.rawRows[i];
            fprintf(fc, ",%s", c < (int)tok.size() ? tok[c].c_str() : "");
        }
        // MSVC 的 %f 遇到 inf 会打出 "-1.#INF00"，CSV 直接不可解析；
        // 不可达点的 clearance 就是 -inf，这里必须显式写成 nan。
        auto num = [](double v, int prec, char* buf, size_t n) {
            if (!std::isfinite(v)) snprintf(buf, n, "nan");
            else snprintf(buf, n, "%.*f", prec, v);
            return buf;
        };
        char b1[32], b2[32], b3[32], b4[32], b5[32], b6[32], b7[32], b8[32], b9[32];
        fprintf(fc, ",%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,", num(q.tiltDeg, 2, b1, sizeof(b1)),
                num(q.clearanceMm, 2, b2, sizeof(b2)), num(q.clearanceToolMm, 2, b9, sizeof(b9)),
                num(q.clearanceLinkMm, 2, b8, sizeof(b8)), segName(q.worstKind, q.worstJ0, q.worstJ1).c_str(),
                segName(q.linkKind, q.linkJ0, q.linkJ1).c_str(),
                num(q.selfGapMm, 2, b3, sizeof(b3)),
                num(q.marginDeg, 2, b4, sizeof(b4)), num(q.standoffErrMm, 6, b5, sizeof(b5)),
                num(q.liftedMm, 1, b6, sizeof(b6)), num(q.cond, 1, b7, sizeof(b7)), q.reverted);
        for (int k = 0; k < 6; ++k) fprintf(fc, "%.3f,", q.theta[k]);
        fprintf(fc, "%s\n", violText(q.violation, vbuf));
    }
    fclose(fc);

    // 基座位姿落盘：对拍 DensoArmKine 时要用同一个 T_base_in_cloud，
    // 否则"两边间隙不一致"会被误判成算法差异。
    {
        std::string bp = std::string(outPath) + ".base.json";
        FILE* fj = fopen(bp.c_str(), "w");
        if (fj) {
            fprintf(fj, "{\n  \"T_base_in_cloud\": [\n");
            for (int row = 0; row < 4; ++row) {
                fprintf(fj, "    [%.9f, %.9f, %.9f, %.9f]%s\n", r.basePose[row * 4 + 0],
                        r.basePose[row * 4 + 1], r.basePose[row * 4 + 2], r.basePose[row * 4 + 3],
                        row < 3 ? "," : "");
            }
            fprintf(fj, "  ],\n  \"auto\": %d,\n  \"tool_len\": %.3f,\n  \"standoff\": %.3f,\n  \"search\": %.3f\n}\n",
                    r.baseAuto, p.toolLen, p.standoff, p.search);
            fclose(fj);
        }
    }

    // ---- 结论
    int p0 = 0;
    for (int b = 0; b < RS_VIOL_COUNT; ++b) {
        const char* nm = robot_safety_rule_name(1u << b);
        // P0：不可达 / 自碰撞 / 连杆撞 / 枪体撞。COLLISION(1) 已被 LINK_HIT/TOOL_HIT 细分，
        // 本层不再置位，所以这里不把它算进 P0（否则 byRule[1] 恒为 0，算了也没意义）。
        bool isP0 = (b == 0 || b == 3 || b == 8 || b == 9);
        if (isP0) p0 += r.byRule[b];
        if (r.byRule[b] > 0 && !quiet) {
            printf("  %-15s %7d  [%s]\n", nm ? nm : "?", r.byRule[b], isP0 ? "P0" : "P1/P2");
        }
    }
    const char* verdict = (p0 > 0) ? "REJECT" : (r.violations > 0 ? "PASS_WITH_WARNINGS" : "PASS");
    printf("结论: %s  (P0=%d)\n", verdict, p0);
    printf("  %d 点 · 可达 %d (%.1f%%) · 倾斜 %d 点 (均值 %.1f° p95 %.1f° 最大 %.1f°)\n", r.count, r.reachable,
           r.count ? 100.0 * r.reachable / r.count : 0.0, r.tiltedPoints, r.tiltMeanDeg, r.tiltP95Deg,
           r.tiltMaxDeg);
    if (std::isfinite(r.minClearanceMm))
        printf("  最小表面距离 %.2fmm · 平滑回退 %d 点 · 喷距偏差 %.6fmm\n", r.minClearanceMm, r.poseReverted,
               r.standoffMaxErrMm);
    else
        printf("  最小表面距离 —（全部点都不可达）· 平滑回退 %d 点\n", r.poseReverted);

    // 枪体 / 连杆分开报：这两类撞车的处置完全不同（枪体撞→改喷距/角度；
    // 连杆撞→改站位或整条轨迹的姿态），混成一个数字现场没法判断。
    {
        char ct[32], cl[32];
        if (std::isfinite(r.minClearanceToolMm)) snprintf(ct, sizeof(ct), "%.2fmm", r.minClearanceToolMm);
        else snprintf(ct, sizeof(ct), "—");
        if (std::isfinite(r.minClearanceLinkMm)) snprintf(cl, sizeof(cl), "%.2fmm", r.minClearanceLinkMm);
        else snprintf(cl, sizeof(cl), "—");
        printf("  枪体最小表面距离 %s · 连杆最小表面距离 %s\n", ct, cl);
    }
    // 最危险连杆段（现场最关心：哪一段连杆几乎蹭到工件）
    {
        double worst = 1e300;
        int wk = -1, wj0 = -1, wj1 = -1, wi = -1;
        for (int i = 0; i < r.count; ++i) {
            const RobotSafetyPoint& q = r.points[i];
            if (!std::isfinite(q.clearanceLinkMm)) continue;
            if (q.clearanceLinkMm < worst) { worst = q.clearanceLinkMm; wk = q.linkKind; wj0 = q.linkJ0; wj1 = q.linkJ1; wi = i; }
        }
        if (wi >= 0)
            printf("  连杆/关节最近：%s  %+.2fmm（idx %d）\n", segName(wk, wj0, wj1).c_str(), worst, wi);
    }
    if (r.baseAuto) printf("  ⚠ 基座位姿为自动假设（未给 -base），可达率与连杆碰撞结论依赖它\n");
    if (r.l0Issues) printf("  ⚠ L0 输入层自检有 %d 项问题（详见 -standoff / 臂展），结论需谨慎\n", r.l0Issues);
    if (r.standoffMedianMm > 0) printf("  反推喷距中位数 %.2fmm（沿枪轴投影）\n", r.standoffMedianMm);
    printf("写出: %s\n", outPath);

    int rc = (p0 > 0) ? 1 : 0;
    robot_safety_free(&r);
    return rc;
}
