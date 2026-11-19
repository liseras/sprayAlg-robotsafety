// ply_io.h - PLY 点云读取（ascii / binary_little_endian，支持任意属性布局）
#pragma once
#include <string>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include "vecmath.h"

struct PointCloud {
    std::vector<Vec3> pts;      // 位置
    std::vector<float> scalar;  // 第一个附加标量属性（如 intensity），无则为空
};

inline bool loadPLY(const std::string& path, PointCloud& cloud, std::string& err) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) { err = "cannot open " + path; return false; }

    std::string format;
    long long vertexCount = 0;
    struct Prop { std::string name, type; };
    std::vector<Prop> props;
    bool inVertex = false;

    char line[1024];
    long long headerEndPos = 0;
    while (fgets(line, sizeof(line), f)) {
        std::string s = line;
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
        if (s.rfind("format", 0) == 0) {
            size_t p1 = s.find(' ');
            size_t p2 = s.find(' ', p1 + 1);
            format = s.substr(p1 + 1, p2 - p1 - 1);
        } else if (s.rfind("element", 0) == 0) {
            inVertex = (s.find("vertex") != std::string::npos);
            if (inVertex) sscanf(s.c_str(), "element vertex %lld", &vertexCount);
        } else if (s.rfind("property", 0) == 0 && inVertex) {
            char type[64], name[256];
            if (sscanf(s.c_str(), "property %63s %255s", type, name) == 2)
                props.push_back({name, type});
        } else if (s.rfind("end_header", 0) == 0) {
            headerEndPos = ftell(f);
            break;
        }
    }
    if (format.empty() || vertexCount <= 0) { fclose(f); err = "invalid PLY header"; return false; }

    int ix = -1, iy = -1, iz = -1, is = -1;
    for (size_t i = 0; i < props.size(); ++i) {
        if (props[i].name == "x") ix = (int)i;
        if (props[i].name == "y") iy = (int)i;
        if (props[i].name == "z") iz = (int)i;
        if (is < 0 && props[i].name != "x" && props[i].name != "y" && props[i].name != "z") is = (int)i;
    }
    if (ix < 0 || iy < 0 || iz < 0) { fclose(f); err = "no x/y/z properties"; return false; }

    auto typeSize = [](const std::string& t) -> int {
        if (t == "char" || t == "int8" || t == "uchar" || t == "uint8") return 1;
        if (t == "short" || t == "int16" || t == "ushort" || t == "uint16") return 2;
        if (t == "int" || t == "int32" || t == "uint" || t == "uint32" || t == "float" || t == "float32") return 4;
        if (t == "double" || t == "float64") return 8;
        return 0;
    };

    cloud.pts.clear(); cloud.scalar.clear();
    cloud.pts.reserve(vertexCount);

    if (format == "ascii") {
        for (long long i = 0; i < vertexCount; ++i) {
            if (!fgets(line, sizeof(line), f)) break;
            double v[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            int nread = 0;
            char* p = line;
            for (int k = 0; k < (int)props.size() && nread < 8; ++k) {
                char* endp = nullptr;
                double val = strtod(p, &endp);
                if (endp == p) break;
                v[nread++] = val;
                p = endp;
            }
            cloud.pts.push_back(Vec3(v[ix], v[iy], v[iz]));
            if (is >= 0 && is < 8) cloud.scalar.push_back((float)v[is]);
        }
    } else if (format == "binary_little_endian" || format == "binary_big_endian") {
        int stride = 0;
        for (auto& pr : props) stride += typeSize(pr.type);
        std::vector<char> rec(stride);
        for (long long i = 0; i < vertexCount; ++i) {
            if (fread(rec.data(), 1, stride, f) != (size_t)stride) break;
            auto readProp = [&](int idx) -> double {
                const char* base = rec.data();
                for (int k = 0; k < idx; ++k) base += typeSize(props[k].type);
                double val = 0;
                if (format == "binary_little_endian") {
                    switch (typeSize(props[idx].type)) {
                        case 1: val = *(int8_t*)base; break;
                        case 2: val = *(int16_t*)base; break;
                        case 4: if (props[idx].type[0] == 'f') val = *(float*)base; else val = *(int32_t*)base; break;
                        case 8: val = *(double*)base; break;
                    }
                }
                return val;
            };
            cloud.pts.push_back(Vec3(readProp(ix), readProp(iy), readProp(iz)));
            if (is >= 0) cloud.scalar.push_back((float)readProp(is));
        }
    } else {
        fclose(f); err = "unsupported format " + format; return false;
    }
    fclose(f);
    return true;
}

// 输出轨迹点 PLY（便于在三维查看器里对照检查）
inline bool saveTrajPLY(const std::string& path, const std::vector<Vec3>& pts) {
    FILE* f = fopen(path.c_str(), "w");
    if (!f) return false;
    fprintf(f, "ply\nformat ascii 1.0\nelement vertex %zu\nproperty float x\nproperty float y\nproperty float z\nend_header\n",
            pts.size());
    for (auto& p : pts) fprintf(f, "%.4f %.4f %.4f\n", p.x, p.y, p.z);
    fclose(f);
    return true;
}
