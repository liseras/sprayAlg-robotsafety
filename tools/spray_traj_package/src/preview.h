// preview.h - 轨迹预览图输出（24 位 BMP，无依赖，便于快速目视检查）
#pragma once
#include <vector>
#include <cstdint>
#include <string>
#include <cstdio>
#include <cstring>
#include <cmath>

struct RGBImage {
    int W = 0, H = 0;
    std::vector<uint8_t> px;  // RGB，自上而下
    RGBImage(int w, int h, uint8_t r, uint8_t g, uint8_t b)
        : W(w), H(h), px((size_t)w * h * 3) {
        for (size_t k = 0; k < px.size(); k += 3) px[k] = r, px[k+1] = g, px[k+2] = b;
    }
    void set(int x, int y, uint8_t r, uint8_t g, uint8_t b) {
        if (x < 0 || x >= W || y < 0 || y >= H) return;
        size_t k = ((size_t)y * W + x) * 3;
        px[k] = r; px[k+1] = g; px[k+2] = b;
    }
    void disc(double cx, double cy, double rad, uint8_t r, uint8_t g, uint8_t b) {
        for (int y = (int)(cy - rad); y <= (int)(cy + rad); ++y)
            for (int x = (int)(cx - rad); x <= (int)(cx + rad); ++x)
                if ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= rad * rad) set(x, y, r, g, b);
    }
    void line(double x0, double y0, double x1, double y1,
              uint8_t r, uint8_t g, uint8_t b) {
        double dx = x1 - x0, dy = y1 - y0;
        int n = std::max(1, (int)std::ceil(std::sqrt(dx * dx + dy * dy)));
        for (int k = 0; k <= n; ++k)
            set((int)(x0 + dx * k / n), (int)(y0 + dy * k / n), r, g, b);
    }
};

inline bool saveBMP(const std::string& path, const RGBImage& img) {
    int rowBytes = (img.W * 3 + 3) / 4 * 4;
    uint32_t dataSize = (uint32_t)rowBytes * img.H;
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return false;
    uint8_t hdr[54] = {0};
    hdr[0]='B'; hdr[1]='M';
    uint32_t fileSize = 54 + dataSize;
    memcpy(hdr + 2, &fileSize, 4);
    uint32_t offset = 54; memcpy(hdr + 10, &offset, 4);
    uint32_t dib = 40; memcpy(hdr + 14, &dib, 4);
    int32_t w = img.W, h = img.H;
    memcpy(hdr + 18, &w, 4); memcpy(hdr + 22, &h, 4);
    uint16_t planes = 1, bpp = 24;
    memcpy(hdr + 26, &planes, 2); memcpy(hdr + 28, &bpp, 2);
    fwrite(hdr, 1, 54, f);
    std::vector<uint8_t> row(rowBytes, 0);
    for (int j = img.H - 1; j >= 0; --j) {   // BMP 自下而上
        memset(row.data(), 0, rowBytes);
        for (int i = 0; i < img.W; ++i) {
            size_t k = ((size_t)j * img.W + i) * 3;
            row[i * 3 + 0] = img.px[k + 2];  // BGR
            row[i * 3 + 1] = img.px[k + 1];
            row[i * 3 + 2] = img.px[k + 0];
        }
        fwrite(row.data(), 1, rowBytes, f);
    }
    fclose(f);
    return true;
}
