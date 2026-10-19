// driver.cpp —— 命令行测试：读 PLY -> 生成轨迹 -> 写 CSV + PLY
#include "spray_traj.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

static bool read_ply(const std::string& path, std::vector<double>& xyz){
    std::ifstream f(path, std::ios::binary);
    if(!f) return false;
    std::string line; int n=0; bool ascii=true;
    while(std::getline(f,line)){
        if(line.rfind("format",0)==0){ ascii = line.find("ascii")!=std::string::npos; }
        if(line.rfind("element vertex",0)==0){ sscanf(line.c_str(),"element vertex %d",&n); }
        if(line.rfind("end_header",0)==0) break;
    }
    if(n<=0) return false;
    xyz.reserve(n*3);
    if(ascii){
        // 面板.ply 每行 4 列: x y z intensity；通用处理：读第一个数是 x，需要知道列数
        // 简化：读到 end_header 后，先探测一行有几列
        std::string fl;
        std::getline(f, fl);
        // 数 token 数
        {
            std::stringstream ss(fl); std::string tok; int cols=0;
            while(ss>>tok) cols++;
            // 回退：无法回退，改为按行处理
        }
        // 重新按行读取（含刚才那行）
        auto parseLine=[&](const std::string& ln){
            std::stringstream ss(ln); std::vector<double> v; double t;
            while(ss>>t) v.push_back(t);
            if(v.size()>=3){ xyz.push_back(v[0]); xyz.push_back(v[1]); xyz.push_back(v[2]); }
        };
        parseLine(fl);
        for(int i=1;i<n;++i){ std::string ln; if(std::getline(f,ln)) parseLine(ln); }
    } else {
        return false;
    }
    return !xyz.empty();
}

int main(int argc, char** argv){
    // 关掉 stdout 全缓冲：算法一旦异常退出，缓冲里的日志会全部丢失，
    // 现场只剩"退出码非0、输出为空"，无法定位。无缓冲保证崩溃前日志已落盘。
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);
    const char* in = argc>1? argv[1] : "panel.ply";
    const char* outcsv = argc>2? argv[2] : "traj_out.csv";
    const char* outply = argc>3? argv[3] : "traj_out.ply";
    std::vector<double> xyz;
    if(!read_ply(in, xyz)){ printf("读取失败(仅支持ascii PLY): %s\n", in); return 1; }
    int n = (int)(xyz.size()/3);
    SprayTrajParams p = spray_traj_default_params();
    // 可选参数：-flat <平面拟合阈值mm>  -curve <曲面升级阈值mm>
    //           -fill <孔洞填补阈值mm>  -dir <0自动/1长边/2短边>
    for(int i=4;i+1<argc;++i){
        if(!strcmp(argv[i],"-flat"))  p.flatnessTol  = atof(argv[i+1]);
        if(!strcmp(argv[i],"-curve")) p.curvatureTol = atof(argv[i+1]);
        if(!strcmp(argv[i],"-fill"))  p.fillHoleSize = atof(argv[i+1]);
        if(!strcmp(argv[i],"-dir"))   p.scanDirMode  = atoi(argv[i+1]);
    }
    printf("参数: flatnessTol=%.1f curvatureTol=%.1f fillHoleSize=%.1f scanDirMode=%d\n",
        p.flatnessTol, p.curvatureTol, p.fillHoleSize, p.scanDirMode);
    SprayTrajResult r = spray_traj_generate(xyz.data(), n, &p);
    if(r.error){ printf("失败: %s\n", r.error); spray_traj_free(&r); return 1; }
    printf("OK: rows=%d points=%d rms=%.3f seg_mean=%.3f seg_max=%.3f\n",
        r.rows, r.count, r.plane_rms, r.seg_len_mean, r.seg_len_max);
    printf("法向=(%.4f,%.4f,%.4f) 内点=%d/%d(%.1f%%) 质心=(%.2f,%.2f,%.2f)\n",
        r.normal_x, r.normal_y, r.normal_z, r.inlier_points,
        r.input_points - r.removed_points,
        (double)r.inlier_points / std::max(1, r.input_points - r.removed_points) * 100.0,
        r.plane_cx, r.plane_cy, r.plane_cz);
    // CSV
    FILE* fc=fopen(outcsv,"w");
    fprintf(fc,"X,Y,Z,I,J,K,speed,spray\n");
    for(int i=0;i<r.count;++i){
        SprayTrajPoint& q=r.points[i];
        fprintf(fc,"%.3f,%.3f,%.3f,%.4f,%.4f,%.4f,%.1f,%d\n",
            q.x,q.y,q.z,q.i,q.j,q.k,q.speed,q.spray);
    }
    fclose(fc);
    // PLY
    FILE* fp=fopen(outply,"w");
    fprintf(fp,"ply\nformat ascii 1.0\ncomment spray trajectory\n");
    fprintf(fp,"element vertex %d\n",r.count);
    fprintf(fp,"property float x\nproperty float y\nproperty float z\n");
    fprintf(fp,"property uchar red\nproperty uchar green\nproperty uchar blue\nend_header\n");
    for(int i=0;i<r.count;++i){ fprintf(fp,"%.3f %.3f %.3f 255 48 32\n",r.points[i].x,r.points[i].y,r.points[i].z); }
    fclose(fp);
    // 重建面片网格导出（供离线检查：每个分割面一个 PLY，序号对应分面顺序）
    for(int i=0;i<r.plane_count && i<8;++i){
        int nv=r.planes[i].polyVertCount, nt=r.planes[i].polyTriCount;
        if(nv<=0||nt<=0) continue;
        char buf[512]; snprintf(buf,sizeof(buf),"%s.plane%d.ply",outply,i);
        FILE* fm=fopen(buf,"w");
        if(!fm) continue;
        fprintf(fm,"ply\nformat ascii 1.0\ncomment reconstructed plane %d\n",i);
        fprintf(fm,"element vertex %d\nproperty float x\nproperty float y\nproperty float z\n",nv);
        fprintf(fm,"element face %d\nproperty list uchar int vertex_index\nend_header\n",nt);
        for(int v=0;v<nv;++v)
            fprintf(fm,"%.3f %.3f %.3f\n", r.planes[i].polyVerts[v*3], r.planes[i].polyVerts[v*3+1], r.planes[i].polyVerts[v*3+2]);
        for(int t=0;t<nt;++t)
            fprintf(fm,"3 %d %d %d\n", r.planes[i].polyTris[t*3], r.planes[i].polyTris[t*3+1], r.planes[i].polyTris[t*3+2]);
        fclose(fm);
    }
    printf("写出: %s / %s\n", outcsv, outply);
    spray_traj_free(&r);
    return 0;
}
