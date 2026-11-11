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
    // 整文件一次读入再解析。
    // 原实现每行建一个 std::stringstream、再 push 进一个临时 vector<double>，
    // 60 万行时要 1s 以上 —— 小件看不出来，大点云上这部分比算法本身还显眼。
    FILE* fp = fopen(path.c_str(), "rb");
    if(!fp) return false;
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    if(sz <= 0){ fclose(fp); return false; }
    std::string buf;
    buf.resize((size_t)sz);
    size_t got = fread(&buf[0], 1, (size_t)sz, fp);
    fclose(fp);
    buf.resize(got);

    // 解析头部：点数 + 编码 + vertex 属性表
    // binary 必须按属性偏移取 x/y/z —— 不能像 ascii 那样"每行前 3 个数"，
    // 现场点云常见 x y z intensity / x y z nx ny nz rgb 两种排布，混着来。
    int n=0; bool ascii=true, binLE=false; size_t dataStart=0;
    struct Prop { std::string name; int bytes; };
    std::vector<Prop> props;
    bool inVertex = false;
    {
        size_t pos = 0;
        while(pos < buf.size()){
            size_t nl = buf.find('\n', pos);
            if(nl == std::string::npos) nl = buf.size();
            std::string line = buf.substr(pos, nl-pos);
            if(!line.empty() && line[line.size()-1]=='\r') line.erase(line.size()-1);
            if(line.rfind("format",0)==0){
                ascii = line.find("ascii")!=std::string::npos;
                binLE = line.find("binary_little_endian")!=std::string::npos;
            } else if(line.rfind("element ",0)==0){
                inVertex = (line.rfind("element vertex",0)==0);
                if(inVertex) sscanf(line.c_str(),"element vertex %d",&n);
            } else if(line.rfind("property",0)==0 && inVertex){
                char type[64]={0}, name[64]={0};
                if(sscanf(line.c_str(),"property %63s %63s",type,name)==2){
                    std::string ts(type);
                    int b=4;
                    if(ts=="char"||ts=="uchar"||ts=="int8"||ts=="uint8") b=1;
                    else if(ts=="short"||ts=="ushort"||ts=="int16"||ts=="uint16") b=2;
                    else if(ts=="double"||ts=="float64") b=8;
                    Prop pr; pr.name = name; pr.bytes = b; props.push_back(pr);
                }
            } else if(line.rfind("end_header",0)==0){
                dataStart = (nl<buf.size())? nl+1 : buf.size(); break;
            }
            pos = nl + 1;
        }
    }
    if(n<=0) return false;

    if(ascii){
        // 数据区：每行取前 3 个数（面板件是 x y z intensity 四列）
        xyz.reserve((size_t)n*3);
        const char* q = buf.data() + dataStart;
        const char* end = buf.data() + buf.size();
        while(q < end && (int)(xyz.size()/3) < n){
            while(q < end && (*q==' '||*q=='\t'||*q=='\r'||*q=='\n')) ++q;
            if(q >= end) break;
            double v[3]; int c=0;
            while(c < 3){
                char* e2 = nullptr;
                double d = std::strtod(q, &e2);
                if(e2 == q) break;
                v[c++] = d; q = e2;
                while(q < end && (*q==' '||*q=='\t')) ++q;
            }
            if(c >= 3){ xyz.push_back(v[0]); xyz.push_back(v[1]); xyz.push_back(v[2]); }
            while(q < end && *q != '\n') ++q;   // 跳到行尾（跳过强度等多余列）
            if(q < end) ++q;
        }
    } else if(binLE){
        // 按属性表算步长与 x/y/z 的字节偏移，逐条记录 memcpy。
        // 只认 4/8 字节的坐标（float/double）；其余宽度按 0 跳过 —— 坐标是整数的点云不存在。
        int stride=0, off[3]={-1,-1,-1}, w[3]={4,4,4};
        for(size_t i=0;i<props.size();++i){
            const std::string& nm = props[i].name;
            if(nm=="x"){ off[0]=stride; w[0]=props[i].bytes; }
            else if(nm=="y"){ off[1]=stride; w[1]=props[i].bytes; }
            else if(nm=="z"){ off[2]=stride; w[2]=props[i].bytes; }
            stride += props[i].bytes;
        }
        if(off[0]<0||off[1]<0||off[2]<0||stride<=0) return false;
        const unsigned char* base = (const unsigned char*)buf.data() + dataStart;
        size_t avail = buf.size() - dataStart;
        xyz.reserve((size_t)n*3);
        for(int i=0;i<n;++i){
            size_t p0 = (size_t)i*(size_t)stride;
            if(p0 + (size_t)stride > avail) break;
            double v[3]={0,0,0};
            for(int k=0;k<3;++k){
                const unsigned char* src = base + p0 + (size_t)off[k];
                if(w[k]==4){ float f; memcpy(&f, src, 4); v[k]=f; }
                else if(w[k]==8){ double d; memcpy(&d, src, 8); v[k]=d; }
            }
            xyz.push_back(v[0]); xyz.push_back(v[1]); xyz.push_back(v[2]);
        }
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
    if(!read_ply(in, xyz)){ printf("读取失败(支持 ascii / binary_little_endian PLY): %s\n", in); return 1; }
    int n = (int)(xyz.size()/3);
    SprayTrajParams p = spray_traj_default_params();
    // 可选参数：-flat <平面拟合阈值mm>  -curve <曲面升级阈值mm>
    //           -fill <孔洞填补阈值mm>  -dir <0自动/1长边/2短边>
    //           -row <行距mm>  -pt <点距mm>  -so <离面距离mm>  -sp <速度mm/s>
    //           -gc <枪体安全半径mm>（= 喷头外径/2 + 安全余量，内壁角落让枪用）
    //           -ee <行端点外扩mm>（行方向；-1=自动 1×行距，0=不外扩）
    //           -ea <跨行外扩mm>（首末行位置；-1=自动 0.25×行距，0=压点云边界）
    //           -st <薄壳双面壁厚mm>（0=关闭；>0 时对"射线两侧都够伸枪"的单层面
    //                补一个反向面，代表这张板的另一面）
    for(int i=4;i+1<argc;++i){
        if(!strcmp(argv[i],"-flat"))  p.flatnessTol  = atof(argv[i+1]);
        if(!strcmp(argv[i],"-curve")) p.curvatureTol = atof(argv[i+1]);
        if(!strcmp(argv[i],"-fill"))  p.fillHoleSize = atof(argv[i+1]);
        if(!strcmp(argv[i],"-dir"))   p.scanDirMode  = atoi(argv[i+1]);
        if(!strcmp(argv[i],"-row"))   p.scanRowSpacing = atof(argv[i+1]);
        if(!strcmp(argv[i],"-pt"))    p.pointSpacing   = atof(argv[i+1]);
        if(!strcmp(argv[i],"-so"))    p.standoff       = atof(argv[i+1]);
        if(!strcmp(argv[i],"-sp"))    p.speed          = atof(argv[i+1]);
        if(!strcmp(argv[i],"-gc"))    p.gunClear       = atof(argv[i+1]);
        if(!strcmp(argv[i],"-ee"))    p.edgeExtend     = atof(argv[i+1]);
        if(!strcmp(argv[i],"-ea"))    p.edgeExtendAcross = atof(argv[i+1]);
        if(!strcmp(argv[i],"-st"))    p.shellThick     = atof(argv[i+1]);
        // -flip 1：把上游所有定向结果强制取反。用于压力测试"方向过滤"这道硬闸——
        // 人为制造一片朝向全错的输入，硬闸必须把能救的翻回来、救不回来的拒喷，
        // 绝不能出现"从背面伸枪"的轨迹。
        if(!strcmp(argv[i],"-flip"))  p.normalFlip     = atoi(argv[i+1]);
        // -up x,y,z：可喷侧参考方向。默认 0,0,-1 —— 相机系里相机只能看到朝向自己那一面，
        // 所以可喷侧恒为朝相机侧（更低 Z），这是事实不是参数，一般不用传。
        // 只有点云已被搬到别的坐标系（如手眼变换后的基座系，+Z 竖直向上）时，才用 0,0,1 覆盖。
        if(!strcmp(argv[i],"-up"))    _putenv_s("SPRAYTRAJ_UP", argv[i+1]);
    }
    printf("参数: flatnessTol=%.1f curvatureTol=%.1f fillHoleSize=%.1f scanDirMode=%d row=%.1f pt=%.1f so=%.1f sp=%.1f gc=%.1f ee=%.1f ea=%.1f st=%.1f\n",
        p.flatnessTol, p.curvatureTol, p.fillHoleSize, p.scanDirMode,
        p.scanRowSpacing, p.pointSpacing, p.standoff, p.speed, p.gunClear, p.edgeExtend, p.edgeExtendAcross, p.shellThick);
    // ABI 自检：C# 侧 Native.cs 必须逐字段一致，改结构体后对一下这三个数
    printf("[abi] sizeof Plane=%zu Point=%zu Result=%zu\n",
        sizeof(SprayTrajPlane), sizeof(SprayTrajPoint), sizeof(SprayTrajResult));
    SprayTrajResult r = spray_traj_generate(xyz.data(), n, &p);
    if(r.error){ printf("失败: %s\n", r.error); spray_traj_free(&r); return 1; }
    printf("OK: rows=%d points=%d rms=%.3f seg_mean=%.3f seg_max=%.3f\n",
        r.rows, r.count, r.plane_rms, r.seg_len_mean, r.seg_len_max);
    printf("法向=(%.4f,%.4f,%.4f) 内点=%d/%d(%.1f%%) 质心=(%.2f,%.2f,%.2f)\n",
        r.normal_x, r.normal_y, r.normal_z, r.inlier_points,
        r.input_points - r.removed_points,
        (double)r.inlier_points / std::max(1, r.input_points - r.removed_points) * 100.0,
        r.plane_cx, r.plane_cy, r.plane_cz);

    // 面级法向清单：nx,ny,nz = 喷距法向（轨迹按它抬升）；n0 = 射线自由空间判出的真实法向。
    // 两者点积为负 → 该面是内壁，喷距方向与真实可喷侧相反。
    // 口径：与 +Z 夹角 ≤ SPRAYTRAJ_NORMAL_WIN_DEG 的面全部列出（含内壁、外壁），
    //       完全朝下的面（如开口盒底板外侧）不列入。
    printf("\n---- 面级法向（+Z = 0°，列出与 +Z 夹角 ≤ %.0f° 的全部面：内壁+外壁）----\n",
           SPRAYTRAJ_NORMAL_WIN_DEG);
    printf("面  点数   喷距法向              真实法向              与+Z   tilt  分类 内壁 依据  自由+/-(mm)  层/壁厚  质心\n");
    int nWin = 0, nSkip = 0, nSingle = 0;
    for(int i=0;i<r.plane_count && i<SPRAYTRAJ_PLANE_SLOTS;++i){
        const SprayTrajPlane& q = r.planes[i];
        if(q.angZ > SPRAYTRAJ_NORMAL_WIN_DEG){ ++nSkip; continue; }   // 朝下的面不进清单
        ++nWin;
        const char* cls = q.side>0 ? "正向" : (q.side<0 ? "反向" : "竖直");
        const char* src = q.sideSrc==1 ? "射线" : (q.sideSrc==3 ? "强制"
                        : (q.sideSrc==4 ? "薄壁" : (q.sideSrc==5 ? "双面"
                        : (q.sideSrc==6 ? "同向" : (q.sideSrc==7 ? "枪位" : "兜底")))));
        const char* inn = q.isInner>0 ? "内 " : (q.isInner<0 ? "未分" : "外 ");
        double flip = q.nx*q.n0x + q.ny*q.n0y + q.nz*q.n0z;
        // 层/壁厚：单层壳 = 这一面在点云里只有一层点、没有壁厚。这种面谈不上
        // "有外壁"，清单里统计成"内壁 N / 外壁 0"会让人以为漏了 N 个外壁。
        // sideSrc==5 是薄板双面补出来的面：它也只有一个层点，但就是"另一面"，
        // 不能算进"单层壳无从区分"。
        const bool dualF = (q.sideSrc == 5);
        if(q.layerCnt == 1 && !dualF) ++nSingle;
        char ly[32];
        if(dualF)                snprintf(ly,sizeof(ly),"双面%.1f",q.layerThick);
        else if(q.layerCnt == 2) snprintf(ly,sizeof(ly),"2层%.1f",q.layerThick);
        else if(q.layerCnt == 1) snprintf(ly,sizeof(ly),"单层%.1f",q.layerThick);
        else                     snprintf(ly,sizeof(ly),"未判");
        printf("%2d %6d  (% .3f,% .3f,% .3f) (% .3f,% .3f,% .3f) %5.1f %6.1f  %s  %s  %s %5.1f/%-5.1f %-8s (%.1f,%.1f,%.1f)%s\n",
            i, q.count, q.nx,q.ny,q.nz, q.n0x,q.n0y,q.n0z, q.angZ, q.tiltDeg,
            cls, inn, src, q.freePlus, q.freeMinus, ly,
            q.cx,q.cy,q.cz, flip<0 ? "  <== 内外反" : "");
    }
    printf("窗口内 %d 面 / 分割共 %d 面", nWin, r.plane_count);
    if(nSkip) printf("，滤除 %d 面（与 +Z 夹角 > %.0f°，朝下）", nSkip, SPRAYTRAJ_NORMAL_WIN_DEG);
    printf("\n");
    if(nWin > 0 && nSingle == nWin)
        printf("  ※ 清单内全部为单层壳：点云在每个面只有一层点、没有壁厚，"
               "数据里不含第二层，内/外壁无从区分\n");

    // 面级清单落盘（Excel 可直接看）
    {
        char buf[512]; snprintf(buf,sizeof(buf),"%s.faces.csv",outcsv);
        FILE* ff=fopen(buf,"w");
        if(ff){
            fprintf(ff,"face,count,nx,ny,nz,n0x,n0y,n0z,angZ,tiltDeg,side,isInner(1内/0外/-1未分),sideSrc(1射线/2兜底/3强制/4薄壁/5双面/6拆层同向/7枪位验收翻转),freePlus,freeMinus,cx,cy,cz,layerCnt(1单层/2双层/0未判定),layerThick,inWin(1=夹角<=%.0f)\n", SPRAYTRAJ_NORMAL_WIN_DEG);
            for(int i=0;i<r.plane_count && i<SPRAYTRAJ_PLANE_SLOTS;++i){
                const SprayTrajPlane& q=r.planes[i];
                fprintf(ff,"%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%.2f,%d,%d,%d,%.2f,%.2f,%.3f,%.3f,%.3f,%d,%.2f,%d\n",
                    i,q.count,q.nx,q.ny,q.nz,q.n0x,q.n0y,q.n0z,q.angZ,q.tiltDeg,
                    q.side,q.isInner,q.sideSrc,q.freePlus,q.freeMinus,q.cx,q.cy,q.cz,
                    q.layerCnt,q.layerThick,
                    q.angZ <= SPRAYTRAJ_NORMAL_WIN_DEG ? 1 : 0);
            }
            fclose(ff);
        }
    }
    // CSV
    FILE* fc=fopen(outcsv,"w");
    if(!fc){ printf("无法写出 CSV: %s\n", outcsv); return 1; }
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
    // 上限必须是槽位数，不能写死 8：薄壁件拆层后面数会超过 8，写死会静默丢掉后几个面。
    for(int i=0;i<r.plane_count && i<SPRAYTRAJ_PLANE_SLOTS;++i){
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
