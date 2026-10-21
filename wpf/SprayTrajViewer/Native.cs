using System;
using System.Runtime.InteropServices;

namespace SprayTrajViewer
{
    [StructLayout(LayoutKind.Sequential)]
    public struct SprayTrajParams
    {
        public double scanRowSpacing;   // 行距 mm
        public double pointSpacing;     // 点距 mm
        public double standoff;         // 离面距离 mm
        public double speed;            // 速度 mm/s
        public double turnRadiusFactor; // 折返半径系数
        public double flatnessTol;      // 平面内点阈值 mm
        public double curvatureTol;     // 曲面升级阈值 mm
        public int    sorK;             // SOR 近邻数
        public double sorZ;             // SOR Z 倍数
        public int    normalFlip;       // 法向翻转：0=自动，1=强制翻转
        public double fillHoleSize;     // 孔洞填补阈值 mm（间隙<此值覆盖重建，行连续跨越）
        public int    scanDirMode;      // 轨迹方向：0=自动，1=沿长边喷，2=沿短边喷
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct SprayTrajPoint
    {
        public double x, y, z;
        public double i, j, k;
        public double speed;
        public int spray;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct SprayTrajPlane
    {
        public double nx, ny, nz;   // 喷距法向（轨迹按它抬升）
        public double cx, cy, cz;
        public int count;
        // 面级法向诊断：与 C++ SprayTrajPlane 逐字段对应，顺序不能改
        public double n0x, n0y, n0z;   // 射线自由空间判出的真实法向（可喷侧）
        public double angZ;            // 与 +Z 夹角 0~180
        public double tiltDeg;         // 仰角 90-angZ，-90~90
        public double freePlus;        // +n 侧自由距离 mm（-1 = 量程内无命中）
        public double freeMinus;       // -n 侧自由距离 mm
        public int side;               // +1 正向 / -1 反向 / 0 竖直
        public int isInner;            // 1 = 内壁
        public int sideSrc;            // 1=射线 2=兜底 3=强制
        public double ux, uy, uz;
        public double vx, vy, vz;
        public double minU, maxU;
        public double minV, maxV;
        // 规则化重建面片：轮廓矢量化 + 矩形吸附 + 三角化后的网格
        public int polyVertCount;                       // 顶点数
        public int polyTriCount;                        // 三角形数
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 1152)]   // 384 顶点 × 3
        public double[] polyVerts;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 2304)]   // 768 三角形 × 3
        public int[] polyTris;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct SprayTrajResult
    {
        public IntPtr points;   // SprayTrajPoint*
        public int count;
        public double plane_rms;
        public double normal_x, normal_y, normal_z;
        public double plane_cx, plane_cy, plane_cz;
        public double seg_len_mean;
        public double seg_len_max;
        public int rows;
        public int input_points;
        public int removed_points;
        public int inlier_points;
        public int plane_count;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 16)]
        public SprayTrajPlane[] planes;
        public IntPtr error;    // const char*
    }

    public static class Native
    {
        const string Dll = "SprayTraj.dll";

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)]
        public static extern SprayTrajParams spray_traj_default_params();

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)]
        public static extern SprayTrajResult spray_traj_generate(
            [In] double[] xyz, int n, ref SprayTrajParams p);

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)]
        public static extern void spray_traj_free(ref SprayTrajResult res);
    }
}
