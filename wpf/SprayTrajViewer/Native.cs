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
        public double nx, ny, nz;
        public double cx, cy, cz;
        public int count;
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
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 8)]
        public SprayTrajPlane[] planes;
        public IntPtr error;    // const char*
    }

    // ---------------------------------------------------------------- 机器人安全层
    //
    // robot_safety.h 的 C ABI 镜像。职责边界：spray_traj 保证"喷得对不对"，
    // robot_safety 保证"机器人能不能安全执行"（可达 / 关节限位 / 枪体与连杆碰撞 / 自碰撞 / 喷角）。
    // 两者编在同一个 DLL 里，互不依赖：不调用 robot_safety_check 时行为与以前完全一致。

    [StructLayout(LayoutKind.Sequential)]
    public struct RobotSafetyParams
    {
        public double standoff;         // 喷距 mm（必须与生成轨迹时一致）
        public double safety;           // 安全间隙 mm
        public double search;           // 判定带 mm（超过只报下界）
        public double toolLen;          // 法兰 → 喷嘴 mm
        public double tiltNonUpDeg;     // α=90°（竖直侧壁）时的倾角；实际 = 它 × α/90°（α 为法向与 +Z 夹角）
        public double upThresholdDeg;   // 已不参与倾角决策（倾角按法向线性给），保留仅为 ABI 稳定
        public double tiltMaxDeg;       // 倾角上限
        public double tiltStepDeg;      // 倾角搜索步长
        public double tiltRateDeg;      // 每点倾角最大变化
        public double maxSprayAngleDeg; // 工艺喷角上限
        public double jointMarginMin;   // 关节限位余量下限
        public double condMax;          // 归一化雅可比条件数上限
        public double transitLiftMax;   // 过渡点最大抬升量
        public int    linkStride;       // 连杆碰撞抽样步长
        public int    includeBase;      // 1=查底座段
        public int    fixTransit;       // 1=允许过渡点抬升
        public int    azFallback;       // 1=允许换倾斜方位兜底
        public int    maxCloudPoints;   // 保留（<=0 = 用全量点云）
    }

    // 违规位（与 RS_VIOL_* 一致）
    public static class RobotSafetyRule
    {
        public const uint IkFail = 1u << 0;
        public const uint Collision = 1u << 1;
        public const uint Clearance = 1u << 2;
        public const uint SelfCollision = 1u << 3;
        public const uint JointLimit = 1u << 4;
        public const uint SprayAngle = 1u << 5;
        public const uint NearSingular = 1u << 6;
        public const uint TiltJump = 1u << 7;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct RobotSafetyPoint
    {
        public double x, y, z;       // 修正后 TCP
        public double i, j, k;       // 修正后枪轴朝向（指向工件）
        public double speed;
        public double tiltDeg;       // 实际喷角（相对**原始**面法向）
        public double clearanceMm;   // 最小表面距离（负 = 穿透）
        public double selfGapMm;     // 自碰撞最小间隙（负 = 干涉）
        public double marginDeg;     // 关节限位最小余量（负 = 超限）
        public double standoffErrMm; // 喷距守恒偏差
        public double liftedMm;      // 过渡点抬升量
        public double cond;          // 归一化雅可比条件数
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 6)]
        public double[] theta;       // 关节角（度）
        public int    spray;
        public int    ikOk;
        public int    reverted;
        public uint   violation;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct RobotSafetyResult
    {
        public IntPtr points;        // RobotSafetyPoint*
        public int    count;
        public int    reachable;
        public int    tiltedPoints;
        public int    poseReverted;
        public int    violations;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 8)]
        public int[]  byRule;
        public double minClearanceMm;
        public double tiltMeanDeg;
        public double tiltP95Deg;
        public double tiltMaxDeg;
        public double standoffMaxErrMm;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 16)]
        public double[] basePose;    // T_base_in_cloud（行主序）
        public int    baseAuto;
        public int    cloudPoints;
        public double cloudSpacingMm;
        public int    l0Issues;
        public double standoffMedianMm;
        public IntPtr error;         // const char*
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

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)]
        public static extern RobotSafetyParams robot_safety_default_params();

        /// <summary>对一条轨迹做可执行性校验 + 倾角修正。</summary>
        /// <param name="pts">轨迹点（spray_traj_generate 的输出）</param>
        /// <param name="cloudXyz">原始点云 xyz，必须与轨迹同坐标系</param>
        /// <param name="basePose16">T_base_in_cloud 4x4 行主序；null = 按工件质心自动假设</param>
        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)]
        public static extern RobotSafetyResult robot_safety_check(
            [In] SprayTrajPoint[] pts, int n,
            [In] double[] cloudXyz, int cloudCount,
            [In] double[]? basePose16,
            ref RobotSafetyParams p);

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)]
        public static extern void robot_safety_free(ref RobotSafetyResult res);

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr robot_safety_rule_name(uint bit);
    }
}
