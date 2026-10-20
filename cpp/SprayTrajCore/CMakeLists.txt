cmake_minimum_required(VERSION 3.16)
project(SprayTraj LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 核心库（动态库，导出 C ABI）
# robot_safety.cpp：机器人可执行性校验（倾角规划 + 枪体/连杆碰撞 + IK 可达）。
# 它与 spray_traj.cpp 各自独立：spray_traj 保证「喷得对不对」，robot_safety 保证
# 「机器人能不能安全执行」。两者编进同一个 DLL，下游按需调用，互不污染。
add_library(SprayTraj SHARED spray_traj.cpp robot_safety.cpp)
target_include_directories(SprayTraj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(SprayTraj PRIVATE SPRAYTRAJ_EXPORTS)
if(MSVC)
    target_compile_options(SprayTraj PRIVATE /utf-8)
endif()
set_target_properties(SprayTraj PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# 命令行驱动（测试用）
add_executable(spray_driver driver.cpp)
target_link_libraries(spray_driver PRIVATE SprayTraj)
if(MSVC)
    target_compile_options(spray_driver PRIVATE /utf-8)
endif()
set_target_properties(spray_driver PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# SprayTrajResult 结构体本身约 145KB（8 个面片各含 384 顶点 + 768 三角的定长数组），
# 驱动把结果放在栈上、DLL 内部也返回一份，默认 1MB 栈在某些机器上会撞到 0xC00000FD
# （STATUS_STACK_OVERFLOW：进程无任何输出直接退出）。这里给测试驱动留 16MB。
if(MSVC)
    target_link_options(spray_driver PRIVATE /STACK:16777216)
endif()

# 安全层命令行驱动（读轨迹 CSV + 点云 PLY -> 校验并输出修正后轨迹）
# 单独一个 exe：spray_driver 的输出决定 9 个非曲面样本的 CSV 基线，安全层不能碰那条链。
add_executable(spray_safety safety_driver.cpp)
target_link_libraries(spray_safety PRIVATE SprayTraj)
if(MSVC)
    target_compile_options(spray_safety PRIVATE /utf-8)
endif()
set_target_properties(spray_safety PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
