cmake_minimum_required(VERSION 3.16)
project(SprayTraj LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 核心库（动态库，导出 C ABI）
add_library(SprayTraj SHARED spray_traj.cpp)
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
