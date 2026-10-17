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
