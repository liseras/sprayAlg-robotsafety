@echo off
setlocal
call "D:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d "%~dp0"
if not exist build_obj mkdir build_obj
cl /nologo /c /O2 /std:c++17 /utf-8 /MD /EHsc /DSPRAYTRAJ_EXPORTS /D_CRT_SECURE_NO_WARNINGS /I. /Fo:build_obj\ /Fd:build_obj\ spray_traj.cpp robot_safety.cpp
if errorlevel 1 exit /b 1
link /nologo /DLL /OUT:build_obj\SprayTraj.dll build_obj\spray_traj.obj build_obj\robot_safety.obj
if errorlevel 1 exit /b 2
cl /nologo /c /O2 /std:c++17 /utf-8 /MD /EHsc /D_CRT_SECURE_NO_WARNINGS /I. /Fo:build_obj\ /Fd:build_obj\ driver.cpp
if errorlevel 1 exit /b 3
link /nologo /OUT:build_obj\spray_driver.exe /STACK:16777216 build_obj\driver.obj build_obj\SprayTraj.lib
if errorlevel 1 exit /b 4
echo BUILD_OK
