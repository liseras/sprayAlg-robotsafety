@echo off
rem build spray_traj.exe with MSVC (VS2022 Community)
call "D:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
if errorlevel 1 ( echo VCVARS_FAIL & exit /b 1 )
cd /d "%~dp0"
if not exist bin mkdir bin
cl /nologo /O2 /std:c++17 /EHsc /utf-8 /W3 /Fe:bin\spray_traj.exe /Fo:bin\ src\main.cpp
if errorlevel 1 ( echo BUILD_FAIL & exit /b 2 )
echo BUILD_OK
