@echo off
set VSROOT=D:\software\VisualStudio\Community
set MSVCVER=14.51.36231
set SDKROOT=D:\Windows Kits\10
set SDKVER=10.0.26100.0
set MSVC=%VSROOT%\VC\Tools\MSVC\%MSVCVER%
set BIN=%MSVC%\bin\Hostx64\x64
set SDKBIN=%SDKROOT%\bin\%SDKVER%\x64
set INCLUDE=%MSVC%\include;%SDKROOT%\Include\%SDKVER%\ucrt;%SDKROOT%\Include\%SDKVER%\um;%SDKROOT%\Include\%SDKVER%\shared;%SDKROOT%\Include\%SDKVER%\winrt
set LIB=%MSVC%\lib\x64;%SDKROOT%\Lib\%SDKVER%\ucrt\x64;%SDKROOT%\Lib\%SDKVER%\um\x64
set PATH=%BIN%;%SDKBIN%;%PATH%
cd /d "%~dp0"
cl /nologo /std:c++17 /EHsc /O2 /MD /utf-8 /DNDEBUG /DSPRAYTRAJ_EXPORTS /I. /c spray_traj.cpp /Fobuild_rls\obj\spray_traj.obj 2>&1
echo CL_EXIT=%ERRORLEVEL%
