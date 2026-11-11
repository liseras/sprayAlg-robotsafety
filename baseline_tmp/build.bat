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
set SRC=%~dp0.
set OUT=%~dp0bin
if not exist "%OUT%" mkdir "%OUT%"
set CFLAGS=/nologo /std:c++17 /EHsc /O2 /MD /utf-8 /DNDEBUG /DSPRAYTRAJ_EXPORTS /I"%SRC%"
"%BIN%\cl.exe" %CFLAGS% /Fe:"%OUT%\spray_driver.exe" /Fo:"%OUT%\\" "%SRC%\spray_traj.cpp" "%SRC%\robot_safety.cpp" "%SRC%\driver.cpp" /link /SUBSYSTEM:CONSOLE /STACK:16777216
if errorlevel 1 exit /b 1
echo BASELINE BUILD OK
