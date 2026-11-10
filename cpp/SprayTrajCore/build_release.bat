@echo off
REM ============================================================
REM Rebuild C++ core with a hand-built MSVC environment.
REM
REM Do NOT use vcvarsall.bat here: it calls reg.exe to locate the VS
REM instance, and reg.exe is blacklisted by the local sandbox, so the
REM env is half-initialised (rc.exe missing) and CMake's compiler
REM self-test always fails. Setting INCLUDE/LIB/PATH directly is
REM equivalent and never touches the registry.
REM
REM NOTE: keep this file pure ASCII. cmd.exe reads .bat as the OEM
REM codepage (GBK here); UTF-8 Chinese comments corrupt the commands.
REM
REM Values measured on this machine:
REM   MSVC 14.51.36231 (ships with VS 2026, v180 toolset)
REM   Windows SDK 10.0.26100.0
REM ============================================================

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
set OUT=%~dp0build_rls
if not exist "%OUT%" mkdir "%OUT%"

REM /utf-8 is mandatory: the sources carry Chinese comments and Chinese
REM string literals. Without it MSVC decodes them as GBK and emits C4819.
set CFLAGS=/nologo /std:c++17 /EHsc /O2 /MD /utf-8 /DNDEBUG /DSPRAYTRAJ_EXPORTS /I"%SRC%"

echo [1/3] SprayTraj.dll
"%BIN%\cl.exe" %CFLAGS% /LD /Fe:"%OUT%\SprayTraj.dll" /Fo:"%OUT%\\" "%SRC%\spray_traj.cpp" "%SRC%\robot_safety.cpp" /link /IMPLIB:"%OUT%\SprayTraj.lib"
if errorlevel 1 goto fail

echo [2/3] spray_safety.exe
"%BIN%\cl.exe" %CFLAGS% /Fe:"%OUT%\spray_safety.exe" /Fo:"%OUT%\\" "%SRC%\safety_driver.cpp" /link /SUBSYSTEM:CONSOLE /LIBPATH:"%OUT%" SprayTraj.lib /STACK:16777216
if errorlevel 1 goto fail

echo [3/3] spray_driver.exe
"%BIN%\cl.exe" %CFLAGS% /Fe:"%OUT%\spray_driver.exe" /Fo:"%OUT%\\" "%SRC%\driver.cpp" /link /SUBSYSTEM:CONSOLE /LIBPATH:"%OUT%" SprayTraj.lib /STACK:16777216
if errorlevel 1 goto fail

echo.
echo BUILD OK: %OUT%
dir /b "%OUT%\*.exe" "%OUT%\*.dll"
exit /b 0

:fail
echo.
echo BUILD FAILED errorlevel=%errorlevel%
exit /b 1
