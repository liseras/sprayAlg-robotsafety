@echo off
REM ============================================================
REM Incremental build: only recompiles changed translation units.
REM Same MSVC environment as build_release.bat (hand-built, no vcvarsall).
REM
REM   build_inc.bat        -> incremental (typical 3~40s)
REM   build_release.bat    -> forced full rebuild (~80s), known-good fallback
REM
REM Keep this file pure ASCII (cmd reads .bat as OEM codepage / GBK).
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

if not exist "%~dp0build_rls" mkdir "%~dp0build_rls"

python "%~dp0mkbuild.py"
if errorlevel 1 goto fail

dir /b "%~dp0build_rls\*.exe" "%~dp0build_rls\*.dll"
exit /b 0

:fail
echo.
echo BUILD FAILED
exit /b 1
