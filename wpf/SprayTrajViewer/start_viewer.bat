@echo off
rem ---------------------------------------------------------------------------
rem Launch SprayTrajViewer without going through Explorer.
rem
rem Why: SmartScreen's app-reputation check only hooks launches started from
rem Explorer. A freshly rebuilt exe has a brand-new hash and no signature, so
rem Windows treats it as "unknown publisher" and blocks it every single build.
rem Starting it from a console skips that hook entirely.
rem
rem For CLI runs (--ply / --shot / --log ...) call the exe directly instead:
rem   bin\Release\net9.0-windows\SprayTrajViewer.exe --ply x.ply --shot y.png
rem ---------------------------------------------------------------------------
setlocal
set "ROOT=%~dp0"
set "EXE=%ROOT%bin\Release\net9.0-windows\SprayTrajViewer.exe"
if not exist "%EXE%" set "EXE=%ROOT%bin\Debug\net9.0-windows\SprayTrajViewer.exe"
if not exist "%EXE%" (
  echo [ERROR] SprayTrajViewer.exe not found.
  echo         Build it first:  dotnet build -c Release
  pause
  exit /b 1
)
start "" "%EXE%" %*
