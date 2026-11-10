$ErrorActionPreference = 'Stop'
$msvc = 'D:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207'
$sdk  = 'C:\Program Files (x86)\Windows Kits\10'
$sdkv = '10.0.26100.0'
$env:PATH = "$msvc\bin\Hostx64\x64;$env:PATH"
$env:INCLUDE = "$msvc\include;$sdk\Include\$sdkv\ucrt;$sdk\Include\$sdkv\um;$sdk\Include\$sdkv\shared"
$env:LIB = "$msvc\lib\x64;$sdk\Lib\$sdkv\ucrt\x64;$sdk\Lib\$sdkv\um\x64"
Set-Location $PSScriptRoot
New-Item -ItemType Directory -Force build_obj | Out-Null

& cl.exe /nologo /c /O2 /std:c++17 /utf-8 /MD /EHsc /DSPRAYTRAJ_EXPORTS /D_CRT_SECURE_NO_WARNINGS /I. /Fo:build_obj\ /Fd:build_obj\ spray_traj.cpp robot_safety.cpp
if ($LASTEXITCODE -ne 0) { Write-Output "CL_CORE_FAIL=$LASTEXITCODE"; exit 1 }

& link.exe /nologo /DLL /OUT:build_obj\SprayTraj.dll build_obj\spray_traj.obj build_obj\robot_safety.obj
if ($LASTEXITCODE -ne 0) { Write-Output "LINK_DLL_FAIL=$LASTEXITCODE"; exit 2 }

& cl.exe /nologo /c /O2 /std:c++17 /utf-8 /MD /EHsc /D_CRT_SECURE_NO_WARNINGS /I. /Fo:build_obj\ /Fd:build_obj\ driver.cpp
if ($LASTEXITCODE -ne 0) { Write-Output "CL_DRIVER_FAIL=$LASTEXITCODE"; exit 3 }

& link.exe /nologo /OUT:build_obj\spray_driver.exe /STACK:16777216 build_obj\driver.obj build_obj\SprayTraj.lib
if ($LASTEXITCODE -ne 0) { Write-Output "LINK_EXE_FAIL=$LASTEXITCODE"; exit 4 }

Write-Output "BUILD_OK"
Get-ChildItem build_obj\SprayTraj.dll, build_obj\spray_driver.exe | Select-Object Name,Length,LastWriteTime | Format-Table -AutoSize
