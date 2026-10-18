@echo off
rem Build scripts\sim.cpp into build\bin\<name>.exe (default name: sim) with MSVC on Windows
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
cd /d %~dp0..
if not exist build\bin mkdir build\bin
if "%1"=="" (set OUT=sim) else (set OUT=%1)
cl /nologo /O2 /EHsc /std:c++17 /Fo:build\bin\%OUT%.obj /Fe:build\bin\%OUT%.exe scripts\sim.cpp
