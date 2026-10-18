@echo off
rem Build scripts/sim.cpp into buildbinsim.exe (MSVC, Windows)
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
cd /d %~dp0..
if not exist build\bin mkdir build\bin
cl /nologo /O2 /EHsc /std:c++17 /Fo:build\bin\ /Fe:build\bin\sim.exe scripts\sim.cpp
