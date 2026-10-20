@echo off
rem Build any cpp with POSIX shim dir build\compat on include path: build_any2.bat <src.cpp> <name> [extra cl flags]
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
cd /d %~dp0..
if not exist build\bin mkdir build\bin
cl /nologo /O2 /EHsc /std:c++17 /Ibuild\compat %3 %4 %5 /Fo:build\bin\%2_ /Fe:build\bin\%2.exe %1
