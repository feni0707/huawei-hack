@echo off
rem Build solutions\<name>.cpp into build\bin\<name>.exe (MSVC /O2)
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
cd /d %~dp0..
if not exist build\bin mkdir build\bin
cl /nologo /O2 /EHsc /std:c++17 /Fo:build\bin\ /Fe:build\bin\%1.exe solutions\%1.cpp
