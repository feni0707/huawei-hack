@echo off
rem Build any cpp: build_any.bat <src.cpp> <name>  -> build\bin\<name>.exe (MSVC /O2)
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
cd /d %~dp0..
cl /nologo /O2 /EHsc /std:c++17 /Fo:build\bin\%2_ /Fe:build\bin\%2.exe %1
