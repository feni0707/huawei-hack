@echo off
rem Build scripts\offline_tool.cpp -> build\bin\offline_tool.exe (MSVC)
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
cd /d %~dp0..
if not exist build\bin mkdir build\bin
cl /nologo /O2 /EHsc /std:c++17 /Fo:build\bin\ot_ /Fe:build\bin\offline_tool.exe scripts\offline_tool.cpp
