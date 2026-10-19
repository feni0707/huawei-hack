@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
cd /d %~dp0..
cl /nologo /O2 /EHsc /std:c++17 /Fo:build\bin\ot10_ /Fe:build\bin\offline_tool10.exe scripts\offline_tool.cpp
