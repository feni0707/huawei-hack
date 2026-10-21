@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
cd /d %~dp0..\..
cl /nologo /O2 /EHsc /std:c++17 /Ibuild\compat /Fo:build\bin\dsearch2_ /Fe:build\bin\dsearch2.exe build\wrap\dsearch2.cpp
