@echo off
rem build_simc.bat <solution.cpp relative to repo> <exe name> [cl flags]  -> build\bin\<name>.exe (feni's simc harness, MSVC)
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
cd /d %~dp0..
if not exist build\bin mkdir build\bin
if not exist build\wrap mkdir build\wrap
set SRC=%1
set SRC=%SRC:\=/%
echo #define SFSRC "%CD:\=/%/%SRC%"> build\wrap\%2.cpp
echo #include "%CD:\=/%/scripts/simc.cpp">> build\wrap\%2.cpp
cl /nologo /O2 /EHsc /std:c++17 /Ibuild\compat %3 %4 %5 /Fo:build\bin\%2_ /Fe:build\bin\%2.exe build\wrap\%2.cpp
