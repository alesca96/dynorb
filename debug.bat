::===========================================================
::Builds the debug_out.exe to debug using Visual Studio Code
::===========================================================
@echo off
if exist .\bin\debug\debug_out.exe (
    del .\bin\debug\debug_out.exe )

:: Compile the source files
g++ ^
-g ^
-I /cygwin64/usr/include ^
-L /cygwin64/usr/lib ^
.\source\test1.cpp ^
-o .\bin\debug\debug_out.exe ^
-lgsl -lgslcblas

