::======================================================
:: Build the out.exe for the release version
::======================================================
@echo off

:: Create the target directory if it doesn't exist
if not exist .\bin\release (
    mkdir .\bin\release
)

if exist .\bin\release\out.exe (
    del .\bin\release\out.exe
)


g++ ^
-O2 ^
-I /cygwin64/usr/include ^
-L /cygwin64/usr/lib ^
.\source\test.cpp ^
-o .\bin\release\out.exe ^
-lgsl -lgslcblas


echo Executing the program: .\bin\release\out.exe
.\bin\release\out.exe

