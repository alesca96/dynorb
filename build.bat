::======================================================
:: Build the out.exe for the release version
::======================================================
@echo off

:: Prompt the user for the source file name
set /p sourcefile=Enter the name of the source file (e.g., test.c): 

:: Create the target directory if it doesn't exist
if not exist .\bin\release (
    mkdir .\bin\release
)

:: Delete the existing out.exe if it exists
if exist .\bin\release\out.exe (
    del .\bin\release\out.exe
)

:: Compile the source file to create out.exe
gcc ^
-O2 ^
-I .\include ^
-I /cygwin64/usr/include ^
-L /cygwin64/usr/lib ^
.\source\%sourcefile% ^
-o .\bin\release\out.exe ^
-lgsl -lgslcblas -lm

:: Execute the program
echo Executing the program: .\bin\release\out.exe
.\bin\release\out.exe


