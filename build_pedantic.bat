::======================================================
:: Build the out.exe for the release version
::======================================================
@echo off

:: Prompt the user for the source file name
set /p sourcefile=Enter the name of the source file (e.g., test.c): 
@REM set sourcefile=ex-01-18b.c

:: Create the target directory if it doesn't exist
if not exist .\bin\release (
    mkdir .\bin\release
)

:: Delete the existing out.exe if it exists
if exist .\bin\release\out.exe (
    del .\bin\release\out.exe
)

@echo on

:: Compile the source file to create out.exe
gcc ^
-Wall ^
-Wextra ^
-Werror ^
-pedantic ^
-std=c18 ^
-O2 ^
-D_GNU_SOURCE ^
-I .\include ^
-I /cygwin64/usr/include ^
-L /cygwin64/usr/lib ^
.\source\%sourcefile% ^
-o .\bin\release\out.exe ^
-lgsl -lgslcblas -lm

@echo off

:: Execute the program
.\bin\release\out.exe
@REM echo Executing the program: .\bin\release\out.exe
@REM .\bin\release\out.exe