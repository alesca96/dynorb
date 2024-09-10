::======================================================
:: Build the test.exe for the test version
::======================================================
@echo off

:: Prompt the user for the source file name
set sourcefile=test.c 

:: Create the target directory if it doesn't exist
if not exist .\bin\test (
    mkdir .\bin\test
)

:: Delete the existing test.exe if it exists
if exist .\bin\test\test.exe (
    del .\bin\test\test.exe
)

@echo on
:: Compile the source file to create test.exe
gcc ^
-Wall ^
-Wextra ^
-Werror ^
-pedantic ^
-std=c18 ^
-O0 ^
-g ^
-I .\include ^
-I /cygwin64/usr/include ^
-L /cygwin64/usr/lib ^
.\source\%sourcefile% ^
-o .\bin\test\test.exe ^
-lgsl -lgslcblas -lm

:: Execute the program
.\bin\test\test.exe


