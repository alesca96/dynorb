::===========================================================
:: Builds the debug_out.exe to debug using Visual Studio Code
::===========================================================
@echo off

:: Create the target directory if it doesn't exist
if not exist .\bin\debug (
    mkdir .\bin\debug
)

:: Delete the old executable if it exists
if exist .\bin\debug\debug_out.exe (
    del .\bin\debug\debug_out.exe
)

:: Filename: 
set filename=test

:: Compile the source files with debug symbols
gcc ^
    -g ^
    -O0 ^
    -Wall ^
    -I "C:\cygwin64\usr\include" ^
    -L "C:\cygwin64\usr\lib" ^
    .\source\%filename%.c ^
    -o .\bin\debug\debug_out.exe ^
    -lgsl -lgslcblas

:: Check for successful compilation
if %errorlevel% neq 0 (
    echo "Build failed!"
    exit /b %errorlevel%
) else (
    echo "Build succeeded!"
)

@REM echo "Starting debugger..."
@REM gdb ^
@REM .\bin\debug\debug_out.exe


