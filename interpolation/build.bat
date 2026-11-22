@echo off
setlocal

rem Build with Intel ifx on Windows.
rem Usage: build.bat [qp^|dp]   (default: qp)

set PRECISION=%1
if "%PRECISION%"=="" set PRECISION=qp

call "D:\Program Files\Intel\oneAPI\compiler\latest\env\vars.bat" intel64 >nul

if not exist build\%PRECISION% mkdir build\%PRECISION%
if not exist bin\%PRECISION% mkdir bin\%PRECISION%

set PPFLAGS=/fpp
if /I "%PRECISION%"=="dp" set PPFLAGS=/fpp /DPRECISION_DP

ifx /O2 %PPFLAGS% /module:build\%PRECISION% ^
  src\interpolation\precision_select.f90 ^
  src\interpolation\scheme_diagnostics.f90 ^
  src\interpolation\scheme_params.f90 ^
  src\interpolation\interp_linear.f90 ^
  src\interpolation\interp_deriv.f90 ^
  src\interpolation\interp_qr.f90 ^
  src\interpolation\interp_rbf.f90 ^
  src\interpolation\interp_traditional.f90 ^
  src\interpolation\interp_framework.f90 ^
  src\interpolation\scheme_dispatch.f90 ^
  src\test\control.f90 ^
  src\test\standard_function.f90 ^
  src\test\diagnostics.f90 ^
  src\test\mesh.f90 ^
  src\test\solver.f90 ^
  src\test\io.f90 ^
  app\main.f90 ^
  -o bin\%PRECISION%\interpolation_test.exe
