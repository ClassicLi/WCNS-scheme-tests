@echo off
setlocal

rem Build the complex-wave test with Intel ifx on Windows.
rem Usage: build.bat [dp^|qp]   (default: dp)
rem The oneAPI path below is a starting point; adjust it if needed.

set PRECISION=%1
if "%PRECISION%"=="" set PRECISION=dp

call "D:\Program Files\Intel\oneAPI\compiler\latest\env\vars.bat" intel64 >nul

if not exist build\%PRECISION% mkdir build\%PRECISION%
if not exist bin\%PRECISION% mkdir bin\%PRECISION%

set PPFLAGS=/fpp
if /I "%PRECISION%"=="dp" set PPFLAGS=/fpp /DPRECISION_DP

ifx /O2 %PPFLAGS% /module:build\%PRECISION% ^
  ..\interpolation\src\interpolation\precision_select.f90 ^
  ..\interpolation\src\interpolation\scheme_diagnostics.f90 ^
  ..\interpolation\src\interpolation\scheme_params.f90 ^
  ..\interpolation\src\interpolation\interp_linear.f90 ^
  ..\interpolation\src\interpolation\interp_deriv.f90 ^
  ..\interpolation\src\interpolation\interp_qr.f90 ^
  ..\interpolation\src\interpolation\interp_rbf.f90 ^
  ..\interpolation\src\interpolation\interp_traditional.f90 ^
  ..\interpolation\src\interpolation\interp_framework.f90 ^
  ..\interpolation\src\interpolation\scheme_dispatch.f90 ^
  src\common\cw_control.f90 ^
  src\common\cw_periodic.f90 ^
  src\common\cw_initial.f90 ^
  src\common\cw_time.f90 ^
  src\common\cw_io.f90 ^
  src\wcns\cw_params.f90 ^
  src\wcns\difference_operators.f90 ^
  src\wcns\cw_wcns_rhs.f90 ^
  app\complex_wave.f90 ^
  -o bin\%PRECISION%\complex_wave.exe
