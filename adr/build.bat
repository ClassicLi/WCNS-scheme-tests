@echo off
setlocal

rem Build adr_weno and adr_wcns with Intel ifx on Windows.
rem Usage: build.bat [qp^|dp]   (default: qp)
rem The oneAPI path below is provided as a starting point; adjust if needed.

set PRECISION=%1
if "%PRECISION%"=="" set PRECISION=qp

call "D:\Program Files\Intel\oneAPI\compiler\latest\env\vars.bat" intel64 >nul

if not exist build\%PRECISION% mkdir build\%PRECISION%
if not exist bin\%PRECISION% mkdir bin\%PRECISION%

set PPFLAGS=/fpp
if /I "%PRECISION%"=="dp" set PPFLAGS=/fpp /DPRECISION_DP

ifx /O2 %PPFLAGS% /module:build\%PRECISION% ^
  ..\interpolation\src\interpolation\precision_select.f90 ^
  src\common\adr_control.f90 ^
  src\common\adr_periodic.f90 ^
  src\common\adr_dft.f90 ^
  src\common\adr_time.f90 ^
  src\common\adr_io.f90 ^
  src\weno\weno_params.f90 ^
  src\weno\weno_schemes.f90 ^
  src\weno\weno_dispatch.f90 ^
  app\adr_weno.f90 ^
  -o bin\%PRECISION%\adr_weno.exe

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
  src\common\adr_control.f90 ^
  src\common\adr_periodic.f90 ^
  src\common\adr_dft.f90 ^
  src\common\adr_time.f90 ^
  src\common\adr_io.f90 ^
  src\wcns\wcns_params.f90 ^
  src\wcns\difference_operators.f90 ^
  src\wcns\wcns_dispatch.f90 ^
  app\adr_wcns.f90 ^
  -o bin\%PRECISION%\adr_wcns.exe
