@ECHO OFF
SETLOCAL
git rev-parse --show-toplevel 1>nul || EXIT /B 1
FOR /F "delims=" %%i IN ('git rev-parse --show-toplevel') DO SET ROOT=%%i
sphinx-build --nitpicky "%ROOT%\source" "%ROOT%\build"
