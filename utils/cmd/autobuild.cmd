@ECHO OFF
SETLOCAL
git rev-parse --show-toplevel >nul 2>&1 || EXIT /B %ERRORLEVEL%
FOR /F "delims=" %%i IN ('git rev-parse --show-toplevel') DO SET ROOT=%%i
sphinx-autobuild --nitpicky --ignore "%ROOT%\source\tags" -- "%ROOT%\source" "%ROOT%\build"
