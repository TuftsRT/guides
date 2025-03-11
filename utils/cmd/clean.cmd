@ECHO OFF
SETLOCAL
git rev-parse --show-toplevel 1>nul || EXIT /B 1
FOR /F "delims=" %%i IN ('git rev-parse --show-toplevel') DO SET ROOT=%%i
IF EXIST "%ROOT%\build" RMDIR /S /Q "%ROOT%\build"
IF EXIST "%ROOT%\jupyter_execute" RMDIR /S /Q "%ROOT%\jupyter_execute"
IF EXIST "%ROOT%\source\tags" RMDIR /S /Q "%ROOT%\source\tags"
