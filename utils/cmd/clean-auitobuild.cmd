@ECHO OFF
CALL "%~dp0clean.cmd"
IF %ERRORLEVEL% NEQ 0 EXIT /B %ERRORLEVEL%
CALL "%~dp0autobuild.cmd"
