@echo off
setlocal

rem Windows launcher for the WSL checkout. Keep MCP stdio attached to stdout.
wsl.exe -d Ubuntu -- bash -lc "exec /home/wangchao/software/AutoZT/bin/autozt %*"
exit /b %ERRORLEVEL%
