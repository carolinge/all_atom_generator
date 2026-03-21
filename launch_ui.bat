@echo off
title all_atom MD Server
call C:\Users\Caroline\miniconda3\Scripts\activate.bat allatom
echo.
echo  all_atom MD Setup Server
echo  ========================
echo  Open browser: http://localhost:5000
echo  To stop: close this window
echo.
start "" "http://localhost:5000"
python "%~dp0app\server.py"
pause
