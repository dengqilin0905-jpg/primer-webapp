@echo off
setlocal enabledelayedexpansion

REM ====== Data (DB) location (survives app updates) ======
set "APPDATA_DIR=%APPDATA%\PrimerWebApp"
if not exist "%APPDATA_DIR%" mkdir "%APPDATA_DIR%"
set "DB_PATH=%APPDATA_DIR%\app.db"

REM ====== Paths ======
set "ROOT=%~dp0"
set "PY=%ROOT%python\python.exe"
set "APP=%ROOT%app"
set "DEPS=%APP%\_deps"

REM ====== Expose DB_PATH to the app ======
set "DB_PATH=%DB_PATH%"

REM ====== Python module search path ======
set "PYTHONPATH=%DEPS%;%APP%"

REM ====== Port ======
set "PORT=18501"

REM ====== Start Streamlit in a new window ======
start "PrimerWebApp" "%PY%" -m streamlit run "%APP%\app.py" --server.address 127.0.0.1 --server.port %PORT% --server.headless true --global.developmentMode=false


REM ====== Open browser ======
timeout /t 2 >nul
start "" "http://127.0.0.1:%PORT%"

echo PrimerWebApp is running at http://127.0.0.1:%PORT%
echo Close the "PrimerWebApp" window to stop it.
