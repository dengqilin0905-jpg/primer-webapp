@echo off
setlocal enabledelayedexpansion

set "ROOT=%~dp0"
set "PY=%ROOT%python\python.exe"
set "APP=%ROOT%app"
set "DEPS=%APP%\_deps"
set "PORT=18501"

echo === Checking files ===
if not exist "%PY%" (echo [ERROR] python.exe not found: %PY% & pause & exit /b 1)
if not exist "%APP%\app.py" (echo [ERROR] app.py not found: %APP%\app.py & pause & exit /b 1)

echo === Fixing python._pth (embeddable isolation) ===
for %%f in ("%ROOT%python\python*._pth") do set "PTH=%%f"
if not defined PTH (echo [ERROR] python._pth not found in %ROOT%python & pause & exit /b 1)

echo Using: %PTH%
findstr /C:"..\app" "%PTH%" >nul || (echo ..\app>>"%PTH%")
findstr /C:"..\app\_deps" "%PTH%" >nul || (echo ..\app\_deps>>"%PTH%")
findstr /C:"import site" "%PTH%" >nul || (echo import site>>"%PTH%")

echo === Show python sys.path and test imports ===
"%PY%" -c "import sys; print('PY OK'); print('\n'.join(sys.path))"
echo.
"%PY%" -c "import streamlit; print('streamlit OK', streamlit.__version__)"
echo.
"%PY%" -c "import pandas as pd; print('pandas OK', pd.__version__)"
echo.

echo === Starting Streamlit on http://127.0.0.1:%PORT% ===
"%PY%" -m streamlit run "%APP%\app.py" --server.address 127.0.0.1 --server.port %PORT% --server.headless true --global.developmentMode=false

echo === Streamlit exited (errors should be above) ===
pause
