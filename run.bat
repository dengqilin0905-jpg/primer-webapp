@echo off
setlocal
if exist "%~dp0.venv\Scripts\primer-webapp.exe" (
  "%~dp0.venv\Scripts\primer-webapp.exe" --server.address 127.0.0.1 --server.port 8501
) else (
  primer-webapp --server.address 127.0.0.1 --server.port 8501
)
