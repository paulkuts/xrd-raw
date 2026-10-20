@echo off
rem XRD converter launcher (portable). Menu is drawn by Python.
chcp 65001 >nul
setlocal
set "PYTHONUTF8=1"
set "PYTHONIOENCODING="
set "HERE=%~dp0"
if not exist "%HERE%python\python.exe" (
  echo [ERROR] Folder "python" not found next to START.bat.
  echo Unpack the whole folder - it must contain python\, tools\, input\.
  pause
  exit /b 1
)
"%HERE%python\python.exe" "%HERE%tools\menu.py" %*
if errorlevel 1 pause
