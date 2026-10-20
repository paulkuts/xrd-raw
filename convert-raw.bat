@echo off
rem Direct conversion without menu: drag and drop a .raw file here.
chcp 65001 >nul
setlocal
set "PYTHONUTF8=1"
set "PYTHONIOENCODING="
if "%~1"=="" (
  echo Drag and drop a .raw file onto this bat, or run START.bat for the menu.
  pause
  exit /b 1
)
"%~dp0python\python.exe" "%~dp0tools\winxrd2xy.py" %*
pause
