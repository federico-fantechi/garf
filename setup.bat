@echo off
REM GARF (binary) - one-time setup for Windows.
REM Creates a shared virtual environment (.venv) here and installs dependencies.
REM Run once (double-click); the run scripts in engine\ and demo\ reuse it.
cd /d "%~dp0"

set "PY="
where py >nul 2>&1 && set "PY=py"
if not defined PY where python >nul 2>&1 && set "PY=python"
if not defined PY (
  echo ERROR: Python was not found. Install Python 3.10+ (see INSTALL.txt) and retry.
  pause
  exit /b 1
)

echo [setup] using %PY%
if not exist ".venv" (
  echo [setup] creating virtual environment .venv ...
  %PY% -m venv .venv
)
call ".venv\Scripts\activate.bat"
python -m pip install --upgrade pip >nul 2>&1
echo [setup] installing dependencies ...
pip install -r requirements.txt
if errorlevel 1 ( echo [setup] ERROR: dependency install failed. & pause & exit /b 1 )
python -c "import pandas, numpy, scipy, sklearn, yaml, tqdm; print('[setup] all dependencies OK')"
echo [setup] Done. Next: open demo\run_demo.bat (test) or engine\run.bat (your data).
pause
