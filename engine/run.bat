@echo off
REM GARF (binary) - run on YOUR data (Windows).
REM 1) put your CSV in engine\input\  (default name: panel.csv)
REM 2) edit engine\config.yaml (column names, T0, k, tuning on/off)
REM 3) double-click this.  Results go to engine\output\.
cd /d "%~dp0"
set "ROOT=%~dp0.."

set "PY="
where py >nul 2>&1 && set "PY=py"
if not defined PY where python >nul 2>&1 && set "PY=python"
if not defined PY (
  echo ERROR: Python was not found. Install Python 3.10+ (see INSTALL.txt) and retry.
  pause
  exit /b 1
)

if not exist "%ROOT%\.venv" (
  echo [run] first run: creating shared virtual environment .venv ...
  %PY% -m venv "%ROOT%\.venv"
)
call "%ROOT%\.venv\Scripts\activate.bat"
python -m pip install --upgrade pip >nul 2>&1
pip install -q -r "%ROOT%\requirements.txt"
if errorlevel 1 ( echo [run] ERROR: dependency install failed. & pause & exit /b 1 )

if not exist input mkdir input
if not exist output mkdir output
python garf_binary.py --config config.yaml
if errorlevel 1 ( echo [run] ERROR: the run failed - see the message above. & pause & exit /b 1 )
echo.
echo [run] Done. Results are in engine\output\.
pause
