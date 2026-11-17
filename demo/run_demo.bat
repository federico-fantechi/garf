@echo off
REM GARF (binary) - DEMO on synthetic data (Windows).
REM Generates demo\input\panel.csv, then runs the engine with demo\config.yaml.
REM Use this to confirm everything works before touching your own data.
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
  echo [demo] first run: creating shared virtual environment .venv ...
  %PY% -m venv "%ROOT%\.venv"
)
call "%ROOT%\.venv\Scripts\activate.bat"
python -m pip install --upgrade pip >nul 2>&1
pip install -q -r "%ROOT%\requirements.txt"
if errorlevel 1 ( echo [demo] ERROR: dependency install failed. & pause & exit /b 1 )

if not exist input mkdir input
if not exist output mkdir output
echo [demo] generating synthetic input\panel.csv ...
python make_example.py
echo [demo] running GARF (binary) ...
python "%ROOT%\engine\garf_binary.py" --config config.yaml
if errorlevel 1 ( echo [demo] ERROR: the run failed - see the message above. & pause & exit /b 1 )
echo.
echo [demo] Done. Look inside demo\output\.
pause
