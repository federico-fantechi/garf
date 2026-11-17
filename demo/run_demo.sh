#!/usr/bin/env bash
# GARF (binary) — DEMO on synthetic data (macOS / Linux).
# Generates demo/input/panel.csv, then runs the engine with demo/config.yaml.
# Use this to confirm everything works before touching your own data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"

if command -v python3 >/dev/null 2>&1; then PY=python3
elif command -v python  >/dev/null 2>&1; then PY=python
else
  echo "ERROR: Python was not found. Install Python 3.10+ (see INSTALL.txt) and retry."
  exit 1
fi

if [ ! -d "$ROOT/.venv" ]; then
  echo "[demo] first run: creating shared virtual environment (.venv) ..."
  "$PY" -m venv "$ROOT/.venv"
fi
# shellcheck disable=SC1091
source "$ROOT/.venv/bin/activate"
python -m pip install --upgrade pip >/dev/null
pip install -q -r "$ROOT/requirements.txt"

cd "$HERE"
mkdir -p input output
echo "[demo] generating synthetic input/panel.csv ..."
python make_example.py
echo "[demo] running GARF (binary) ..."
python "$ROOT/engine/garf_binary.py" --config config.yaml
echo "[demo] Done. Look inside demo/output/."
