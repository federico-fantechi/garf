#!/usr/bin/env bash
# GARF (binary) — run on YOUR data (macOS / Linux).
# 1) put your CSV in engine/input/  (default name: panel.csv)
# 2) edit engine/config.yaml (column names, T0, k, tuning on/off)
# 3) run this.  Results go to engine/output/.
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
  echo "[run] first run: creating shared virtual environment (.venv) ..."
  "$PY" -m venv "$ROOT/.venv"
fi
# shellcheck disable=SC1091
source "$ROOT/.venv/bin/activate"
python -m pip install --upgrade pip >/dev/null
pip install -q -r "$ROOT/requirements.txt"

cd "$HERE"
mkdir -p input output
python garf_binary.py --config config.yaml
