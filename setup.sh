#!/usr/bin/env bash
# GARF (binary) — one-time setup for macOS / Linux.
# Creates a shared virtual environment (.venv) in this folder and installs the
# dependencies. Run this once; the run scripts in engine/ and demo/ reuse it.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

if command -v python3 >/dev/null 2>&1; then PY=python3
elif command -v python  >/dev/null 2>&1; then PY=python
else
  echo "ERROR: Python was not found. Install Python 3.10+ (see INSTALL.txt) and retry."
  exit 1
fi

echo "[setup] using $($PY --version)"
if [ ! -d ".venv" ]; then
  echo "[setup] creating virtual environment (.venv) ..."
  "$PY" -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install --upgrade pip >/dev/null
echo "[setup] installing dependencies ..."
pip install -r requirements.txt
python -c "import pandas, numpy, scipy, sklearn, yaml, tqdm; print('[setup] all dependencies OK')"
echo "[setup] Done. Next:  cd demo && ./run_demo.sh   (test)   or   cd engine && ./run.sh   (your data)"
