#!/usr/bin/env bash
# One-command reproducibility: setup + full pipeline + verification.
# Usage:  bash run_all.sh
# Runtime: ~20 minutes (the three simulation benchmarks dominate).
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -d .venv ]; then
  echo "Creating virtual environment..."
  python3 -m venv .venv
fi
.venv/bin/pip install -q -r requirements.txt

echo "=== 01 real-data baseline ==="
.venv/bin/python analysis/01_real_data_baseline.py

echo "=== 02 simulation study (~8 min) ==="
.venv/bin/python analysis/02_simulation_study.py

echo "=== 03 real-data validation ==="
.venv/bin/python analysis/03_real_data_validation.py

echo "=== 04 structured GxE robustness (~7 min) ==="
.venv/bin/python analysis/04_robustness_structured.py

echo "=== 05 unbalanced-data robustness (~4 min) ==="
.venv/bin/python analysis/05_robustness_unbalanced.py

echo "=== 07 factor-analytic rank analysis (~1.5 min) ==="
.venv/bin/python analysis/07_robustness_factor_analytic.py

echo "=== 06 summary figure ==="
.venv/bin/python analysis/06_summary_figure.py

echo "=== publication figures (300 DPI) ==="
.venv/bin/python analysis/figures_publication.py

echo "=== editable Word manuscript ==="
.venv/bin/python analysis/make_docx.py

echo "=== verification: manuscript tables vs outputs ==="
.venv/bin/python analysis/check_all_tables.py

echo "ALL DONE"
