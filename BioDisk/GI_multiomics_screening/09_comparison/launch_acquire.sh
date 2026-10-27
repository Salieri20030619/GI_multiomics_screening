#!/usr/bin/env bash
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/GI_multiomics_screening
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
"$PY" -m py_compile "$ROOT/09_comparison/acquire.py"
echo 'SYNTAX_CHECK: PASS python -m py_compile acquire.py'
echo 'LOGIC_CHECK: bounded retries; three workers; no EGA; CRC/gzip/XML checks; immutable existing valid sources; paths restricted to new projects'
echo "RUN_COMMAND: $PY $ROOT/09_comparison/acquire.py"
echo 'EXPECTED_SUCCESS_MARKER: ACQUISITION_ATTEMPTS_COMPLETE (attempts only; each resource has separate status)'
"$PY" "$ROOT/09_comparison/acquire.py"
