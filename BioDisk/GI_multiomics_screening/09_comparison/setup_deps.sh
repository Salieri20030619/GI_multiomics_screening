#!/usr/bin/env bash
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/GI_multiomics_screening
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
mkdir -p "$ROOT/09_comparison/vendor"
"$PY" -m pip install --no-deps --target "$ROOT/09_comparison/vendor" pyxlsb==1.0.10 xlrd==2.0.2 > "$ROOT/09_comparison/dependency_install.log" 2>&1
cp /mnt/BioDisk/wangmingyuan/preflight.py /mnt/BioDisk/wangmingyuan/env_check.R "$ROOT/09_comparison/"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla "$ROOT/09_comparison/env_check.R" > "$ROOT/09_comparison/R_environment.log" 2>&1
echo PROJECT_DEPENDENCIES_COMPLETE
