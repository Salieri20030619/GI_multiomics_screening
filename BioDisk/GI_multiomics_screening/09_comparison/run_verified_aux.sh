#!/usr/bin/env bash
set -euo pipefail
R=/mnt/BioDisk/wangmingyuan/GI_multiomics_screening
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
RS=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
"$RS" --vanilla -e "invisible(parse(file='$R/09_comparison/kinase_screen.R'))"
"$PY" -m py_compile "$R/09_comparison/verify_statistics.py"
echo SYNTAX_CHECK:PASS
echo LOGIC_CHECK:PSP_human_exact_isoform_coordinates_5sites_3genes_patient_units_and_independent_Scipy_verification
echo "RUN_COMMAND: $RS --vanilla $R/09_comparison/kinase_screen.R"
echo EXPECTED_SUCCESS_MARKER:KINASE_SCREEN_COMPLETE
"$RS" --vanilla "$R/09_comparison/kinase_screen.R"
echo "RUN_COMMAND: $PY $R/09_comparison/verify_statistics.py"
echo EXPECTED_SUCCESS_MARKER:STATISTICAL_VERIFICATION_COMPLETE
"$PY" "$R/09_comparison/verify_statistics.py"
