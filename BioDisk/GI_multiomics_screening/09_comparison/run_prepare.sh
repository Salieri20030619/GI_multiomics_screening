#!/usr/bin/env bash
set -euo pipefail
R=/mnt/BioDisk/wangmingyuan/GI_multiomics_screening
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
RS=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
"$PY" -m py_compile "$R/09_comparison/prepare_analysis.py"
"$RS" --vanilla -e "invisible(parse(file='$R/09_comparison/map_reference.R'))"
echo SYNTAX_CHECK:PASS
echo LOGIC_CHECK:unique_patients_labels_explicit_zeros_preserved_exact_mapping_scale_amendment_before_tests
echo "RUN_COMMAND: $PY $R/09_comparison/prepare_analysis.py"
echo EXPECTED_SUCCESS_MARKER:ANALYSIS_INPUT_PREPARATION_COMPLETE
"$PY" "$R/09_comparison/prepare_analysis.py"
echo "RUN_COMMAND: $RS --vanilla $R/09_comparison/map_reference.R"
echo EXPECTED_SUCCESS_MARKER:REFERENCE_MAPPING_COMPLETE
"$RS" --vanilla "$R/09_comparison/map_reference.R"
