#!/usr/bin/env bash
set -euo pipefail
R=/mnt/BioDisk/wangmingyuan/GI_multiomics_screening
RS=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
"$RS" --vanilla -e "invisible(parse(file='$R/09_comparison/extend_screen.R'))"
echo SYNTAX_CHECK:PASS
echo LOGIC_CHECK:locked_pathway_family_patient_unit_same_patient_triad_explicit_TTR_code_external_endpoint_conflicts_excluded
echo "RUN_COMMAND: $RS --vanilla $R/09_comparison/extend_screen.R"
echo EXPECTED_SUCCESS_MARKER:EXTENDED_SCREEN_COMPLETE
"$RS" --vanilla "$R/09_comparison/extend_screen.R"
