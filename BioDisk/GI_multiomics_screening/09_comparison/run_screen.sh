#!/usr/bin/env bash
set -euo pipefail
R=/mnt/BioDisk/wangmingyuan/GI_multiomics_screening
RS=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
"$RS" --vanilla -e "invisible(parse(file='$R/09_comparison/screen.R'))"
echo SYNTAX_CHECK:PASS
echo LOGIC_CHECK:patient_pairs_independent_B_fits_NA_preserved_distinct_FDR_families_fixed_validation_candidates_Z_scale_not_subtracted
echo "RUN_COMMAND: $RS --vanilla $R/09_comparison/screen.R"
echo EXPECTED_SUCCESS_MARKER:MAIN_SCREEN_COMPLETE
"$RS" --vanilla "$R/09_comparison/screen.R"
