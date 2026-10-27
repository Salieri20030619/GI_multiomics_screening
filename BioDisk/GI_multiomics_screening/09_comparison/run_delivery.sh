#!/usr/bin/env bash
set -euo pipefail
R=/mnt/BioDisk/wangmingyuan/GI_multiomics_screening
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
RS=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
"$PY" -m py_compile "$R/09_comparison/build_summary.py" "$R/09_comparison/build_reports.py" "$R/09_comparison/package_delivery.py"
"$RS" --vanilla "$R/09_comparison/check_all_syntax.R" > "$R/09_comparison/all_R_syntax.log"
echo SYNTAX_CHECK:PASS
echo LOGIC_CHECK:actual_result_values_separate_unmet_stages_patient_metadata_no_raw_matrices_complete_statistics_30MiB_gate
"$PY" "$R/09_comparison/build_summary.py" > "$R/09_comparison/summary.log"
echo "RUN_COMMAND: $PY $R/09_comparison/build_reports.py"
echo EXPECTED_SUCCESS_MARKER:REPORTS_COMPLETE
"$PY" "$R/09_comparison/build_reports.py"
echo "RUN_COMMAND: $PY $R/09_comparison/package_delivery.py"
echo EXPECTED_SUCCESS_MARKER:DELIVERY_COMPLETE
"$PY" "$R/09_comparison/package_delivery.py"
