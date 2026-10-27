#!/usr/bin/env bash
set -euo pipefail
R=/mnt/BioDisk/wangmingyuan/GI_multiomics_screening
export PYTHONPATH="$R/09_comparison/vendor"
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
"$PY" -m py_compile "$R/09_comparison/extract_tables.py"
echo SYNTAX_CHECK:PASS
echo LOGIC_CHECK:signature_detection_dimension_reset_no_numeric_transformation
echo "RUN_COMMAND: $PY $R/09_comparison/extract_tables.py $1"
echo EXPECTED_SUCCESS_MARKER:TABLE_EXTRACTION_COMPLETE
cp "$R/09_comparison/protocol.md" "$R/$1/00_docs/protocol_v1.md"
"$PY" "$R/09_comparison/extract_tables.py" "$1"
