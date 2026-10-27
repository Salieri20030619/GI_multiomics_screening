from pathlib import Path
import json,openpyxl,zipfile,pandas as pd
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening')
for p in R.glob('*/04_qc/source_workbooks.json'):
 for s in json.loads(p.read_text()):
  if any(x in s['sheet'] for x in ['Somatic','Clinical factor','mutated','GSEA','IPA','Glyco']):continue
  print(Path(s['source']).name,s['sheet'],s['shape'])
  for row in s['preview'][:3]:print(row)
for p in (R/'A_CRLM_phosphoregulation/02_raw').glob('Cell2020*.xlsx'):
 if p.name=='Cell2020_mmc2.xlsx':continue
 w=openpyxl.load_workbook(p,read_only=True,data_only=True)
 for s in w:
  print(p.name,s.title,s.max_row,s.max_column)
  print(list(s.values)[:2] if s.max_row and s.max_row<100 else list(next(s.values)))
print('INPUT_SUMMARY_COMPLETE')
