from pathlib import Path
import pandas as pd,json
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation'
for s in json.loads((A/'04_qc/source_workbooks.json').read_text()):
 if any(x in s['source'] for x in ['s005.xls','s006.xls']) or s['sheet']=='5 Total Cohort Abundances':
  print(Path(s['source']).name,s['sheet'],s['shape']);print([r[:15] for r in s['preview'][:6]])
for p in (R/'B_ICC_lymphatic_invasion/00_docs').glob('GSE*.brief.soft.txt'):
 print(p.name)
 lines=[s for s in p.read_text().splitlines() if any(x in s for x in ['!Sample_title','!Sample_source_name','!Sample_characteristics','!Sample_supplementary_file'])]
 print('\n'.join(lines[:100]))
