from pathlib import Path
import pandas as pd,numpy as np,openpyxl,zipfile,json,re,hashlib,sys
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening')
sys.path.insert(0,str(R/'09_comparison/vendor'))
from pyxlsb import open_workbook
import xlrd
project=sys.argv[1] if len(sys.argv)>1 else 'B_ICC_lymphatic_invasion'
P=R/project;manifest=[]
for p in sorted((P/'02_raw').rglob('*')):
 if p.suffix not in ['.xlsx','.xls']:continue
 if p.name.startswith('Cell2020_') and p.name!='Cell2020_mmc2.xlsx':continue
 dest=P/'03_processed/source_tables'/p.stem;dest.mkdir(parents=True,exist_ok=True)
 marker=dest/'extraction.json'
 if marker.exists():
  manifest.extend(json.loads(marker.read_text()));continue
 records=[]
 def save(s,rows,engine):
  # Preserve preamble, header and all source cell values; do not guess header here.
  frame=pd.DataFrame(rows);fname=re.sub(r'[^A-Za-z0-9._-]','_',s)+'.tsv.gz'
  frame.to_csv(dest/fname,sep='\t',index=False,header=False)
  rec={'source':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'sheet':s,'engine':engine,'shape':list(frame.shape),'table':str(dest/fname),'preview':frame.head(7).fillna('').astype(str).values.tolist()}
  records.append(rec);print(p.name,s,frame.shape,flush=True)
 if p.read_bytes()[:4]==b'PK\x03\x04':
  with zipfile.ZipFile(p) as z:isbin='xl/workbook.bin' in z.namelist();assert z.testzip() is None
  if isbin:
   with open_workbook(str(p)) as w:
    for s in w.sheets:
     with w.get_sheet(s) as sh:save(s,[[c.v for c in row] for row in sh.rows()], 'pyxlsb')
  else:
   with p.open('rb') as f:
    w=openpyxl.load_workbook(f,read_only=True,data_only=True)
    for s in w:
     s.reset_dimensions();save(s.title,list(s.values),'openpyxl_reset_dimensions')
 else:
  w=xlrd.open_workbook(str(p))
  for s in w.sheets():save(s.name,[s.row_values(i) for i in range(s.nrows)],'xlrd_OLE')
 marker.write_text(json.dumps(records,indent=2,ensure_ascii=False));manifest.extend(records)
(P/'04_qc/source_workbooks.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
print('TABLE_EXTRACTION_COMPLETE',project,flush=True)
