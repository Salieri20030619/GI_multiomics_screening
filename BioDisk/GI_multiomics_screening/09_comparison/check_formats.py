from pathlib import Path
import zipfile,json,sys
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening')
for p in R.glob('A*/02_raw/A_supp/*.xls'):
 print(p.name,p.stat().st_size,p.read_bytes()[:8],flush=True)
 if zipfile.is_zipfile(p):
  z=zipfile.ZipFile(p);print(z.namelist()[:20],flush=True)
  if '[Content_Types].xml' in z.namelist():print(z.read('[Content_Types].xml')[:1300],flush=True)
