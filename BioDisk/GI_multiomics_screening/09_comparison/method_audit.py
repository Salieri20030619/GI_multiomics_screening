from pathlib import Path
import xml.etree.ElementTree as E,json,openpyxl,subprocess
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening')
for p in R.glob('*/00_docs/*full.xml'):
 t=E.parse(p);print('\nARTICLE',p.name)
 for x in t.findall('.//sec'):
  title=x.find('title');s=''.join(title.itertext()) if title is not None else ''
  if any(v in s.lower() for v in ['statistic','data analy','quantifi','mass spectrom','availability','proteomic analy']):print(s, ''.join(x.itertext()))
 for x in t.findall('.//supplementary-material'):
  print('ATTACH',x.get('id'),[(y.tag,y.attrib) for y in x.iter() if y.tag in ['media','ext-link']])
for p in R.glob('*/02_raw/*.xlsx'):
 if not any(s in p.name for s in ['MOESM3_','MOESM11_']):continue
 with p.open('rb') as f:
  w=openpyxl.load_workbook(f,read_only=True,data_only=True)
  for s in w:
   s.reset_dimensions();rows=list(s.values)
   import pandas as pd
   d=pd.DataFrame(rows[1:],columns=rows[0]);out=p.parent.parent/'01_metadata'/p.with_suffix('.tsv').name;d.to_csv(out,sep='\t',index=False)
   print('CLIN',p.name,'COLUMNS',d.columns.tolist())
   for c in d.columns:print(c,d[c].value_counts(dropna=False).to_dict() if d[c].nunique()<20 else list(d[c].head(5)))
print('METHOD_AUDIT_COMPLETE')
