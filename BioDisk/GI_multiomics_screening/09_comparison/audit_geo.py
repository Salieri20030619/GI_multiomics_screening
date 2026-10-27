from pathlib import Path
import re,pandas as pd,numpy as np,json,gzip
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening/B_ICC_lymphatic_invasion')
for g in ['GSE107943','GSE138709']:
 p=R/f'00_docs/{g}.brief.soft.txt'
 if not p.exists():continue
 samples=[]
 for block in p.read_text().split('^SAMPLE = ')[1:]:
  d={'accession':block.splitlines()[0]}
  for line in block.splitlines():
   if line.startswith('!Sample_') and ' = ' in line:
    k,v=line.split(' = ',1);k=k[8:];d[k]=d.get(k,'')+(' | ' if k in d else '')+v
  samples.append(d)
 md=pd.DataFrame(samples);md.to_csv(R/f'01_metadata/{g}_sample_registry.tsv',sep='\t',index=False)
 print(g,'n_GSM',len(md));print(md[[c for c in ['accession','title','source_name_ch1','characteristics_ch1'] if c in md]].to_string(index=False))
 if g=='GSE107943':
  print('phenotype_labels_present',[(c,md[c].str.contains('surviv|invas|recurr|lymph|follow',case=False,na=False).sum()) for c in md if 'character' in c])
p=R/'02_raw/GSE107943_rawread.txt.gz'
if p.exists():
 d=pd.read_csv(p,sep='\t');print('RNA_SHAPE',d.shape,'COLS',d.columns.tolist());print(d.head(2).iloc[:,:10].to_string(index=False))
 report={'rows':len(d),'columns':d.columns.tolist(),'numeric_nonnegative':bool((d.select_dtypes('number')>=0).all().all())}
 (R/'04_qc/GSE107943_matrix_audit.json').write_text(json.dumps(report,indent=2))
print('GEO_SOURCE_AUDIT_COMPLETE')
