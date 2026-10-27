from pathlib import Path
import pandas as pd,numpy as np,json
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening')
for project in ['A_CRLM_phosphoregulation','B_ICC_lymphatic_invasion']:
 f=R/project/'04_qc/source_workbooks.json'
 if not f.exists():continue
 for s in json.loads(f.read_text()):
  if project.startswith('A') and not (s['sheet'] in ['S2B','S2C','Table S1A','Table S1B']):continue
  if project.startswith('B') and not 'MOESM' in s['source']:continue
  d=pd.read_csv(s['table'],sep='\t',header=None,keep_default_na=False)
  print('\nSOURCE',Path(s['source']).name,'SHEET',s['sheet'],'SHAPE',d.shape)
  print('HEADER',d.iloc[0].tolist())
  print('ROW1',d.iloc[1,:15].tolist() if len(d)>1 else [])
  print('ROW2',d.iloc[2,:15].tolist() if len(d)>2 else [])
  if len(d)>500:
   n=d.iloc[1:].apply(pd.to_numeric,errors='coerce');n=n.loc[:,n.notna().mean()>.1]
   q=pd.DataFrame({'column':n.columns,'n':n.notna().sum(),'mean':n.mean(),'sd':n.std(),'median':n.median(),'min':n.min(),'max':n.max(),'zero':(n==0).sum()})
   print(q.iloc[:12].to_string(index=False));print('MEAN_SD_RANGE',q[['mean','sd','median']].agg(['min','max']).to_string())
   q.to_csv(R/project/'04_qc'/(Path(s['source']).stem+'_'+s['sheet'].replace(' ','_')+'_scale.tsv'),sep='\t',index=False)
  else:
   print(d.iloc[:6].to_string(index=False,header=False))
print('SCALE_AUDIT_COMPLETE')
