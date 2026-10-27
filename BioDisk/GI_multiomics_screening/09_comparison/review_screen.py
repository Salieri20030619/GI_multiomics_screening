from pathlib import Path
import pandas as pd,numpy as np,json
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation';B=R/'B_ICC_lymphatic_invasion'
def rd(p):return pd.read_csv(p,sep='\t')
for P in [A,B]:
 print('\nPROJECT',P.name)
 for f in sorted((P/'06_tables').glob('*')):
  d=rd(f)
  if 'FDR' not in d:continue
  print(f.name,'nrows',len(d),'minFDR',d.FDR.min(),'FDR05',(d.FDR<.05).sum())
  cols=[c for c in ['feature','gene','site','effect','P','FDR','n','n_pos','n_neg','family'] if c in d]
  print(d.sort_values('FDR').head(6)[cols].to_string(index=False))
 for f in sorted((P/'04_qc').glob('*sample_qc.tsv')):
  d=rd(f);print('QC',f.name);print(d.sort_values('n_observed').head(4).to_string(index=False))
d=rd(B/'04_qc/discovery_event_audit.tsv');print('EVENT_CROSSTAB',pd.crosstab(d.censored_raw,d.has_rec_days).to_dict());print('EVENT_EXCEPTION',d.loc[(d.censored_raw==1)&~d.event_equals_rec|(d.censored_raw==0)&~d.event_equals_fu].to_string(index=False))
ex=rd(A/'01_metadata/external_specimens.tsv');print('EXTERNAL',ex.groupby('Batch').agg(specimens=('Sample','size'),patients=('Patient_Nr','nunique')));print('RECURRENCE',ex.Recurrence.value_counts(dropna=False).to_dict())
print('REVIEW_SUMMARY_COMPLETE')
