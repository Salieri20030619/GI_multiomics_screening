from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,zipfile,gzip,xml.etree.ElementTree as ET,re
from scipy import stats
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation';B=R/'B_ICC_lymphatic_invasion'
checks=[]
def ck(name,ok,detail=''):
 checks.append({'check':name,'pass':bool(ok),'detail':str(detail)})
def rd(p):return pd.read_csv(p,sep='\t')
for P in [A,B]:
 for f in sorted((P/'06_tables').glob('*.tsv*')):
  d=rd(f)
  for col in ['P','FDR','P_unmoderated','P_Welch']:
   if col in d:ck(str(f.relative_to(R))+':'+col+'_bounds',d[col].dropna().between(0,1).all())
  if set(['effect','CI_low','CI_high']).issubset(d):
   z=d.dropna(subset=['effect','CI_low','CI_high']);ck(str(f.relative_to(R))+':CI',((z.CI_low<=z.effect+1e-10)&(z.CI_high>=z.effect-1e-10)).all())
  if 'family' in d and 'FDR' in d and 'P' in d:
   for fam,z in d.groupby(['cohort','family'] if 'cohort' in d.columns else 'family'):
    if fam=='descriptive' or (isinstance(fam,tuple) and fam[-1]=='descriptive'):continue
    good=z.P.notna()
    if good.sum():
     # R p.adjust includes NA entries in family length only where present? All tested main families should be complete.
     p=z.loc[good,'P'].values;want=stats.false_discovery_control(p,method='bh');obs=z.loc[good,'FDR'].values
     ck(str(f.relative_to(R))+':BH:'+str(fam),np.allclose(want,obs,rtol=1e-6,atol=1e-10,equal_nan=True),f'n={len(p)}')
for kind in ['protein','site']:
 x=rd(A/f'03_processed/{kind}_paired_delta.tsv.gz').set_index('feature');st=rd(A/f'06_tables/{kind}_paired_complete.tsv.gz').set_index('feature');a=x.mean(axis=1)
 ck('A_'+kind+'_all_mean_effects',np.allclose(a,st.loc[a.index,'effect'],equal_nan=True,atol=1e-10))
 for feat in st.sort_values('P').head(5).index:
  v=x.loc[feat].dropna().values;pv=stats.ttest_1samp(v,0).pvalue;ck('A_'+kind+'_'+feat+'_paired_t',np.isclose(pv,st.loc[feat,'P_unmoderated'],rtol=1e-6,atol=1e-12))
for co in ['discovery','validation']:
 x=rd(B/f'03_processed/{co}_log2.tsv.gz').set_index('feature');m=rd(B/f'01_metadata/{co}_samples.tsv');m=m[(m.Group=='Tumor')&m.L.notna()];s=rd(B/f'06_tables/{co}_protein_complete.tsv.gz').set_index('feature')
 pos=x[m.loc[m.L==1,'sample']];neg=x[m.loc[m.L==0,'sample']];effects=pos.mean(axis=1)-neg.mean(axis=1)
 ck('B_'+co+'_all_effects',np.allclose(effects,s.loc[effects.index,'effect'],atol=1e-10,equal_nan=True))
 ck('B_'+co+'_unique_patients',m.patient.is_unique)
 for feat in s.sort_values('P').head(5).index:
  pv=stats.ttest_ind(pos.loc[feat].dropna(),neg.loc[feat].dropna(),equal_var=False).pvalue;ck('B_'+co+'_'+feat+'_Welch',np.isclose(pv,s.loc[feat,'P_Welch'],rtol=1e-6,atol=1e-12))
# Audit original file bytes and report final availability independently of failed earlier routes.
manifest=[]
for P in [A,B]:
 for p in sorted((P/'02_raw').rglob('*')):
  if not p.is_file() or p.name.endswith('.part') or p.suffix not in ['.xlsx','.xls','.zip','.gz','.docx']:continue
  sig=p.read_bytes()[:8];integrity='signature_and_sha256'
  if sig[:4]==b'PK\x03\x04':
   with zipfile.ZipFile(p) as z:ck(str(p.relative_to(R))+':CRC',z.testzip() is None)
   integrity='ZIP_CRC_and_sha256'
  elif p.suffix=='.gz':
   with gzip.open(p,'rb') as f:
    while f.read(1024*1024):pass
   integrity='gzip_CRC_and_sha256'
  manifest.append({'project':P.name,'file':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'integrity':integrity,'signature_hex':sig.hex()})
pd.DataFrame(manifest).to_csv(R/'09_comparison/final_raw_integrity.tsv',sep='\t',index=False)
for P in [A,B]:pd.DataFrame([m for m in manifest if m['project']==P.name]).to_csv(P/'01_metadata/final_raw_integrity.tsv',sep='\t',index=False)
# Exact format/pair counts, independent checks not hardcoded selection rules.
mapping=rd(A/'01_metadata/patient_sample_omics.tsv');ck('A_protein_samples101',mapping.protein_in_matrix.sum()==101);ck('A_site_samples99',mapping.site_in_matrix.sum()==99)
ss=rd(A/'04_qc/site_sample_qc.tsv');ck('A_S2C_zscore_confirmed',np.allclose(ss['mean'],0,atol=1e-12)&np.allclose(ss.sd,1,atol=1e-12))
cm=rd(B/'04_qc/discovery_event_audit.tsv');ck('B_event1_means_recorded_recurrence',np.array_equal(cm.censored_raw==1,cm.has_rec_days))
pd.DataFrame(checks).to_csv(R/'09_comparison/statistical_integrity_checks.tsv',sep='\t',index=False)
summary={'n_checks':len(checks),'passed':sum(x['pass'] for x in checks),'failed':[x for x in checks if not x['pass']]};(R/'09_comparison/statistical_verification.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if summary['failed']:raise SystemExit('STATISTICAL_VERIFICATION_FAILED')
print('STATISTICAL_VERIFICATION_COMPLETE')
