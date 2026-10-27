from pathlib import Path
import pandas as pd,numpy as np,json,re,xml.etree.ElementTree as ET
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation';B=R/'B_ICC_lymphatic_invasion'
def rd(p):return pd.read_csv(p,sep='\t')
summary={'A':{},'B':{}}
for lab,P in [('A',A),('B',B)]:
 rows=[]
 for p in sorted((P/'06_tables').glob('*.tsv*')):
  d=rd(p)
  if 'FDR' not in d or 'P' not in d:continue
  group=['cohort','mode'] if 'cohort' in d and 'mode' in d else ['cohort','family'] if 'cohort' in d and 'family' in d else ['family'] if 'family' in d else ['cohort'] if 'cohort' in d else None
  for f,v in d.groupby(group) if group else [('all',d)]:
   rows.append({'table':p.name,'family':str(f),'n_rows':len(v),'n_tested':int(v.P.notna().sum()),'nominal_P05':int((v.P<.05).sum()),'FDR05':int((v.FDR<.05).sum()),'minimum_FDR':v.FDR.min()})
 st=pd.DataFrame(rows);st.to_csv(P/'06_tables/statistical_family_summary.tsv',sep='\t',index=False);summary[lab]['families']=rows
sm=rd(A/'01_metadata/site_protein_mapping.tsv');cond=rd(A/'06_tables/site_conditional_diagnostic.tsv.gz');summary['A']['exact_matches']=int((sm.mapping_detail=='exact_unique_accession_gene_site').sum());summary['A']['conditional_primary_n']=int((cond.family=='primary_diagnostic').sum())
ext=rd(A/'01_metadata/external_patient_metadata.tsv');summary['A']['external_patient_n']=len(ext);summary['A']['external_survival_n']=int((ext.OS_event.notna()&ext.OS_years.notna()&(ext.OS_years>0)).sum());summary['A']['external_events']=int(ext.loc[ext.OS_event.notna()&ext.OS_years.notna()&(ext.OS_years>0),'OS_event'].sum())
summary['A']['external_survival_by_center']=ext.loc[ext.OS_event.notna()&ext.OS_years.notna()].groupby('Batch').size().to_dict()
conf=rd(A/'04_qc/external_patient_field_conflicts.tsv');summary['A']['external_conflict_patients']=conf.patient.nunique();summary['A']['survival_time_conflicts']=conf.loc[conf.field=='survival time (Years)','patient'].tolist()
for co in ['discovery','validation']:
 m=rd(B/f'01_metadata/{co}_samples.tsv');t=m[m.Group=='Tumor'];summary['B'][co]={'n_tumor':len(t),'L':t.L.value_counts(dropna=False).to_dict(),'patients':t.patient.nunique()}
 bv=rd(B/f'06_tables/{co}_protein_complete.tsv.gz');summary['B'][co]['unmapped_proteins']=int(bv.gene.isna().sum())
 if co=='validation':
  for col in ['Lymphangiosis','Haemangiosis','Perineural.Sheath.Invasion','fustat','X5yrsurv']:summary['B'][co][col]=t[col].value_counts(dropna=False).to_dict()
  pd.crosstab(t.fustat,t.X5yrsurv,dropna=False).to_csv(B/'04_qc/validation_survival_code_crosstab.tsv',sep='\t')
 geo=rd(B/'01_metadata/GSE107943_sample_registry.tsv')
records=[]
for row in geo.itertuples():
 rec={'accession':row.accession,'sample':row.title}
 for s in str(row.characteristics_ch1).split(' | '):
  if ': ' in s:k,v=s.split(': ',1);rec[k]=v
 records.append(rec)
gm=pd.DataFrame(records);gm.to_csv(B/'01_metadata/GSE107943_clinical_labels.tsv',sep='\t',index=False)
summary['B']['GSE107943_available_fields']=gm.columns.tolist()
files=rd(A/'01_metadata/iprox_public_file_registry.tsv');summary['A']['iprox_n_raw']=int((files.fileType=='RAW').sum());summary['A']['iprox_raw_reported_size_sum']=float(pd.to_numeric(files.loc[files.fileType=='RAW','fileSize']).sum());summary['A']['iprox_size_unit']='provider_fileSize_unconfirmed'
# Parse public XML representations, retaining byte-integrity caveat.
cfg=[]
for p in (A/'00_docs/iprox_lookup').glob('*.xml'):
 try:
  root=ET.parse(p).getroot();experiments=[x.text for x in root.findall('.//experiments/string')];filepaths=[x.text for x in root.findall('.//filePaths/string')]
  txt=p.read_text(errors='replace');patient_ids=sorted(set(re.findall(r'\bP\d+_(?:LM|T|N)\b',txt)))
  cfg.append({'file':p.name,'root':root.tag,'n_raw_paths':len(filepaths),'experiments':sorted(set(experiments)),'patient_tissue_labels':patient_ids,'integrity':'endpoint_XML_representation_not_byte_identical_to_uploaded_file'})
 except Exception as e:cfg.append({'file':p.name,'error':str(e),'integrity':'not_verified_XML'})
(A/'04_qc/maxquant_parameter_audit.json').write_text(json.dumps(cfg,indent=2));summary['A']['maxquant_configs']=cfg
def clean(v):
 if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
 if isinstance(v,list):return [clean(x) for x in v]
 if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
 return v
summary=clean(summary)
(R/'09_comparison/result_summary.json').write_text(json.dumps(summary,indent=2,default=str,allow_nan=False));print(json.dumps(summary,indent=2,default=str)[:13000]);print('RESULT_SUMMARY_COMPLETE')
