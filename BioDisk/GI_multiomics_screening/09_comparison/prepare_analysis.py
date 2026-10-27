from pathlib import Path
import pandas as pd,numpy as np,json,re,hashlib
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation';B=R/'B_ICC_lymphatic_invasion'
def table(P,source,sheet=None):
 m=json.loads((P/'04_qc/source_workbooks.json').read_text());m=[s for s in m if source in s['source'] and (sheet is None or s['sheet']==sheet)];assert len(m)==1,(source,sheet,len(m))
 return pd.read_csv(m[0]['table'],sep='\t',header=0,keep_default_na=False)
def numeric(d):return d.apply(pd.to_numeric,errors='coerce')
def save(d,p,index=True):d.to_csv(p,sep='\t',index=index,index_label='feature' if index else None)
acc=set()
for cohort,mi,ci in [('discovery',4,3),('validation',7,11)]:
 d=table(B,f'MOESM{mi}_');md=table(B,f'MOESM{ci}_')
 assert md.Sample.is_unique;assert set(d.columns[1:])==set(md.Sample)
 md=md.set_index('Sample').loc[d.columns[1:]].rename_axis('Sample').reset_index();md['cohort']=cohort
 md['patient']=cohort+'_'+md.Patient_No.astype(str);md['sample']=md.Sample
 field='Lymphangioinvasion' if cohort=='discovery' else 'Lymphangiosis'
 md['L']=md[field].map({'Yes':1,'No':0,'L1':1,'L0':0})
 md['sex']=md.Sex.str.lower().map({'m':1,'male':1,'f':0,'female':0})
 age={'<50':45,'50-59':55,'60-69':65,'70-79':75,'>79':85};md['age_ordinal']=md.age_group.map(age)
 if cohort=='discovery':md['duct']=md['Duct Type'].map({'SDT':0,'LDT':1})
 else:md['grade']=md.Grading.map({'G1':1,'G2':2,'G3':3})
 tumor=md.Group=='Tumor';assert md.loc[tumor,'patient'].is_unique
 save(md,B/f'01_metadata/{cohort}_samples.tsv',False)
 x=numeric(d.set_index('Protein'));assert x.index.is_unique
 save(x,B/f'03_processed/{cohort}_log2.tsv.gz');acc.update(x.index)
 qc=pd.DataFrame({'sample':x.columns,'n_observed':x.notna().sum(),'missing_fraction':x.isna().mean(),'median':x.median(),'mean':x.mean(),'sd':x.std()}).reset_index(drop=True).merge(md[['sample','patient','Group','L']])
 save(qc,B/f'04_qc/{cohort}_sample_qc.tsv',False)
 print(cohort,'tumor',int(tumor.sum()),'L counts',md.loc[tumor,'L'].value_counts(dropna=False).to_dict(), 'proteins',len(x))
 if cohort=='discovery':
  t=md.loc[tumor].copy();t['rec_days']=pd.to_numeric(t.daystorecurrence,errors='coerce');t['fu_days']=pd.to_numeric(t.daystolastfu,errors='coerce');t['event_days']=pd.to_numeric(t.daystoevent,errors='coerce');t['censored_raw']=pd.to_numeric(t.censored,errors='coerce')
  t['has_rec_days']=t.rec_days.notna();t['event_equals_rec']=np.isclose(t.event_days,t.rec_days,equal_nan=False);t['event_equals_fu']=np.isclose(t.event_days,t.fu_days,equal_nan=False)
  save(t[['sample','patient','daystorecurrence','daystolastfu','censored_raw','daystoevent','Alive','cluster','has_rec_days','event_equals_rec','event_equals_fu']],B/'04_qc/discovery_event_audit.tsv',False)

p=table(A,'s009.xls','S2B');s=table(A,'s009.xls','S2C');clin=table(A,'s006.xls','Table S1A');mapping=table(A,'s006.xls','Table S1B')
mapping.Patient=mapping.Patient.replace('',np.nan).ffill();mapping['sample']=mapping.Patient+'_'+mapping.Sample
mapping['protein_in_matrix']=mapping['sample'].isin(p.columns);mapping['site_in_matrix']=mapping['sample'].isin(s.columns)
save(mapping,A/'01_metadata/patient_sample_omics.tsv',False);save(clin,A/'01_metadata/patients_clinical.tsv',False)
pm=p.iloc[:,:2].copy();pm.columns=['gene','accession'];pm['feature']=['PR'+str(i+1) for i in range(len(pm))];save(pm,A/'01_metadata/protein_features.tsv',False)
x=numeric(p.iloc[:,2:]);x.index=pm.feature
sm=s.iloc[:,:2].copy();sm.columns=['site','gene'];sm['feature']=['PH'+str(i+1) for i in range(len(sm))];sm['accession']=sm.site.str.rsplit('_',n=1).str[0];sm['residue_position']=sm.site.str.rsplit('_',n=1).str[1]
sm['site_syntax_ok']=sm.residue_position.str.fullmatch('[STY][0-9]+').fillna(False)
single=pm.gene.str.fullmatch('[A-Za-z0-9_-]+').fillna(False)&pm.accession.str.fullmatch('[A-Z0-9]+(?:-[0-9]+)?').fillna(False)
keys=pm.loc[single].groupby(['accession','gene']).feature.agg(list).to_dict()
matches=[]
for row in sm.itertuples():
 v=keys.get((row.accession,row.gene),[]);reason='exact_unique' if len(v)==1 and row.site_syntax_ok else ('site_ambiguous' if not row.site_syntax_ok else 'no_exact_unique_accession_gene')
 matches.append((v[0] if reason=='exact_unique' else '',reason))
sm['protein_feature']=[v[0] for v in matches];sm['match_status']=[v[1] for v in matches];sm['absolute_protein_adjustment_status']='not_estimable_site_sample_z_scores'
save(sm,A/'01_metadata/site_protein_mapping.tsv',False)
y=numeric(s.iloc[:,2:]);y.index=sm.feature
save(x,A/'03_processed/protein_log2.tsv.gz');save(y,A/'03_processed/site_z.tsv.gz');acc.update(pm.accession)
for name,v in [('protein',x),('site',y)]:
 q=pd.DataFrame({'sample':v.columns,'patient':[c.split('_')[0] for c in v.columns],'tissue':[c.split('_')[1] for c in v.columns],'n_observed':v.notna().sum().values,'missing_fraction':v.isna().mean().values,'mean':v.mean().values,'sd':v.std().values,'median':v.median().values})
 save(q,A/f'04_qc/{name}_sample_qc.tsv',False)
 pats=sorted(set(c.split('_')[0] for c in v.columns if c.endswith('_T'))&set(c.split('_')[0] for c in v.columns if c.endswith('_LM')),key=lambda z:int(z[1:]))
 delta=v[[a+'_LM' for a in pats]].copy();delta.columns=pats;t=v[[a+'_T' for a in pats]].copy();t.columns=pats;delta=delta-t
 save(delta,A/f'03_processed/{name}_paired_delta.tsv.gz')
 print('A',name,v.shape,'pairs',len(pats),'effective pairs',delta.notna().sum(axis=1).describe().to_dict())
print('A exact site matches',sm.match_status.value_counts().to_dict())
# Existing external cohort data retained at specimen level; patient aggregation occurs downstream.
em=table(A,'external_supp/mmc4.xlsx','1 Sample Overview');save(em,A/'01_metadata/external_specimens.tsv',False)
ext=table(A,'external_supp/mmc4.xlsx','5 Total Cohort Abundances');save(ext,A/'03_processed/external_author_abundances.tsv.gz',False)
print('EXTERNAL COLUMNS',ext.columns[:8].tolist());print('EXTERNAL patients',em.Patient_Nr.nunique(),'specimens',len(em))
(R/'09_comparison/accessions.txt').write_text('\n'.join(sorted(acc)))
amend='''# 尺度审计后的预设方案修订（主比较统计前）
S2C 99个样本均值0、样本SD1；S2B SD为0.35–0.81。补充方法明确磷酸化在批次去除后逐样本Z标准化。
因此不执行不可解释的 Z−log2 相减。位点原始配对差异单位为样本Z，不能报告为log2FC。
对精确匹配位点可另做患者内 siteZ差值~对应蛋白log2差值的条件线性模型，截距代表蛋白差值为0时的标准化位点差；仅作诊断性关联，不能替代用户要求的磷酸化量校正或称占有率。
严格量校正统计以不可估计记录；需要公开未Z标准化 reporter/IRS比值与批次映射。此缺口直接影响A立项优先级。
B矩阵为作者log2及中位数中心化的DIA-NN表达量，保留数值0及缺失NA，不再次log转换。统计阈值不变。
'''
(A/'00_docs/protocol_scale_amendment.md').write_text(amend)
print('ANALYSIS_INPUT_PREPARATION_COMPLETE')
