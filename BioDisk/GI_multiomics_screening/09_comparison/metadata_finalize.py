from pathlib import Path
import pandas as pd,numpy as np,json,re,hashlib,requests,concurrent.futures
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation';B=R/'B_ICC_lymphatic_invasion'
def rd(p):return pd.read_csv(p,sep='\t')
em=rd(A/'01_metadata/external_specimens.tsv');print('EXTERNAL STATUS',em['survival status'].value_counts(dropna=False).to_dict());print('BYBATCH',pd.crosstab(em.Batch,em['survival status']).to_dict())
sm=rd(A/'01_metadata/site_protein_mapping.tsv');pm=rd(A/'01_metadata/protein_features.tsv')
accs=set(pm.accession);base=set(pm.accession.str.replace(r'-\d+$','',regex=True));counts=pm.groupby(['accession','gene']).size().to_dict()
def why(r):
 if not isinstance(r.gene,str) or not re.fullmatch('[A-Za-z0-9_-]+',r.gene):return 'ambiguous_or_missing_gene'
 if not re.fullmatch('[A-Z0-9]+(?:-\d+)?',str(r.accession)):return 'ambiguous_accession'
 if not re.fullmatch('[STY]\d+',str(r.residue_position)):return 'ambiguous_residue_position'
 n=counts.get((r.accession,r.gene),0)
 if n>1:return 'multiple_totalprotein_rows_same_key'
 if n==1:return 'exact_unique_accession_gene_site'
 if r.accession in accs:return 'accession_present_gene_label_mismatch'
 if re.sub(r'-\d+$','',r.accession) in base:return 'isoform_mismatch_sequence_alignment_required'
 return 'totalprotein_accession_not_measured'
sm['mapping_detail']=sm.apply(why,axis=1);sm.to_csv(A/'01_metadata/site_protein_mapping.tsv',sep='\t',index=False)
summary=sm.mapping_detail.value_counts().rename_axis('reason').reset_index(name='n_sites');summary.to_csv(A/'04_qc/site_mapping_exclusion_summary.tsv',sep='\t',index=False);print(summary.to_string(index=False))
g=rd(B/'01_metadata/GSE107943_sample_registry.tsv')
print('GEO107943 CHAR',g[['title','characteristics_ch1']].head(6).to_string(index=False))
# Parse matrix sample IDs separately from GEO label records and preserve the mapping assumptions.
cols=pd.read_csv(B/'02_raw/GSE107943_rawread.txt.gz',sep='\t',nrows=0).columns
rows=[]
for c in cols:
 m=re.fullmatch(r'(\d+)([NT])',c)
 if m:rows.append({'sample':c,'patient':m[1],'tissue':'tumor' if m[2]=='T' else 'normal'})
mm=pd.DataFrame(rows);mm['paired']=mm.patient.isin(mm.groupby('patient').tissue.nunique().loc[lambda x:x==2].index);mm.to_csv(B/'01_metadata/GSE107943_matrix_sample_mapping.tsv',sep='\t',index=False);print('RNA_MAPPING',mm.tissue.value_counts().to_dict(),'pairs',mm.loc[mm.paired,'patient'].nunique())
files=json.loads((A/'00_docs/iprox_lookup/subproject.json').read_text())['data']['subdatafilesInfo']
xmls=[f for f in files if f['fileName'].endswith('.xml')]
def get(f):
 u='https://www.iprox.cn/FAF016Controller/readXml.jsonp?fileId='+f['fileId'];p=A/'00_docs/iprox_lookup'/f['fileName']
 try:
  r=requests.get(u,timeout=(15,45));r.raise_for_status();p.write_bytes(r.content)
  print('CONFIG',p.name,len(r.content),'SHA1_MATCH',hashlib.sha1(r.content).hexdigest()==f['sha1'],flush=True)
 except Exception as e:print('CONFIG_FAIL',f['fileName'],repr(e),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(get,xmls))
print('METADATA_FINALIZATION_COMPLETE')
