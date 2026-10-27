from pathlib import Path
import pandas as pd,json,collections,hashlib,zipfile,gzip,numpy as np
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation';B=R/'B_ICC_lymphatic_invasion'
p=A/'00_docs/iprox_lookup/subproject.json';d=json.loads(p.read_text())['data'];print('IPROX KEYS',list(d));files=d.get('subdatafilesInfo',[]);f=pd.DataFrame(files);f.to_csv(A/'01_metadata/iprox_public_file_registry.tsv',sep='\t',index=False)
print('IPROX FILE_COUNT',len(f),'TYPES',f.fileType.value_counts(dropna=False).to_dict())
print('NONRAW',f.loc[~f.fileName.str.lower().str.endswith('.raw'),['fileName','fileType','fileSize','sha1']].to_string(index=False))
for k,v in d.items():
 if k not in ['subdatafilesInfo'] and isinstance(v,(list,dict)):print(k,str(v)[:2400])
for name in ['kinase_patient_scores_complete.tsv','composition_paired_tests.tsv','external_pathway_OS_complete.tsv']:
 p=A/'06_tables'/name
 if p.exists():
  z=pd.read_csv(p,sep='\t');print(name);print(z.sort_values('FDR').head(8).to_string(index=False))
for p in [B/'06_tables/pathway_shared_gene_sensitivity_triad.tsv',B/'06_tables/discovery_pathway_TTR_association.tsv',A/'04_qc/external_patient_field_conflicts.tsv']:
 if p.exists():print(p.name);print(pd.read_csv(p,sep='\t').head(24).to_string(index=False))
print('EVIDENCE_AUDIT_INPUTS_COMPLETE')
