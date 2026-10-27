from pathlib import Path
import requests,concurrent.futures,json,hashlib,gzip
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation';B=R/'B_ICC_lymphatic_invasion'
tasks=[(A/'00_docs/iprox_proxi.json','https://www.iprox.cn/proxi/datasets/IPX0007391001'),(A/'00_docs/iprox_project.html','https://www.iprox.cn/page/project.html?id=IPX0007391001'),(A/'00_docs/NODE_OEP000729.html','https://www.biosino.org/node/project/detail/OEP000729'),(R/'09_comparison/reference/KSEAdata.RData','https://raw.githubusercontent.com/cran/KSEAapp/master/data/KSData.RData')]
for g in ['GSE107943','GSE138709','GSE26566']:
 tasks.append((B/f'00_docs/{g}.brief.soft.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}&targ=gsm&form=text&view=brief'))
def get(t):
 p,u=t;rec={'url':u,'file':str(p)}
 try:
  q=requests.get(u,timeout=(15,75));q.raise_for_status();p.write_bytes(q.content);rec.update(status='downloaded',bytes=len(q.content),sha256=hashlib.sha256(q.content).hexdigest())
 except Exception as e:rec.update(status='failed',error=str(e))
 print(rec,flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:res=list(pool.map(get,tasks))
(R/'09_comparison/resource_followup.json').write_text(json.dumps(res,indent=2))
print('RESOURCE_FOLLOWUP_ATTEMPTS_COMPLETE')
