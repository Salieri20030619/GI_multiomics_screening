from pathlib import Path
import requests,json,concurrent.futures,hashlib,zipfile,xml.etree.ElementTree as E
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation';B=R/'B_ICC_lymphatic_invasion'
tasks=[(B/'00_docs/ICC_github_tree.json','https://api.github.com/repos/SchillingLabProteomics/ICC/git/trees/main?recursive=1'),(A/'00_docs/Cell2020_publisher.html','https://www.cell.com/cancer-cell/fulltext/S1535-6108(20)30413-X')]
tasks += [(A/f'02_raw/Cell2020_mmc{i}.xlsx',f'https://ars.els-cdn.com/content/image/1-s2.0-S153561082030413X-mmc{i}.xlsx') for i in range(2,9)]
tasks += [(A/'00_docs/Cell2020_mmc1.pdf','https://ars.els-cdn.com/content/image/1-s2.0-S153561082030413X-mmc1.pdf')]
tasks += [(A/'00_docs/A_oa.xml','https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC12822478'),(A/'00_docs/external_oa.xml','https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC12335997')]
def get(t):
 p,u=t;rec={'url':u,'file':str(p)}
 try:
  if p.exists():rec['status']='existing';return rec
  r=requests.get(u,timeout=(15,40));rec['http']=r.status_code;r.raise_for_status()
  content=r.content
  if p.suffix=='.xlsx':
   import io
   with zipfile.ZipFile(io.BytesIO(content)) as z:assert z.testzip() is None
  if p.suffix=='.pdf':assert content.startswith(b'%PDF')
  p.write_bytes(content);rec.update(status='verified' if p.suffix in ['.pdf','.xlsx'] else 'downloaded',bytes=len(content),sha256=hashlib.sha256(content).hexdigest())
 except Exception as e:rec['error']=str(e);rec['status']='failed'
 print(json.dumps(rec),flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:res=list(pool.map(get,tasks))
(R/'09_comparison/extra_manifest.json').write_text(json.dumps(res,indent=2))
tree=B/'00_docs/ICC_github_tree.json'
if tree.exists():
 js=json.loads(tree.read_text());print('REPO_PATHS',[x['path'] for x in js.get('tree',[])],flush=True)
 for x in js.get('tree',[]):
  if x['path'].endswith(('.R','.Rmd','.r','.md','.txt')) and x.get('size',0)<200000:
   p=B/'00_docs/author_code'/x['path'];p.parent.mkdir(parents=True,exist_ok=True)
   get((p,'https://raw.githubusercontent.com/SchillingLabProteomics/ICC/main/'+x['path']))
print('EXTRA_LOOKUP_ATTEMPTS_COMPLETE')
