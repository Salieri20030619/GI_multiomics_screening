from pathlib import Path
import requests,zipfile,json,hashlib,concurrent.futures
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation'
def get(name):
 p=A/'02_raw'/name
 for base in ['https://pmc.ncbi.nlm.nih.gov/articles/instance/12822478/bin/','https://pmc.ncbi.nlm.nih.gov/articles/PMC12822478/bin/']:
  url=base+name
  try:
   r=requests.get(url,timeout=(15,60));r.raise_for_status()
   b=r.content
   if not (b.startswith(b'PK') or b.startswith(bytes.fromhex('d0cf11e0'))):raise ValueError('not workbook/docx')
   p.write_bytes(b)
   if zipfile.is_zipfile(p):
    with zipfile.ZipFile(p) as z:assert z.testzip() is None
   rec={'url':url,'file':str(p),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'status':'verified'}
   (A/'logs'/(name+'.direct.json')).write_text(json.dumps(rec,indent=2));print(rec,flush=True);return
  except Exception as e:print(url,repr(e),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(get,['ADVS-13-e11744-s009.xls','ADVS-13-e11744-s006.xls','ADVS-13-e11744-s004.docx']))
print('A_DIRECT_ATTEMPTS_COMPLETE')
