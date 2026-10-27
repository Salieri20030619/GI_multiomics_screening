from pathlib import Path
import requests,concurrent.futures,json,re
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');D=R/'A_CRLM_phosphoregulation/00_docs/iprox_lookup';D.mkdir(exist_ok=True)
paths=['app/js/common.js','app/js/MSV02201.js','app/js/MSV02202.js','app/js/MSV02203.js','page/ProjectFileList.html?projectId=IPX0007391001']
def get(s):
 u='https://www.iprox.cn/'+s
 try:
  r=requests.get(u,timeout=(15,45));r.raise_for_status();name=s.split('/')[-1].replace('?','_');(D/name).write_text(r.text)
  print(u,r.status_code,len(r.content),flush=True)
  for line in r.text.splitlines():
   if ('url:' in line or 'jsonServer' in line or 'FileList' in line):print(line[:400],flush=True)
 except Exception as e:print(u,repr(e),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(get,paths))
print('IPROX_FRONTEND_LOOKUP_COMPLETE')
