from pathlib import Path
import requests,re,json
D=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening/A_CRLM_phosphoregulation/00_docs/iprox_lookup')
for f in ['MSV02201.js','MSV02202.js']:
 text=(D/f).read_text()
 for pat in ['findBySubProjectId','findFilesBySubProjectID','findByProjectId']:
  for m in re.finditer(pat,text):print(f,text[max(0,m.start()-180):m.start()+900],flush=True)
r=requests.get('https://www.iprox.cn/app/js/ProjectFileList.js',timeout=(15,45));r.raise_for_status();(D/'ProjectFileList.js').write_text(r.text)
for m in re.finditer('url:',r.text):print(r.text[max(0,m.start()-80):m.start()+650],flush=True)
for key in ['subProjectId','subprojectId']:
 u='https://www.iprox.cn/PMD009Controller/findBySubProjectId.jsonp'
 try:
  r=requests.get(u,params={key:'IPX0007391001'},timeout=(15,45));r.raise_for_status();(D/(key+'_info.json')).write_text(r.text);print(key,r.text[:3500],flush=True)
 except Exception as e:print(repr(e),flush=True)
print('IPROX_API_LOOKUP_COMPLETE')
