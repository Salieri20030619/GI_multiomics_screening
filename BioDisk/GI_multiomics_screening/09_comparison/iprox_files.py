from pathlib import Path
import requests,json
D=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening/A_CRLM_phosphoregulation/00_docs/iprox_lookup')
S=requests.Session()
tasks=[('subproject','POST','PMD009Controller/findBySubProjectId.jsonp',{'subProjectId':'IPX0007391001','language':'en'}),('files','GET','PMD009Controller/findFilesBySubProjectID.jsonp',{'subProjectId':'IPX0007391001','pageNum':1,'pageSize':10000}),('tree','GET','projectFileList/selectProjectTree.jsonp',{'projectId':'IPX0007391000','url':''})]
for n,method,u,args in tasks:
 try:
  r=S.request(method,'https://www.iprox.cn/'+u,data=args if method=='POST' else None,params=args if method=='GET' else None,timeout=(15,60));r.raise_for_status();(D/(n+'.json')).write_text(r.text)
  print(n,r.status_code,len(r.content),r.text[:2200],flush=True)
 except Exception as e:print(n,repr(e),flush=True)
print('IPROX_FILE_METADATA_ATTEMPTS_COMPLETE')
