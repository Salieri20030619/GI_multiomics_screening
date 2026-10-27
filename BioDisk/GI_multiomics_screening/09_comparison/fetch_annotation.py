from pathlib import Path
import requests,json,hashlib,concurrent.futures,xml.etree.ElementTree as ET,zipfile
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');D=R/'09_comparison/reference';D.mkdir(exist_ok=True)
tasks=[('hallmark.gmt','https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt'),('KSEAdata.RData','https://raw.githubusercontent.com/cran/KSEAapp/master/data/KSData.RData')]
def get(t):
 n,u=t;p=D/n;r={'file':str(p),'url':u}
 try:
  q=requests.get(u,timeout=(20,90));q.raise_for_status();p.write_bytes(q.content)
  r.update(status='downloaded',sha256=hashlib.sha256(q.content).hexdigest(),bytes=len(q.content))
 except Exception as e:r.update(status='failed',error=str(e))
 print(r,flush=True);return r
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:res=list(ex.map(get,tasks))
(D/'manifest.json').write_text(json.dumps(res,indent=2))
for p in R.glob('A*/02_raw/A_supp/*.docx'):
 with zipfile.ZipFile(p) as z:
  e=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
  paras=[''.join(t.itertext()) for t in e.findall('.//w:t',ns)]
  (R/'A_CRLM_phosphoregulation/00_docs/supp_methods.txt').write_text('\n'.join(paras))
print('REFERENCE_ACQUISITION_ATTEMPTS_COMPLETE')
