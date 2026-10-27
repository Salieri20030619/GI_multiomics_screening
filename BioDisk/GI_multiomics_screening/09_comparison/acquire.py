from pathlib import Path
import subprocess, concurrent.futures, json, hashlib, zipfile, tarfile, gzip, datetime, os
ROOT=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening')
A=ROOT/'A_CRLM_phosphoregulation';B=ROOT/'B_ICC_lymphatic_invasion'
resources=[
 (A,'00_docs/A_full.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12822478/fullTextXML'),
 (A,'02_raw/A_supp.zip','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12822478/supplementaryFiles'),
 (A,'00_docs/external_full.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12335997/fullTextXML'),
 (A,'02_raw/external_supp.zip','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12335997/supplementaryFiles'),
 (B,'00_docs/B_full.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13013968/fullTextXML'),
]
prefix='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-70817-1/MediaObjects/'
for i in [1,3,4,7,11]:
 name=f'41467_2026_70817_MOESM{i}_ESM.'+('pdf' if i==1 else 'xlsx')
 resources.append((B,('00_docs/' if i==1 else '02_raw/')+name,prefix+name))
for g in ['GSE107943','GSE138709','GSE26566']:
 resources.append((B,f'00_docs/{g}.soft.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}&targ=all&form=text&view=full'))
resources.append((B,'02_raw/GSE107943_rawread.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107943/suppl/GSE107943_rawread.txt.gz'))

def check(p):
 if p.stat().st_size<100:raise ValueError('too small')
 if p.suffix in ['.zip','.xlsx']:
  with zipfile.ZipFile(p) as z:
   bad=z.testzip()
   if bad:raise ValueError('CRC fail '+bad)
 elif p.suffix=='.gz':
  with gzip.open(p,'rb') as f:
   while f.read(1024*1024):pass
 elif p.suffix=='.xml':
  from lxml import etree
  tree=etree.parse(str(p))
  if tree.getroot().tag!='article':raise ValueError('not article XML')
 elif p.suffix=='.pdf' and p.read_bytes()[:4]!=b'%PDF':raise ValueError('not PDF')
 return hashlib.sha256(p.read_bytes()).hexdigest()

def get(r):
 project,name,url=r;p=project/name;p.parent.mkdir(exist_ok=True,parents=True)
 rec={'file':str(p),'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  if p.exists():
   rec['sha256']=check(p);rec['status']='reused_verified'
  else:
   part=p.with_name(p.name+'.part')
   cmd=['curl','-fL','--retry','2','--retry-delay','4','--connect-timeout','30','--max-time','420','-C','-','-o',str(part),url]
   proc=subprocess.run(cmd,capture_output=True,text=True)
   rec['command']=cmd;rec['exit_status']=proc.returncode
   (project/'logs'/ (p.name+'.download.log')).write_text(proc.stderr)
   if proc.returncode:raise RuntimeError(f'curl exit {proc.returncode}: '+proc.stderr[-600:])
   part.rename(p)
   rec['sha256']=check(p);rec['status']='downloaded_verified'
  rec['bytes']=p.stat().st_size
  if p.suffix=='.zip':
   dest=p.parent/p.stem;dest.mkdir(exist_ok=True)
   with zipfile.ZipFile(p) as z:
    for n in z.namelist():
     if not (dest/n).resolve().is_relative_to(dest.resolve()):raise ValueError('unsafe archive path')
    z.extractall(dest)
 except Exception as e:rec['status']='failed';rec['error']=str(e)
 print(json.dumps(rec),flush=True)
 (project/'logs'/(p.name+'.download.json')).write_text(json.dumps(rec,indent=2))
 return rec
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:result=list(ex.map(get,resources))
 for project in [A,B]:
  selected=[r for r in result if r['file'].startswith(str(project))]
  (project/'01_metadata/download_manifest.json').write_text(json.dumps(selected,indent=2))
 print('ACQUISITION_ATTEMPTS_COMPLETE',flush=True)
