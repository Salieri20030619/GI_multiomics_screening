from pathlib import Path
import sys, importlib, json, subprocess, datetime
root=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening')
for project in ['A_CRLM_phosphoregulation','B_ICC_lymphatic_invasion']:
    for d in ['00_docs','01_metadata','02_raw','03_processed','04_qc','05_analysis','06_tables','07_figures','08_reports','scripts','logs']:
        (root/project/d).mkdir(parents=True,exist_ok=True)
for d in ['09_comparison','10_delivery']:(root/d).mkdir(exist_ok=True)
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.executable,'version':sys.version,'packages':{}}
for p in ['numpy','pandas','scipy','statsmodels','matplotlib','seaborn','openpyxl','pyxlsb','xlrd','requests','lxml','scanpy','anndata','lifelines','pypdf']:
    try:
        m=importlib.import_module(p);out['packages'][p]=getattr(m,'__version__','installed')
    except Exception as e:out['packages'][p]=str(e)
for cmd in [['df','-h',str(root)],['free','-h'],['lscpu']]:
    out[' '.join(cmd)]=subprocess.run(cmd,capture_output=True,text=True).stdout
(root/'09_comparison/preflight.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2));print('PREFLIGHT_COMPLETE')
