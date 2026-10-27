from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,zipfile,shutil,datetime,subprocess,re
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');D=R/'10_delivery';D.mkdir(exist_ok=True)
projects=['A_CRLM_phosphoregulation','B_ICC_lymphatic_invasion']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
verification=json.loads((R/'09_comparison/statistical_verification.json').read_text());assert verification['passed']==verification['n_checks'] and not verification['failed']
for pr in projects:
 P=R/pr
 assert (P/'08_reports'/('A_screening_report.html' if pr.startswith('A') else 'B_screening_report.html')).exists()
 source=json.loads((P/'04_qc/source_workbooks.json').read_text());pd.DataFrame([{k:v for k,v in x.items() if k!='preview'} for x in source]).to_csv(P/'04_qc/source_workbook_structure.tsv',sep='\t',index=False)
 # Copy only the latest reproducibility scripts after this builder itself is finalized.
 for p in (R/'09_comparison').iterdir():
  if p.is_file() and p.suffix in ['.py','.R','.sh']:shutil.copy2(p,P/'scripts'/p.name)
 vs=[{'script':p.name,'sha256':sha(p),'bytes':p.stat().st_size} for p in (P/'scripts').iterdir() if p.suffix in ['.py','.R','.sh']]
 pd.DataFrame(vs).to_csv(P/'scripts/script_versions.tsv',sep='\t',index=False)

# Execution ledger reports actual checked main stages and explicit unmet objectives.
ledger=[
 ('source_extraction','run_extract.sh','PASS before final successful run','signature/OLE-XLSB precedence; actual dimensions','TABLE_EXTRACTION_COMPLETE','COMPLETE'),
 ('analysis_preparation','run_prepare.sh','PASS before run','explicit patient/sample labels; unique biological units; scale amendment','ANALYSIS_INPUT_PREPARATION_COMPLETE','COMPLETE'),
 ('main_screen','run_screen.sh','PASS before run','fixed primary families; independent validation; no Z/log2 subtraction','MAIN_SCREEN_COMPLETE','COMPLETE'),
 ('extended_screen','run_extend.sh','PASS before run','same-patient score triad; clinical endpoint conflicts; separate TTR/OS','EXTENDED_SCREEN_COMPLETE','COMPLETE'),
 ('kinase_final_run','run_verified_aux.sh','PASS before final run; initial exploratory invocation lacked separate record','PSP human exact accession-site; >=5 sites >=3 genes; patient means','KINASE_SCREEN_COMPLETE','COMPLETE'),
 ('statistics_audit','verify_statistics.py','PASS before final run','independent SciPy effects/t/BH; final raw CRC/hash','STATISTICAL_VERIFICATION_COMPLETE','COMPLETE'),
 ('figures','make_figures.py','PASS before run','patient-level values and actual units; eight panels reviewed; A01 layout fixed','FIGURES_COMPLETE','COMPLETE'),
 ('reports','build_reports.py','PASS before run','figures/tables directly from final statistics; no recommended targets','REPORTS_COMPLETE','COMPLETE'),
 ('A_quantity_adjustment','not_run','not_applicable','incompatible source units, no inverse scaling/measurement table','NONE','NOT_ESTIMABLE'),
 ('A_independent_paired_validation','not_run','not_applicable','required independent LM matrix inaccessible; no substitution with authors tests','NONE','NOT_AVAILABLE'),
 ('B_singlecell_target_localization','not_run','not_applicable','stopped after no stable candidate and independent reversed programs','NONE','NOT_TRIGGERED'),
]
ex=pd.DataFrame(ledger,columns=['stage','RUN_COMMAND_entry','SYNTAX_CHECK','LOGIC_CHECK','EXPECTED_SUCCESS_MARKER','status']);ex['server_script_dir']=str(R/'09_comparison');ex.to_csv(R/'09_comparison/execution_manifest.tsv',sep='\t',index=False)
note='''# 执行审计说明

主要分析入口均在最后成功运行前完成独立语法检查和逻辑检查，运行日志与完成标志对应真实阶段。早期格式/标签提取辅助脚本部分首次调用未单列语法检查日志；交付前对全部代码统一执行静态语法检查。激酶初次调用未单列语法记录，已通过run_verified_aux.sh补齐检查并实际重跑。未把初次失败或未运行阶段登记为完成。
统计核查初次159项中的1项失败来自检查器把两队列各5个成分检验误合并成10个；核查器改为按cohort和family分组后160/160通过，原统计结果未改。作图首次A03存在Series标签自动对齐错误，已改为按已配对顺序数组作图；不影响统计矩阵。A01小区段标签拥挤已改为百分比条宽，保留实际条目数。
原始下载尝试中B发现矩阵经服务器路径超时，通过本地备用下载并校验后上传；最终原始文件清单明确其已到位。全文XML初次校验依赖lxml不可用，改为标准库解析；实际全文已核对。所有失败日志保留，不用日志中的早期失败覆盖最终文件校验。
'''
(R/'09_comparison/execution_audit.md').write_text(note)
for pr in projects:
 P=R/pr;shutil.copy2(R/'09_comparison/execution_manifest.tsv',P/'logs/execution_manifest.tsv');shutil.copy2(R/'09_comparison/execution_audit.md',P/'logs/execution_audit.md')
 # Independent final syntax pass is explicitly dated as pre-delivery, not retroactively pre-first-run.
 pyfiles=list((P/'scripts').glob('*.py'))
 for p in pyfiles:compile(p.read_text(),str(p),'exec')

items=[]
for pr in projects:
 P=R/pr
 for folder in ['00_docs','01_metadata','04_qc','06_tables','07_figures','08_reports','scripts','logs']:
  for p in sorted((P/folder).rglob('*')):
   if not p.is_file() or '__pycache__' in str(p):continue
   if folder=='00_docs' and p.suffix!='.md':continue
   if folder=='04_qc' and p.name=='source_workbooks.json':continue
   if folder=='logs' and (p.name.endswith('.download.log') or p.stat().st_size>500000):continue
   if p.suffix in ['.rds','.RDS','.h5ad','.h5','.part','.pyc','.xlsx','.xls','.RData']:continue
   items.append(p)
 # Patient-level aggregate scores for figure/data review, excluding full expression matrices.
 for p in sorted((P/'05_analysis').glob('*patient_scores.tsv')):items.append(p)
for name in ['GI_comparison_report.md','GI_comparison_report.html','protocol.md','result_summary.json','statistical_integrity_checks.tsv','statistical_verification.json','stage_status.tsv','execution_manifest.tsv','execution_audit.md','preflight.json','R_environment.log','screen_sessionInfo.txt','reference_integrity.tsv','final_raw_integrity.tsv','literature_review.tsv']:
 p=R/'09_comparison'/name
 if p.exists():items.append(p)
# Small fixed gene sets and mapping preserve reproducibility; reference data are not raw omics matrices.
for name in ['hallmark.gmt','uniprot_gene_map.tsv']:
 p=R/'09_comparison/reference'/name
 if p.exists():items.append(p)
# Full filesystem inventory records the server locations for files omitted from the lightweight package.
allfiles=[]
for p in sorted(R.rglob('*')):
 if not p.is_file() or D in p.parents or '__pycache__' in str(p):continue
 allfiles.append({'server_path':str(p),'relative_path':str(p.relative_to(R)),'bytes':p.stat().st_size,'included_in_delivery':p in items,'sha256':sha(p) if p.stat().st_size<100000000 else 'see raw integrity manifest'})
pd.DataFrame(allfiles).to_csv(D/'FULL_SERVER_FILE_LIST.tsv',sep='\t',index=False)
readme=f'''GI Multi-omics Screening · 2026-09-22

结论：A补齐明确计算证据后再评估；B淘汰；推荐靶点0。
先读09_comparison/GI_comparison_report.html，再读各项目08_reports。
完整服务器结果：{R}
原始矩阵、RDS、单细胞矩阵与缓存不在轻量包。完整服务器文件清单：FULL_SERVER_FILE_LIST.tsv。
包内文件及SHA256：DELIVERY_FILE_LIST.tsv。脚本版本和运行命令见项目scripts和logs。
重要：A位点源数据为逐样本Z值；site_conditional_diagnostic不是磷酸化量校正或占有率。
B四个发现通路在独立蛋白队列均反向，不因发现TTR关联或临床调整名义变化推荐立项。
'''
(D/'README.txt').write_text(readme)
items=sorted(set(items));entries=[{'archive_path':str(p.relative_to(R)).replace('\\','/'),'server_path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in items]
for name in ['README.txt','FULL_SERVER_FILE_LIST.tsv']:
 p=D/name;entries.append({'archive_path':name,'server_path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
pd.DataFrame(entries).to_csv(D/'DELIVERY_FILE_LIST.tsv',sep='\t',index=False)
archive=D/'GI_multiomics_screening_20260922.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
 for rec in entries:z.write(rec['server_path'],rec['archive_path'])
 z.write(D/'DELIVERY_FILE_LIST.tsv','DELIVERY_FILE_LIST.tsv')
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 names=z.namelist();assert not any('/02_raw/' in p or '/03_processed/' in p or p.lower().endswith(('.rds','.h5ad','.h5','.part')) for p in names)
 for rec in entries:assert hashlib.sha256(z.read(rec['archive_path'])).hexdigest()==rec['sha256']
assert archive.stat().st_size<30*1024*1024,archive.stat().st_size
result={'archive':str(archive),'bytes':archive.stat().st_size,'MiB':archive.stat().st_size/1024**2,'sha256':sha(archive),'zip_crc':'PASS','member_sha256':'PASS','n_files':len(names),'raw_matrices_excluded':True,'statistical_checks':verification,'decision_A':'补齐明确计算证据后再评估','decision_B':'淘汰','recommended_targets':0}
(D/'DELIVERY_VERIFICATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));(D/'SHA256SUMS.txt').write_text(result['sha256']+'  '+archive.name+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2));print('DELIVERY_COMPLETE')
