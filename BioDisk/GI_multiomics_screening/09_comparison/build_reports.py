from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,html,re,shutil,base64,datetime
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation';B=R/'B_ICC_lymphatic_invasion'
def rd(p):return pd.read_csv(p,sep='\t')
def fmt(v):
 if pd.isna(v):return 'NA'
 if isinstance(v,(float,np.floating)):return f'{v:.3g}'
 return str(v)
def table(d):
 d=d.copy();cols=list(d.columns);rows=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
 for vals in d.itertuples(index=False,name=None):rows.append('| '+' | '.join(fmt(v).replace('|','/') for v in vals)+' |')
 return '\n'.join(rows)
def inline(s):
 s=html.escape(s);s=re.sub(r'\*\*(.*?)\*\*',r'<strong>\1</strong>',s);s=re.sub(r'`([^`]+)`',r'<code>\1</code>',s)
 s=re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',r'<a href="\2">\1</a>',s);return s
def html_report(md,path):
 parts=[];lines=md.splitlines();i=0
 while i<len(lines):
  l=lines[i]
  if not l.strip():i+=1;continue
  if l.startswith('#'):
   n=len(l)-len(l.lstrip('#'));parts.append(f'<h{min(n,4)}>'+inline(l[n:].strip())+f'</h{min(n,4)}>');i+=1;continue
  if l.startswith('| '):
   group=[]
   while i<len(lines) and lines[i].startswith('| '):group.append(lines[i]);i+=1
   parts.append('<div class="table"><table>')
   for j,row in enumerate(group):
    if j==1:continue
    tag='th' if j==0 else 'td';parts.append('<tr>'+''.join(f'<{tag}>'+inline(x.strip())+f'</{tag}>' for x in row.strip('|').split('|'))+'</tr>')
   parts.append('</table></div>');continue
  if l.startswith('!['):
   m=re.match(r'!\[(.*?)\]\((.*?)\)',l)
   if m:
    p=(path.parent/m[2]).resolve();data=base64.b64encode(p.read_bytes()).decode();parts.append(f'<figure><img alt="{html.escape(m[1])}" src="data:image/png;base64,{data}"></figure>')
   i+=1;continue
  if l.startswith('- '):parts.append('<p class="item">• '+inline(l[2:])+'</p>');i+=1;continue
  parts.append('<p>'+inline(l)+'</p>');i+=1
 css='body{font-family:"Microsoft YaHei","Noto Sans CJK SC",sans-serif;max-width:1100px;margin:35px auto;padding:0 30px;color:#26343e;line-height:1.75}h1{font-size:28px;color:#164960}h2{margin-top:32px;color:#164960;border-bottom:1px solid #c7d6de}h3{color:#425d6b}table{border-collapse:collapse;width:100%;font-size:13px}th{background:#e8f0f5;text-align:left}td,th{padding:7px 9px;border-bottom:1px solid #d7e0e5;vertical-align:top}.table{overflow-x:auto}figure{margin:22px 0}img{max-width:100%;height:auto}code{font-size:12px;overflow-wrap:anywhere;background:#f1f4f6}a{color:#19668e}.item{margin:5px 0}@media print{body{margin:0;max-width:none;font-size:11pt}h2,h3{break-after:avoid}figure,table{break-inside:avoid}}'
 path.with_suffix('.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>GI公开多组学初筛</title><style>'+css+'</style><body>'+''.join(parts)+'</body></html>',encoding='utf-8')
def write_report(P,name,md):
 p=P/name;p.write_text(md,encoding='utf-8');html_report(md,p)

s=json.loads((R/'09_comparison/result_summary.json').read_text());lit=rd(R/'09_comparison/literature_review.tsv')
for label,P in [('A',A),('B',B)]:lit[lit.project==label].to_csv(P/'06_tables/literature_neighborhood_audit.tsv',sep='\t',index=False)
ap=rd(A/'06_tables/protein_paired_complete.tsv.gz');asp=rd(A/'06_tables/site_paired_complete.tsv.gz');cond=rd(A/'06_tables/site_conditional_diagnostic.tsv.gz');ev=rd(A/'06_tables/cross_cohort_site_evidence_complete.tsv.gz')
top=ev[(ev.family=='primary_diagnostic')].sort_values('FDR_conditional').drop_duplicates('gene').head(5)
top.to_csv(A/'06_tables/diagnostic_examples_not_recommended.tsv',sep='\t',index=False)
counts=[]
for label,d in [('总蛋白',ap),('Class I位点',asp)]:
 counts.append({'层级':label,'总条目':len(d),'≥20对':int((d.n>=20).sum()),'主覆盖率':f'{(d.n>=20).mean()*100:.1f}%','10–19对':int(d.n.between(10,19).sum()),'<10对':int((d.n<10).sum()),'主家族FDR<0.05':int(((d.family=='primary')&(d.FDR<.05)).sum())})
tp=top[['gene','site','n','effect_conditional','FDR_conditional','protein_effect','composition_FDR']].rename(columns={'gene':'基因','site':'位点','n':'有效配对','effect_conditional':'条件ΔZ','FDR_conditional':'诊断FDR','protein_effect':'蛋白Δlog2','composition_FDR':'成分模型FDR'})
sites=rd(A/'04_qc/site_mapping_exclusion_summary.tsv').rename(columns={'reason':'匹配记录','n_sites':'位点数'})
adoc=f'''# A：结直肠癌肝转移磷酸化调控公开初筛

结论日期：2026-09-22。服务器：sgs10 / osv01。项目位置：`{A}`。

## 当前不推荐立项；补齐明确计算证据后再评估

原发灶与肝转移灶存在稳定的患者内蛋白和标准化位点变化，但现有公开附件不能完成本题要求的“总蛋白量校正后磷酸化改变”。本轮推荐靶点为 **0个**。如仅决定哪个方向值得补计算证据，A优先于B；这不等于推荐客户进入湿实验立项。

决定性依据是 S2C 已逐样本Z标准化，而 S2B 仍为log2比值。99个磷酸化样本的均值为0、样本标准差为1，总蛋白样本标准差为0.351–0.812。直接相减会混合单位。本轮没有生成伪“校正后log2FC”，所有位点均在量校正状态表中明确登记为不可估计。

## 方法

### 数据与患者单位

发现研究为 [Advanced Science 原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC12822478/)。读取原始临床表S1A/S1B、患者级S2B/S2C及补充方法。s009.xls文件签名为ZIP，内部含 `xl/workbook.bin`，以pyxlsb读取。其他旧版OLE工作簿优先识别文件头，以xlrd读取；不把其内嵌主题ZIP误认成XLSX。保存未经重命名的源表与派生映射。

主比较固定为同一患者LM−T。总蛋白101列，包括N33、T34、LM34；位点99列，包括N/T/LM各33。蛋白保留34对，磷酸化保留33对。P34缺少N蛋白列及全部磷酸化列。正常结直肠仅用于背景核对，不把N−LM作为转移机制。

作者方法采用TMT10-plex，IRS在126通道，其余9通道为每组3份组织；位点定位概率阈值≥0.75。总蛋白作样本中位数归一化、IRS比值及log2转换；磷酸化以相应蛋白样本中位数归一化，作IRS/log2转换及批次处理后逐样本Z标准化。本轮使用附件实际数值，保留NA，不把缺失设0，不作分组填补。

### 统计与敏感性

主家族要求每条目≥20对有效观测；10–19对单独探索；不足10对只描述。以患者内差值拟合limma截距模型，使用trend经验贝叶斯方差，报告实际单位效应、95%CI、双侧P及分家族BH-FDR。未收缩配对t检验与留一患者效应作为稳定性记录。

严格按accession、基因和残基编号匹配总蛋白。仅在匹配唯一时，额外拟合“患者内位点ΔZ ~ 患者内蛋白Δlog2”的条件模型；其截距为蛋白差值为0时的标准化位点差，列为诊断性关联。量校正目标与诊断模型分别保存。

通路采用预先固定的MSigDB Hallmark 2025.1.Hs；至少10个可测基因，患者分数至少覆盖50%的可测成员。先在基因内汇总位点，避免多位点基因过度加权。激酶仅采用KSEAapp内人类PhosphoSitePlus实验支持映射，精确匹配accession与位点，至少5个位点和3个底物基因；按底物基因均等汇总后在患者层面检验。通路、激酶与蛋白/位点分别BH校正。

固定肝细胞、血液、免疫、间质、上皮蛋白代理集，完成患者内成分变化及条件模型敏感性。另检查样本中位数中心化口径。主分析不按结果剔除患者。

外部蛋白资源为 [2025年三中心CRLM研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC12335997/)。152份标本来自111名患者：UKSH42/42、UKE50/18、Sydney60/51（标本/患者）。在患者内取多病灶中位数后分析预定义蛋白程序，逐字段登记冲突；复发字段保留为采样相关复发记录，不作随访复发终点。外部连续程序与OS采用Cox模型，提供每SD的HR、95%CI、P/FDR及比例风险检验。

## 结果

### 蛋白配对变化充分，位点主覆盖率仅24.4%

{table(pd.DataFrame(counts))}

完整统计同时保留正向、负向、未显著及不可检验条目。≥20对门槛使蛋白保留76.4%，位点保留24.4%；不能以16,300个位点总数代表主比较具有足够配对的信息量。

![A01](../07_figures/A01_scale_coverage.png)

图A01：Prot=总蛋白；Phos=磷酸化；N=正常结直肠；T=原发灶；LM=肝转移。每点为一份实际样本，右图标明各配对数区间的条目数。SD=标准差。

![A02](../07_figures/A02_paired_effects.png)

图A02：仅显示主家族；蛋白有效配对20–34，位点20–33。Δlog2与ΔZ是不同单位。虚线为FDR=0.05；位点图不表示总蛋白校正后效应。

### 肝细胞相关成分差异是核心解释因素

肝细胞代理分数LM−T为1.059，95%CI 0.843–1.275，P=1.24×10⁻¹³，FDR=6.22×10⁻¹³，31/34名患者升高。APOE、PLIN2、PON1、ADH4等位于蛋白变化前列，蛋白通路包括异生物代谢、胆汁酸代谢和凝血程序。该结果确认组织背景与肝相关程序发生变化，不能直接归因于癌细胞内的转移驱动机制。

蛋白Hallmark有11/49条达到FDR<0.05；标准化位点Hallmark有6/32条达到FDR<0.05。完整结果见相应通路表。成分模型与原主集使用相同患者；其条件效应与未校正效应并列保存。

### 精确位点匹配与诊断关联仍未形成量校正证据链

{table(sites)}

7,130个位点具有唯一精确匹配，其中1,753个位点满足条件模型≥20对的主诊断门槛。按诊断FDR排序选择5个不同基因展示，目的仅为核查，不构成靶点推荐：

{table(tp)}

![A03](../07_figures/A03_patient_examples.png)

图A03：各连线连接同一患者。APOE为蛋白log2比值，POF1B S413与FBN1 S2702为源位点Z值；例图展示确实存在的配对方向，不证明总蛋白量独立调控。POF1B与POU5F1B、FBN1与FN1为不同基因。

### CHEK1底物分数降低，但不能推出干预靶点

9个激酶通过覆盖门槛。CHEK1映射5个位点、3个底物基因；标准化底物分数ΔZ=−0.259，95%CI −0.394至−0.124，P=0.000452，FDR=0.00407；26/33名患者下降，3个底物基因方向均下降。其余8个激酶未通过FDR。CHEK1没有可用的唯一对应总蛋白结果，也没有独立位点复现；不据此提出CHEK1抑制策略。

![A04](../07_figures/A04_composition_kinase.png)

图A04：Hep=肝细胞蛋白代理，非病理肿瘤纯度；subs=底物。CHEK1各点为一名患者的底物均值差，横轴为患者排列，非时间。

### 外部队列尚未补齐本题的独立验证

三中心外部队列只有LM，能够检验蛋白程序的患者异质性和临床关联，不能检验T−LM差异或位点。23名患者至少一个字段存在跨标本冲突，10名UKE患者的生存时间不一致。结合缺失和冲突处理，OS分析实际为Sydney的51名患者、21个死亡事件；完整预定义通路OS家族未达到FDR<0.05。该结果仅对应可用生存子集。

[2020年Cancer Cell研究](https://doi.org/10.1016/j.ccell.2020.08.002) 的8个公开附件已定位，取得146列蛋白和146列磷酸化定量表，以及患者临床映射；公开附件未给出本轮所需的LM患者级量化列。NODE OEP000729访问返回403。作者差异统计已登记，但未替代配对矩阵重分析。

### 公开原始文件不能在本轮轻量范围恢复关键输入

iProX IPX0007391001实际公开清单450项：444个RAW、5个参数XML、1个FASTA；未列出proteinGroups或Phospho(STY)Sites等未Z标准化定量表。参数XML网页表示可解析出TMT实验批次编号，但没有P患者—组织—通道对应标签。参数表示与登记SHA1不完全一致，不能当作字节校验通过的原始上传文件。注册原始RAW不等于量化结果已到位。

## 文献邻近性核查

FBN1的CAF相关性、CLDN3的上皮屏障/STAT3/Wnt机制、C5AR1的CRC功能及PLEC的骨架作用均有既有研究；这些不能验证本轮相应位点。CLDN3不同研究的蛋白方向存在矛盾，已保留。POF1B的上皮黏附功能有文献背景，但本轮没有S413的独立功能证明。完整检索词、原始研究链接与逐项判断见 `literature_neighborhood_audit.tsv`。

SHMT1–甲酸–AMPK、PIM–NDRG1 S330、EIF4B S422，以及原文已有标志物不列为新发现。作者S8A的DepMap均值仅作已发表CRC依赖背景；未在该筛选名单中不等于CRISPR阴性，也不是位点功能证明。

## 讨论与决策边界

当前数据确立了患者内蛋白和标准化位点差异；关键缺口在于其能否形成总蛋白量独立、细胞来源明确、独立复现的磷酸化调控机制。样本Z变换缺少逆变换参数、患者到TMT通道的完整映射未到位、独立配对位点矩阵未取得，这三项阻止本轮靶点推荐。成分代理不能代替病理纯度或已分析的单细胞定位。

本轮判定：**补齐明确计算证据后再评估**。需要的具体计算输入是未经逐样本Z变换的蛋白/位点IRS量化、对应批次与患者通道映射，以及可独立重分析的T−LM蛋白/位点矩阵；然后按本轮固定标准完成量校正、成分敏感性和复现。此次不开展大规模RAW重搜库，不建议客户用实验判断弱假说是否成立。没有进入实验承接阶段的推荐候选。

## 交付与复现

原始数据：`{A}/02_raw`；清洗矩阵：`{A}/03_processed`；完整汇总统计：`{A}/06_tables`；脚本与日志：`{A}/scripts`、`{A}/logs`。轻量包不包含原矩阵及RDS。表中NA按具体状态解释，不能改为0。数据清单、校验、方案冻结、尺度修订和统计家族均随包交付。
'''
write_report(A,'08_reports/A_screening_report.md',adoc)

dp=rd(B/'06_tables/discovery_protein_complete.tsv.gz');vp=rd(B/'06_tables/validation_protein_complete.tsv.gz');path=rd(B/'06_tables/pathway_cross_cohort_complete.tsv');p=path[path.selected_D][['feature','effect_D','FDR_D','effect_V','P_V','V_FDR_selected']].copy();p.feature=p.feature.str.replace('HALLMARK_','',regex=False);p.columns=['固定通路','发现效应','发现FDR','验证效应','验证P','验证候选家族FDR']
clin=rd(B/'06_tables/discovery_clinical_sensitivity.tsv.gz');z=clin[(clin.family=='primary')&(clin.FDR<.05)].merge(vp[['feature','effect','P','FDR']],on='feature',how='left',suffixes=('_D_clinical','_V_primary'));z['decision']='not_primary_discovery_candidate';z.to_csv(B/'06_tables/clinical_sensitivity_signal_audit.tsv',sep='\t',index=False)
bdoc=f'''# B：肝内胆管癌淋巴管侵犯公开初筛

结论日期：2026-09-22。服务器：sgs10 / osv01。项目位置：`{B}`。

## 当前数据不支持稳定独立复现，淘汰本轮方向

发现队列40例L+与40例L−的主蛋白筛查没有FDR<0.05条目；独立验证队列19例L+与41例L−同样没有。发现端4个显著Hallmark通路在独立验证中**全部方向相反**，且同一可测基因集合的复核结论一致。本轮推荐靶点 **0个**，不追加模型或用实验救回假说。

## 方法

### 数据与主终点

来源为 [Nature Communications原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC13013968/)，使用公开MOESM3/4/7/11，不使用EGA。工作表重新检测维度，临床表实际分别149和123份组织记录；矩阵分别10,443×149、8,526×123。每队列单独处理，使用作者已完成log2与中位数中心化的DIA-NN矩阵，不再次log转换或填补。

发现队列肿瘤80、癌旁69；独立验证肿瘤62、癌旁61。发现Lymphangioinvasion为Yes40/No40；验证Lymphangiosis为L1 19/L0 41/未知2，未知者不进入主表型比较。血管侵犯V1 15/V0 45/未知2，神经侵犯Pn1 20/Pn0 40/未知2，均与L分别保存。所有推断以独立患者为单位。

### 预先固定的统计与验证

蛋白主检验要求每组至少10个有效观测且全部肿瘤样本≥70%有值，分别以limma拟合L+−L−，trend经验贝叶斯方差，报告log2效应、95%CI、P和BH-FDR；低覆盖但每组≥10的条目为独立探索家族，不进入主筛选。主筛选阈值FDR<0.05，没有按结果变更。

发现候选表先锁定，验证仅检验固定候选。本轮主蛋白候选为空，另保留验证全蛋白组完整统计。Hallmark固定50集合，按预设可测覆盖构造患者分数；发现端通过的4集合组成独立验证候选家族，验证不参与基因集选择。评分仅用于预定义程序描述，不是预后风险模型。

临床敏感性包含性别、年龄组顺序和发现导管亚型/验证分级；成分敏感性包含肝细胞、血液、免疫、间质和上皮代理。未机械纳入包含侵犯定义或其下游结果的综合分期。对共同可测基因评分同时保存原主集、固定患者原评分和相同患者新评分三种比较。

发现队列TTR使用源 `daystoevent` 和已核实的事件编码：`censored=1`为复发事件，与56例已记录复发时间一一对应，`0`为24例未记录复发；作者代码直接以该值作为Surv事件。两例患者源事件时间2213天短于末次随访，保留源TTR删失时间并单列审计。连续程序TTR关联采用Cox模型，每SD报告HR、CI、P/FDR及比例风险检验。

## 结果

### 两个独立蛋白队列均未产生主筛选候选

| 队列 | 实际主终点患者 | 矩阵蛋白总数 | 主检验蛋白 | 最小主FDR | 主FDR<0.05 |
| --- | --- | --- | --- | --- | --- |
| 发现 | L+40 / L−40 | 10,443 | 6,944 | 0.0572 | 0 |
| 验证 | L+19 / L−41 | 8,526 | 4,321 | 0.466 | 0 |

![B01](../07_figures/B01_protein_screen.png)

图B01：Disc=发现，Valid=验证，L=淋巴管侵犯。每个蛋白使用实际有效观测；标题为队列上限样本量，具体每组n见完整表。虚线FDR=0.05。

### 四个发现通路在独立队列全部反向

{table(p)}

发现端E2F/G2M升高、Angio/Coag降低；验证端四者相反。所有预选通路的验证FDR≥0.346。共测基因敏感性保留相同80/60名患者，方向仍然相反，排除了“仅因两个队列评分基因数量不同”解释该现象。

![B02](../07_figures/B02_pathway_replication.png)

图B02：E2F=E2F targets，G2M=G2M checkpoint，Angio=angiogenesis，Coag=coagulation。横轴为患者评分L+−L−，误差线95%CI；全部为预先锁定的发现集合。

![B03](../07_figures/B03_patient_pathways.png)

图B03：每点一名患者，黑线为中位数；显示E2F和Angio的患者分布及实际分组n，不将差异推断成侵犯机制。

### 临床调整不能补回独立复现

临床协变量模型在发现端出现12个FDR<0.05的蛋白，包含SMC2、SMC4和MCM家族；这些均未通过固定主筛选，未提升为候选。成分敏感性未得到FDR<0.05蛋白，验证临床与成分模型也未得到主家族显著蛋白。12项临床敏感性信号及对应验证结果完整保留，避免选择性展示主模型阴性。

发现队列最少观测样本KB107仅519个蛋白，属于癌旁，不进入L主比较；对应肿瘤KB34有5,580个蛋白，保留在主集。没有因观察到差异结果而更换患者。

![B04](../07_figures/B04_effect_qc.png)

图B04：左图4,289个两个主家族共同蛋白的独立效应；右图每点为一份肿瘤样本。D=发现，V=验证，Prot=总蛋白；分组n分别40/40和41/19。

### 发现端复发关联成立，但不是独立侵犯验证

80名患者、56个复发事件中，Angio分数每SD的TTR HR=0.522（95%CI 0.394–0.693，FDR=2.65×10⁻⁵），Coag HR=0.648（0.481–0.874，FDR=0.00887），E2F HR=1.363（1.052–1.766，FDR=0.0254）；G2M FDR=0.0557。该关联对应发现队列的复发终点，不能替代已经失败的独立L表型验证，也不能与验证队列OS混为同一终点。

### 转录与单细胞资源已核对，未作为救回证据

GSE107943计数矩阵已下载并通过校验，57,773个条目，30肿瘤、27正常，27对配对患者。GEO标签包含复发、无病时间（月）、死亡、OS时间（月）、血管侵犯及病理类型，未提供淋巴管侵犯。不能把血管侵犯或肿瘤表达当L表型验证。由于主蛋白候选为0，未追加候选转录差异或筛选新的临床终点。

GSE138709公开样本登记包括5份肿瘤文库，其中ICC_24_Tumor1/2具有相同患者标签；登记可识别4个肿瘤患者ID和3个癌旁患者ID，总计5个患者ID。文库数不能直接算5名独立肿瘤患者。尚未分析其单细胞表达矩阵及细胞注释，未作细胞来源已验证的表述。GSE26566不纳入本轮后续分析。

## 文献与创新性

原研究已报道ECM/增殖亚型、其侵犯及生存关联，并完成EIF4A1靶点验证。本轮E2F/G2M与部分基质相关程序与既有主线邻近，且独立侵犯方向不复现。未重做原分类器、未将EIF4A1或临床调整后的SMC2包装成新靶点。文献核查表与全部阴性统计随包交付。

## 讨论与最终判定

两队列临床L字段在语义上均指病理淋巴管侵犯，V/Pn/N分开处理；当前附件没有提供两中心逐病例统一复阅或染色判定规则，不能宣称病理判读已完全标准化。两队列仪器、覆盖度、组织取材不同，不能据观察到的反向结果断言统一生物学效应不存在；但按本轮“公开数据承担独立复现”的选题标准，当前证据不支持推进。

验证队列futime的单位由论文确认为月，但公开临床表与取得的代码未找到可独立核实fustat 0/1的明确数据字典。本轮不根据字段名称猜测事件编码，故未提交验证OS推断；原始值和交叉表完整保留。它不改变独立L表型未通过的结论。

本轮判定：**淘汰**。不通过更换阈值、分组、通路或预后模型追加阳性；不建议客户先做实验决定本假说是否成立。本轮没有可以让湿实验承接的稳定候选机制。

## 交付与复现

原始数据、处理矩阵和内部检查点分别位于 `{B}/02_raw`、`{B}/03_processed`、`{B}/05_analysis`。轻量包包括完整统计、患者元数据、图和脚本，不包括原始矩阵、RDS及单细胞大表。运行方法见 `scripts/RUNBOOK.md`。
'''
write_report(B,'08_reports/B_screening_report.md',bdoc)

comparison=f'''# 消化系统肿瘤公开多组学初筛对比

完成日期：2026-09-22。服务器：sgs10（osv01）。根目录：`{R}`。

## 本轮没有可推荐立项的候选

**A列为“补齐明确计算证据后再评估”，B淘汰。推荐靶点总数0。** 若只决定补计算证据的次序，优先A；当前两方向均不建议进入客户湿实验立项。

| 证据项目 | A：CRLM磷酸化 | B：ICC淋巴管侵犯 |
| --- | --- | --- |
| 主设计 | 34对蛋白、33对位点的患者内LM−T | 发现40/40；独立验证19/41 |
| 主蛋白筛查 | 6,181项；1,434项FDR<0.05 | 发现6,944项、验证4,321项；均0项 |
| 主位点筛查 | 3,972项；480项FDR<0.05，单位ΔZ | 不适用 |
| 核心证据障碍 | 位点逐样本Z值与总蛋白log2单位不兼容；不能量校正 | 4个发现通路独立队列全部反向 |
| 混杂与覆盖 | 肝细胞代理31/34升高；位点主覆盖24.4%；TMT患者通道映射未到位 | 临床与成分敏感性不能补回独立复现；共同基因评分仍反向 |
| 独立验证范围 | 111患者外部LM蛋白，不能验证配对差异或位点；可用OS仅Sydney51患者 | 独立蛋白验证已完成并未通过 |
| 决策 | 补齐明确计算证据后再评估 | 淘汰 |
| 推荐靶点 | 0 | 0 |

## 生信已完成的证据链

A已完成原文与格式核对、源文件校验、患者和组学交集、缺失/尺度审计、蛋白与标准化位点配对差异、各自通路、9个合格激酶的患者级多底物检验、精确位点映射和条件诊断、成分及中心化敏感性、外部多标本患者汇总与可用OS分析、功能背景和文献邻近性核查。原始量校正和独立配对位点验证明确未完成，未借诊断模型替代。

B已完成两队列独立主比较、完整蛋白及通路统计、4个锁定通路的独立验证、临床与成分敏感性、共同基因评分的三种患者固定比较、发现TTR关联，以及公开转录/单细胞标签与终点核对。四通路发现关联是真实计算结果；独立反向同样完整保留，不用复发关联掩盖验证失败。

全部主效应、关键未收缩检验、BH家族和源附件经过独立Python/SciPy复核，160项检查通过；R使用R_4.5.2、Python使用batch_bench。未升级既有环境、未修改胰腺癌项目、未使用EGA、未联系作者或提交受控申请。

## A进一步评估所需的明确计算输入

- 未逐样本Z变换的蛋白/位点IRS量化，以及对应的患者—组织—TMT通道/批次映射；现有公开清单仅有RAW/参数/FASTA，参数未含患者标签。
- 独立T−LM患者级蛋白及位点矩阵；2020附件取得的是146列量化表，所需LM列未到位，NODE访问403。
- 基于上述输入完成固定的蛋白量校正、异构体/位点匹配及独立复现，再判断细胞来源与机制。

这些是再次计算评估的前置条件，不是向客户推荐实验去“救回”假说。本轮没有推荐候选，故不提供泛化实验方案。

## 交付阅读顺序

先读本对比报告，再读两个项目 `08_reports` 的中文报告。完整阳性/阴性统计在各自 `06_tables`，患者与源标签在 `01_metadata`，关键审计在 `04_qc`，脚本和运行命令在 `scripts`。量校正不可估计状态单独成表，不会被当作无显著结果或成功完成。

完整服务器结果保留在 `{R}`。ZIP不含原始矩阵、细胞矩阵、RDS/H5AD或缓存；完整文件清单、哈希和服务器路径随包提供。
'''
write_report(R,'09_comparison/GI_comparison_report.md',comparison)

dictionary='''# 汇总表字段与证据层级

- feature：源条目唯一编号；A的PR/PH编号对应各自完整特征映射，B以源蛋白accession为编号。
- n / n_pos / n_neg：A的n为有效患者配对数，n_pos/n_neg为差值正/负人数；B的n_pos/n_neg分别为L+/L−有效观测人数。不得跨设计混用含义。
- effect：实际模型系数。A protein为LM−T的log2比值变化；A site为LM−T的样本Z变化；B protein为L+−L− log2表达差；通路为基因Z均值评分差。
- CI_low / CI_high：95%置信区间。主模型使用limma后验方差；条件诊断与未收缩敏感性采用相应普通模型。
- P / FDR：双侧P与指定家族BH-FDR。A protein、site、conditional、pathway、kinase分别校正；主家族与10–19对探索分开。B两个队列与蛋白/通路分开；锁定候选验证只在相应候选家族内校正。
- *_conditional：标准化位点对蛋白变化条件关联，不是用户目标的磷酸化量校正；不可报告为占有率或校正后的log2FC。
- P_unmoderated / P_Welch：未收缩检验敏感性；其名义P不替代主FDR。
- LOO_min_effect / LOO_max_effect：逐一移除患者后的平均差值范围，不是独立验证。
- fixed_patient_effect：协变量完整病例对应同一患者子集的原模型系数。共测评分triad表包含原患者/原评分、固定患者/原评分、固定患者/新评分。
- V_FDR_all：验证全蛋白或全通路家族；V_FDR_candidates / V_FDR_selected：冻结发现候选后的验证家族。发现样本内再检验属于选择后统计，不能称全局显著。
- HR：连续分数每SD的风险比；TTR与OS分表、时间单位明确。NA可能表示未测、样本不足、映射歧义或不可估计，结合status/reason解释。
- author_DepMap_mean：原研究S8A已筛选CRC依赖基因的平均分数。缺失不代表基因无依赖，不是独立CRISPR重分析或位点功能。
- family=descriptive：未做推断检验；primary_diagnostic与primary统计目标不同。recommended候选为0。
'''
for P in [A,B]:(P/'06_tables/TABLE_DICTIONARY.md').write_text(dictionary)

# Version all local analytical code, retaining failed attempt logs for audit.
script_files=[p for p in (R/'09_comparison').iterdir() if p.is_file() and p.suffix in ['.py','.R','.sh']]
runbook=f'''# 运行与阶段

根目录：{R}
Python：/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
Rscript：/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla
隔离读取依赖：{R}/09_comparison/vendor（pyxlsb 1.0.10、xlrd 2.0.2；未更改Conda）

按依赖顺序：launch_acquire.sh → run_extract.sh <项目名> → run_prepare.sh → run_screen.sh → run_extend.sh → run_verified_aux.sh → make_figures.py → build_reports.py → package_delivery.py。
在线资源若已完整，复用SHA256/格式校验后的原文件。原始尝试失败与本地备用下载记录并存；最终有效性以final_raw_integrity.tsv和source_registry_final.tsv为准。
主入口位于09_comparison；两个项目scripts目录保留同版本副本，脚本使用指定服务器绝对路径。

主分析前独立语法检查：Python -m py_compile；Rscript --vanilla -e "invisible(parse(file='script.R'))"（仅远程语法探针，分析均保存.R文件）。分析逻辑检查与实际运行命令见各run_*.sh、logs和execution_manifest.tsv。
完成标志：TABLE_EXTRACTION_COMPLETE、ANALYSIS_INPUT_PREPARATION_COMPLETE、MAIN_SCREEN_COMPLETE、EXTENDED_SCREEN_COMPLETE、KINASE_SCREEN_COMPLETE、STATISTICAL_VERIFICATION_COMPLETE、FIGURES_COMPLETE、REPORTS_COMPLETE、DELIVERY_COMPLETE。每个标志只对应其阶段。
未完成目标：A位点量校正及独立配对位点验证；B验证OS事件编码未充分核实。未触发阶段：因无稳定候选停止后不做单细胞矩阵重分析、候选转录推断或新增模型。没有将这些目标写成成功标志。

统计代码版本：GI-screen-v1-20260922；SHA256在script_versions.tsv。R会话信息及环境审计另附。仅重跑所需阶段，勿无故扩大分析。
'''
versions=[]
for p in script_files:
 versions.append({'script':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
 for P in [A,B]:shutil.copy2(p,P/'scripts'/p.name)
for P in [A,B]:(P/'scripts/RUNBOOK.md').write_text(runbook);pd.DataFrame(versions).to_csv(P/'scripts/script_versions.tsv',sep='\t',index=False)
logs=['preflight.json','R_environment.log','screen_sessionInfo.txt','statistical_integrity_checks.tsv','statistical_verification.json','statistics_final_check.log','prepare.log','screen.log','extend.log','kinase.log','verified_aux.log','acquisition.log','resource_followup.json','extra_manifest.json','metadata_finalize.log','geo_audit.log','figures.log']
for P in [A,B]:
 for name in logs:
  q=R/'09_comparison'/name
  if q.exists() and q.stat().st_size<2000000:shutil.copy2(q,P/'logs'/name)
# Final source registry distinguishes acquisition attempts from usable files.
for P in [A,B]:
 reg=rd(P/'01_metadata/final_raw_integrity.tsv');urls=[]
 for path in reg.file:
  n=Path(path).name
  if 'A_supp' in path: u='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12822478/supplementaryFiles'
  elif 'external_supp' in path:u='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12335997/supplementaryFiles'
  elif 'MOESM' in n:u='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-70817-1/MediaObjects/'+n
  elif n.startswith('Cell2020_'):u='https://ars.els-cdn.com/content/image/1-s2.0-S153561082030413X-'+n.replace('Cell2020_','')
  elif n.startswith('GSE107943'):u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107943/suppl/'+n
  else:u='see download manifest'
  urls.append(u)
 reg['source_url']=urls;reg['final_status']='available_integrity_verified';reg.to_csv(P/'01_metadata/source_registry_final.tsv',sep='\t',index=False)
ref=R/'09_comparison/reference';refs=[]
for p in ref.iterdir():
 if p.is_file():refs.append({'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
pd.DataFrame(refs).to_csv(R/'09_comparison/reference_integrity.tsv',sep='\t',index=False)
for P in [A,B]:shutil.copy2(R/'09_comparison/reference_integrity.tsv',P/'01_metadata/reference_integrity.tsv')
stages=[('core_data','COMPLETE','source_registry_final.tsv'),('A_paired_source_scales','COMPLETE','MAIN_SCREEN_COMPLETE'),('A_quantity_adjustment','NOT_ESTIMABLE','site_quantity_adjusted_status.tsv.gz'),('A_independent_paired_site','NOT_AVAILABLE','NODE403; no LM columns in obtained supplement'),('A_external_protein','COMPLETE','EXTENDED_SCREEN_COMPLETE'),('B_discovery_and_validation','COMPLETE_NO_REPLICATED_CANDIDATE','MAIN_SCREEN_COMPLETE'),('B_pathway_adjudication','COMPLETE_ALL_FOUR_REVERSED','EXTENDED_SCREEN_COMPLETE'),('B_single_cell','NOT_TRIGGERED_AFTER_STOP_RULE','no stable candidate'),('statistics_verification','COMPLETE','160/160 checks'),('figures','COMPLETE','8 PNG + 8 SVG; visual QA'),('reports','COMPLETE','REPORTS_COMPLETE')]
pd.DataFrame(stages,columns=['stage','status','evidence']).to_csv(R/'09_comparison/stage_status.tsv',sep='\t',index=False)
print('REPORTS_COMPLETE')
