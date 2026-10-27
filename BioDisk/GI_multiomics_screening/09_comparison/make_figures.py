from pathlib import Path
import pandas as pd,numpy as np,matplotlib,sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening');A=R/'A_CRLM_phosphoregulation';B=R/'B_ICC_lymphatic_invasion'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','savefig.dpi':180})
C={'T':'#2878A6','LM':'#CE6248','D':'#2878A6','V':'#CE6248','N':'#9AABB6'};rng=np.random.default_rng(20260922)
def rd(p):return pd.read_csv(p,sep='\t')
def save(fig,P,name):
 fig.savefig(P/'07_figures'/f'{name}.png',bbox_inches='tight');fig.savefig(P/'07_figures'/f'{name}.svg',bbox_inches='tight');plt.close(fig)
def pair(ax,x,y,unit,title):
 x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float);good=np.isfinite(x)&np.isfinite(y);x=x[good];y=y[good]
 for a,b in zip(x,y):ax.plot([0,1],[a,b],color='#CDD5DC',lw=.65,zorder=1)
 ax.scatter(np.zeros(len(x)),x,color=C['T'],s=17,zorder=2);ax.scatter(np.ones(len(y)),y,color=C['LM'],s=17,zorder=2)
 ax.set(xticks=[0,1],xticklabels=['T','LM'],ylabel=unit,title=f'{title}\nn={len(x)} pairs',xlim=(-.35,1.35))
def dots(ax,v,g,title,unit='Score'):
 for i,label in enumerate([0,1]):
  vals=np.array(v)[np.array(g)==label];vals=vals[np.isfinite(vals)]
  ax.scatter(i+rng.normal(0,.065,len(vals)),vals,s=18,color=[C['T'],C['LM']][i],alpha=.8)
  if len(vals):ax.plot([i-.17,i+.17],[np.median(vals)]*2,color='black',lw=2)
 ax.set(xticks=[0,1],xticklabels=[f'L−\nn={sum((np.array(g)==0)&np.isfinite(v))}',f'L+\nn={sum((np.array(g)==1)&np.isfinite(v))}'],ylabel=unit,title=title,xlim=(-.4,1.4))
def volcano(ax,d,title,unit):
 d=d[d.family=='primary'].dropna(subset=['effect','FDR']);sig=d.FDR<.05
 ax.scatter(d.effect,-np.log10(d.FDR.clip(lower=1e-300)),s=6,c=np.where(sig,C['LM'],'#AAB5BD'),alpha=.65,rasterized=True)
 ax.axhline(-np.log10(.05),ls='--',color='#4A5863',lw=.8);ax.axvline(0,color='#D0D5DA',lw=.7)
 ax.set(xlabel=unit,ylabel='−log10 FDR',title=title+'\n'+f'{len(d):,} tests; {sig.sum():,} hits')

fig,axs=plt.subplots(1,3,figsize=(11,3.5),layout='constrained')
for name,label in [('protein','Prot'),('site','Phos')]:
 q=rd(A/f'04_qc/{name}_sample_qc.tsv');j=0 if name=='protein' else 1
 for k,tissue in enumerate(['N','T','LM']):
  z=q[q.tissue==tissue];axs[j].scatter(np.arange(len(z)),z.sd,s=14,color=C[tissue],label=f'{tissue} n={len(z)}')
 axs[j].set(title=f'{label}: sample SD',xlabel='Sample rank',ylabel='SD');axs[j].legend(frameon=False,fontsize=8)
stats=[rd(A/'06_tables/protein_paired_complete.tsv.gz'),rd(A/'06_tables/site_paired_complete.tsv.gz')]
for i,d in enumerate(stats):
 n=np.array([(d.n>=20).sum(),d.n.between(10,19).sum(),(d.n<10).sum()]);start=0
 for val,co in zip(n,['#2878A6','#E4B15D','#B9C3CB']):
  width=val/n.sum()*100;axs[2].barh(i,width,left=start,color=co);axs[2].text(start+width/2,i,str(val),ha='center',va='center',fontsize=8);start+=width
axs[2].set(yticks=[0,1],yticklabels=['Prot','Phos'],xlabel='Features (%)',title='Pair count bins',xlim=(0,100));axs[2].legend(['≥20','10–19','<10'],frameon=False,fontsize=8,loc='center',ncol=3)
save(fig,A,'A01_scale_coverage')
if '--only-a01' in sys.argv:print('A01_VISUAL_FIX_COMPLETE');raise SystemExit(0)
fig,axs=plt.subplots(1,2,figsize=(9,3.8),layout='constrained');volcano(axs[0],stats[0],'Prot: ≤34 pairs','Δ log2');volcano(axs[1],stats[1],'Phos: ≤33 pairs','Δ Z');save(fig,A,'A02_paired_effects')
px=rd(A/'03_processed/protein_log2.tsv.gz').set_index('feature');sx=rd(A/'03_processed/site_z.tsv.gz').set_index('feature');pm=rd(A/'01_metadata/protein_features.tsv');sm=rd(A/'01_metadata/site_protein_mapping.tsv')
fig,axs=plt.subplots(1,3,figsize=(10,4),layout='constrained')
for ax,name,tab,ann,col,unit in [(axs[0],'APOE',px,pm,'gene','log2 ratio'),(axs[1],'Q8WVV4-1_S413',sx,sm,'site','Z score'),(axs[2],'P35555_S2702',sx,sm,'site','Z score')]:
 f=ann.loc[ann[col]==name,'feature'].iloc[0];pats=sorted(set(c.split('_')[0] for c in tab.columns if c.endswith('_T'))&set(c.split('_')[0] for c in tab.columns if c.endswith('_LM')),key=lambda s:int(s[1:]));title=name if col=='gene' else ('POF1B S413' if 'Q8W' in name else 'FBN1 S2702');pair(ax,tab.loc[f,[p+'_T' for p in pats]],tab.loc[f,[p+'_LM' for p in pats]],unit,title)
save(fig,A,'A03_patient_examples')
cs=rd(A/'05_analysis/composition_patient_scores.tsv').set_index('set');ks=rd(A/'05_analysis/kinase_patient_scores.tsv').set_index('kinase');fig,axs=plt.subplots(1,2,figsize=(8,3.8),layout='constrained');pats=[c[:-3] for c in cs.columns if c.endswith('_LM')];pair(axs[0],cs.loc['Hep',[p+'_T' for p in pats]],cs.loc['Hep',[p+'_LM' for p in pats]],'Score','Hep proxy')
v=ks.loc['CHEK1'].dropna();axs[1].scatter(range(1,len(v)+1),v,color=np.where(v>0,C['LM'],C['T']),s=25);axs[1].axhline(0,color='black',lw=.8);axs[1].set(xlabel='Patient rank',ylabel='Mean Δ Z',title=f'CHEK1 subs\nn={len(v)} pairs; 3 genes');save(fig,A,'A04_composition_kinase')

ds=rd(B/'06_tables/discovery_protein_complete.tsv.gz');vs=rd(B/'06_tables/validation_protein_complete.tsv.gz');fig,axs=plt.subplots(1,2,figsize=(9,3.8),layout='constrained');volcano(axs[0],ds,'Disc: L+40 / L−40','Δ log2');volcano(axs[1],vs,'Valid: L+19 / L−41','Δ log2');save(fig,B,'B01_protein_screen')
pp=rd(B/'06_tables/pathway_cross_cohort_complete.tsv');pp=pp[pp.selected_D].copy();names={'HALLMARK_ANGIOGENESIS':'Angio','HALLMARK_COAGULATION':'Coag','HALLMARK_E2F_TARGETS':'E2F','HALLMARK_G2M_CHECKPOINT':'G2M'}
fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
for j,row in enumerate(pp.itertuples()):
 for c,off in [('D',-.13),('V',.13)]:
  e=getattr(row,'effect_'+c);lo=getattr(row,'CI_low_'+c);hi=getattr(row,'CI_high_'+c);ax.errorbar(e,j+off,xerr=[[e-lo],[hi-e]],fmt='o',color=C[c],capsize=3,label={'D':'Disc n=80','V':'Valid n=60'}[c] if j==0 else None)
ax.axvline(0,color='#5F6970',lw=.8);ax.set(yticks=range(len(pp)),yticklabels=[names[n] for n in pp.feature],xlabel='L+ − L− score',title='Locked paths');ax.legend(frameon=False);save(fig,B,'B02_pathway_replication')
fig,axs=plt.subplots(2,2,figsize=(8,7),layout='constrained')
for i,co in enumerate(['discovery','validation']):
 m=rd(B/f'01_metadata/{co}_samples.tsv');m=m[(m.Group=='Tumor')&m.L.notna()];sc=rd(B/f'05_analysis/{co}_pathway_patient_scores.tsv').set_index('set')
 for j,setname in enumerate(['HALLMARK_E2F_TARGETS','HALLMARK_ANGIOGENESIS']):dots(axs[i,j],sc.loc[setname,m['sample']].values,m.L.values,('Disc' if i==0 else 'Valid')+': '+names[setname])
save(fig,B,'B03_patient_pathways')
z=ds[ds.family=='primary'][['feature','effect']].merge(vs[vs.family=='primary'][['feature','effect']],on='feature',suffixes=('_D','_V'));fig,axs=plt.subplots(1,2,figsize=(9,4),layout='constrained');axs[0].scatter(z.effect_D,z.effect_V,s=5,color='#5B7E93',alpha=.5,rasterized=True);axs[0].axhline(0,color='#BCC5CB',lw=.8);axs[0].axvline(0,color='#BCC5CB',lw=.8);axs[0].set(xlabel='Disc Δ log2',ylabel='Valid Δ log2',title=f'Shared proteins\nn={len(z):,}')
for i,co in enumerate(['discovery','validation']):
 q=rd(B/f'04_qc/{co}_sample_qc.tsv');q=q[(q.Group=='Tumor')&q.L.notna()]
 for j,l in enumerate([0,1]):
  v=q[q.L==l].n_observed;axs[1].scatter(i*3+j+rng.normal(0,.055,len(v)),v,s=13,color=[C['T'],C['LM']][j],alpha=.7)
axs[1].set(xticks=[0,1,3,4],xticklabels=['D L−','D L+','V L−','V L+'],ylabel='Prot count',title='Tumor sample QC');save(fig,B,'B04_effect_qc')
print('FIGURES_COMPLETE')
