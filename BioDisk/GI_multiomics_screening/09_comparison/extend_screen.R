# Load only shared I/O/statistical helpers; do not rerun completed main comparisons.
code<-readLines('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening/09_comparison/screen.R')
eval(parse(text=code[seq_len(grep('^# A source comparison',code)-1)]))
suppressPackageStartupMessages(library(survival))
coxone <- function(time,event,value,strata_var=NULL) {
 d<-data.frame(time=time,event=event,value=value);if(!is.null(strata_var))d$batch<-strata_var
 d<-d[complete.cases(d)&d$time>0,,drop=FALSE];out<-data.frame(n=nrow(d),events=sum(d$event),HR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,PH_P=NA_real_)
 if(nrow(d)<20||sum(d$event)<10||sd(d$value)==0)return(out)
 d$value<-as.numeric(scale(d$value));f<-if(is.null(strata_var))Surv(time,event)~value else Surv(time,event)~value+strata(batch)
 fit<-try(coxph(f,data=d),silent=TRUE);if(inherits(fit,'try-error'))return(out)
 z<-summary(fit);out$HR<-z$conf.int[1,1];out$CI_low<-z$conf.int[1,3];out$CI_high<-z$conf.int[1,4];out$P<-z$coefficients[1,5]
 ph<-try(cox.zph(fit),silent=TRUE);if(!inherits(ph,'try-error'))out$PH_P<-ph$table[1,3];out
}

# B: explicitly validate the pathway family selected in discovery, then close route.
dp<-rd(file.path(B,'06_tables/discovery_pathways.tsv'));vp<-rd(file.path(B,'06_tables/validation_pathways.tsv'))
selected<-dp$feature[!is.na(dp$FDR)&dp$FDR<.05&dp$family=='primary'];wr(data.frame(set=selected),file.path(B,'01_metadata/discovery_pathways_locked.tsv'))
pjoin<-merge(dp,vp,by='feature',all.x=TRUE,suffixes=c('_D','_V'));pjoin$selected_D<-pjoin$feature%in%selected;pjoin$V_FDR_selected<-NA_real_
pp<-pjoin$P_V[pjoin$selected_D];pp[is.na(pp)]<-1;pjoin$V_FDR_selected[pjoin$selected_D]<-p.adjust(pp,'BH');pjoin$same_direction<-sign(pjoin$effect_D)==sign(pjoin$effect_V)
wr(pjoin,file.path(B,'06_tables/pathway_cross_cohort_complete.tsv'))
cat('B PRESELECTED PATHWAYS\n');print(as.data.frame(pjoin[pjoin$selected_D,c('feature','effect_D','FDR_D','effect_V','P_V','FDR_V','V_FDR_selected','same_direction')]))
obj<-readRDS(file.path(B,'05_analysis/screen_checkpoint.rds'));sens<-list();survres<-list();component_tests<-list()
for(cohort in c('discovery','validation')) {
 md<-obj$metadata[[cohort]];sc<-obj$pathway_scores[[cohort]];comp<-rd(file.path(B,paste0('05_analysis/',cohort,'_composition_patient_scores.tsv')));comp<-comp[match(md$sample,comp$sample),];stopifnot(identical(comp$sample,md$sample))
 covars<-setdiff(names(comp),'sample');mc<-cbind(md,comp[,-1,drop=FALSE]);cc<-complete.cases(mc[,covars,drop=FALSE]);mc<-mc[cc,];ss<-sc[,cc,drop=FALSE];des<-model.matrix(reformulate(c('L',covars)),mc)
 st<-stats_fit(ss,des,2,'B',mc$L);wr(st,file.path(B,paste0('06_tables/',cohort,'_pathway_composition_sensitivity.tsv')))
 component_tests[[cohort]]<-stats_fit(t(as.matrix(comp[,-1])),model.matrix(~L,md),2,'B',md$L)
 if(cohort=='discovery') {
  # Source 'censored'=1 is a recurrence event, verified against all 56 recurrence dates and author code.
  tm<-as.numeric(md$daystoevent);ev<-as.numeric(md$censored)
  for(set in intersect(selected,rownames(sc))) {sr<-coxone(tm,ev,sc[set,]);sr$set<-set;sr$cohort<-cohort;sr$endpoint<-'TTR_days';survres[[paste(cohort,set)]]<-sr}
 }
 # No OS inference for validation until independent endpoint coding evidence is available.
}
wr(rbindlist(component_tests,idcol='cohort'),file.path(B,'06_tables/composition_group_tests.tsv'))
if(length(survres)) {z<-rbindlist(survres);z$FDR<-p.adjust(z$P,'BH');wr(z,file.path(B,'06_tables/discovery_pathway_TTR_association.tsv'))}
# Same-gene scoring sensitivity, with the three required patient-composition comparisons.
gx<-lapply(c('discovery','validation'),function(co) {st<-obj$stats[[co]];collapse(obj$data[[co]][st$family=='primary',,drop=FALSE],obj$genes[[co]][st$family=='primary'])});names(gx)<-c('discovery','validation')
common<-intersect(rownames(gx$discovery),rownames(gx$validation));triad<-list();coverage<-list()
for(set in selected) {
 gg<-intersect(sets[[set]],common);coverage[[set]]<-data.frame(set=set,n_shared=length(gg),genes=paste(gg,collapse=';'))
 if(length(gg)<10)next
 for(co in names(gx)) {
  ns<-scores(gx[[co]][gg,,drop=FALSE],setNames(list(gg),set))$matrix
  os<-obj$pathway_scores[[co]];if(!set%in%rownames(os))next
  md<-obj$metadata[[co]];ok<-is.finite(ns[1,])&is.finite(os[set,]);selmd<-md[ok,];
  run<-function(v,g){if(min(table(g))<10)return(data.frame(effect=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_));q<-t.test(v[g==1],v[g==0]);data.frame(effect=mean(v[g==1])-mean(v[g==0]),CI_low=q$conf.int[1],CI_high=q$conf.int[2],P=q$p.value)}
  for(mode in c('original_patients_original_score','fixed_patients_original_score','fixed_patients_shared_score')) {
   if(mode=='original_patients_original_score'){v<-os[set,];g<-md$L}else{v<-if(mode=='fixed_patients_original_score')os[set,ok] else ns[1,ok];g<-selmd$L}
   v0<-is.finite(v);r<-run(v[v0],g[v0]);r$set<-set;r$cohort<-co;r$mode<-mode;r$n_pos<-sum(v0&g==1);r$n_neg<-sum(v0&g==0);triad[[paste(set,co,mode)]]<-r
  }
 }
}
wr(rbindlist(coverage),file.path(B,'04_qc/selected_pathway_shared_genes.tsv'))
if(length(triad)) {tr<-rbindlist(triad);tr[,FDR:=p.adjust(P,'BH'),by=.(cohort,mode)];wr(tr,file.path(B,'06_tables/pathway_shared_gene_sensitivity_triad.tsv'))}

# A: tissue composition and sample normalization sensitivities for actual measured endpoints.
cs<-mat(file.path(A,'05_analysis/composition_patient_scores.tsv'));ann<-list(protein=pm<-rd(file.path(A,'01_metadata/protein_features.tsv')),site=sm<-rd(file.path(A,'01_metadata/site_protein_mapping.tsv')))
for(kind in c('protein','site')) {
 delta<-mat(file.path(A,paste0('03_processed/',kind,'_paired_delta.tsv.gz')));pats<-colnames(delta)
 cd<-t(cs[,paste0(pats,'_LM'),drop=FALSE]-cs[,paste0(pats,'_T'),drop=FALSE]);des<-cbind(Intercept=1,cd);colnames(des)<-c('Intercept',rownames(cs));stopifnot(qr(des)$rank==ncol(des))
 st<-stats_fit(delta,des,1,'A');st$unit<-if(kind=='protein')'log2_conditional' else 'site_z_conditional';st$gene<-ann[[kind]]$gene[match(st$feature,ann[[kind]]$feature)];wr(st,file.path(A,paste0('06_tables/',kind,'_composition_sensitivity.tsv.gz')))
 wr(data.frame(patient=pats,des,check.names=FALSE),file.path(A,paste0('04_qc/',kind,'_composition_design.tsv')))
 full<-mat(file.path(A,if(kind=='protein')'03_processed/protein_log2.tsv.gz' else '03_processed/site_z.tsv.gz'));full<-sweep(full,2,apply(full,2,median,na.rm=TRUE),'-');dd<-full[,paste0(pats,'_LM'),drop=FALSE]-full[,paste0(pats,'_T'),drop=FALSE];colnames(dd)<-pats
 zz<-stats_fit(dd,matrix(1,ncol(dd),1),1,'A');zz$gene<-ann[[kind]]$gene[match(zz$feature,ann[[kind]]$feature)];wr(zz,file.path(A,paste0('06_tables/',kind,'_median_sensitivity.tsv.gz')))
}
cd<-t(cs[,grep('_LM$',colnames(cs)),drop=FALSE]-cs[,sub('_LM$','_T',grep('_LM$',colnames(cs),value=TRUE)),drop=FALSE]);cmp<-stats_fit(t(cd),matrix(1,nrow(cd),1),1,'A');wr(cmp,file.path(A,'06_tables/composition_paired_tests.tsv'))

# A external cohort: parse actual header, aggregate repeated lesions before all tests.
source<-file.path(A,'03_processed/source_tables/mmc4/5_Total_Cohort_Abundances.tsv.gz')
raw<-as.data.frame(fread(source,skip=2,na.strings=c('NA','NaN',''),check.names=FALSE));md<-rd(file.path(A,'01_metadata/external_specimens.tsv'))
scs<-intersect(md$Sample,names(raw));stopifnot(length(scs)==152);m<-as.matrix(raw[,scs,drop=FALSE]);storage.mode(m)<-'double';rownames(m)<-raw[[1]]
g<-raw[[2]];g[is.na(g)]<-gene(rownames(m))[is.na(g)];gx<-collapse(m,g)
ids<-unique(md$Patient_Nr);pg<-sapply(ids,function(id)apply(gx[,md$Sample[md$Patient_Nr==id],drop=FALSE],1,median,na.rm=TRUE));pg[!is.finite(pg)]<-NA
meta<-list();conflicts<-list()
fields<-c('Batch','Sex','Occurrence Time','Pre-OP treatment','survival status','survival time (Years)','Recurrence','Age at Diagnosis CRC-LM')
for(id in ids) {z<-md[md$Patient_Nr==id,];r<-list(patient=id,n_specimens=nrow(z));for(f in fields){vv<-unique(na.omit(z[[f]]));r[[f]]<-if(length(vv)==1)as.character(vv) else NA_character_;if(length(vv)>1)conflicts[[paste(id,f)]]<-data.frame(patient=id,field=f,values=paste(vv,collapse=';'))};meta[[id]]<-as.data.frame(r,check.names=FALSE)}
meta<-rbindlist(meta);meta$OS_event<-ifelse(meta[['survival status']]=='dead',1,ifelse(meta[['survival status']]=='alive',0,NA));meta$OS_years<-as.numeric(meta[['survival time (Years)']]);wr(meta,file.path(A,'01_metadata/external_patient_metadata.tsv'))
if(length(conflicts))wr(rbindlist(conflicts),file.path(A,'04_qc/external_patient_field_conflicts.tsv'))
es<-scores(pg[rowMeans(!is.na(pg))>=.7,,drop=FALSE],sets);ec<-scores(pg,panels,2);wr(es$coverage,file.path(A,'04_qc/external_pathway_coverage.tsv'));wr(data.frame(set=rownames(es$matrix),es$matrix,check.names=FALSE),file.path(A,'05_analysis/external_pathway_patient_scores.tsv'))
wr(data.frame(gene=rownames(pg),pg,check.names=FALSE),file.path(A,'03_processed/external_patient_gene_log2.tsv.gz'))
out<-list()
for(set in rownames(es$matrix)) {v<-es$matrix[set,meta$patient];r<-coxone(meta$OS_years,meta$OS_event,v,meta$Batch);r$set<-set;r$endpoint<-'OS_years';out[[set]]<-r}
z<-rbindlist(out);z$FDR<-p.adjust(z$P,'BH');wr(z,file.path(A,'06_tables/external_pathway_OS_complete.tsv'))
# Source dependency list is secondary annotation; absence is not a negative CRISPR result.
deps<-rd(file.path(A,'03_processed/source_tables/ADVS-13-e11744-s005/S8A.tsv.gz'))
cond<-rd(file.path(A,'06_tables/site_conditional_diagnostic.tsv.gz'));base<-rd(file.path(A,'06_tables/site_paired_complete.tsv.gz'));pr<-rd(file.path(A,'06_tables/protein_paired_complete.tsv.gz'));sen<-rd(file.path(A,'06_tables/site_composition_sensitivity.tsv.gz'))
ev<-merge(cond,base[,c('feature','effect','FDR','LOO_min_effect','LOO_max_effect')],by='feature',suffixes=c('_conditional','_site'))
ev$protein_effect<-pr$effect[match(ev$protein_feature,pr$feature)];ev$protein_FDR<-pr$FDR[match(ev$protein_feature,pr$feature)];ev$composition_effect<-sen$effect[match(ev$feature,sen$feature)];ev$composition_FDR<-sen$FDR[match(ev$feature,sen$feature)]
ev$external_gene_measured<-ev$gene%in%rownames(pg);ev$author_DepMap_mean<-as.numeric(deps[[2]][match(ev$gene,deps[[1]])]);ev$dependency_status<-ifelse(is.na(ev$author_DepMap_mean),'not_in_author_selected_list_not_a_negative_test','author_selected_CRC_dependency_not_site_specific')
ev$independent_paired_site_validation<-'unavailable';ev$protein_quantity_adjustment<-'not_estimable';ev$recommendation<-'not_recommended_current_evidence'
wr(ev,file.path(A,'06_tables/cross_cohort_site_evidence_complete.tsv.gz'))
cat('EXTENDED_SCREEN_COMPLETE\n')
