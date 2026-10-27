suppressPackageStartupMessages({library(limma);library(data.table)})
R <- '/mnt/BioDisk/wangmingyuan/GI_multiomics_screening'
A <- file.path(R,'A_CRLM_phosphoregulation'); B <- file.path(R,'B_ICC_lymphatic_invasion')
rd <- function(p) as.data.frame(fread(p,na.strings=c('NA',''),check.names=FALSE))
mat <- function(p) { d<-rd(p); rn<-d[[1]]; d<-as.matrix(d[,-1,drop=FALSE]);storage.mode(d)<-'double';rownames(d)<-rn;d }
wr <- function(d,p) fwrite(as.data.frame(d),p,sep='\t',na='NA',quote=FALSE)
map <- rd(file.path(R,'09_comparison/reference/uniprot_gene_map.tsv'))
maps <- split(map$SYMBOL,map$UNIPROT)
gene <- function(x) vapply(sub('-[0-9]+$','',x),function(a) {v<-unique(na.omit(maps[[a]]));if(length(v)==1)v else NA_character_},character(1))
gmt <- strsplit(readLines(file.path(R,'09_comparison/reference/hallmark.gmt')),'\t');sets<-setNames(lapply(gmt,function(v)unique(v[-c(1,2)])),vapply(gmt,`[`,character(1),1))
panels <- list(Hep=c('ALB','TTR','APOA1','APOA2','APOB','ASGR1','ASGR2','CPS1','ARG1','ALDOB','HPD'),Blood=c('HBA1','HBA2','HBB','HBD','AHSP'),Immune=c('PTPRC','LST1','TYROBP','FCER1G','CD3D','CD3E','CD74'),Stroma=c('COL1A1','COL1A2','COL3A1','DCN','LUM','COL6A1','COL6A2','COL6A3'),Epi=c('EPCAM','KRT8','KRT18','KRT19'))
collapse <- function(x,g) {ok<-!is.na(g)&grepl('^[A-Za-z0-9_-]+$',g);x<-x[ok,,drop=FALSE];g<-g[ok]; v<-lapply(split(seq_along(g),g),function(i)apply(x[i,,drop=FALSE],2,median,na.rm=TRUE));d<-do.call(rbind,v);d[!is.finite(d)]<-NA;d }
scores <- function(x,gs,min_n=10) {
 z<-(x-rowMeans(x,na.rm=TRUE))/apply(x,1,sd,na.rm=TRUE);z[!is.finite(z)]<-NA
 ans<-list();coverage<-list()
 for(n in names(gs)) {ii<-intersect(rownames(z),gs[[n]]);coverage[[n]]<-data.frame(set=n,reference_n=length(gs[[n]]),measured_n=length(ii),genes=paste(ii,collapse=';'))
  if(length(ii)<min_n)next
  v<-z[ii,,drop=FALSE];a<-colMeans(v,na.rm=TRUE);a[colSums(!is.na(v))<ceiling(length(ii)*.5)]<-NA;ans[[n]]<-a
 }
 list(matrix=do.call(rbind,ans),coverage=do.call(rbind,coverage))
}
stats_fit <- function(x,design,coef,kind='B',group=NULL,main_coverage=.7) {
 n<-rowSums(!is.na(x));res<-data.frame(feature=rownames(x),n=n,effect=NA_real_,SE=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,FDR=NA_real_,family='descriptive',df=NA_real_)
 if(kind=='A') {res$effect<-rowMeans(x,na.rm=TRUE);res$n_pos<-rowSums(x>0,na.rm=TRUE);res$n_neg<-rowSums(x<0,na.rm=TRUE);res$family[n>=20]<-'primary';res$family[n>=10&n<20]<-'exploratory_10_19';ok<-n>=10}
 else {res$n_pos<-rowSums(!is.na(x[,group==1,drop=FALSE]));res$n_neg<-rowSums(!is.na(x[,group==0,drop=FALSE]));res$effect<-rowMeans(x[,group==1,drop=FALSE],na.rm=TRUE)-rowMeans(x[,group==0,drop=FALSE],na.rm=TRUE);ok<-pmin(res$n_pos,res$n_neg)>=10;res$family[ok]<-'exploratory_low_coverage';res$family[ok & n/ncol(x)>=main_coverage]<-'primary'}
 if(sum(ok)>1) {
  fit<-lmFit(x[ok,,drop=FALSE],design);fit<-eBayes(fit,trend=TRUE)
  se<-fit$stdev.unscaled[,coef]*sqrt(fit$s2.post);ci<-qt(.975,fit$df.total)*se
  res$effect[ok]<-fit$coefficients[,coef];res$SE[ok]<-se;res$CI_low[ok]<-res$effect[ok]-ci;res$CI_high[ok]<-res$effect[ok]+ci;res$P[ok]<-fit$p.value[,coef];res$df[ok]<-fit$df.total
  for(f in unique(res$family[ok])){ii<-which(res$family==f);res$FDR[ii]<-p.adjust(res$P[ii],'BH')}
 }
 res
}

# A source comparison in its actual units.
pm<-rd(file.path(A,'01_metadata/protein_features.tsv'));sm<-rd(file.path(A,'01_metadata/site_protein_mapping.tsv'))
for(kind in c('protein','site')) {
 d<-mat(file.path(A,paste0('03_processed/',kind,'_paired_delta.tsv.gz')));design<-matrix(1,ncol(d),1,dimnames=list(colnames(d),'Intercept'))
 st<-stats_fit(d,design,1,'A');ann<-if(kind=='protein')pm else sm
 st<-merge(ann,st,by='feature',all.x=TRUE,sort=FALSE);st$unit<-if(kind=='protein')'log2_ratio' else 'sample_z_score'
 st$P_unmoderated<-vapply(st$feature,function(id){v<-na.omit(d[id,]);if(length(v)>=10&&sd(v)>0)t.test(v)$p.value else NA_real_},numeric(1))
 st$LOO_min_effect<-vapply(st$feature,function(id){v<-na.omit(d[id,]);if(length(v)>1)min((sum(v)-v)/(length(v)-1)) else NA_real_},numeric(1))
 st$LOO_max_effect<-vapply(st$feature,function(id){v<-na.omit(d[id,]);if(length(v)>1)max((sum(v)-v)/(length(v)-1)) else NA_real_},numeric(1))
 wr(st,file.path(A,paste0('06_tables/',kind,'_paired_complete.tsv.gz')))
 cat('A',kind,'primary',sum(st$family=='primary'),'FDR05',sum(st$family=='primary'&st$FDR<.05,na.rm=TRUE),'\n')
 full<-mat(file.path(A,if(kind=='protein')'03_processed/protein_log2.tsv.gz' else '03_processed/site_z.tsv.gz'))
 pats<-colnames(d);full<-full[,c(paste0(pats,'_T'),paste0(pats,'_LM')),drop=FALSE]
 g<-ann$gene[match(rownames(full),ann$feature)];gx<-collapse(full[rowSums(!is.na(d))>=20,,drop=FALSE],g[rowSums(!is.na(d))>=20])
 sc<-scores(gx,sets);wr(sc$coverage,file.path(A,paste0('04_qc/',kind,'_pathway_coverage.tsv')))
 if(!is.null(sc$matrix)) {ss<-sc$matrix;dd<-ss[,paste0(pats,'_LM'),drop=FALSE]-ss[,paste0(pats,'_T'),drop=FALSE];colnames(dd)<-pats;ps<-stats_fit(dd,matrix(1,ncol(dd),1),1,'A');ps$unit<-'gene_z_mean';wr(ps,file.path(A,paste0('06_tables/',kind,'_pathways.tsv')));wr(data.frame(set=rownames(ss),ss,check.names=FALSE),file.path(A,paste0('05_analysis/',kind,'_pathway_patient_scores.tsv')))}
 if(kind=='protein') {
  cs<-scores(collapse(full,g),panels,2);wr(cs$coverage,file.path(A,'04_qc/composition_coverage.tsv'));wr(data.frame(set=rownames(cs$matrix),cs$matrix,check.names=FALSE),file.path(A,'05_analysis/composition_patient_scores.tsv'))
 }
}
# Explicitly conditional diagnostic, not a stoichiometric/quantity correction.
pd<-mat(file.path(A,'03_processed/protein_paired_delta.tsv.gz'));sd<-mat(file.path(A,'03_processed/site_paired_delta.tsv.gz'))
cond<-lapply(seq_len(nrow(sm)),function(i) {
 r<-sm[i,];out<-data.frame(feature=r$feature,site=r$site,gene=r$gene,protein_feature=r$protein_feature,n=0,effect=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,protein_slope=NA_real_,family='not_mapped',unit='site_z_conditional')
 if(is.na(r$protein_feature)||!nzchar(r$protein_feature))return(out)
 yp<-sd[r$feature,];xp<-pd[r$protein_feature,colnames(sd)];ok<-is.finite(yp)&is.finite(xp);out$n<-sum(ok)
 out$family<-if(out$n>=20)'primary_diagnostic' else if(out$n>=10)'exploratory_diagnostic' else 'descriptive'
 if(out$n>=10&&sd(xp[ok])>0) {fit<-lm(yp[ok]~xp[ok]);z<-summary(fit)$coefficients;out$effect<-z[1,1];out$P<-z[1,4];ci<-confint(fit)[1,];out$CI_low<-ci[1];out$CI_high<-ci[2];out$protein_slope<-z[2,1]}
 out
})
cond<-rbindlist(cond);cond[,FDR:=p.adjust(P,'BH'),by=family];wr(cond,file.path(A,'06_tables/site_conditional_diagnostic.tsv.gz'))
wr(data.frame(feature=sm$feature,site=sm$site,gene=sm$gene,status='NOT_ESTIMABLE',reason='Samplewise phosphosite Z transformation; inverse scaling parameters unavailable'),file.path(A,'06_tables/site_quantity_adjusted_status.tsv.gz'))

# B discovery and independent validation: independent fits and fixed discovery selection.
bst<-list();bdata<-list();bmeta<-list();bgene<-list();bscores<-list()
for(cohort in c('discovery','validation')) {
 x<-mat(file.path(B,paste0('03_processed/',cohort,'_log2.tsv.gz')));md<-rd(file.path(B,paste0('01_metadata/',cohort,'_samples.tsv')))
 md<-md[md$Group=='Tumor'&!is.na(md$L),];x<-x[,md$sample,drop=FALSE];stopifnot(!anyDuplicated(md$patient),identical(colnames(x),md$sample))
 g<-gene(rownames(x));design<-model.matrix(~L,md);st<-stats_fit(x,design,2,'B',md$L);st$gene<-g;st$cohort<-cohort;st$unit<-'log2'
 st$P_Welch<-vapply(seq_len(nrow(x)),function(i){a<-na.omit(x[i,md$L==1]);b<-na.omit(x[i,md$L==0]);if(min(length(a),length(b))>=10&&sd(c(a,b))>0)t.test(a,b)$p.value else NA_real_},numeric(1))
 st$P_missingness<-vapply(seq_len(nrow(x)),function(i)fisher.test(table(factor(md$L,levels=0:1),factor(is.na(x[i,]),levels=c(FALSE,TRUE))))$p.value,numeric(1));st$FDR_missingness<-p.adjust(st$P_missingness,'BH')
 wr(st,file.path(B,paste0('06_tables/',cohort,'_protein_complete.tsv.gz')));bst[[cohort]]<-st;bdata[[cohort]]<-x;bmeta[[cohort]]<-md;bgene[[cohort]]<-g
 cat('B',cohort,'groups',table(md$L),'primary',sum(st$family=='primary'),'FDR05',sum(st$family=='primary'&st$FDR<.05,na.rm=TRUE),'\n')
 if(cohort=='discovery') {sel<-st[st$family=='primary'&!is.na(st$FDR)&st$FDR<.05&!is.na(st$gene),];wr(sel,file.path(B,'01_metadata/discovery_candidates_locked.tsv'))}
 gx<-collapse(x[st$family=='primary',,drop=FALSE],g[st$family=='primary']);sc<-scores(gx,sets);bscores[[cohort]]<-sc$matrix
 wr(sc$coverage,file.path(B,paste0('04_qc/',cohort,'_pathway_coverage.tsv')))
 ps<-stats_fit(sc$matrix,design,2,'B',md$L);ps$unit<-'gene_z_mean';wr(ps,file.path(B,paste0('06_tables/',cohort,'_pathways.tsv')))
 wr(data.frame(set=rownames(sc$matrix),sc$matrix,check.names=FALSE),file.path(B,paste0('05_analysis/',cohort,'_pathway_patient_scores.tsv')))
 cs<-scores(collapse(x,g),panels,2);wr(cs$coverage,file.path(B,paste0('04_qc/',cohort,'_composition_coverage.tsv')))
 wr(data.frame(sample=colnames(cs$matrix),t(cs$matrix),check.names=FALSE),file.path(B,paste0('05_analysis/',cohort,'_composition_patient_scores.tsv')))
 comp<-as.data.frame(t(cs$matrix));md<-cbind(md,comp[md$sample,,drop=FALSE]);
 # Full-cohort sensitivity retains patients; no sample filtering based on outcomes.
 for(sens in c('clinical','composition')) {
  covars<-if(sens=='clinical')c('sex','age_ordinal',if(cohort=='discovery')'duct' else 'grade') else intersect(c('Hep','Blood','Immune','Stroma','Epi'),names(md))
  cc<-complete.cases(md[,c('L',covars),drop=FALSE]);mc<-md[cc,,drop=FALSE]; xx<-x[,cc,drop=FALSE]
  covars<-covars[vapply(mc[,covars,drop=FALSE],function(v)sd(v)>1e-10,logical(1))]
  design2<-model.matrix(reformulate(c('L',covars)),mc);stopifnot(qr(design2)$rank==ncol(design2))
  base_fixed<-stats_fit(xx,model.matrix(~L,mc),2,'B',mc$L);zz<-stats_fit(xx,design2,2,'B',mc$L);zz$gene<-g;zz$fixed_patient_effect<-base_fixed$effect;zz$fixed_patient_P<-base_fixed$P
  wr(zz,file.path(B,paste0('06_tables/',cohort,'_',sens,'_sensitivity.tsv.gz')))
  wr(data.frame(sample=mc$sample,patient=mc$patient,design2,check.names=FALSE),file.path(B,paste0('04_qc/',cohort,'_',sens,'_design.tsv')))
 }
}
v<-bst$validation;ii<-match(sel$feature,v$feature);cv<-cbind(sel[,c('feature','gene','effect','CI_low','CI_high','P','FDR','n_pos','n_neg')],v[ii,c('effect','CI_low','CI_high','P','FDR','n_pos','n_neg','family')]);names(cv)<-c('feature','gene','D_effect','D_CI_low','D_CI_high','D_P','D_FDR','D_n_pos','D_n_neg','V_effect','V_CI_low','V_CI_high','V_P','V_FDR_all','V_n_pos','V_n_neg','V_family')
# Use original locked discovery family size, including unavailable validation candidates as P=1.
pv<-cv$V_P;pv[is.na(pv)]<-1;cv$V_FDR_candidates<-p.adjust(pv,'BH');cv$same_direction<-sign(cv$D_effect)==sign(cv$V_effect);cv$replicated<-!is.na(cv$V_P)&cv$V_FDR_candidates<.05&cv$same_direction
wr(cv,file.path(B,'06_tables/cross_cohort_candidate_evidence.tsv'));cat('B LOCKED',nrow(cv),'VALIDATION_FDR05_SAME',sum(cv$replicated,na.rm=TRUE),'\n')
saveRDS(list(stats=bst,data=bdata,metadata=bmeta,genes=bgene,pathway_scores=bscores),file.path(B,'05_analysis/screen_checkpoint.rds'))
writeLines(capture.output(sessionInfo()),file.path(R,'09_comparison/screen_sessionInfo.txt'))
cat('MAIN_SCREEN_COMPLETE\n')
