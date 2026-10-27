code<-readLines('/mnt/BioDisk/wangmingyuan/GI_multiomics_screening/09_comparison/screen.R');eval(parse(text=code[seq_len(grep('^# A source comparison',code)-1)]))
k<-rd(file.path(R,'09_comparison/reference/KSData.tsv'));k<-k[k$KIN_ORGANISM=='human'&k$SUB_ORGANISM=='human'&k$Source=='PhosphoSitePlus',]
sm<-rd(file.path(A,'01_metadata/site_protein_mapping.tsv'));d<-mat(file.path(A,'03_processed/site_paired_delta.tsv.gz'));x<-mat(file.path(A,'03_processed/site_z.tsv.gz'))
pm<-rd(file.path(A,'06_tables/protein_paired_complete.tsv.gz'));sp<-rd(file.path(A,'06_tables/site_paired_complete.tsv.gz'))
sm$n_pairs<-rowSums(!is.na(d))[match(sm$feature,rownames(d))]
# Exact accession and residue number; isoforms are not silently mapped to canonical coordinates.
key<-paste(sm$accession,sm$residue_position,sep='_');k$key<-paste(k$SUB_ACC_ID,k$SUB_MOD_RSD,sep='_')
hits<-merge(k,sm[,c('feature','site','gene','accession','residue_position','n_pairs')],by.x='key',by.y='site');hits<-unique(hits[,c('GENE','feature','gene','key','n_pairs')]);wr(hits,file.path(A,'01_metadata/kinase_substrate_exact_mapping.tsv'))
out<-list();cov<-list();ss<-list()
for(kin in unique(k$GENE)) {
 h<-hits[hits$GENE==kin&hits$n_pairs>=20,];ids<-unique(h$feature);ng<-length(unique(h$gene));ok<-length(ids)>=5&&ng>=3
 cov[[kin]]<-data.frame(kinase=kin,n_sites=length(ids),n_substrate_genes=ng,status=if(ok)'tested' else 'insufficient_5_sites_3_genes',substrates=paste(unique(h$gene),collapse=';'))
 if(!ok)next
 # Average within substrate gene first; then equal-weight the multiple substrate proteins per patient.
 z<-collapse(d[ids,,drop=FALSE],h$gene[match(ids,h$feature)]);v<-colMeans(z,na.rm=TRUE);v[colSums(!is.na(z))<ceiling(nrow(z)*.5)]<-NA
 score<-matrix(v,nrow=1,dimnames=list(kin,names(v)));n<-sum(!is.na(v));q<-if(n>=10)t.test(v) else NULL
 r<-data.frame(kinase=kin,n=n,n_sites=length(ids),n_substrates=ng,effect=mean(v,na.rm=TRUE),CI_low=if(!is.null(q))q$conf.int[1] else NA,CI_high=if(!is.null(q))q$conf.int[2] else NA,P=if(!is.null(q))q$p.value else NA,patient_positive=sum(v>0,na.rm=TRUE),patient_negative=sum(v<0,na.rm=TRUE),substrate_positive=sum(rowMeans(z,na.rm=TRUE)>0),substrate_negative=sum(rowMeans(z,na.rm=TRUE)<0),unit='site_z_delta',quantity_adjustment='not_estimable')
 pv<-pm[pm$gene==kin,];r$kinase_protein_effect<-if(nrow(pv)==1)pv$effect else NA_real_;r$kinase_protein_FDR<-if(nrow(pv)==1)pv$FDR else NA_real_;out[[kin]]<-r;ss[[kin]]<-v
}
coverage<-rbindlist(cov);wr(coverage,file.path(A,'04_qc/kinase_coverage_all.tsv'))
if(length(out)) {z<-rbindlist(out);z$FDR<-p.adjust(z$P,'BH');wr(z,file.path(A,'06_tables/kinase_patient_scores_complete.tsv'));mm<-do.call(rbind,ss);wr(data.frame(kinase=rownames(mm),mm,check.names=FALSE),file.path(A,'05_analysis/kinase_patient_scores.tsv'));cat('KINASES',nrow(z),'FDR05',sum(z$FDR<.05,na.rm=TRUE),'\n');print(as.data.frame(z[order(z$P),][1:min(8,nrow(z)),]))}
cat('KINASE_SCREEN_COMPLETE\n')
