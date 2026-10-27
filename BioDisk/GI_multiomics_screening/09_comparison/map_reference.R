suppressPackageStartupMessages({library(org.Hs.eg.db);library(AnnotationDbi)})
root <- '/mnt/BioDisk/wangmingyuan/GI_multiomics_screening/09_comparison'
ids <- readLines(file.path(root,'accessions.txt'))
ids <- unique(sub('-[0-9]+$','',ids))
valid <- intersect(ids,keys(org.Hs.eg.db,keytype='UNIPROT'))
m <- select(org.Hs.eg.db,keys=valid,keytype='UNIPROT',columns=c('SYMBOL','ENTREZID'))
write.table(m,file.path(root,'reference/uniprot_gene_map.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
if(file.exists(file.path(root,'reference/KSEAdata.RData'))) {
  objs <- load(file.path(root,'reference/KSEAdata.RData'))
  for(n in objs) { d <- get(n); if(is.data.frame(d)) write.table(d,file.path(root,paste0('reference/',n,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA') }
}
cat('REFERENCE_MAPPING_COMPLETE\n')
