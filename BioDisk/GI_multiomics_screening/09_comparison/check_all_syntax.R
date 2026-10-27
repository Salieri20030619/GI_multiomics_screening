root <- '/mnt/BioDisk/wangmingyuan/GI_multiomics_screening/09_comparison'
files <- list.files(root,pattern='[.]R$',full.names=TRUE)
for (f in files) {invisible(parse(file=f));cat('SYNTAX_CHECK PASS',basename(f),'\n')}
cat('ALL_R_SYNTAX_COMPLETE\n')
