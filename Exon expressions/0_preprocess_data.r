##############################################################################################
# Purpose: Download reference; Make reference index for alignment; Preprocess ref for counting
# This script has to run once only
##############################################################################################

rm(list=ls())

dir.create('../data')

## Download reference...
# Download GTF file from Ensembl
xx <- "Homo_sapiens.GRCh38.105.gtf.gz"
system(paste0("wget -O ../data/",xx," http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/",xx))
system(paste0("gunzip ../data/",xx))

## Preprocess ref for counting
system('sh ./preprocess_exons.sh')

