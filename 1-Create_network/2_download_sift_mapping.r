##############################################################################################
# Purpose: download all xml files of the cross mapping from SIFTS
##############################################################################################

rm(list=ls())
library(data.table)
library(gtools)
# turning warnings into errors
options(warn=2)


#################################################
uniprot_all <- fread('../data/processed_complex_final.txt', sep='\t', header=TRUE)
allpdbs <- unique(tolower(uniprot_all$ID))
######################################################
## check which SIFTS data are already present
# check if the folder exists
store_dir <- '../data/SIFTS'

if(dir.exists(store_dir)){
	allfiles <- list.files(store_dir)
	allpresent <- unlist(lapply(strsplit(allfiles, '[.]'), '[[', 1))
	todownload <- setdiff(allpdbs, allpresent)
}else{
  	dir.create(store_dir)
	todownload <- allpdbs
}

## download SIFT data
allsifts <- substr(todownload, 2,3)

for(k in 1:length(todownload)){

	output_name <- paste0(todownload[k],'.xml.gz')
	query <- paste0('https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/',allsifts[k],'/',output_name)
	cmd1 <- paste0('wget -O ',store_dir,'/',output_name,' ',query)
	cmd2 <- paste0("gunzip --force ", store_dir,'/',output_name)
	system(cmd1)
	system(cmd2)

}
