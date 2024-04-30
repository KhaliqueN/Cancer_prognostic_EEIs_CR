##############################################################################################
# Purpose: run eppic with the given pdb ids
#
######################################################################################################

rm(list=ls())

outputDir <- '../data/EPPIC'
if(!dir.exists(outputDir)){
	dir.create(outputDir)
}

allpdbs <- data.table::fread('../data/uniprot_pdb_Ensembl_finalized.txt', header=TRUE)[[1]]


##-- which ids are already present
allids <- list.files(outputDir, pattern='*.scores')
toprocess <- setdiff(allpdbs, unlist(lapply(strsplit(allids, '[.]'), '[[',1)))

##---
## EPPIC not giving output for the following eight PDBIDS: c("1n6j", "6l9z", "2as5", "2j6f", "5jcz", "2j6o", "7aew", "6xdl")

for(k in 1:length(toprocess)){

	cmd <- paste0('eppic -i ',toprocess[k],' -o ',outputDir,' -s -a 30')
	system(cmd)

}



