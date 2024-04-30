##############################################################################################
# Purpose: parse the PISA xml files to save a csv file per complex
##############################################################################################

rm(list=ls())

pisa_pdbs <- list.files('../data/PISA_files', full.names=TRUE, recursive=FALSE)
pisabase <- basename(pisa_pdbs)

store_dir <- '../data/PISA_files_parsed'
if(!dir.exists(store_dir)){
	dir.create(store_dir)
}

pisa_parsed <- list.files('../data/PISA_files_parsed', full.names=TRUE, recursive=FALSE)

pisads <- setdiff(pisabase, basename(pisa_parsed))

wh <- which(pisabase %in% pisads)
pisa_pdbs <- pisa_pdbs[wh]

# ##--- check which pdb folders have no interaction table -------when dir is created but the files not processed--------------------------
# allp <- list.files('../data/PISA_files_parsed', full.names=TRUE)
# start <- 1
# end <- length(allp)
# not_processed <- c()

# for(k in start:end){

# 	temp <- paste0(allp[k],'/interaction_table.csv')

# 	if(!file.exists(temp)){
# 		not_processed <- c(not_processed, allp[k])
# 	}
	
# }
# left_files <- unlist(lapply(strsplit(basename(not_processed),'[.]'), '[[', 1))
# wh <- which(pisabase %in% left_files)
# pisa_pdbs <- pisa_pdbs[wh]
# ##----------------------------------------------------------------------------------------

for(k in 1:length(pisa_pdbs)){

	store_dir1 <- paste0(store_dir,'/',basename(pisa_pdbs[k]))
	if(!dir.exists(store_dir1)){
		dir.create(store_dir1)
	}

	tempf <- paste0(store_dir1,'/interaction_table.csv')
	if(!file.exists(tempf)){

		temp <- paste0(pisa_pdbs[k],'/interfacetable.xml')
		tmpcmd <- paste('python Pisa_table_parser.py',temp,store_dir1,'interaction_table')
		system(tmpcmd)

		residue_files <- gtools::mixedsort(list.files(pisa_pdbs[k], pattern='^residue', full.names=TRUE))

		for(j in 1:length(residue_files)){
			tmpcmd <- paste('python Pisa_residue_parser.py',residue_files[j],store_dir1,paste0('residue_',j))
			system(tmpcmd)
		}
	}
	cat('PDB', k, 'of',length(pisa_pdbs),'done\n')
}

