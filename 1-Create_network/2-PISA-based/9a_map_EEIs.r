##############################################################################################
# Purpose: build EEI network using the PISA files
## Residues with non-zero Solvation energies are interface residues
## First, select whether interface is biological/crystal (based on "P-value" < delta)
## Then, select the residues of the two chains that are part of the interface (Solvation energy != 0)
######################################################################################################

rm(list=ls())
library(data.table)
library(ggplot2)
source("eein_cancer_util.r")

saveDirectory <- '../data/PISA_networks'
dir.create(saveDirectory, recursive=TRUE)
pval_thres <- c(0.5)

allpdbs <- data.table::fread('../data/uniprot_pdb_Ensembl_finalized.txt', header=TRUE)[[1]]

pisafiles <- list.files('../data/PISA_files', full.names=TRUE)

##-- files not present on PISA
not_pisa <- setdiff(allpdbs, basename(pisafiles))
##-- 26 files (PDB IDs)

allmaps <- list.files('../data/uniprot_EnsemblExonPDB_map', full.names=TRUE)
input_dir <- list.files('../data/PISA_files_parsed', full.names=TRUE)

##-- files present on PISA but the residues files could not be downloaded-----
##-- 20 files (PDB IDs) ------
start <- 1
end <- length(pisafiles)
not_processed <- c()

for(k in start:end){

	temp <- list.files(pisafiles[k], pattern='^residue')

	if(length(temp) == 0){
		not_processed <- c(not_processed, pisafiles[k])
	}
	
}

##----------------------------------------------------------------------------
##--- So out of 2,064 PDB IDs, 2,064-26-20 = 2,018 PDB IDs remain ----
missed_pdbs <- union(basename(not_processed), not_pisa)
processed_pdbs <- setdiff(allpdbs, missed_pdbs)

for(pp in 1:length(pval_thres)){

	EXON1 <- c()
	EXON2 <- c()
	PROTEIN1 <- c()
	PROTEIN2 <- c()
	INTFAA1 <- c()
	INTFAA2 <- c()
	FREE_ENERGY <- c()
	BURIED_AREA <- c()
	HYDROGEN <- c()
	DISULPHIDE <- c()
	SALTBRIDGE <- c()
	COVALENT <- c()
	BURIED_AREA_ABS <- c()
    ASA_AREA_ABS <- c()
    PDBID <- c()
	not_mapped_PDB <- c()

	for(k in 1:length(input_dir)){

		##-- check whether the PDB IDs is processed or not
		if(basename(input_dir[k]) %in% processed_pdbs){

			##--- select the map files relevant for this PDB id
			maps <- allmaps[which(allmaps %like% basename(input_dir[k]))]

			if(length(maps) != 0){

				relevant_chains <- trimws(unlist(lapply(strsplit(unlist(lapply(strsplit(basename(maps),
				 '[_]'), '[[', 3)), '[.]'),'[[',1)))

				##-- evaluate whether the interface is biological or not --> and select the corresponding chains
				int_tab <- data.table::fread(paste0(input_dir[k],'/interaction_table.csv'))
				colnames(int_tab) <- paste0('X',seq(1,length(int_tab)))
				int_tab_filt <- int_tab[int_tab$X9 < pval_thres[pp], ]

				if(nrow(int_tab_filt) != 0){ ## if the at least one entry remains in the interface file
					temp_call <- select_edges(int_tab_filt, relevant_chains, input_dir[k], maps)
					EXON1 <- c(EXON1,  temp_call[[3]])
					EXON2 <- c(EXON2, temp_call[[4]])
					PROTEIN1 <- c(PROTEIN1, temp_call[[1]])
					PROTEIN2 <- c(PROTEIN2, temp_call[[2]])
					INTFAA1 <- c(INTFAA1, temp_call[[5]])
					INTFAA2 <- c(INTFAA2, temp_call[[6]])
					FREE_ENERGY <- c(FREE_ENERGY, temp_call[[7]])
					BURIED_AREA <- c(BURIED_AREA, temp_call[[8]])
					HYDROGEN <- c(HYDROGEN, temp_call[[9]])
					DISULPHIDE <- c(DISULPHIDE, temp_call[[10]])
					SALTBRIDGE <- c(SALTBRIDGE, temp_call[[11]])
					COVALENT <- c(COVALENT, temp_call[[12]])
					BURIED_AREA_ABS <- c(BURIED_AREA_ABS, temp_call[[13]])
				    ASA_AREA_ABS <- c(ASA_AREA_ABS, temp_call[[14]])
				    PDBID <- c(PDBID, rep(basename(input_dir[k]), length(temp_call[[1]])))
				}
			}else{
				not_mapped_PDB <- c(not_mapped_PDB, basename(input_dir[k]))
			}
		}

	}

	allData <- data.frame(exon1=EXON1, exon2=EXON2, AA1=INTFAA1, AA2=INTFAA2, protein1=PROTEIN1, protein2=PROTEIN2, 
	FreeEnergy=FREE_ENERGY, BuriedArea=BURIED_AREA, Hydrogen=HYDROGEN, Disulphide=DISULPHIDE, Saltbridge=SALTBRIDGE,
	Covalent=COVALENT, BuriedAreaAbs=BURIED_AREA_ABS, SolAccAreaAbs=ASA_AREA_ABS, PDBID=PDBID)

	data.table::fwrite(allData, paste0(saveDirectory,'/PISA_EEIN_',pval_thres[pp],'.txt'), 
		row.names=FALSE, quote=FALSE, sep='\t')
	cat('Threshold', pp, 'of',length(pval_thres),'done\n')

}

