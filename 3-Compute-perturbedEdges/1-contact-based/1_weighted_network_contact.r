##############################################################################################
# Purpose: build a weighted network for each cancer type 
##############################################################################################

rm(list=ls())
library(data.table)
library(GenomicDataCommons)
source("eein_cancer_util.r")

cancer_type <- c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA')
cpm_threshold <- c(0.05, 0.1, 0.2, 0.5, 1, 2)

##-- read all expr data -- to speeed up ---
all_expr_data <- list()
for(j in 1:length(cancer_type)){
	all_expr_data[[j]] <- data.table::fread(paste0('../data/normalized_exons_sv/', cancer_type[j],'/normalized_all.txt'), header=TRUE)
}

##-- read all network data -- to speeed up ---
all_net_data <- data.table::fread(paste0('../data/CONTACT_networks/CONTACT_net_6_1.txt'))

##-- read all manifest data -- to speeed up ---
all_manifest <- list()
for(j in 1:length(cancer_type)){
	all_manifest[[j]] <- data.table::fread(paste0('../data/',cancer_type[j],'_manifest_final.txt'))
}

out_dir <- paste0('../data/CONTACT_weighted_networks/thres_6_1')
dir.create(out_dir, recursive=TRUE)

net_file <- all_net_data[,c(3,4)]
nodes <- union(net_file[[1]], net_file[[2]])

for(thr in 1:length(cpm_threshold)){

	out_dir1 <- paste0(out_dir,'/threshold_',cpm_threshold[thr])
	if(!dir.exists(out_dir1)){
		dir.create(out_dir1)
	}

	##-- for plot -----
	control_edges <- c()
	condition_edges <- c()
	intersect_edges <- c()
	pcancer <- c()
	##-----------------

	for(cn in 1:length(cancer_type)){

		all_survival <- get_survival(cancer_type[cn], 0)
		tempids <- all_manifest[[cn]]
		tempf <- as.data.frame(all_expr_data[[cn]])
		tempf <- tempf[tempf$EXON %in% nodes, ]
		sample_names_all <- unlist(lapply(strsplit(colnames(tempf), '[_]'), '[[', 1))
		whc <- which(sample_names_all %in% union(tempids$nid,'EXON'))
		tempf <- tempf[,whc]
		exons <- tempf$EXON ## this order is constant
		tempf <- tempf[,-length(tempf)]
		##-- only keep patients with control and condition ----------------------------
		temp_brk <- unlist(lapply(strsplit(colnames(tempf), '[_]'), '[[', 2))
		sample_names <- unlist(lapply(strsplit(colnames(tempf), '[_]'), '[[', 1))
		condition_pos_all <- which( (temp_brk == '1') | (temp_brk == '01') )
		control_pos <- which(temp_brk == '11')
		keep_samples <- sample_names[control_pos]
		condition_pos <- intersect(which(sample_names %in% keep_samples), condition_pos_all)
		tempf_control <- tempf[,control_pos]
		tempf_condition <- tempf[,condition_pos]
		tempf_control[tempf_control < cpm_threshold[thr] ] <- 0 
		tempf_condition[tempf_condition < cpm_threshold[thr] ] <- 0 

		##----for real expression--------------------------------
		out_real <- get_perturbed(tempf_control, tempf_condition, exons, net_file)
		td <- out_real[[1]]
		td$key <- paste0(td$from,'_',td$to)
		td <- td[order(td$key), ]
		data.table::fwrite(td, paste0(out_dir1,'/',cancer_type[cn],'_gained.txt'), row.names=FALSE, quote=FALSE, sep='\t')

		td <- out_real[[2]]
		td$key <- paste0(td$from,'_',td$to)
		td <- td[order(td$key), ]
		data.table::fwrite(td, paste0(out_dir1,'/',cancer_type[cn],'_lost.txt'), row.names=FALSE, quote=FALSE, sep='\t')

	}
}







