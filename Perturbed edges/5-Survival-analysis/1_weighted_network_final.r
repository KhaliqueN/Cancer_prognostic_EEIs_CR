##############################################################################################--
# Purpose: build a weighted network for each cancer type, for each of the three networks 
########################################################--######################################--

args=commandArgs(TRUE)

library(data.table)
library(GenomicDataCommons)

cancer_type <- c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA')
cpm_threshold <- c(0.5)
level <- args[1]

##-- read all expr data -- to speeed up
all_expr_data <- list()
for(j in 1:length(cancer_type)){
	all_expr_data[[j]] <- data.table::fread(paste0('../data/normalized_exons_sv/', cancer_type[j],'/normalized_all.txt'), header=TRUE)
}

##-- read all network data -- to speeed up -----
allnets <- gtools::mixedsort(list.files(paste0('../data/Final_networks/',level),full.names=TRUE))
all_net_data <- list()
for(j in 1:length(allnets)){
	temp <- data.table::fread(allnets[j], header=FALSE)
	colnames(temp) <- c('exon1','exon2')
	all_net_data[[j]] <- temp
}

##-- read all manifest data -- to speeed up
all_manifest <- list()
for(j in 1:length(cancer_type)){
	all_manifest[[j]] <- data.table::fread(paste0('../data/Manifests/',cancer_type[j],'_manifest_final.txt'))
}

for(qq in 1:length(allnets)){

		out_dir <- paste0('../data/Final_weighted_networks/',level,'/',strsplit(basename(allnets[qq]),'[.]')[[1]][1])
		dir.create(out_dir, recursive=TRUE)

		net_file <- all_net_data[[qq]][,c(1,2)]
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
			union_edges <- c()
			pcancer <- c()
			gained_edges <- c()
			lost_edges <- c()
			gained_edges_num <- c()
			lost_edges_num <- c()
			##-----------------

			for(cn in 1:length(cancer_type)){

				ctype <- cancer_type[cn]
				gain_net <- net_file
				loss_net <- net_file

				tempf <- all_expr_data[[cn]]
				tempf <- tempf[tempf$EXON %in% nodes, ]

				exons <- tempf$EXON
				tempf1 <- tempf[,EXON:=NULL]
				tempf1[tempf1 < cpm_threshold[thr] ] <- 0 ## <-- variable
				tempf1 <- as.data.frame(tempf1)

				##---filter sample names to only those that have both control and condition
				tempids <- all_manifest[[cn]]
				sample_names_all <- unlist(lapply(strsplit(colnames(tempf1), '[_]'), '[[', 1))
				whc <- which(sample_names_all %in% tempids$nid)
				tempf1 <- tempf1[,whc]

				temp_brk <- unlist(lapply(strsplit(colnames(tempf1), '[_]'), '[[', 2))
				sample_names <- unlist(lapply(strsplit(colnames(tempf1), '[_]'), '[[', 1))
				condition_pos <- which( (temp_brk == '1') | (temp_brk == '01') )
				control_pos <- which(temp_brk == '11')

				control_graph <- list()
				condition_graph <- list()
				intersect_graph <- list()
				union_graph <- list()
				lost_graph <- list()
				gained_graph <- list()

				for(k in 1:length(control_pos)){

					control_file <- tempf1[,control_pos[k]]
					control_wh <- which(control_file != 0)

					control_nodes <- exons[control_wh]
					wh1 <- which(net_file$exon1 %in% control_nodes)
					wh2 <- which(net_file$exon2 %in% control_nodes)
					wh <- intersect(wh1, wh2)
					control_net <- net_file[wh, ]
					control_ig <- igraph::graph_from_data_frame(control_net, directed=FALSE)
					control_graph[[k]] <- control_ig

					condition_file <- tempf1[,condition_pos[k]]
					condition_wh <- which(condition_file != 0)
					condition_nodes <- exons[condition_wh]
					wh1 <- which(net_file$exon1 %in% condition_nodes)
					wh2 <- which(net_file$exon2 %in% condition_nodes)
					wh <- intersect(wh1, wh2)
					condition_net <- net_file[wh, ]
					condition_ig <- igraph::graph_from_data_frame(condition_net, directed=FALSE)
					condition_graph[[k]] <- condition_ig

					intersect_graph[[k]] <- igraph::intersection(control_ig, condition_ig)
					union_graph[[k]] <- igraph::union(control_ig, condition_ig)

					lle <- igraph::difference(control_ig, condition_ig, byname=TRUE)
					lost_graph[[k]] <- lle
					lled <- as.data.frame(igraph::get.edgelist(lle))
					colnames(lled) <- c('exon1', 'exon2')
					loss_net <- rbind(loss_net, lled)

					gge <- igraph::difference(condition_ig, control_ig, byname=TRUE)
					gained_graph[[k]] <- gge
					gged <- as.data.frame(igraph::get.edgelist(gge))
					colnames(gged) <- c('exon1', 'exon2')
					gain_net <- rbind(gain_net, gged)

				}

				##--- for plot -------
				control_edges <- c(control_edges, unlist(lapply(control_graph, function(x) igraph::ecount(x))))
				condition_edges <- c(condition_edges, unlist(lapply(condition_graph, function(x) igraph::ecount(x))))
				intersect_edges <- c(intersect_edges, unlist(lapply(intersect_graph, function(x) igraph::ecount(x))))
				union_edges  <- c(union_edges, unlist(lapply(union_graph, function(x) igraph::ecount(x))))
				pcancer <- c(pcancer, rep(cancer_type[cn], length(control_pos)))
				g_jac <- unlist(lapply(gained_graph, function(x) igraph::ecount(x)))/unlist(lapply(control_graph, function(x) igraph::ecount(x)))
				l_jac <- unlist(lapply(lost_graph, function(x) igraph::ecount(x)))/unlist(lapply(control_graph, function(x) igraph::ecount(x)))
				gained_edges <- c(gained_edges, g_jac)
				lost_edges <- c(lost_edges, l_jac)

				gained_edges_num <- c(gained_edges_num, unlist(lapply(gained_graph, function(x) igraph::ecount(x))))
				lost_edges_num <- c(lost_edges_num, unlist(lapply(lost_graph, function(x) igraph::ecount(x))))
				##-----------------

				all_gain_net <- igraph::graph_from_data_frame(gain_net, directed=FALSE)
				igraph::E(all_gain_net)$weight <- 1
				all_gain_net1 <- igraph::simplify(all_gain_net, edge.attr.comb=list(weight="sum"))
				all_gained <- igraph::as_data_frame(all_gain_net1)
				all_gained <- all_gained[order(-all_gained$weight), ]
				all_gained$weight <- all_gained$weight-1 ##-- subtracting 1 to denote that edge was not present in any patient
				all_gained$patient <- all_gained$weight
				all_gained$weight <- all_gained$weight/length(control_pos) ##-- fraction of patients in which the edge is gained/lost

				all_loss_net <- igraph::graph_from_data_frame(loss_net, directed=FALSE)
				igraph::E(all_loss_net)$weight <- 1
				all_loss_net1 <- igraph::simplify(all_loss_net, edge.attr.comb=list(weight="sum"))
				all_lost <- igraph::as_data_frame(all_loss_net1)
				all_lost <- all_lost[order(-all_lost$weight), ]
				all_lost$weight <- all_lost$weight-1 ##-- subtracting 1 to denote that edge was not present in any patient
				all_lost$patient <- all_lost$weight
				all_lost$weight <- all_lost$weight/length(control_pos) ##-- fraction of patients in which the edge is gained/lost

				data.table::fwrite(all_gained, paste0(out_dir1,'/',ctype,'_gained.txt'), row.names=FALSE, quote=FALSE, sep='\t')
				data.table::fwrite(all_lost, paste0(out_dir1,'/',ctype,'_lost.txt'), row.names=FALSE, quote=FALSE, sep='\t')
			}

			##-- save list object for later use
			slist <- list(pcancer, gained_edges, lost_edges)
			saveRDS(slist, paste0('../data/gl_edges_',level,'_',qq,'_',cpm_threshold[thr],'.RDS'))

			slist2 <- list(pcancer, gained_edges_num, lost_edges_num)
			saveRDS(slist2, paste0('../data/gl_edges_',level,'_',qq,'_',cpm_threshold[thr],'_num.RDS'))


		}

		cat('Network', qq, 'of', length(allnets),'done\n')

}





