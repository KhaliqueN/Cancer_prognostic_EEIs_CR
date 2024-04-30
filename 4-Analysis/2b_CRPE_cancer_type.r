##############################################################################################--
# Purpose: CRPE correlations with control, condition, perturbed, and number of patient samples
##############################################################################################--

rm(list=ls())
library(dplyr)
library(data.table)
library(ggplot2)
library(GenomicDataCommons)
source("eein_cancer_util.r")

cancer_type <- c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA')
cpm_threshold <- 0.5
pval_thres <- 0.05

get_num_edges <- function(tempf_control, tempf_condition, exons, net_file){

	 control_graph <- list()
    condition_graph <- list()
    gained_graph <- list()
    lost_graph <- list()
    perturbed_edges <-  list()
    perturbed_graph <-  list()

    gain_net <- net_file
    loss_net <- net_file

    control_names <- unlist(lapply(strsplit(colnames(tempf_control), '[_]'), '[[', 1)) 
    condition_names <- unlist(lapply(strsplit(colnames(tempf_condition), '[_]'), '[[', 1)) 

    for(k in 1:length(tempf_control)){

        control_file <- tempf_control[,k]
        control_wh <- which(control_file != 0)
        control_nodes <- exons[control_wh]
        wh1 <- which(net_file$exon1 %in% control_nodes)
        wh2 <- which(net_file$exon2 %in% control_nodes)
        wh <- intersect(wh1, wh2)
        control_net <- net_file[wh, ]
        control_ig <- igraph::graph_from_data_frame(control_net, directed=FALSE)
        control_graph[[k]] <- control_ig

        whc <- which(condition_names == control_names[k])
        condition_file <- tempf_condition[,whc]
        condition_wh <- which(condition_file != 0)
        condition_nodes <- exons[condition_wh]
        wh1 <- which(net_file$exon1 %in% condition_nodes)
        wh2 <- which(net_file$exon2 %in% condition_nodes)
        wh <- intersect(wh1, wh2)
        condition_net <- net_file[wh, ]
        condition_ig <- igraph::graph_from_data_frame(condition_net, directed=FALSE)
        condition_graph[[k]] <- condition_ig

        lle <- igraph::difference(control_ig, condition_ig, byname=TRUE)
        lost_graph[[k]] <- lle
        if(igraph::ecount(lle) != 0){
            lled <- as.data.frame(igraph::get.edgelist(lle))
            colnames(lled) <- c('exon1', 'exon2')
            loss_net <- rbind(loss_net, lled)
        }

        gge <- igraph::difference(condition_ig, control_ig, byname=TRUE)
        gained_graph[[k]] <- gge
        if(igraph::ecount(lle) != 0){
            gged <- as.data.frame(igraph::get.edgelist(gge))
            colnames(gged) <- c('exon1', 'exon2')
            gain_net <- rbind(gain_net, gged)
        }

        perturbed_edges[[k]] <- igraph::ecount(igraph::union(lle, gge))
        perturbed_graph[[k]] <- igraph::union(lle, gge)
    }

    return(list(control_graph, condition_graph, perturbed_graph))
}

get_pert_edges <- function(tempf_control, tempf_condition, exons, net_file, allsurv){

	 control_graph <- list()
	 condition_graph <- list()
	 gained_graph <- list()
	 lost_graph <- list()
	 perturbed_graph <-  list()
	 cancer_relevant_edges <-  list()
	 cancer_relevant_graph <-  list()

	 gain_net <- net_file
	 loss_net <- net_file

	 control_names <- unlist(lapply(strsplit(colnames(tempf_control), '[_]'), '[[', 1)) 
	 condition_names <- unlist(lapply(strsplit(colnames(tempf_condition), '[_]'), '[[', 1)) 

	 for(k in 1:length(tempf_control)){

	     control_file <- tempf_control[,k]
	     control_wh <- which(control_file != 0)
	     control_nodes <- exons[control_wh]
	     wh1 <- which(net_file$exon1 %in% control_nodes)
	     wh2 <- which(net_file$exon2 %in% control_nodes)
	     wh <- intersect(wh1, wh2)
	     control_net <- net_file[wh, ]
	     control_ig <- igraph::graph_from_data_frame(control_net, directed=FALSE)
	     control_graph[[k]] <- igraph::ecount(control_ig)

	     whc <- which(condition_names == control_names[k])
	     condition_file <- tempf_condition[,whc]
	     condition_wh <- which(condition_file != 0)
	     condition_nodes <- exons[condition_wh]
	     wh1 <- which(net_file$exon1 %in% condition_nodes)
	     wh2 <- which(net_file$exon2 %in% condition_nodes)
	     wh <- intersect(wh1, wh2)
	     condition_net <- net_file[wh, ]
	     condition_ig <- igraph::graph_from_data_frame(condition_net, directed=FALSE)
	     condition_graph[[k]] <- igraph::ecount(condition_ig)

	     lle <- igraph::difference(control_ig, condition_ig, byname=TRUE)
	     lost_graph[[k]] <- lle
	     if(igraph::ecount(lle) != 0){
	         lled <- as.data.frame(igraph::get.edgelist(lle))
	         colnames(lled) <- c('exon1', 'exon2')
	         loss_net <- rbind(loss_net, lled)
	     }

	     gge <- igraph::difference(condition_ig, control_ig, byname=TRUE)
	     gained_graph[[k]] <- gge
	     if(igraph::ecount(lle) != 0){
	         gged <- as.data.frame(igraph::get.edgelist(gge))
	         colnames(gged) <- c('exon1', 'exon2')
	         gain_net <- rbind(gain_net, gged)
	     }

	     pert_graph <- igraph::union(lle, gge)
	     perturbed_graph[[k]] <- pert_graph
		  cancer_relevant_edges[[k]] <- igraph::ecount(igraph::intersection(pert_graph, allsurv))
		  cancer_relevant_graph[[k]] <- igraph::intersection(pert_graph, allsurv)

	 }

	 return(list(control_graph, condition_graph, cancer_relevant_edges, perturbed_graph, cancer_relevant_graph))
}


###--- NETLOW ---------------------------------------------------------------
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


contact_net <- data.table::fread(paste0('../data/CONTACT_networks/CONTACT_net_6_1.txt'))
contact_net <- igraph::graph_from_data_frame(contact_net[,c(3,4)], directed=FALSE)
pisa_net <- data.table::fread(paste0('../data/PISA_networks_filtered/PISA_EEIN_0.5.txt'))
pisa_net <- igraph::graph_from_data_frame(pisa_net[,c(1,2)], directed=FALSE)
eppic_net <- data.table::fread('../data/EPPIC_EEIN_filtered.txt')
eppic_net <- igraph::graph_from_data_frame(eppic_net[,c(1,2)], directed=FALSE)
net_file <- igraph::as_data_frame(igraph::union(igraph::union(contact_net, pisa_net), eppic_net))
colnames(net_file) <- c('exon1', 'exon2')
nodes <- union(net_file[[1]], net_file[[2]])

##-- for plot -----
control_edges <- c()
condition_edges <- c()
perturbed_edges <- c()
survival_edges <- c()
pcancer <- c()
thresholds <- c()
lsamples <- c()
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
	tempf_control[tempf_control < cpm_threshold ] <- 0 
	tempf_condition[tempf_condition < cpm_threshold ] <- 0 

	lsamples <- c(lsamples, length(all_survival$submitter_id))
	##----for real expression--------------------------------
	out_real <- get_num_edges(tempf_control, tempf_condition, exons, net_file)
	temp_graph <- out_real[[1]][[1]]
	for(i in 2:length(out_real[[1]])){
		temp_graph <- igraph::union(temp_graph, out_real[[1]][[i]])
	}
	control_edges <- c(control_edges, igraph::ecount(temp_graph))

	temp_graph <- out_real[[2]][[1]]
	for(i in 2:length(out_real[[2]])){
		temp_graph <- igraph::union(temp_graph, out_real[[2]][[i]])
	}
	condition_edges <- c(condition_edges, igraph::ecount(temp_graph))

	temp_graph <- out_real[[3]][[1]]
	for(i in 2:length(out_real[[3]])){
		temp_graph <- igraph::union(temp_graph, out_real[[3]][[i]])
	}
	perturbed_edges <- c(perturbed_edges, igraph::ecount(temp_graph))

	pcancer <- c(pcancer, cancer_type[cn])
	thresholds <- c(thresholds, cpm_threshold)

	##--- cancer-relevant edges ---
	pisa_g <- data.table::fread(paste0('../results/PISA_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Gained_Surv_.txt'))
	pisa_g <- pisa_g[pisa_g$pval <= pval_thres, ][,c(1,2)]
	pisa_l <- data.table::fread(paste0('../results/PISA_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Lost_Surv_.txt'))
	pisa_l <- pisa_l[pisa_l$pval <= pval_thres, ][,c(1,2)]
	pisa <- igraph::graph_from_data_frame(rbind(pisa_g,pisa_l), directed=FALSE)

	eppic_g <- data.table::fread(paste0('../results/EPPIC_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Gained_Surv_.txt'))
	eppic_g <- eppic_g[eppic_g$pval <= pval_thres, ][,c(1,2)]
	eppic_l <- data.table::fread(paste0('../results/EPPIC_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Lost_Surv_.txt'))
	eppic_l <- eppic_l[eppic_l$pval <= pval_thres, ][,c(1,2)]
	eppic <- igraph::graph_from_data_frame(rbind(eppic_g,eppic_l), directed=FALSE)

	contact_g <- data.table::fread(paste0('../results/CONTACT_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Gained_Surv_.txt'))
	contact_g <- contact_g[contact_g$pval <= pval_thres, ][,c(1,2)]
	contact_l <- data.table::fread(paste0('../results/CONTACT_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Lost_Surv_.txt'))
	contact_l <- contact_l[contact_l$pval <= pval_thres, ][,c(1,2)]
	contact <- igraph::graph_from_data_frame(rbind(contact_g,contact_l), directed=FALSE)

	allsurv <- igraph::union(igraph::union(pisa, eppic), contact)
	out_surv_real <- get_pert_edges(tempf_control, tempf_condition, exons, net_file, allsurv)

	survival_graph <- out_surv_real[[5]][[1]]
	for(i in 2:length(out_surv_real[[5]])){
		survival_graph <- igraph::union(survival_graph, out_surv_real[[5]][[i]])
	}

	survival_edges <- c(survival_edges, igraph::ecount(survival_graph))

}


pdata <- data.frame(Cancer=pcancer, Control=control_edges, Condition=condition_edges, 
	Perturbed=perturbed_edges, Survival=survival_edges, Threshold=thresholds, Samples=lsamples)


###----- correlation with number of edges ------------------------------------------------
coral <- cor.test(x=pdata$Control, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Control, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of control edges") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=16000,y=6000, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETLOW_correlation_control.png"),width=4.8, height=2.5, dpi=400)


coral <- cor.test(x=pdata$Condition, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Condition, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of condition edges") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=19000,y=6000, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETLOW_correlation_condition.png"),width=4.8, height=2.5, dpi=400)


coral <- cor.test(x=pdata$Perturbed, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Perturbed, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of perturbed edges") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=16000,y=6000, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETLOW_correlation_perturbed.png"),width=4.8, height=2.5, dpi=400)



coral <- cor.test(x=pdata$Samples, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Samples, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of patient samples") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=600,y=6000, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETLOW_correlation_samples.png"),width=4.8, height=2.5, dpi=400)



###----------------------------------------------------------------



###--- NETMEDIUM --------------
contact_net <- data.table::fread(paste0('../data/CONTACT_networks/CONTACT_net_6_1.txt'))
contact_net <- igraph::graph_from_data_frame(contact_net[,c(3,4)], directed=FALSE)
pisa_net <- data.table::fread(paste0('../data/PISA_networks_filtered/PISA_EEIN_0.5.txt'))
pisa_net <- igraph::graph_from_data_frame(pisa_net[,c(1,2)], directed=FALSE)
eppic_net <- data.table::fread('../data/EPPIC_EEIN_filtered.txt')
eppic_net <- igraph::graph_from_data_frame(eppic_net[,c(1,2)], directed=FALSE)

net_file <- igraph::union(igraph::union(igraph::intersection(pisa_net,eppic_net,keep.all.vertices=FALSE),
		igraph::intersection(contact_net,eppic_net,keep.all.vertices=FALSE)), 
		igraph::intersection(contact_net,pisa_net,keep.all.vertices=FALSE))

net_file <- igraph::as_data_frame(net_file)
colnames(net_file) <- c('exon1', 'exon2')
nodes <- union(net_file[[1]], net_file[[2]])

##-- for plot -----
control_edges <- c()
condition_edges <- c()
perturbed_edges <- c()
survival_edges <- c()
pcancer <- c()
thresholds <- c()
lsamples <- c()
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
	tempf_control[tempf_control < cpm_threshold ] <- 0 
	tempf_condition[tempf_condition < cpm_threshold ] <- 0 

	lsamples <- c(lsamples, length(all_survival$submitter_id))
	##----for real expression--------------------------------------------------
	out_real <- get_num_edges(tempf_control, tempf_condition, exons, net_file)
	temp_graph <- out_real[[1]][[1]]
	for(i in 2:length(out_real[[1]])){
		temp_graph <- igraph::union(temp_graph, out_real[[1]][[i]])
	}
	control_edges <- c(control_edges, igraph::ecount(temp_graph))

	temp_graph <- out_real[[2]][[1]]
	for(i in 2:length(out_real[[2]])){
		temp_graph <- igraph::union(temp_graph, out_real[[2]][[i]])
	}
	condition_edges <- c(condition_edges, igraph::ecount(temp_graph))

	temp_graph <- out_real[[3]][[1]]
	for(i in 2:length(out_real[[3]])){
		temp_graph <- igraph::union(temp_graph, out_real[[3]][[i]])
	}
	perturbed_edges <- c(perturbed_edges, igraph::ecount(temp_graph))

	pcancer <- c(pcancer, cancer_type[cn])
	thresholds <- c(thresholds, cpm_threshold)

	##--- cancer-relevant edges ---
	pisa_g <- data.table::fread(paste0('../results/PISA_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Gained_Surv_.txt'))
	pisa_g <- pisa_g[pisa_g$pval <= pval_thres, ][,c(1,2)]
	pisa_l <- data.table::fread(paste0('../results/PISA_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Lost_Surv_.txt'))
	pisa_l <- pisa_l[pisa_l$pval <= pval_thres, ][,c(1,2)]
	pisa <- igraph::graph_from_data_frame(rbind(pisa_g,pisa_l), directed=FALSE)

	eppic_g <- data.table::fread(paste0('../results/EPPIC_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Gained_Surv_.txt'))
	eppic_g <- eppic_g[eppic_g$pval <= pval_thres, ][,c(1,2)]
	eppic_l <- data.table::fread(paste0('../results/EPPIC_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Lost_Surv_.txt'))
	eppic_l <- eppic_l[eppic_l$pval <= pval_thres, ][,c(1,2)]
	eppic <- igraph::graph_from_data_frame(rbind(eppic_g,eppic_l), directed=FALSE)

	contact_g <- data.table::fread(paste0('../results/CONTACT_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Gained_Surv_.txt'))
	contact_g <- contact_g[contact_g$pval <= pval_thres, ][,c(1,2)]
	contact_l <- data.table::fread(paste0('../results/CONTACT_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Lost_Surv_.txt'))
	contact_l <- contact_l[contact_l$pval <= pval_thres, ][,c(1,2)]
	contact <- igraph::graph_from_data_frame(rbind(contact_g,contact_l), directed=FALSE)

	allsurv <- igraph::union(igraph::union(pisa, eppic), contact)
	out_surv_real <- get_pert_edges(tempf_control, tempf_condition, exons, net_file, allsurv)

	survival_graph <- out_surv_real[[5]][[1]]
	for(i in 2:length(out_surv_real[[5]])){
		survival_graph <- igraph::union(survival_graph, out_surv_real[[5]][[i]])
	}

	survival_edges <- c(survival_edges, igraph::ecount(survival_graph))

}



pdata <- data.frame(Cancer=pcancer, Control=control_edges, Condition=condition_edges, 
	Perturbed=perturbed_edges, Survival=survival_edges, Threshold=thresholds, Samples=lsamples)


###----- correlation with number of edges ------------------------------------------------
coral <- cor.test(x=pdata$Control, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Control, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of control edges") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=16000,y=6000, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETMEDIUM_correlation_control.png"),width=4.8, height=2.5, dpi=400)


coral <- cor.test(x=pdata$Condition, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Condition, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of condition edges") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=19000,y=6000, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETMEDIUM_correlation_condition.png"),width=4.8, height=2.5, dpi=400)


coral <- cor.test(x=pdata$Perturbed, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Perturbed, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of perturbed edges") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=16000,y=6000, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETMEDIUM_correlation_perturbed.png"),width=4.8, height=2.5, dpi=400)



coral <- cor.test(x=pdata$Samples, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Samples, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of patient samples") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=600,y=6000, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETMEDIUM_correlation_samples.png"),width=4.8, height=2.5, dpi=400)





###--- NETHIGH -------------------------------------------------------------------------------------
contact_net <- data.table::fread(paste0('../data/CONTACT_networks/CONTACT_net_6_1.txt'))
contact_net <- igraph::graph_from_data_frame(contact_net[,c(3,4)], directed=FALSE)
pisa_net <- data.table::fread(paste0('../data/PISA_networks_filtered/PISA_EEIN_0.5.txt'))
pisa_net <- igraph::graph_from_data_frame(pisa_net[,c(1,2)], directed=FALSE)
eppic_net <- data.table::fread('../data/EPPIC_EEIN_filtered.txt')
eppic_net <- igraph::graph_from_data_frame(eppic_net[,c(1,2)], directed=FALSE)

net_file <- igraph::intersection(igraph::intersection(pisa_net,eppic_net,keep.all.vertices=FALSE),
	contact_net,keep.all.vertices=FALSE) 

net_file <- igraph::as_data_frame(net_file)
colnames(net_file) <- c('exon1', 'exon2')
nodes <- union(net_file[[1]], net_file[[2]])

##-- for plot -----
control_edges <- c()
condition_edges <- c()
perturbed_edges <- c()
survival_edges <- c()
pcancer <- c()
thresholds <- c()
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
	tempf_control[tempf_control < cpm_threshold ] <- 0 
	tempf_condition[tempf_condition < cpm_threshold ] <- 0 

	lsamples <- c(lsamples, length(all_survival$submitter_id))
	##----for real expression--------------------------------
	out_real <- get_num_edges(tempf_control, tempf_condition, exons, net_file)
	temp_graph <- out_real[[1]][[1]]
	for(i in 2:length(out_real[[1]])){
		temp_graph <- igraph::union(temp_graph, out_real[[1]][[i]])
	}
	control_edges <- c(control_edges, igraph::ecount(temp_graph))

	temp_graph <- out_real[[2]][[1]]
	for(i in 2:length(out_real[[2]])){
		temp_graph <- igraph::union(temp_graph, out_real[[2]][[i]])
	}
	condition_edges <- c(condition_edges, igraph::ecount(temp_graph))

	temp_graph <- out_real[[3]][[1]]
	for(i in 2:length(out_real[[3]])){
		temp_graph <- igraph::union(temp_graph, out_real[[3]][[i]])
	}
	perturbed_edges <- c(perturbed_edges, igraph::ecount(temp_graph))

	pcancer <- c(pcancer, cancer_type[cn])
	thresholds <- c(thresholds, cpm_threshold)

	##--- cancer-relevant edges ---
	pisa_g <- data.table::fread(paste0('../results/PISA_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Gained_Surv_.txt'))
	pisa_g <- pisa_g[pisa_g$pval <= pval_thres, ][,c(1,2)]
	pisa_l <- data.table::fread(paste0('../results/PISA_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Lost_Surv_.txt'))
	pisa_l <- pisa_l[pisa_l$pval <= pval_thres, ][,c(1,2)]
	pisa <- igraph::graph_from_data_frame(rbind(pisa_g,pisa_l), directed=FALSE)

	eppic_g <- data.table::fread(paste0('../results/EPPIC_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Gained_Surv_.txt'))
	eppic_g <- eppic_g[eppic_g$pval <= pval_thres, ][,c(1,2)]
	eppic_l <- data.table::fread(paste0('../results/EPPIC_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Lost_Surv_.txt'))
	eppic_l <- eppic_l[eppic_l$pval <= pval_thres, ][,c(1,2)]
	eppic <- igraph::graph_from_data_frame(rbind(eppic_g,eppic_l), directed=FALSE)

	contact_g <- data.table::fread(paste0('../results/CONTACT_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Gained_Surv_.txt'))
	contact_g <- contact_g[contact_g$pval <= pval_thres, ][,c(1,2)]
	contact_l <- data.table::fread(paste0('../results/CONTACT_survival/threshold_',cpm_threshold,'/',cancer_type[cn],'_Lost_Surv_.txt'))
	contact_l <- contact_l[contact_l$pval <= pval_thres, ][,c(1,2)]
	contact <- igraph::graph_from_data_frame(rbind(contact_g,contact_l), directed=FALSE)

	allsurv <- igraph::union(igraph::union(pisa, eppic), contact)
	out_surv_real <- get_pert_edges(tempf_control, tempf_condition, exons, net_file, allsurv)

	survival_graph <- out_surv_real[[5]][[1]]
	for(i in 2:length(out_surv_real[[5]])){
		survival_graph <- igraph::union(survival_graph, out_surv_real[[5]][[i]])
	}

	survival_edges <- c(survival_edges, igraph::ecount(survival_graph))

}


pdata <- data.frame(Cancer=pcancer, Control=control_edges, Condition=condition_edges, 
	Perturbed=perturbed_edges, Survival=survival_edges, Threshold=thresholds, Samples=lsamples)


###----- correlation with number of edges ------------------------------------------------
coral <- cor.test(x=pdata$Control, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Control, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of control edges") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=3500,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETHIGH_correlation_control.png"),width=4.8, height=2.5, dpi=400)


coral <- cor.test(x=pdata$Condition, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Condition, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of condition edges") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=4500,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETHIGH_correlation_condition.png"),width=4.8, height=2.5, dpi=400)


coral <- cor.test(x=pdata$Perturbed, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Perturbed, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of perturbed edges") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=4000,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETHIGH_correlation_perturbed.png"),width=4.8, height=2.5, dpi=400)



coral <- cor.test(x=pdata$Samples, y=pdata$Survival, method = 'spearman')
p <- ggplot(pdata, aes(Samples, Survival, color=Cancer)) + 
geom_point(size=2)+
theme(legend.text=element_text(size=12))
basesize <- 12
p <- p + theme_bw(base_size = basesize * 0.8) +
scale_x_continuous(name="# of patient samples") + 
scale_y_continuous(name="# of cancer relevant \nperturbed edges (CRPEs)") +
geom_text(aes(x=600,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,color='black',show.legend = FALSE)+
scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
guides(color=guide_legend(title="Cancer type",ncol=3))
ggsave(p,filename=paste0("../results/Final_results/NETHIGH_correlation_samples.png"),width=4.8, height=2.5, dpi=400)



