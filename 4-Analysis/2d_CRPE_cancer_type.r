 ##############################################################################################
# Purpose: Figure for biological properties 
##############################################################################################

rm(list=ls())
library(data.table)
library('biomaRt')
library(Rcpp)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(GenomicDataCommons)
source("eein_cancer_util.r")

cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))
out_dir <- '../results/Final_results'
cpm_threshold <- 0.5
pval_thres <- 0.05
allnets <- gtools::mixedsort(list.files('../data/Final_networks/level1',full.names=TRUE))
net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')
# allnets <- gtools::mixedsort(list.files('../data/Final_networks/level1',full.names=TRUE))

gnet <- data.table::fread(paste0('../results/PISA_survival/PISA_net_final_',cpm_threshold,'.txt'), header=FALSE)
gnet <- mapProtein(gnet[[1]], gnet[[2]], data.table::fread('../data/PISA_networks_filtered/PISA_EEIN_0.5.txt'))
gnet <- gnet[,c(1,2,5,6)]

enet <- data.table::fread(paste0('../results/EPPIC_survival/EPPIC_net_final_',cpm_threshold,'.txt'), header=FALSE)
enet <- mapProtein(enet[[1]], enet[[2]], data.table::fread('../data/EPPIC_EEIN_filtered.txt'))
enet <- enet[,c(1,2,5,6)]

anet <- data.table::fread(paste0('../results/CONTACT_survival/CONTACT_net_final_',cpm_threshold,'.txt'), header=FALSE)
xxp <- data.table::fread('../data/CONTACT_networks/CONTACT_net_6_1.txt')
xxp$p1 <- unlist(lapply(strsplit(xxp[[1]],'[_]'),'[[',1))
xxp$p2 <- unlist(lapply(strsplit(xxp[[2]],'[_]'),'[[',1))
xxp <- xxp[, -c(1,2)]
anet <- mapProtein(anet[[1]], anet[[2]], xxp)
anet <- anet[, c(1,2,10,11)]
colnames(anet) <- c('exon1','exon2','protein1','protein2')
aq <- rbind(gnet, anet)
unet <- rbind(aq, enet)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


for(qq in 1:length(allnets)){

    in_dir <- paste0('../results/Final_survival/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/')
    net_file <- data.table::fread(allnets[qq], header=FALSE)

    all_pert <- list()

    for(k in 1:length(cancer_type)){
        c_type <- cancer_type[k]
        gained <- data.table::fread(paste0(in_dir,c_type,'_Gained_Surv_.txt'))
        gained_ig <- igraph::graph_from_data_frame(gained[gained$pval <= pval_thres, ][,c(1,2)], directed=FALSE)
        lost <- data.table::fread(paste0(in_dir,c_type,'_Lost_Surv_.txt'))
        lost_ig <- igraph::graph_from_data_frame(lost[lost$pval <= pval_thres, ][,c(1,2)], directed=FALSE)
        all_pert[[k]] <- igraph::union(gained_ig, lost_ig)
    }

   
    ###------------------------------------------------------------------------------

    ###--- KEGG term computations ----------------------------------------------------------------------
    net_file1 <- mapProtein(net_file[[1]], net_file[[2]], unet)
    all_kegg <- data.table::fread('../data/all_human_pathways.txt', sep='\t')
    colnames(all_kegg) <- c('ensembl_gene_id','uniprot_id','entrez_id','go_id','pathways_name')
    all_kegg <- all_kegg[all_kegg$uniprot_id %in% union(net_file1[[3]], net_file1[[4]]), ]
    # ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

    attrs <- c('ensembl_gene_id', 'ensembl_exon_id')
    alld <- getBM(attributes=attrs, filters='ensembl_gene_id', values=all_kegg[[1]], mart=ensembl) ##--> CHECK THI AND GO BGSIZE
    path_name <- rep('',length(alld[[1]]))
    path_id <- rep('',length(alld[[1]]))
    for(k in 1:length(all_kegg[[1]])){
        wh <- which(alld$ensembl_gene_id == all_kegg[[1]][k])
        if(length(wh > 0)){
            path_name[wh] <- all_kegg$pathways_name[k]
            path_id[wh] <- all_kegg$go_id[k]
        }
    }

    alld$pathways_name <- path_name
    alld$path_id <- path_id
    alldd <- alld[alld$path_id != '', ]
    ## only keep the exons in my background network -------------------------------------
    alldd <- alldd[alldd$ensembl_exon_id %in% union(net_file[[1]], net_file[[2]]), ] 

    all_path <- plyr::count(alldd$path_id)
    wh_path <- unique(all_path[which(all_path$freq > 2), ]$x)
    allkegg <- alldd[alldd$path_id %in% wh_path, ]
    colnames(allkegg) <- c('ensembl_gene_id','ensembl_exon_id','name_1006','go_id')
    kg_bgSize <- length(unique(allkegg$ensembl_exon_id))

    allKE1 <- goEnrich(allkegg, kg_bgSize, all_pert, cancer_type)
    allKE1$qval <- p.adjust(allKE1$pval, 'fdr')
    allKE <- allKE1[allKE1$qval <= pval_thres, ]

    #-- save excel file ---
    wb1 <- openxlsx::createWorkbook(paste0(out_dir,'/Supplementary_Table_S3_',net_type[qq],'.xlsx'))
    openxlsx::saveWorkbook(wb1, paste0(out_dir,'/Supplementary_Table_S3_',net_type[qq],'.xlsx'), overwrite=T)

    for(k in 1:length(cancer_type)){

        temp <- allKE[allKE$temp_cancer == cancer_type[k], ]
        temp1a <- temp[order(temp$qval), ]

        wb <- openxlsx::loadWorkbook(paste0(out_dir,'/Supplementary_Table_S3_',net_type[qq],'.xlsx'))
        openxlsx::addWorksheet(wb, sheetName = cancer_type[k])
        openxlsx::writeData(wb, sheet = cancer_type[k], temp1a)
        openxlsx::saveWorkbook(wb, paste0(out_dir,'/Supplementary_Table_S3_',net_type[qq],'.xlsx'), overwrite = T)

    }

    

    ####----- Pairwise overlaps of enriched KEGG terms ---------------------------------------------------------------------------------
    c1 <- c()
    c2 <- c()
    value <- c()
    pvalue <- c()
    pvalueu <- c()
    loop1 <- length(cancer_type)-1
    loop2 <- length(cancer_type)

    for(k in 1:loop1){

        tempx1 <- allKE[allKE$temp_cancer == cancer_type[k], ]
        m  <- k+1
        for(j in m:loop2){

            # if(k != j){
            tempx2 <- allKE[allKE$temp_cancer == cancer_type[j], ]
            e1 <- unique(tempx1$GOterm)
            e2 <- unique(tempx2$GOterm)

            ovr <- intersect(e1, e2)
            ur <- union(e1, e2)
            pv <- phyper(length(ovr)-1, length(e1), kg_bgSize-length(e1), length(e2), lower.tail=FALSE)
            pvu <- phyper(length(ovr), length(e1), kg_bgSize-length(e1), length(e2), lower.tail=TRUE)

            if(length(e1) < length(e2)){
                tempo <- length(ovr)/length(e1)
            }else{
                tempo <- length(ovr)/length(e2)
            }
            value <- c(value, length(ovr))

            pvalue <- c(pvalue, pv)
            pvalueu <- c(pvalueu, pvu)
            c1 <- c(c1, paste0(cancer_type[k],' (',length(e1),')'))
            c2 <- c(c2, paste0(cancer_type[j],' (',length(e2),')'))
        # }
        }
    }

    qval <- signif(p.adjust(pvalue,'fdr'),2)
    qvalu <- signif(p.adjust(pvalueu,'fdr'),2)
    qvalx <- rep('Overlap as \nexpected by chance',length(pvalue))
    for(i in 1:length(pvalue)){
        if(qval[i] <= pval_thres){
            qvalx[i] <- 'Higher overlap than\nexpected by chance'
        }
        if(qvalu[i] <= pval_thres){
            qvalx[i] <- 'Lower overlap than\nexpected by chance'
        }
    }

    pdata <- data.frame(c1=c1, c2=c2, val=value, pval=qvalx)

    # pdata$pvalue <- cut(pdata$pval, breaks = c(0,0.01,0.05,1), include.lowest=TRUE)

    cols <- rev(c('#377eb8','#e41a1c'))#rev(brewer.pal(3,"Spectral"))

    p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pval),colour = "white")+
      theme(legend.text=element_text(size=8))+scale_fill_manual(values=cols,drop=FALSE)
    basesize <- 8
    p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type", y = "Cancer type") +
      scale_y_discrete() +
      scale_x_discrete()+coord_fixed()+
      guides(fill=guide_legend(title="Category"))+
      geom_text(data=pdata,aes(y=c2,label=val),size=2.5)+
      theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0.5, colour = "black"),
        axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0.5, colour = "black"))#+
      # guides(fill='none')
    ggsave(p,filename=paste0(out_dir,'/KEGG_overlap_',net_type[qq],'.png'),width=4.5, height=3, dpi=300)



    ###---- unique KEGG terms -------
    unique_KE <- list()
    uKEs <- c()
    for(k in 1:length(cancer_type)){
        counter <- k
        temp_uni <- unique(allKE[allKE$temp_cancer == cancer_type[k], ]$GOtermn)

        for(j in 1:length(cancer_type)){
            if(counter != j){
                temp_uni <- setdiff(temp_uni, allKE[allKE$temp_cancer == cancer_type[j], ]$GOtermn)
            }
        }
        unique_KE[[k]] <- temp_uni
        uKEs <- c(uKEs, length(temp_uni))
    }

    tempg <- data.frame(cancer=cancer_type, count=uKEs)
    p <- ggplot(tempg, aes(cancer, count)) + 
    geom_bar(stat="identity",position=position_dodge())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of uniquely enriched \nKEGG pathways", limits=c(0,(max(tempg$count))+2)) +
    geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=0, size=3)+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_unique_KEGG.png"),width=3.5, height=3, dpi=400)


    ##--- common KEGG---
    common_KE <- plyr::count(allKE$GOtermn)
    common_KE <- common_KE[order(-common_KE$freq),]
    common_KE <- common_KE[common_KE$freq > 1, ]

    ## q-values of common KEs
    tempq <- allKE[allKE$GOtermn %in% common_KE[[1]], ][,c(4,5)]

    tempqq <- tempq%>%
    group_by(GOtermn)%>% 
    summarise(Mean=mean(qval), Max=max(qval), Min=min(qval), Median=median(qval), Std=sd(qval))

    tempqq <- as.data.frame(tempqq)
    tempqq <- tempqq[order(tempqq$Median), ]


    ## q-values of unique KEs
    tempq <- allKE[allKE$GOtermn %in% unlist(unique_KE), ][,c(1,4,5)]
    tempqq <- tempq%>%
    group_by(temp_cancer)%>% 
    summarise(Min=min(qval))

