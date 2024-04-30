##############################################################################################
# Purpose: KEGG enrichment of biomarkers
##############################################################################################

rm(list=ls())
library(data.table)
library(ggplot2)
source("eein_cancer_util.r")
library(GenomicDataCommons)
library('biomaRt')

cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))
out_dir <- '../results/Final_results'
dir.create(out_dir, recursive=TRUE)
cpm_threshold <- 0.5
pval_thres <- 0.05
allnets <- gtools::mixedsort(list.files('../data/Final_networks/level1',full.names=TRUE))
net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')

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


for(qq in 3:length(allnets)){

    biocrpes <- readRDS(paste0(net_type[qq],"_biocrpe.Rds"))
    ##--- all bio crpes
    all_bio <- biocrpes[[1]]
    for(j in 2:length(biocrpes)){
        all_bio <- igraph::union(all_bio, biocrpes[[j]])
    }

    net_file <- data.table::fread(allnets[qq], header=FALSE)

    ###--- KEGG term computations ----------------------------------------------------------------------
    net_file1 <- mapProtein(net_file[[1]], net_file[[2]], unet)
    all_kegg <- data.table::fread('../data/all_human_pathways.txt', sep='\t')
    colnames(all_kegg) <- c('ensembl_gene_id','uniprot_id','entrez_id','go_id','pathways_name')
    all_kegg <- all_kegg[all_kegg$uniprot_id %in% union(net_file1[[3]], net_file1[[4]]), ]
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

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
    allKE1 <- goEnrich_uniset(allkegg, kg_bgSize, all_bio)
    allKE1$qval <- p.adjust(allKE1$pval, 'fdr')
    allKE <- allKE1[allKE1$qval <= pval_thres, ]




    ##----- GO-term enrichment of gained/lost edges ---------
    p1 <- net_file[[1]]
    p2 <- net_file[[2]]
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    attrs <- c('ensembl_exon_id', 'go_id', 'namespace_1003', 'name_1006')
    alld <- getBM(attributes=attrs, filters='ensembl_exon_id', values=union(p1,p2), mart=ensembl)
    allgo <- plyr::count(alld$go_id)
    wh_go <- unique(allgo[which(allgo$freq > 2), ]$x)
    alldd <- alld[alld$go_id %in% wh_go, ]

    ## all GO terms to consider
    all_go_bp <- alldd[alldd$namespace_1003 == 'biological_process', ]
    all_go_mf <- alldd[alldd$namespace_1003 == 'molecular_function', ]
    all_go_cc <- alldd[alldd$namespace_1003 == 'cellular_component', ]
    bp_bgSize <- length(unique(all_go_bp$ensembl_exon_id))
    mf_bgSize <- length(unique(all_go_mf$ensembl_exon_id))
    cc_bgSize <- length(unique(all_go_cc$ensembl_exon_id))

    goenrich_bp <- goEnrich_uniset(all_go_bp, bp_bgSize, all_bio)## BP
    goenrich_bp$flag <- rep('BP', length(goenrich_bp[[1]]))## BP
    goenrich_mf <- goEnrich_uniset(all_go_mf, mf_bgSize, all_bio)## MF
    goenrich_mf$flag <- rep('MF', length(goenrich_mf[[1]]))## MF
    goenrich_cc <- goEnrich_uniset(all_go_cc, cc_bgSize, all_bio)## CC
    goenrich_cc$flag <- rep('CC', length(goenrich_cc[[1]]))## CC
    allGO1 <- rbind(goenrich_bp, goenrich_mf)
    allGO <- rbind(allGO1, goenrich_cc)
    allGO$qval <- p.adjust(allGO$pval, 'fdr')
    allGO <- allGO[allGO$qval <= pval_thres, ]

    #-- save excel file ---
    # openxlsx::saveWorkbook(wb1, paste0(out_dir,'/GO_BioCRPE_',net_type[qq],'.xlsx'), overwrite=T)

    temp <- allGO
    temp1a <- temp[temp$flag == 'BP', ]
    temp1a <- temp1a[order(temp1a$qval), ]
    wb1 <- openxlsx::createWorkbook(paste0(out_dir,'/GO_BioCRPE_',net_type[qq],'.xlsx'))
    # wb <- openxlsx::loadWorkbook(paste0(out_dir,'/GO_BioCRPE_',net_type[qq],'.xlsx'))
    openxlsx::addWorksheet(wb1, sheetName = 'BP')
    openxlsx::writeData(wb1, sheet = 'BP', temp1a)
    openxlsx::saveWorkbook(wb1, paste0(out_dir,'/GO_BioCRPE_',net_type[qq],'.xlsx'), overwrite = T)

    temp1b <- temp[temp$flag == 'MF', ]
    temp1b <- temp1b[order(temp1b$qval), ]
    wb <- openxlsx::loadWorkbook(paste0(out_dir,'/GO_BioCRPE_',net_type[qq],'.xlsx'))
    openxlsx::addWorksheet(wb, sheetName = 'MF')
    openxlsx::writeData(wb, sheet = 'MF', temp1b)
    openxlsx::saveWorkbook(wb, paste0(out_dir,'/GO_BioCRPE_',net_type[qq],'.xlsx'), overwrite = T)

    temp1c <- temp[temp$flag == 'CC', ]
    temp1c <- temp1c[order(temp1c$qval), ]
    wb <- openxlsx::loadWorkbook(paste0(out_dir,'/GO_BioCRPE_',net_type[qq],'.xlsx'))
    openxlsx::addWorksheet(wb, sheetName = 'CC')
    openxlsx::writeData(wb, sheet = 'CC', temp1c)
    openxlsx::saveWorkbook(wb, paste0(out_dir,'/GO_BioCRPE_',net_type[qq],'.xlsx'), overwrite = T)


}


