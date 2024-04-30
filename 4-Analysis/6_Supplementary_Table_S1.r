##############################################################################################
# Purpose: Save excel files with different properties 
##############################################################################################

rm(list=ls())
library(data.table)
library('biomaRt')
library(Rcpp)
library(RColorBrewer)
library(ggplot2)
library(GenomicDataCommons)
source("eein_cancer_util.r")

cppFunction("List getpropa(CharacterVector ex1, CharacterVector ex2, CharacterVector exon1, CharacterVector exon2, NumericVector prop){

    int loop1 = ex1.size();
    int loop2 = exon1.size();
    NumericVector tprop;
    int demo=100000;

    for(int k=0; k<loop1; k++){

        int flag = 0;
        for(int j=0; j<loop2; j++){

            if((ex1[k] == exon1[j]) & (ex2[k] == exon2[j])){
                tprop.push_back(prop[j]);
                flag = 1;
                break;
            }

            if((ex1[k] == exon2[j]) & (ex2[k] == exon1[j])){
                tprop.push_back(prop[j]);
                flag = 1;
                break;
            }
        }

        if(flag == 0){
            tprop.push_back(demo);
        }

    }

    List L = List::create(tprop);
    return L;
  
}")


cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 
    'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))

net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')

outdir <- '../results/Final_results'
allnets <- gtools::mixedsort(list.files('../data/Final_networks/level1',full.names=TRUE))
cpm_threshold <- 0.5

gnet <- data.table::fread(paste0('../data/PISA_survival/PISA_net_final_',cpm_threshold,'.txt'), header=FALSE)
gnet <- mapProtein(gnet[[1]], gnet[[2]], data.table::fread('../data/PISA_networks_filtered/PISA_EEIN_0.5.txt'))
gnetx <- gnet[,c(1,2,5,6)]

enet <- data.table::fread(paste0('../data/EPPIC_survival/EPPIC_net_final_',cpm_threshold,'.txt'), header=FALSE)
enet <- mapProtein(enet[[1]], enet[[2]], data.table::fread('../data/EPPIC_EEIN_filtered.txt'))
enetx <- enet[,c(1,2,5,6)]

anet <- data.table::fread(paste0('../data/CONTACT_survival/CONTACT_net_final_',cpm_threshold,'.txt'), header=FALSE)
xxp <- data.table::fread('../data/CONTACT_networks/CONTACT_net_6_1.txt')
xxp$p1 <- unlist(lapply(strsplit(xxp[[1]],'[_]'),'[[',1))
xxp$p2 <- unlist(lapply(strsplit(xxp[[2]],'[_]'),'[[',1))
xxp <- xxp[, -c(1,2)]
anet <- mapProtein(anet[[1]], anet[[2]], xxp)
anetx <- anet[, c(1,2,10,11)]
colnames(anetx) <- c('exon1','exon2','protein1','protein2')

aq <- rbind(gnetx, anetx)
unet <- rbind(aq, enetx)

##-------- for survival analysis edges -------------------------
for(qq in 1:length(allnets)){

    indir <- paste0('../data/Final_survival/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold)
    in_dir <- paste0('../data/Final_weighted_networks/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1])

    #-- save excel file ---
    wb1 <- openxlsx::createWorkbook(paste0(outdir,'/Supplementary_Table_1_',strsplit(basename(net_type[qq]),
        '[.]')[[1]][1],'.xlsx'))
    openxlsx::saveWorkbook(wb1, paste0(outdir,'/Supplementary_Table_1_',strsplit(basename(net_type[qq]),
        '[.]')[[1]][1],'.xlsx'), overwrite=T)

    for(k in 1:length(cancer_type)){

        c_type <- cancer_type[k]
        
        gained <- data.table::fread(paste0(indir,'/',c_type,'_Gained_Surv_.txt'))
        gained <- gained[gained$pval <= 0.05, ]
        lost <- data.table::fread(paste0(indir,'/',c_type,'_Lost_Surv_.txt'))
        lost <- lost[lost$pval <= 0.05, ]

        only_gained <- igraph::as_data_frame(igraph::difference(igraph::graph_from_data_frame(gained[,c(1,2)], directed=FALSE), igraph::graph_from_data_frame(lost[,c(1,2)], directed=FALSE)))
        only_lost <- igraph::as_data_frame(igraph::difference(igraph::graph_from_data_frame(lost[,c(1,2)], directed=FALSE), igraph::graph_from_data_frame(gained[,c(1,2)], directed=FALSE)))
        gained_lost <- igraph::as_data_frame(igraph::intersection(igraph::graph_from_data_frame(lost[,c(1,2)], directed=FALSE), igraph::graph_from_data_frame(gained[,c(1,2)], directed=FALSE)))

        og_pv <- getpropa(only_gained[[1]], only_gained[[2]], gained[[1]], gained[[2]], gained[[3]])
        only_gained$pvalue <- og_pv[[1]]
        og_pt <- getpropa(only_gained[[1]], only_gained[[2]], gained[[1]], gained[[2]], gained[[7]])
        only_gained$gained_in_patient <- og_pt[[1]]
        only_gained$lost_in_patient <- rep(0,length(only_gained[[1]]))

        og_pv <- getpropa(only_lost[[1]], only_lost[[2]], lost[[1]], lost[[2]], lost[[3]])
        only_lost$pvalue <- og_pv[[1]]
        og_pt <- getpropa(only_lost[[1]], only_lost[[2]], lost[[1]], lost[[2]], lost[[7]])
        only_lost$gained_in_patient <- rep(0,length(only_lost[[1]]))
        only_lost$lost_in_patient <- og_pt[[1]]

        og_pv <- getpropa(gained_lost[[1]], gained_lost[[2]], lost[[1]], lost[[2]], lost[[3]])
        gained_lost$pvalue <- og_pv[[1]]
        og_pt <- getpropa(gained_lost[[1]], gained_lost[[2]], gained[[1]], gained[[2]], gained[[7]])
        gained_lost$gained_in_patient <- og_pt[[1]]
        og_pt <- getpropa(gained_lost[[1]], gained_lost[[2]], lost[[1]], lost[[2]], lost[[7]])
        gained_lost$lost_in_patient <- og_pt[[1]]

        ##------ e) 
        t2 <- rbind(only_gained, only_lost)
        t3 <- rbind(t2, gained_lost)
        t5 <- t3[order(t3$pvalue), ]
        xx <- mapProtein(t5[[1]], t5[[2]], unet)
        xx$pvalue <- t5$pvalue
        xx$gained_in_patient <- t5$gained_in_patient
        xx$lost_in_patient <- t5$lost_in_patient

        wb <- openxlsx::loadWorkbook(paste0(outdir,'/Supplementary_Table_1_',strsplit(basename(net_type[qq]),
            '[.]')[[1]][1],'.xlsx'))
        openxlsx::addWorksheet(wb, sheetName = c_type)
        openxlsx::writeData(wb, sheet = c_type, xx)
        openxlsx::saveWorkbook(wb, paste0(outdir,'/Supplementary_Table_1_',strsplit(basename(net_type[qq]),
            '[.]')[[1]][1],'.xlsx'), overwrite = T)
        
    }
}
