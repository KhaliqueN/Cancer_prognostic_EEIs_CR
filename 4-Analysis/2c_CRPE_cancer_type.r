##############################################################################################
# Purpose: Get perturbed EEIs and PPIs
# This is to test for why there is significantly high numbers of cancer relevant edges in KIRC
## Are thre any collection of PPIs that over-represnets the cancer relevant edges?
## Splicing events from --> Kahles et al. Comprehensive Analysis of Alternative Splicing Across Tumors from 8,705 Patients, Cancer Cell, Volume 34, Issue 2, 2018
## Mutations from -->  Pan-cancer analysis of whole genomes (https://www.nature.com/articles/s41586-020-1969-6#MOESM3) supplementary table S1
## 
################################################################################################

rm(list=ls())
library(data.table)
library(ggplot2)
library(RColorBrewer)
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

temp_bar_limitl <- 50
temp_bar_limit <- 500
temp_y_limit <- 1400
temp_x_limit <- 5000

for(qq in 1:length(allnets)){

    if(qq==1){
        temp_bar_limit <- 1000
        temp_y_limit <- 4000
        temp_bar_limitl <- 500
    }else if(qq==2){
        temp_y_limit <- 3000
        temp_bar_limitl <- 300
    }

    ppi_count <- c()
    eei_count <- c()
    ppi_max <- c()
    perturbed_graph <- list()

    xx <- c() ## number of eeis
    ff <- c() ## number of ppis
    ss <- c() ## Fraction of all PPIs
    Cancer <- c()
    eeis <- list()
    ppis <- list()
    patient_count <- as.data.frame(matrix(nrow=0,ncol=2))
    mdn <- c()
    mxp <- c()


    net_file <- igraph::graph_from_data_frame(data.table::fread(allnets[qq]), directed=FALSE)
    net_filed <- igraph::as_data_frame(net_file)
    netx <- mapProtein(net_filed[[1]], net_filed[[2]], unet)

    for(gg in 1:length(cancer_type)){

        c_type <- cancer_type[gg]
        ##--- get survival results ------------
        gained_surv <- data.table::fread(paste0('../results/Final_survival/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
        c_type,'_Gained_Surv_.txt'))
        gained_surv$qvalue <- p.adjust(gained_surv$pval, 'fdr')
        gained_surv1 <- gained_surv[gained_surv$pval <= pval_thres, ]
        gs <- gained_surv1[,c(1,2,6)]
        # mdn <- c(mdn, gained_surv1[[7]][1]/gained_surv1[[6]][1])
        # gs <- igraph::graph_from_data_frame(gained_surv1[,c(1,2,6)], directed=FALSE)
        ##-------------------------------------
        
        ##--- get survival results ------------
        lost_surv <- data.table::fread(paste0('../results/Final_survival/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
        c_type,'_Lost_Surv_.txt'))
        lost_surv1 <- lost_surv[lost_surv$pval <= pval_thres, ]
        ls <- lost_surv1[,c(1,2,6)]
        # ls <- igraph::graph_from_data_frame(lost_surv1[,c(1,2,6)], directed=FALSE)
        ##-------------------------------------
        gs_ls_u <- rbind(gs,ls)
        gs_ls_u <- igraph::graph_from_data_frame(gs_ls_u, directed=FALSE)
        gs_ls_u <- igraph::simplify(gs_ls_u, edge.attr.comb=list(weight="sum"))
        # gs_ls_u <- igraph::simplify(gs_ls_u, edge.attr.comb=list(patient="sum"))

        # igraph::union(ls, gs)
        eeis[[gg]] <- gs_ls_u
        tempg <- igraph::as_data_frame(gs_ls_u)
        mdn <- c(mdn, median(tempg[[3]]*(gained_surv1[[7]][1]/gained_surv1[[6]][1])))
        mxp <- c(mxp, max(100*tempg[[3]]))

        ## store patient counts of perturbed edges
        patient_count <- rbind(patient_count,data.frame(cancer=rep(cancer_type[gg],length(tempg[[3]])), 100*tempg[[3]]))
        # patient_count <- rbind(patient_count,data.frame(cancer=rep(cancer_type[gg],length(tempg[[3]])), tempg[[3]]))

        p1 <- c()
        p2 <- c()

        for(i in 1:length(tempg[[1]])){
            wh1 <- which(unet$exon1 == tempg[[1]][i])
            wh2 <- which(unet$exon2 == tempg[[2]][i])
            wha <- intersect(wh1, wh2)
            if(length(wha) > 0){
                tempx <- unet[wha, ]
                p1 <- c(p1, tempx$protein1[1])
                p2 <- c(p2, tempx$protein2[1])
                next
            }
            wh1 <- which(unet$exon2 == tempg[[1]][i])
            wh2 <- which(unet$exon1 == tempg[[2]][i])
            whb <- intersect(wh1, wh2)
            tempx <- unet[whb, ]
            p2 <- c(p2, tempx$protein1[1])
            p1 <- c(p1, tempx$protein2[1])
        }

        tempg$protein1 <- p1
        tempg$protein2 <- p2

        ppig <- igraph::graph_from_data_frame(tempg[,c(3,4)], directed=FALSE)
        ppi_graph <- igraph::simplify(ppig)
        ppis[[gg]] <- ppi_graph
        ppi_count <- c(ppi_count, igraph::ecount(ppi_graph))

        igraph::E(ppig)$weight <- 1
        ppig1 <- igraph::simplify(ppig, edge.attr.comb=list(weight="sum"))
        ppig2 <- igraph::as_data_frame(ppig1)
        ppig2 <- ppig2[order(-ppig2$weight), ]
        ppi_max <- c(ppi_max, ppig2$weight[1])

        eei_graph <- igraph::simplify(igraph::graph_from_data_frame(tempg[,c(1,2)], directed=FALSE))
        eei_count <- c(eei_count, igraph::ecount(eei_graph))

        ## store for plot -----
        ppig2_count <- plyr::count(ppig2$weight)
        ppig2_count$pert <- ppig2_count$freq/sum(ppig2_count$freq)
        xx <- c(xx, ppig2_count$x)
        ff <- c(ff, ppig2_count$freq)
        ss <- c(ss,ppig2_count$pert)
        Cancer <- c(Cancer, rep(cancer_type[gg], length(ppig2_count[[1]])))
        perturbed_graph[[gg]] <- gs_ls_u

    }

    ##-- plot -----------------------------------
    # pdata <- data.frame(Cancer=Cancer, eei=xx, ppi=ff, Frac=ss)
    # frac <- c()
    # for(gg in 1:length(cancer_type)){
    #     tempx <- pdata[pdata$Cancer == cancer_type[gg], ]
    #     tempx <- tempx[tempx$eei == 1, ]
    #     frac <- c(frac, sum(tempx$Frac)*100)
    # }

    # pdata1 <- data.frame(Cancer=cancer_type, Frac=frac)
    # p <- ggplot(pdata1, aes(Cancer, Frac, fill=Cancer)) + 
    # geom_bar(stat="identity",position=position_dodge())+
    # theme(legend.text=element_text(size=12))
    # basesize <- 12
    # p <- p + theme_bw(base_size = basesize * 0.8) +
    # # geom_errorbar(aes(ymin=Perturbed-std, ymax=Perturbed+std), size=0.3, width=.05, position=position_dodge(0.1)) + 
    # # geom_hline(color='black', linetype='dashed', size=0.25, yintercept=overall_avg)+
    # scale_x_discrete(name="Cancer type") + 
    # scale_y_continuous(name="Percentage of PPIs that contain \nonly one cancer relevant EEI",breaks = seq(0, 100, by = 25)) +
    # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    # theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    # axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    # strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    # guides(fill=guide_legend(title="Cancer type",ncol=2))
    # ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_eei_per_ppi.png"),width=7, height=3.5, dpi=400)

    # ##--- correlation plot -----------
    # pdata2 <- data.frame(Cancer=cancer_type, Frac=frac, EEI=eei_count)

    # p <- ggplot(pdata2, aes(Frac, EEI, color=Cancer)) + 
    # geom_point(size=2)+
    # theme(legend.text=element_text(size=12))
    # basesize <- 12
    # p <- p + theme_bw(base_size = basesize * 0.8) +
    # scale_x_continuous(name="Percentage of PPIs that contain \nonly one cancer relevant EEI") + 
    # scale_y_continuous(name="Total # of cancer relevant edges") +
    # scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    # theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    # axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    # strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    # guides(color=guide_legend(title="Cancer type",ncol=2))
    # ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_eei_per_ppi_cor.png"),width=7, height=3.5, dpi=400)

    ##----- number of cancer-relevant edges plot -----------
    pdata <- data.frame(Cancer=cancer_type, count=eei_count)
    # p <- ggplot(pdata, aes(Cancer, count, fill=Cancer)) + 
    p <- ggplot(pdata, aes(Cancer, count)) + 
    geom_bar(stat="identity",position=position_dodge())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,(max(pdata$count))+temp_bar_limit)) +
    geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
    # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_CRPE.png"),width=3.5, height=3, dpi=400)


    ### number of patients to which the perturbed edges belong -------
    pdatax <- patient_count
    pmdn <- data.frame(Cancer=cancer_type,count=signif(mdn,1),pos=mxp+8)
    colnames(pdatax) <- c('Cancer', 'count')
    p <- ggplot(pdatax, aes(Cancer, count, fill=Cancer)) + 
    geom_boxplot(alpha=0.9)+
    # geom_point(aes(fill = Cancer), shape = 21, position = position_jitterdodge()) +
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="% of patients \nwith paired samples", limits=c(0,110)) +
    geom_text(data=pmdn, aes(y=pos, x=Cancer,label=count), position=position_dodge(width=0.9),hjust=0.5, vjust=0.5, angle=0, size=3)+
    scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme_bw()+theme(axis.text.x = element_text(size = 8, angle = 60, vjust=1, hjust=1, colour = "black"),
        axis.text.y = element_text(size = 8, angle = 0, colour = "black"))+#,
        # panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill='none')#+guides(size=guide_legend(title="# of patients",ncol=1))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_pateint_violin.png"),width=4.5, height=2.5, dpi=400)

#breaks = seq(0, 100, by = 20)

    ### ---------------------- correlation with known number of cancer type-specific splicing events ------------------
    pdata <- pdata[order(pdata$Cancer), ]

    # exon_skip <- data.table::fread('../data/exon_skip.csv')
    # exon_skip <- exon_skip[exon_skip$V1 %in% cancer_type, ]
    # exon_skip <- exon_skip[order(exon_skip$V1), ]
    # exon_skip$V2 <- floor(exon_skip$V2)

    # mutual_exon <- data.table::fread('../data/mutual_exon.csv')
    # mutual_exon <- mutual_exon[mutual_exon$V1 %in% cancer_type, ]
    # mutual_exon <- mutual_exon[order(mutual_exon$V1), ]
    # mutual_exon$V2 <- floor(mutual_exon$V2)

    # all_events <- data.frame(V1=exon_skip[[1]], V2=exon_skip[[2]]+mutual_exon[[2]])

    # pdata$skip_events <- exon_skip$V2
    # pdata$mutual_events <- mutual_exon$V2
    # pdata$all_events <- all_events$V2

    # coral <- cor.test(x=pdata$count, y=pdata$all_events, method = 'spearman')
    # p <- ggplot(pdata, aes(all_events, count, color=Cancer)) + 
    # geom_point(size=2)+
    # theme(legend.text=element_text(size=12))
    # basesize <- 12
    # p <- p + theme_bw(base_size = basesize * 0.8) +
    # scale_x_continuous(name="# of exon skipping or \nmutually-exclusive exons splicing events", limits=c(0,max(pdata$all_events))) + 
    # geom_text(aes(x=6000,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5)+
    # scale_y_continuous(name="# of cancer-relevant perturbed \nedges (CRPEs)", limits=c(0,max(pdata$count))) +
    # scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    # theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    # axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    # strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    # guides(color=guide_legend(title="Cancer type",ncol=2))
    # ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_allEvents.png"),width=4.5, height=3, dpi=400)

    ###------------------------------------------------------------------------------------------------------------------

    ###----- AS events and survival ----
    splicing <- readxl::read_excel('../data/oncosplicing.xlsx', 1)
    splicing <- as.data.frame(splicing)
    colnames(splicing) <- as.vector(splicing[1,])
    splicing <- splicing[-1,]
    splicing <- splicing[1:31,]
    splicing1 <- splicing[,c(1,5,8)]
    splicing1 <- splicing1[splicing1$`TCGA cancer` %in% cancer_type, ]
    splicing1 <- splicing1[order(splicing1$`TCGA cancer`), ]
    pdata$overall_splicing <- as.numeric(splicing1[[2]])
    pdata$splicing <- as.numeric(splicing1[[3]])
    wh <- which(!is.na(pdata$splicing))
    pdatax <- pdata[wh,]

    coral <- cor.test(x=pdatax$count, y=pdatax$splicing, method = 'spearman')
    p <- ggplot(pdatax, aes(splicing, count, color=Cancer)) + 
    geom_point(size=2)+
    theme(legend.text=element_text(size=12))
    basesize <- 10
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_continuous(name="# of survival associated \nalternative splicing events", limits=c(0,max(pdatax$splicing))) + 
    geom_text(aes(x=8500,y=temp_y_limit, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=3,show.legend = FALSE)+
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,max(pdatax$count))) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','black','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(color=guide_legend(title="Cancer type",ncol=3))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_splicingSurvival.png"),width=4.5, height=2.5, dpi=400)

    coral <- cor.test(x=pdata$count, y=pdata$overall_splicing, method = 'spearman')
    p <- ggplot(pdata, aes(overall_splicing, count, color=Cancer)) + 
    geom_point(size=2)+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_continuous(name="# of alternative splicing\n events", limits=c(0,max(pdata$overall_splicing))) + 
    geom_text(aes(x=45000,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5,show.legend = FALSE)+
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,max(pdata$count))) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(color=guide_legend(title="Cancer type",ncol=2))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_splicing.png"),width=4, height=2.5, dpi=400)



    # ### ---------------------- correlation with known number of mutations ------------------
    # mutations <- readxl::read_excel('../data/mutations.xlsx', 1)
    # mutations <- as.data.frame(mutations)
    # cols_in <- c('project_code','Coding.SNVs','Non.coding.SNVs','Coding.MNVs','Non.coding.MNVs','Coding.Indels','Non.coding.Indels','all.SNVs','all.MNVs','all.Indels')
    # mutations1 <- mutations[,which(colnames(mutations) %in% cols_in)]
    # mutations1$project_code <- unlist(lapply(strsplit(mutations1$project_code, '[-]'), '[[', 1))
    # mutations2 <- mutations1[mutations1$project_code %in% cancer_type, ]
    # pdatax <- pdata[pdata$Cancer %in% unique(mutations2$project_code), ]
    # mutations3 <- aggregate(.~project_code, mutations2, median)
    # mutations3 <- mutations3[order(mutations3$project_code), ]
    # pdatax$mutations <- mutations3$`Coding.SNVs`+mutations3$`Coding.Indels`

    # coral <- cor.test(x=pdatax$count, y=pdatax$mutations, method = 'spearman')
    # p <- ggplot(pdatax, aes(mutations, count, color=Cancer)) + 
    # geom_point(size=2)+
    # theme(legend.text=element_text(size=12))
    # basesize <- 12
    # p <- p + theme_bw(base_size = basesize * 0.8) +
    # scale_x_continuous(name="Median # of coding region-related \nsingle nucleotide variations \nor indels per patient", limits=c(0,max(pdatax$mutations))) + 
    # geom_text(aes(x=100,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=2.5)+
    # scale_y_continuous(name="# of cancer-relevant perturbed \n edges (CRPEs)", limits=c(0,max(pdatax$count))) +
    # scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    # theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    # axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    # strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    # guides(color='none')
    # ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_mutations.png"),width=3, height=3, dpi=400)
    # ###----------------------------------------------------------------

    # ### ---------------------- correlation with known degree of heterogeneity ------------------
    # mutations <- data.table::fread('../data/heterogeniety.csv')
    # mutations <- mutations[order(mutations$V1), ]
    # pdata$mutations <- mutations$V2

    # coral <- cor.test(x=pdata$count, y=pdata$mutations, method = 'spearman')
    # p <- ggplot(pdata, aes(mutations, count, color=Cancer)) + 
    # geom_point(size=2)+
    # theme(legend.text=element_text(size=12))
    # basesize <- 12
    # p <- p + theme_bw(base_size = basesize * 0.8) +
    # scale_x_continuous(name="Median # of clones", limits=c(0,max(pdata$mutations))) + 
    # scale_y_continuous(name="# of cancer-relevant \n perturbed edges", limits=c(0,max(pdata$count))) +
    # scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    # theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    # axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    # strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    # guides(color='none')
    # ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_heterogeneity.png"),width=3, height=3, dpi=400)

    # ###----------------------------------------------------------------

    # ### ---------------------- correlation with known numbers of genes significantly correlated with survival ------------------
    # positions <- c(2,3,5,6,9,10,11,12,13,14,19,22,25,26)
    # tcan <- c('BLCA', 'BRCA', 'COAD', 'ESCA', 'HNSC', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'PRAD', 'STAD', 'THCA', 'UCEC')
    # len <- c()
    # for(pp in 1:length(positions)){
    #     survival <- readxl::read_excel('../data/cancer_hallmark_genes.xlsx', positions[pp])
    #     survival <- as.data.frame(survival)
    #     colnames(survival) <- as.vector(survival[1,])
    #     survival <- survival[-1,]
    #     survival$pvalue <- as.numeric(survival$pvalue)
    #     survivalx <- survival[survival$pvalue <= 0.05, ]
    #     len <- c(len, length(survivalx[[1]]))
    # }

    # survival_data <- data.frame(Cancer=tcan, count=len)

    # p <- ggplot(survival_data, aes(Cancer, count)) + 
    # geom_bar(stat="identity",position=position_dodge())+
    # theme(legend.text=element_text(size=12))
    # basesize <- 12
    # p <- p + theme_bw(base_size = basesize * 0.8) +
    # scale_x_discrete(name="Cancer type") + 
    # scale_y_continuous(name="# of cancer hallmark genes (671)\n correlating with survival", limits=c(0,(max(survival_data$count))+100)) +
    # geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
    # # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    # theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    # axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    # strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    # guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
    # ggsave(p,filename=paste0("../results/Final_results/Survival_hallmark_genes.png"),width=3.5, height=3, dpi=400)


    ####--------------------------------------------------------------------------------------------------------------------------------

    # ###----- Log-rank science ----
    # positions <- c(1,3,5,6,7,10,11,12,13,14,19,22,25,26)
    # tcan <- c('BRCA','COAD', 'UCEC', 'LIHC', 'HNSC', 'ESCA', 'HNSC', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'PRAD', 'STAD', 'THCA', 'UCEC')
    # len <- c()
    # for(pp in 1:length(positions)){
    #     survival <- readxl::read_excel('../data/logrank_science.xlsx', positions[pp])
    #     survival <- as.data.frame(survival)
    #     colnames(survival) <- as.vector(survival[1,])
    #     survival <- survival[-1,]
    #     survival$pvalue <- as.numeric(survival$pvalue)
    #     survivalx <- survival[survival$pvalue <= 0.05, ]
    #     len <- c(len, length(survivalx[[1]]))
    # }




    ###---------------------------------------------------------------






    ### ---------------------- correlation with known numbers of genes significantly correlated with survival ------------------
    survival <- data.table::fread('../data/survival_analysis_all.txt', sep='\t')
    survs <- survival[,-1]

    # signif_cn1 <- 2*(sapply(survs[[7]], function(x) pnorm(q=abs(x), lower.tail=FALSE) ))
    # signif_cn2 <- p.adjust(signif_cn1, 'fdr')
    # signif_cn2 <- apply(signif_cn1, 2, function(x) p.adjust(x, 'fdr'))
    # signif_can <- apply(signif_cn2, 2, function(x) length(which(x <= 0.05)))

    signif_can <- sapply(survs, function(x) length(which(x > 1.96 | x < -1.96))) ## 3.1 in case of the 0.001 p-value
    survival_data <- data.frame(Cancer=names(signif_can), count=signif_can)
    survival_data <- survival_data[order(survival_data$Cancer), ]
    pdata$genes <- survival_data$count


    p <- ggplot(survival_data, aes(Cancer, count)) + 
    geom_bar(stat="identity",position=position_dodge())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of genes", limits=c(0,max(survival_data$count)+1500)) +
    geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
    # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))
    # +guide_legend(title="Cancer type",ncol=2)
    ggsave(p,filename=paste0("../results/Final_results/Survival_genes.png"),width=3.5, height=3, dpi=400)



    coral1 <- cor.test(x=pdata$count, y=pdata$genes, method = 'spearman')
    # coral1 <- cor.test(x=pdata$count, y=pdata$genes, method = 'pearson')

    p <- ggplot(pdata, aes(genes, count, color=Cancer)) + 
    geom_point(size=2)+
    theme(legend.text=element_text(size=10))
    basesize <- 10
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_continuous(name="# of genes significantly \nassociated with survival", limits=c(0,max(pdata$genes))) + 
    geom_text(aes(x=5000,y=temp_y_limit, label=paste0('Spearman correlation: ',signif(coral1$estimate[[1]],3),'\np-value: ',signif(coral1$p.value,3))), size=3,show.legend = FALSE)+
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,max(pdata$count))) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(color=guide_legend(title="Cancer type",ncol=3))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_geneSurvival.png"),width=4.5, height=2.5, dpi=400)

    
    ##-- mirna---
    survival <- data.table::fread('../data/survival_analysis_all_mirna.txt', sep='\t')
    survs <- survival[,-1]
    signif_can <- sapply(survs, function(x) length(which(x > 1.96 | x < -1.96)))
    survival_data <- data.frame(Cancer=names(signif_can), count=signif_can)
    survival_data <- survival_data[order(survival_data$Cancer), ]
    pdata$mirna <- survival_data$count

    coral2 <- cor.test(x=pdata$count, y=pdata$mirna, method = 'spearman')
    # coral2 <- cor.test(x=pdata$count, y=pdata$mirna, method = 'pearson')

    p <- ggplot(pdata, aes(mirna, count, color=Cancer)) + 
    geom_point(size=2)+
    theme(legend.text=element_text(size=10))
    basesize <- 10
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_continuous(name="# of miRNAs significantly \nassociated with survival", limits=c(0,max(pdata$mirna))) + 
    geom_text(aes(x=150,y=1400, label=paste0('Spearman correlation: ',signif(coral2$estimate[[1]],3),'\np-value: ',signif(coral2$p.value,3))), size=3,show.legend = FALSE)+
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,max(pdata$count))) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(color='none')
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_miRNASurvival.png"),width=2.5, height=2.5, dpi=400)


    ##-- proteins ---
    survival <- data.table::fread('../data/survival_analysis_all_protein.txt', sep='\t')
    survs <- survival[,-1]
    signif_can <- sapply(survs, function(x) length(which(x > 1.96 | x < -1.96)))
    survival_data <- data.frame(Cancer=names(signif_can), count=signif_can)
    survival_data <- survival_data[order(survival_data$Cancer), ]
    pdata$protein <- survival_data$count

    coral3 <- cor.test(x=pdata$count, y=pdata$protein, method = 'spearman')
    p <- ggplot(pdata, aes(protein, count, color=Cancer)) + 
    geom_point(size=2)+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_continuous(name="# of proteins significantly \nassociated with survival", limits=c(0,max(pdata$protein))) + 
    geom_text(aes(x=100,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=3,show.legend = FALSE)+
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,max(pdata$count))) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(color=guide_legend(title="Cancer type",ncol=4))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_proteinSurvival.png"),width=5.5, height=2.5, dpi=400)




    ##-- mutations ---
    survival <- data.table::fread('../data/survival_analysis_all_Mutation.txt', sep='\t')
    survs <- survival[,-1]
    signif_can <- sapply(survs, function(x) length(which(x > 1.96 | x < -1.96)))
    survival_data <- data.frame(Cancer=names(signif_can), count=signif_can)
    survival_data <- survival_data[order(survival_data$Cancer), ]
    pdata$mutation <- survival_data$count

    coral4 <- cor.test(x=pdata$count, y=pdata$mutation, method = 'spearman')
    p <- ggplot(pdata, aes(mutation, count, color=Cancer)) + 
    geom_point(size=2)+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_continuous(name="# of mutations significantly \nassociated with survival", limits=c(0,max(pdata$mutation))) + 
    geom_text(aes(x=800,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=3,show.legend = FALSE)+
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,max(pdata$count))) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(color=guide_legend(title="Cancer type",ncol=4))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_mutationSurvival.png"),width=5.5, height=2.5, dpi=400)


    ##-- methylations ---
    survival <- data.table::fread('../data/survival_analysis_all_Methylation.txt', sep='\t')
    survs <- survival[,-1]
    signif_can <- sapply(survs, function(x) length(which(x > 1.96 | x < -1.96)))
    survival_data <- data.frame(Cancer=names(signif_can), count=signif_can)
    survival_data <- survival_data[order(survival_data$Cancer), ]
    pdata$methylation <- survival_data$count

    coral5 <- cor.test(x=pdata$count, y=pdata$methylation, method = 'spearman')
    p <- ggplot(pdata, aes(methylation, count, color=Cancer)) + 
    geom_point(size=2)+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_continuous(name="# of methylations significantly \nassociated with survival", limits=c(0,max(pdata$methylation))) + 
    geom_text(aes(x=800,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=3,show.legend = FALSE)+
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,max(pdata$count))) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(color=guide_legend(title="Cancer type",ncol=4))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_methylationSurvival.png"),width=5.5, height=2.5, dpi=400)


    ##-- CNA ---
    survival <- data.table::fread('../data/survival_analysis_all_CNA.txt', sep='\t')
    survs <- survival[,-1]
    signif_can <- sapply(survs, function(x) length(which(x > 1.96 | x < -1.96)))
    survival_data <- data.frame(Cancer=names(signif_can), count=signif_can)
    survival_data <- survival_data[order(survival_data$Cancer), ]
    pdata$CNA <- survival_data$count

    coral6 <- cor.test(x=pdata$count, y=pdata$CNA, method = 'spearman')
    p <- ggplot(pdata, aes(CNA, count, color=Cancer)) + 
    geom_point(size=2)+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_continuous(name="# of Copy Number Abberations (CNAs) \nsignificantly associated with survival", limits=c(0,max(pdata$CNA))) + 
    geom_text(aes(x=800,y=1400, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=3,show.legend = FALSE)+
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,max(pdata$count))) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(color=guide_legend(title="Cancer type",ncol=4))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_per_cancer_cnaSurvival.png"),width=5.5, height=2.5, dpi=400)



    # ###----------------------------------------------------------------
    # ###----- Unique CRPEs ------------------------------------------
    unique_eeis <- list()
    ueeis <- c()
    for(k in 1:length(cancer_type)){
        counter <- k
        temp_uni <- eeis[[k]]

        for(j in 1:length(cancer_type)){
            if(counter != j){
                temp_uni <- igraph::difference(temp_uni, eeis[[j]])
            }
        }
        unique_eeis[[k]] <- temp_uni
        ueeis <- c(ueeis, igraph::ecount(temp_uni))
    }

    ##----- number of cancer-relevant edges plot -----------
    pdata$ueeis <- ueeis
    # p <- ggplot(pdata, aes(Cancer, count, fill=Cancer)) + 
    p <- ggplot(pdata, aes(Cancer, ueeis)) + 
    geom_bar(stat="identity",position=position_dodge())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of unique cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,(max(pdata$ueeis))+temp_bar_limitl)) +
    geom_text(aes(label=ueeis), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_unique_CRPE.png"),width=3.5, height=3, dpi=400)

    ##--- write unique EEIs in a file ---
    alld <- data.frame(matrix(ncol=3, nrow=0))
    for(k in 1:length(cancer_type)){
        temp <- igraph::as_data_frame(unique_eeis[[k]])
        temp <- temp[,-3]
        temp$cancer <- rep(cancer_type[k], length(temp[[1]]))
        alld <- rbind(alld, temp)
    }

    data.table::fwrite(alld, '../data/unique_edges_nethigh.txt', quote=FALSE, col.names=FALSE, row.names=FALSE)



    # ## which PPIs are unique ------
    # unique_ppis <- list()
    # for(k in 1:length(cancer_type)){
    #     counter <- k
    #     temp_uni <- ppis[[k]]

    #     for(j in 1:length(cancer_type)){
    #         if(counter != j){
    #             temp_uni <- igraph::difference(temp_uni, ppis[[j]])
    #         }
    #     }
    #     unique_ppis[[k]] <- temp_uni
    # }


    # for(k in 1:length(cancer_type)){
    #     print(igraph::ecount(unique_ppis[[k]]))
    # }












    ####-------------------------------------------------------------------

    ##--- figure for the overlap of CRPEs ---------------------------
    c1 <- c()
    c2 <- c()
    oval <- c()
    pval <- c()
    pvalu <- c()
    loop1 <- length(cancer_type)-1
    loop2 <- length(cancer_type)

    for(ov1 in 1:loop1){
        m <- ov1+1
        for(ov2 in m:loop2){

            # if(ov1 != ov2){
            c1 <- c(c1, paste0(cancer_type[ov1],' (', eei_count[ov1],')'))
            c2 <- c(c2, paste0(cancer_type[ov2],' (', eei_count[ov2],')'))
            e1 <- igraph::ecount(perturbed_graph[[ov1]])
            e2 <- igraph::ecount(perturbed_graph[[ov2]])
            inter_e <- igraph::intersection(perturbed_graph[[ov1]], perturbed_graph[[ov2]])
            if(e1 < e2){
                tempo <- igraph::ecount(inter_e)/e1
            }else{
                tempo <- igraph::ecount(inter_e)/e2
            }
            oval <- c(oval, igraph::ecount(inter_e))
            # oval <- c(oval, tempo)

            hyp <- phyper(igraph::ecount(inter_e)-1,e1,igraph::ecount(net_file)-e1,e2,lower.tail = FALSE)
            hyp <- signif(hyp, digits=3)
            pval <- c(pval, hyp)

            ## under representation
            hyp <- phyper(igraph::ecount(inter_e),e1,igraph::ecount(net_file)-e1,e2,lower.tail = TRUE)
            hyp <- signif(hyp, digits=3)
            pvalu <- c(pvalu, hyp)
        # }
        }
    }

    qval <- p.adjust(pval, 'fdr')
    qvalu <- p.adjust(pvalu, 'fdr')
    qvalx <- rep('Overlap as \nexpected by chance',length(pval))
    for(i in 1:length(pval)){
        if(qval[i] <= pval_thres){
            qvalx[i] <- 'Higher overlap than\nexpected by chance'
        }
        if(qvalu[i] <= pval_thres){
            qvalx[i] <- 'Lower overlap than\nexpected by chance'
        }
    }

    pdata <- data.frame(c1=c1, c2=c2, val=oval, pval=qvalx)


    # pdata <- data.frame(c1=c1, c2=c2, val=signif(oval*100,2), qval=signif(qval, 2))
    # pdata <- data.frame(c1=c1, c2=c2, val=oval, qval=signif(qval, 2))
    # pdata$pvalue <- cut(pdata$qval, breaks = c(0,0.01, 0.05, 1), include.lowest=TRUE)

    cols <- rev(c('#377eb8','#4daf4a','#e41a1c'))#rev(brewer.pal(3,"Spectral"))
    p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pval),colour = "white")+
      theme(legend.text=element_text(size=8))+scale_fill_manual(values=cols,drop=FALSE)
    basesize <- 8
    p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type", y = "Cancer type") +
      scale_y_discrete() +
      scale_x_discrete()+coord_fixed()+
      guides(fill=guide_legend(title="Overlap category"))+
      geom_text(data=pdata,aes(y=c2,label=val),size=1.8)+
      theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0, colour = "black"),
        axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0, colour = "black"))
    ggsave(p,filename=paste0(out_dir,'/CRPE_edge_overlap_smaller_',net_type[qq],'.png'),width=4.5, height=3, dpi=300)


    # pdata <- data.frame(c1=c1, c2=c2, val=oval, qval=signif(qvalu, 2))

    # pdata$pvalue <- cut(pdata$qval, breaks = c(0,0.01, 0.05, 1), include.lowest=TRUE)

    # cols <- c('#e41a1c','#4daf4a','#377eb8')#rev(brewer.pal(3,"Spectral"))
    # p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pvalue),colour = "white")+
    #   theme(legend.text=element_text(size=10))+scale_fill_manual(values=cols,drop=FALSE)
    # basesize <- 10
    # p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type", y = "Cancer type") +
    #   scale_y_discrete() +
    #   scale_x_discrete()+coord_fixed()+
    #   guides(fill=guide_legend(title="Q-value"))+
    #   geom_text(data=pdata,aes(y=c2,label=val),size=2.5)+
    #   theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0.5, colour = "black"),
    #     axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0.5, colour = "black"))
    # ggsave(p,filename=paste0(out_dir,'/CRPE_edge_overlap_smaller_',net_type[qq],'_under.png'),width=5, height=5, dpi=300)

    ##-------

    ###-----= shared CRPEs --------
    c1 <- c()
    c2 <- c()
    for(i in 1:length(cancer_type)){
        tempgg <- as.data.frame(igraph::get.edgelist(eeis[[i]]))
        c1 <- c(c1, tempgg[[1]])
        c2 <- c(c2, tempgg[[2]])
    }

    shared_eeis <- data.frame(c1,c2)
    shared_eeis <- igraph::graph_from_data_frame(shared_eeis, directed=FALSE)
    igraph::E(shared_eeis)$weight <- 1
    shared_eeis <- igraph::simplify(shared_eeis, edge.attr.comb=list(weight="sum"))
    shared_eeis <- igraph::as_data_frame(shared_eeis)
    shared_eeis <- shared_eeis[order(-shared_eeis$weight), ]

    sharedx <- shared_eeis[shared_eeis$weight >= 5, ]
    sin <- c()

    for(i in 1:length(sharedx[[1]])){
        wh1 <- which(unet[[1]] == sharedx[[1]][i])
        wh2 <- which(unet[[2]] == sharedx[[2]][i])
        wha <- intersect(wh1, wh2)

        wh1 <- which(unet[[1]] == sharedx[[2]][i])
        wh2 <- which(unet[[2]] == sharedx[[1]][i])
        whb <- intersect(wh1, wh2)

        wh <- union(wha, whb)
        sin <- c(sin, wh[1])
    }

    sharedy <- unet[sin, ]
    uppi <- igraph::simplify(igraph::graph_from_data_frame(sharedy[,c(3,4)], directed=FALSE))


    ## number of edges shared 
    tempg <- plyr::count(shared_eeis$weight)
    tempg$log <- log2(tempg$freq)

    p <- ggplot(tempg, aes(x, log)) + 
    geom_point(size=3, alpha=0.4)+#geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_y_continuous(name="# of cancer-relevant perturbed \nedges (CRPEs) in log2") + 
    scale_x_continuous(name="# of cancer types") +
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(shape='none')
    ggsave(p,filename=paste0(out_dir, "/", net_type[qq],'_shared_CRPEs.png'),width=3, height=3, dpi=400)

}




