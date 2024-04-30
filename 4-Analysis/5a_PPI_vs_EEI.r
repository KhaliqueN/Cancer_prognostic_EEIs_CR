##############################################################################################
# Purpose: Figures for importance of EEI vs PPI
# Highlighting CRPEs
##############################################################################################

rm(list=ls())
library(data.table)
library('biomaRt')
library(Rcpp)
library(RColorBrewer)
library(ggplot2)
library(GenomicDataCommons)
library(dplyr)
source("eein_cancer_util.r")

cppFunction("List getprop(CharacterVector ex1, CharacterVector ex2, CharacterVector exon1, CharacterVector exon2, NumericVector prop){

    int loop1 = ex1.size();
    int loop2 = exon1.size();
    NumericVector tprop;

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
            tprop.push_back(0);
            //Rcout << k << std::endl;
            //break;
        }

    }

    List L = List::create(tprop);
    return L;
  
}")


cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))
net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')

outdir <- '../results/Final_results'
allnets <- gtools::mixedsort(list.files('../data/Final_networks/level1',full.names=TRUE))
cpm_threshold <- 0.5

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
pval_thres <- 0.05

alldd <- data.table::fread('../data/gene_map.txt')


save_impPPIs <- function(tempcc){

    gene1 <- c()
    gene2 <- c()
    for(j in 1:length(tempcc[[1]])){
        p1 <- alldd[alldd$uniprotswissprot == tempcc$protein1[j],]
        if(nrow(p1) != 0){
            gene1 <- c(gene1, p1$hgnc_symbol[1])
        }else{
            gene1 <- c(gene1, 'Symbol')
        }

        p2 <- alldd[alldd$uniprotswissprot == tempcc$protein2[j],]
        if(nrow(p2) != 0){
            gene2 <- c(gene2, p2$hgnc_symbol[1])
        }else{
            gene2 <- c(gene2, 'Symbol')
        }
    }

    tempcc$gene1 <- gene1
    tempcc$gene2 <- gene2
    return(tempcc)

}


for(qq in 1:length(allnets)){

    in_dir <- paste0('../data/Final_weighted_networks/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1])
    net_file <- data.table::fread(allnets[qq], header=FALSE)
    net_file <- mapProtein(net_file[[1]], net_file[[2]], unet)
    fdata <- data.frame(matrix(ncol=8, nrow=0))
    qdata <- data.frame(matrix(ncol=3, nrow=0))

    wb1 <- openxlsx::createWorkbook(paste0(outdir,'/Supplementary_Table_2_',net_type[qq],'.xlsx'))
    mdn <- c()

    for(k in 1:length(cancer_type)){
        c_type <- cancer_type[k]
        indir <- paste0(in_dir,'/threshold_',cpm_threshold)
        ## ----- CRPEs --------
        gained_surv <- data.table::fread(paste0('../results/Final_survival/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
        c_type,'_Gained_Surv_.txt'))
        gained_surv1 <- gained_surv[gained_surv$pval <= pval_thres, ]
        gs <- igraph::graph_from_data_frame(gained_surv1[,c(1,2)], directed=FALSE)
        lost_surv <- data.table::fread(paste0('../results/Final_survival/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
        c_type,'_Lost_Surv_.txt'))
        lost_surv1 <- lost_surv[lost_surv$pval <= pval_thres, ]
        ls <- igraph::graph_from_data_frame(lost_surv1[,c(1,2)], directed=FALSE)
        pert <- igraph::union(gs, ls)
        allpert <- igraph::as_data_frame(pert)
        fpert <- mapProtein(allpert[[1]], allpert[[2]], unet)
        ### unique PPIs 
        tempgp <- igraph::graph_from_data_frame(fpert[,c(3,4)], directed=FALSE)
        tempgp <- igraph::simplify(tempgp)
        tempgp <- igraph::as_data_frame(tempgp)
        tempdata <- data.frame(matrix(ncol=5, nrow=0))
        counterg <- 0
        countg <- c()
        countgp <- c()
        patg <- c()
        ppic1 <- c()
        ppic2 <- c()
        tot_eeis <- c()
        pert_eeis <- c()
        patientg <- c()
        crpes <- c()
        temp_mdn <- c()
        num_patients <- gained_surv[[7]][1]/gained_surv[[6]][1]
        save_ex <- data.frame(matrix(nrow=0, ncol=5))
        ###--- both gained and lost together
        tempgl <- tempgp

        ##--- for each PPI select the EEIs 
        for(j in 1:length(tempgl[[1]])){

            wh1 <- which(net_file$protein1 == tempgl[[1]][j])
            wh2 <- which(net_file$protein2 == tempgl[[2]][j])
            wha <- intersect(wh1, wh2)
            wh1 <- which(net_file$protein1 == tempgl[[2]][j])
            wh2 <- which(net_file$protein2 == tempgl[[1]][j])
            whb <- intersect(wh1, wh2)
            wh <- union(wha, whb)
            tempxg <- net_file[wh, ] ## Among all EEIs for this PPI which of them are not perturbed 

            ## get crpe
            wh1 <- which(fpert$protein1 == tempgl[[1]][j])
            wh2 <- which(fpert$protein2 == tempgl[[2]][j])
            wha <- intersect(wh1, wh2)
            wh1 <- which(fpert$protein1 == tempgl[[2]][j])
            wh2 <- which(fpert$protein2 == tempgl[[1]][j])
            whb <- intersect(wh1, wh2)
            wh <- union(wha, whb)
            temp_crpe <- fpert[wh, ] ## Among all EEIs for this PPI which of them are not perturbed 
            tempcr_gr <- igraph::graph_from_data_frame(temp_crpe[,c(1,2)], directed=FALSE) ## make a graph of all eeis perturbed in at least one patient

            
            tempa <- getprop(tempxg[[1]], tempxg[[2]], gained_surv[[1]], gained_surv[[2]], gained_surv[[7]])
            tempb <- getprop(tempxg[[1]], tempxg[[2]], lost_surv[[1]], lost_surv[[2]], lost_surv[[7]])

            tempxg$patient <- tempa[[1]]+tempb[[1]] ## perturbed

            ##-- where are the crpes--
            crpe_pos <- c()
            for(i in 1:length(temp_crpe[[1]])){
                wh1 <- which(tempxg$exon1 == temp_crpe[[1]][i])
                wh2 <- which(tempxg$exon2 == temp_crpe[[2]][i])
                wha <- intersect(wh1, wh2)
                wh1 <- which(tempxg$exon2 == temp_crpe[[1]][i])
                wh2 <- which(tempxg$exon1 == temp_crpe[[2]][i])
                whb <- intersect(wh1, wh2)
                wh <- union(wha, whb)
                crpe_pos <- c(crpe_pos, wh)
            }
            

            # tempxg_gr <- igraph::graph_from_data_frame(tempxg[tempxg$patient != 0][,c(1,2)], directed=FALSE) ## make a graph of all eeis perturbed in at least one patient
            # diff1 <- igraph::difference(tempcr_gr, tempxg_gr)
            # diff2 <- igraph::difference(tempxg_gr,tempcr_gr)

            tempyg <- tempxg[tempxg$patient == 0, ]

            if((nrow(tempyg) != 0)){# & igraph::ecount(diff1) == 0 & igraph::ecount(diff2) == 0){ ## meaning that there are some EEIs that are not perturbed in any patient

                counterg  <- counterg+1 ## counting the number of PPIs with at least one non-perturbed EEI
                countg <- c(countg, (length(tempyg[[1]])/length(tempxg[[1]]))*100) ## % of non-perturbed
                tempzg <- tempxg[tempxg$patient != 0, ]
                countgp <- c(countgp, (length(tempzg[[1]])/length(tempxg[[1]]))*100) ## % of perturbed
                patg <- c(patg, mean(tempzg$patient))#/num_patients)*100) ## mean number of patients in which the edges are perturbed
                ppic1 <- c(ppic1, tempzg$protein1[1])
                ppic2 <- c(ppic2, tempzg$protein2[1])
                tot_eeis <- c(tot_eeis, length(tempxg[[1]]))
                pert_eeis <- c(pert_eeis, length(tempzg[[1]]))
                # crpes <- crpes+length(tempag[[1]])
                tempag <- tempxg[crpe_pos,]
                save_ex <- rbind(save_ex, save_impPPIs(tempxg))
                patientg <- c(patientg,(tempag$patient/num_patients)*100)
                temp_mdn <- c(temp_mdn, tempag$patient)
            }
            tempxg$flag <- rep('Perturbed', length(tempxg[[1]]))
            tempdata <- rbind(tempdata, tempxg)
        }

        ##-- store the fractions ---
        if(length(ppic1) != 0){

        tdatag <- data.frame(Protein1=ppic1, Protein2=ppic2, Perturbed=countgp, 
            Number_of_PPIs=rep(length(tempgp[[1]]), length(countgp)), 
            cancer=cancer_type[k], Patients=patg, Non_perturbed=countg, 
            Number_of_desired_PPIs=rep(counterg, length(countgp)), Total_EEIs=tot_eeis, Perturbed_EEIs=pert_eeis)
        fdata <- rbind(fdata, tdatag)

        tdatagg <- data.frame(Cancer=rep(cancer_type[k], length(patientg)), Patients=patientg, 
            Number_of_PPIs=rep(counterg, length(patientg)), Number_of_EEIs=rep(length(patientg), length(patientg)))
        qdata <- rbind(qdata, tdatagg)

        mdn <- c(mdn, median(temp_mdn))

        }

        ## save excel sheet ----
        openxlsx::addWorksheet(wb1, sheetName = cancer_type[k])
        openxlsx::writeData(wb1, sheet = cancer_type[k], save_ex)
        openxlsx::saveWorkbook(wb1, paste0(outdir,'/Supplementary_Table_2_',net_type[qq],'.xlsx'), overwrite = T)
    }

    # #temp
    # temp <- fdata[fdata$cancer == 'LIHC' & fdata$Total_EEIs > 20, ]
    # temp <- fdata[fdata$cancer == 'THCA' & fdata$Total_EEIs > 20, ]
    ufl <- 400
    if(qq==2){
        ufl <- 250
    }else if(qq==3){
        ufl <- 20
    }

    ##--- number of interesting PPIs ------
    pdatax <- unique(fdata[,c(5,8)])
    # p <- ggplot(pdata, aes(Cancer, count, fill=Cancer)) + 
    p <- ggplot(pdatax, aes(cancer, Number_of_desired_PPIs)) + 
    geom_bar(stat="identity",position=position_dodge())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of PPIs", limits=c(0,(max(pdatax$Number_of_desired_PPIs))+5)) +
    geom_text(aes(label=Number_of_desired_PPIs), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
    # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_intPPIs.png"),width=3.5, height=3, dpi=400)


    ##--- number of interesting EEIs ------
    pdatax <- unique(qdata[,c(1,4)])
    # p <- ggplot(pdata, aes(Cancer, count, fill=Cancer)) + 
    p <- ggplot(pdatax, aes(Cancer, Number_of_EEIs)) + 
    geom_bar(stat="identity",position=position_dodge())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of AS-selective cancer \nrelevant perturbed egdes (CRPEs)", limits=c(0,(max(pdatax$Number_of_EEIs))+ufl)) +
    geom_text(aes(label=Number_of_EEIs), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
    # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
    ggsave(p,filename=paste0("../results/Final_results/",net_type[qq],"_intEEIs.png"),width=3.5, height=3, dpi=400)



    # ##--- distribution of total number of EEIs in interesting PPIs and affected EEIs ------
    # pdatax <- fdata
    # pdatax$`# of patients` <- cut(fdata$Patients, breaks = c(1, 5, 10, 20, 30, 45), include.lowest=TRUE)
    # ppx <- ggplot(data = pdatax, aes(x=Total_EEIs, y=Perturbed_EEIs, size = `# of patients`)) + 
    # geom_point(alpha=0.4)+
    # xlab("Total number of EEIs in a PPI")+
    # ylab("# of EEIs in a PPI that are perturbed")+
    # scale_x_continuous()+
    # scale_y_continuous()+
    # theme_bw()+theme(axis.text.x = element_text(size = 8, angle = 60, vjust=1, hjust=1, colour = "black"),
    #     axis.text.y = element_text(size = 8, angle = 0, colour = "black"),
    #     panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    # guides(fill='none')#+guides(size=guide_legend(title="# of patients",ncol=1))
    # # guides(color=guide_legend(title="Cancer type",ncol=3))
    # # guides(color=guide_legend(title="Cancer type", ncol=2,bycol=TRUE), shape=guide_legend(title='Edge type'))
    # ggsave(ppx,filename=paste0(outdir,"/",net_type[qq],"_intPPIs_dist.png"),width=4, height=2.5, dpi=400)

    # if(qq == 3){
    #     col_val <- c('#a6cee3','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#e31a1c','black','#9e0142','#053061')
    # }else{
    #     col_val <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061')
    # }

    col_val <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061')

    ##--- distribution of total number of EEIs in interesting PPIs and affected EEIs ------
    pdatax <- fdata
    pdatax$`# of patients` <- cut(fdata$Patients, breaks = c(1, 5, 10, 20, 30, 45), include.lowest=TRUE)
    ppx <- ggplot(data = pdatax, aes(x=Total_EEIs, y=Perturbed_EEIs, color=cancer, group=cancer, size = `# of patients`)) + 
    geom_point(alpha=0.9)+
    xlab("Total number of EEIs in a PPI")+
    ylab("# of perturbed edges in a PPI")+
    scale_x_continuous(limits=c(0, max(pdatax$Total_EEIs)))+
    scale_y_continuous(limits=c(0, max(pdatax$Perturbed_EEIs)))+
    scale_color_manual(values=col_val)+
    # scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    # geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.01, size=0.1) +
    theme_bw()+theme(axis.text.x = element_text(size = 8, angle = 60, vjust=1, hjust=1, colour = "black"),
        axis.text.y = element_text(size = 8, angle = 0, colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    # guides(fill='none')#+guides(size=guide_legend(title="# of patients",ncol=1))
    guides(color=guide_legend(title="Cancer type",ncol=3))
    # guides(color=guide_legend(title="Cancer type", ncol=2,bycol=TRUE), shape=guide_legend(title='Edge type'))
    ggsave(ppx,filename=paste0(outdir,"/",net_type[qq],"_intPPIs_dist_p.png"),width=7, height=3, dpi=400)


    # pdatax <- fdata
    # pdatax$`# of patients` <- cut(fdata$Patients, breaks = c(1, 5, 10, 20, 30, 45), include.lowest=TRUE)
    # ppx <- ggplot(data = pdatax, aes(x=Total_EEIs, y=Perturbed_EEIs, color=cancer, group=cancer, size = `# of patients`)) + 
    # geom_point(alpha=0.4)+
    # xlab("Total number of EEIs in a PPI")+
    # ylab("# of non-perturbed EEIs in a PPI")+
    # scale_x_continuous()+
    # scale_y_continuous()+
    # scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    # # geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.01, size=0.1) +
    # theme_bw()+theme(axis.text.x = element_text(size = 8, angle = 60, vjust=1, hjust=1, colour = "black"),
    #     axis.text.y = element_text(size = 8, angle = 0, colour = "black"),
    #     panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    # # guides(fill='none')#+guides(size=guide_legend(title="# of patients",ncol=1))
    # guides(color=guide_legend(title="Cancer type",ncol=3))
    # # guides(color=guide_legend(title="Cancer type", ncol=2,bycol=TRUE), shape=guide_legend(title='Edge type'))
    # ggsave(ppx,filename=paste0(outdir,"/",net_type[qq],"_intPPIs_dist_np.png"),width=7, height=4, dpi=400)


    ##--- boxplot showing the fraction of PPIs in which at least one EEI is not perturbed 
    td <- qdata[,c(1,2)]
    mxp <-  td %>% group_by(Cancer) %>% summarise(mv = max(Patients))
    # mdn <-  td %>% group_by(Cancer) %>% summarise(mv = median(Patients))

    pmdn <- data.frame(Cancer=cancer_type,count=mdn,pos=mxp$mv+10)
    ppx <- ggplot(data = qdata, aes(y=Patients, x=Cancer, fill=Cancer, group=Cancer)) + geom_boxplot(alpha=0.9)+
    # geom_point(aes(fill = Cancer), shape = 21, position = position_jitterdodge()) +
    xlab("Cancer type")+
    ylab("% of patients \nwith paired samples")+
    scale_x_discrete()+
    scale_y_continuous(limits=c(0,100))+
    scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    geom_text(data=pmdn, aes(y=pos, x=Cancer,label=count), position=position_dodge(width=0.9),hjust=0.5, vjust=0.5, angle=60, size=3)+
    # geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.01, size=0.1) +
    theme_bw()+theme(axis.text.x = element_text(size = 8, angle = 60, vjust=1, hjust=1, colour = "black"),
        axis.text.y = element_text(size = 8, angle = 0, colour = "black"))+#,
        # panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill='none')#+guides(size=guide_legend(title="# of patients",ncol=1))

    # guides(color=guide_legend(title="Cancer type", ncol=2,bycol=TRUE), shape=guide_legend(title='Edge type'))
    ggsave(ppx,filename=paste0(outdir,"/PPI-vs-EEI_",net_type[qq],".png"),width=4.5, height=2.5, dpi=400)


    # ##--- boxplot showing the fraction of PPIs in which at least one EEI is not perturbed 
    # fdata$`# of patients` <- cut(fdata$Patients, breaks = c(1, 5, 10, 20, 30, 45), include.lowest=TRUE)
    # ppx <- ggplot(data = fdata, aes(y=Perturbed_EEIs, x=cancer, fill=cancer, group=cancer)) + geom_boxplot(alpha=0.4,outlier.shape=NA)+
    # geom_point(aes(fill = cancer, size = `# of patients`), shape = 21, position = position_jitterdodge()) +
    # xlab("Cancer type")+
    # ylab("% of CRPEs in a PPI with \nat least one non-perturbed EEI")+
    # scale_x_discrete()+
    # scale_y_continuous(limits=c(0, 100))+
    # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    # # geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.01, size=0.1) +
    # theme_bw()+theme(axis.text.x = element_text(size = 8, angle = 60, vjust=1, hjust=1, colour = "black"),
    #     axis.text.y = element_text(size = 8, angle = 0, colour = "black"),
    #     panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    # guides(fill='none')+guides(size=guide_legend(title="# of patients",ncol=1))

    # # guides(color=guide_legend(title="Cancer type", ncol=2,bycol=TRUE), shape=guide_legend(title='Edge type'))
    # ggsave(ppx,filename=paste0(outdir,"/PPI-vs-EEI_",net_type[qq],".png"),width=9, height=4, dpi=400)


    # ##--- select the PPIs with less than 25% perturbed in more than 5 patients
    # pdatax <- fdata[fdata$Perturbed_EEIs < 3 & fdata$Patients > 5, ]
    # pdatax$`# of patients` <- cut(pdatax$Patients, breaks = c(1, 5, 10, 20, 30, 45), include.lowest=TRUE)
    # ppx <- ggplot(data = pdatax, aes(x=Total_EEIs, y=Perturbed_EEIs, color=cancer, group=cancer, size = `# of patients`)) + 
    # geom_point(alpha=0.99)+
    # xlab("Total # of EEIs in a PPI")+
    # ylab("# of cancer-relevant perturbed \nedges (CRPEs) in a PPI")+
    # scale_x_continuous(limits=c(0,45))+
    # scale_y_continuous(limits=c(0,45))+
    # scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'),
    #     limits=cancer_type)+
    # # geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.01, size=0.1) +
    # theme_bw()+theme(axis.text.x = element_text(size = 8, angle = 60, vjust=1, hjust=1, colour = "black"),
    #     axis.text.y = element_text(size = 8, angle = 0, colour = "black"),
    #     panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    # # guides(fill='none')#+guides(size=guide_legend(title="# of patients",ncol=1))
    # guides(color=guide_legend(title="Cancer type",ncol=3))
    # # guides(color=guide_legend(title="Cancer type", ncol=2,bycol=TRUE), shape=guide_legend(title='Edge type'))
    # ggsave(ppx,filename=paste0(outdir,"/",net_type[qq],"_intPPIs_dist_s.png"),width=7, height=4, dpi=400)



}


## number of such PPIs and CRPEs ----




# #### survival differences example between patients of THCA 
# ##

# temp_expr <- data.table::fread(paste0('../data/normalized_exons_sv/THCA/normalized_all.txt'), header=TRUE)

# exon1 <- 'ENSE00003604293'
# exon2 <- 'ENSE00001454729'
# wh <- which(temp_expr$EXON %in% c(exon1,exon2))
# temp_exprq <- temp_expr[wh,]

# all_survival <- get_survival('THCA', 0)
# tempids <- data.table::fread(paste0('../data/THCA_manifest_final.txt'))
# tempf <- as.data.frame(temp_exprq)
# sample_names_all <- unlist(lapply(strsplit(colnames(tempf), '[_]'), '[[', 1))
# whc <- which(sample_names_all %in% union(tempids$nid,'EXON'))
# tempf <- tempf[,whc]
# exons <- tempf$EXON ## this order is constant
# tempf <- tempf[,-length(tempf)]
# ##-- only keep patients with control and condition ----------------------------
# temp_brk <- unlist(lapply(strsplit(colnames(tempf), '[_]'), '[[', 2))
# sample_names <- unlist(lapply(strsplit(colnames(tempf), '[_]'), '[[', 1))
# condition_pos_all <- which( (temp_brk == '1') | (temp_brk == '01') )
# control_pos <- which(temp_brk == '11')
# keep_samples <- sample_names[control_pos]
# condition_pos <- intersect(which(sample_names %in% keep_samples), condition_pos_all)
# tempf_control <- tempf[,control_pos]
# tempf_condition <- tempf[,condition_pos]
# tempf_control[tempf_control < 0.5 ] <- 0 
# tempf_condition[tempf_condition < 0.5 ] <- 0

# condition_vec <- rep(1,length(tempf_condition))
# for(k in 1:length(tempf_condition)){
#     tx <- tempf_condition[,k]
#     if(tx[[1]] == 0 | tx[[2]] == 0){
#         condition_vec[k] <- 0
#     }
# }

# control_vec <- rep(1,length(tempf_condition))
# for(k in 1:length(tempf_condition)){
#     tx <- tempf_control[,k]
#     if(tx[[1]] == 0 | tx[[2]] == 0){
#         control_vec[k] <- 0
#     }
# }

# pert <- xor(control_vec, condition_vec)

# pert_patients <- unlist(lapply(strsplit(colnames(tempf_control)[which(pert == TRUE)],'[_]'),'[[',1))
# npert_patients <- unlist(lapply(strsplit(colnames(tempf_control)[which(pert == FALSE)],'[_]'),'[[',1))

# pert_surv <- all_survival[all_survival$submitter_id %in% pert_patients, ]$overall_survival
# npert_surv <- all_survival[all_survival$submitter_id %in% npert_patients, ]$overall_survival

# pv <- wilcox.test(pert_surv, npert_surv)

# # plot
# pdata1 <- data.frame(ID=rep('Perturbed', length(pert_surv)), Survival=pert_surv)
# pdata2 <- data.frame(ID=rep('Not perturbed', length(npert_surv)), Survival=npert_surv)
# pdata <- rbind(pdata1, pdata2)

# ppx <- ggplot(data = pdata, aes(x=ID, y=Survival, color=ID)) + 
# geom_boxplot()+
# xlab("")+
# ylab("Overall survival")+
# scale_x_discrete()+
# scale_y_continuous()+
# theme_bw()+theme(axis.text.x = element_text(size = 8, angle = 60, vjust=1, hjust=1, colour = "black"),
#     axis.text.y = element_text(size = 8, angle = 0, colour = "black"),
#     panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
# guides(color='none')
# ggsave(ppx,filename=paste0(outdir,"/THCA_ENSE00003604293_ENSE00001454729.png"),width=2, height=4, dpi=400)

# # EEI between \nENSE00003604293\n and ENSE00001454729

