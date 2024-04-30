##############################################################################################
# Purpose: numbers of gained and lost CRPEs
# for each category
## categories are:
# a: Non-perturbed
# b: Only gained
# c: Only lost
# d: Gained/lost
# e: Patient-specific gain
# f: Patient-specific loss
# These categories partitions the edges of the entire network into six groups
## Visualize the numbers of edges that significantly correlate with survival
# in these categories for each cancer type
##############################################################################################

rm(list=ls())
library(data.table)
library(ggplot2)
library(RColorBrewer)

cancer_type <- c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA')
out_dir <- '../results/Final_results'
cpm_threshold <- 0.5
pval_thres <- 0.05
allnets <- gtools::mixedsort(list.files('../data/Final_networks/level1',full.names=TRUE))
net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')

a <- c()
b <- c()
c <- c()
d <- c()
e <- c()
f <- c()

x <- c()
y <- c()
cancer <- c()
netl <- c()

mut_var <- c()
mut_strict <- c()

for(qq in 1:length(allnets)){

	in_dir <- paste0('../data/Final_weighted_networks/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1])
	net_file <- data.table::fread(allnets[qq], header=FALSE)
    perturbed_graph <- list()

	for(gg in 1:length(cancer_type)){
		c_type <- cancer_type[gg]
        indir <- paste0(in_dir,'/threshold_',cpm_threshold)
        gained <- data.table::fread(paste0(indir,'/',c_type,'_gained.txt'))
        lost <- data.table::fread(paste0(indir,'/',c_type,'_lost.txt'))
        ##-- a) Non-perturbed edges --------------
        nonp1 <- gained[gained$patient == 0, ]
        nonp1 <- igraph::graph_from_data_frame(nonp1[,c(1,2)], directed=FALSE)
        nonp2 <- lost[lost$patient == 0, ]
        nonp2 <- igraph::graph_from_data_frame(nonp2[,c(1,2)], directed=FALSE)
        nonp <- igraph::intersection(nonp1, nonp2)
        ##-- e) f) Patient-spcific gained and lost ----
        ps1 <- gained[gained$patient == 1, ]
        ps1 <- igraph::graph_from_data_frame(ps1[,c(1,2)], directed=FALSE)
        lost_ig <- igraph::graph_from_data_frame(lost[lost$patient > 0, ][,c(1,2)], directed=FALSE)
        ps1 <- igraph::difference(ps1, lost_ig) ## e) patient-specific gain
        ps2 <- lost[lost$patient == 1, ]
        ps2 <- igraph::graph_from_data_frame(ps2[,c(1,2)], directed=FALSE)
        gained_ig <- igraph::graph_from_data_frame(gained[gained$patient > 0, ][,c(1,2)], directed=FALSE)
        ps2 <- igraph::difference(ps2, gained_ig) ## f) patient-specific loss
        ##-- b) c) only gained/ only lost
        og <- gained[gained$patient > 1, ]
        og <- igraph::graph_from_data_frame(og[,c(1,2)], directed=FALSE)
        og <- igraph::difference(og, lost_ig) ## b) only gained --. gained in at least 2, but never lost
		ol <- lost[lost$patient > 1, ]
        ol <- igraph::graph_from_data_frame(ol[,c(1,2)], directed=FALSE)
        ol <- igraph::difference(ol, gained_ig)
        ##-- d) gained/lost
        gl <- igraph::intersection(lost_ig, gained_ig) ## d) gained/lost edges
        xx <- igraph::ecount(gl)+igraph::ecount(ol)+igraph::ecount(og)+igraph::ecount(ps1)+igraph::ecount(ps2)+igraph::ecount(nonp)
        if(xx != length(net_file[[1]])){
        	cat('Error')
        	break
        } 
        ##--- get survival results ------------
        gained_surv <- data.table::fread(paste0('../results/Final_survival/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
            c_type,'_Gained_Surv_.txt'))
        gained_surv1 <- gained_surv[gained_surv$pval <= pval_thres, ]
        # gained_surv1 <- gained_surv1[gained_surv1$frac1 <= pval_thres, ]
        # gained_surv1 <- gained_surv1[gained_surv1$frac2 <= pval_thres, ]
        gs <- igraph::graph_from_data_frame(gained_surv1[,c(1,2)], directed=FALSE)
        gs_og <- igraph::intersection(gs, og) ## only gained
        gs_ps1 <- igraph::intersection(gs, ps1) ## patient-specific
        ##-------------------------------------
        ##--- get survival results ------------
        lost_surv <- data.table::fread(paste0('../results/Final_survival/level1/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
            c_type,'_Lost_Surv_.txt'))
        lost_surv1 <- lost_surv[lost_surv$pval <= pval_thres, ]
        # lost_surv1 <- lost_surv1[lost_surv1$frac1 <= pval_thres, ]
        # lost_surv1 <- lost_surv1[lost_surv1$frac2 <= pval_thres, ]
        ls <- igraph::graph_from_data_frame(lost_surv1[,c(1,2)], directed=FALSE)
        ls_ol <- igraph::intersection(ls, ol) ## only lost
        ls_ps2 <- igraph::intersection(ls, ps2) ## patient-specific
        ##-------------------------------------
        gs_ls_u <- igraph::union(ls, gs)
        gs_gl <- igraph::intersection(gs_ls_u, gl) ## gained/lost
        ##-- store -------------------------------
        a <- c(a, igraph::ecount(gs_og)/igraph::ecount(gs_ls_u))
        c <- c(c, igraph::ecount(gs_gl)/igraph::ecount(gs_ls_u))
        d <- c(d, igraph::ecount(gs_ps1)/igraph::ecount(gs_ls_u))
        x <- c(x, igraph::ecount(gs_ls_u))
        b <- c(b, igraph::ecount(ls_ol)/igraph::ecount(gs_ls_u))
        e <- c(e, igraph::ecount(ls_ps2)/igraph::ecount(gs_ls_u))
        y <- c(y, igraph::ecount(gs_ls_u))
        cancer <- c(cancer, cancer_type[gg])
        netl <- c(netl, strsplit(basename(allnets[qq]),'[.]')[[1]][1])
        tempss <- igraph::union(gs_og, gs_gl)
        tempzz <- igraph::union(tempss, gs_ps1)
        tempjj <- igraph::union(tempzz, ls_ol)
        perturbed_graph[[gg]] <- igraph::union(tempjj, ls_ps2)

        mut_var <- c(mut_var, igraph::ecount(gs_gl))
        mut_strict <- c(mut_strict, igraph::ecount(igraph::union(gs_og, ls_ol)))

	}

    all_pert <- perturbed_graph[[1]]
    for(cn in 2:length(cancer_type)){
        all_pert <- igraph::union(all_pert, perturbed_graph[[cn]])
    }

    # print(igraph::ecount(all_pert))

    # print(igraph::ecount(all_pert)/length(net_file[[1]]))


}


pData <- data.frame(A=a*100, B=b*100, C=c*100, D=d*100, E=e*100, X=x, Cancer=cancer, Network=netl)

lbls <- c('atleast_1', 'atleast_2', 'atleast_3')

for(k in 1:length(lbls)){

    pData1 <- pData[pData$Network == lbls[k], ]
    if(k == 3){print(pData1)}
    pData1 <- reshape2::melt(pData1,id=c('Cancer','Network','X'))
    pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C", "D", "E","X"))  # setting level order
    basesize <- 8
    ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
    geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
    scale_y_continuous(limits=c(0,110), breaks = seq(0, 100, by = 10))+
    xlab("Cancer type")+ylab("% of cancer-relevant \nperturbed edges (CRPEs)")+
    scale_fill_manual(labels=c("A" = "Strictly gained", 
        "B"="Strictly lost", "C"="Gained or lost","D"="Patient-specific gained","E"="Patient-specific lost"), 
    values=c('#e41a1c','#006837','#fdc086','#ff7f00','#7fc97f'))+
    theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
        axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
    guides(fill=guide_legend(title="Edge category"))
    ggsave(ppx,filename=paste0(out_dir, "/CRPE_categories_",lbls[k],".png"),width=7, height=3, dpi=400)

}


# ##--- correlation plots ----------------

# # pdatax <- data.frame(gl=mut_var, st=mut_strict, Cancer=cancer, Network=netl)
# pdatax <- data.frame(gl=c*100, st=(a+b)*100, Cancer=cancer, Network=netl)

# for(k in 1:length(lbls)){

#     pdata <- pdatax[pdatax$Network==lbls[k], ]
#     pdata <- pdata[order(pdata$Cancer), ]

#     ### ---------------------- correlation with known number of mutations ------------------
#     mutations <- data.table::fread('../data/mutations.csv')
#     mutations <- mutations[order(mutations$V1), ]
#     pdata$mutations <- mutations$V2

#     # coral <- cor.test(x=pdata$count, y=pdata$mutations, method = 'spearman')
#     p <- ggplot(pdata, aes(mutations, gl, color=Cancer)) + 
#     geom_point(size=2)+
#     theme(legend.text=element_text(size=12))
#     basesize <- 12
#     p <- p + theme_bw(base_size = basesize * 0.8) +
#     scale_x_continuous(name="Median # of somatic mutations \n per-million base", limits=c(0,max(pdata$mutations))) + 
#     scale_y_continuous(name="% of CRPEs in the \ngained or lost category", limits=c(0,max(pdata$gl))) +
#     scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
#     theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
#     axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
#     strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
#     guides(color=guide_legend(title="Cancer type",ncol=2))
#     ggsave(p,filename=paste0("../results/Final_results/Mutations_",net_type[qq],"_gl.png"),width=4.5, height=3, dpi=400)

#     p <- ggplot(pdata, aes(mutations, st, color=Cancer)) + 
#     geom_point(size=2)+
#     theme(legend.text=element_text(size=12))
#     basesize <- 12
#     p <- p + theme_bw(base_size = basesize * 0.8) +
#     scale_x_continuous(name="Median # of somatic mutations \n per-million base", limits=c(0,max(pdata$mutations))) + 
#     scale_y_continuous(name="% of CRPEs in the strictly gained + \nstrictly lost categories", limits=c(0,max(pdata$st))) +
#     scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
#     theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
#     axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
#     strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
#     guides(color=guide_legend(title="Cancer type",ncol=2))
#     ggsave(p,filename=paste0("../results/Final_results/Mutations_",net_type[qq],"_strict.png"),width=4.5, height=3, dpi=400)
#     ###----------------------------------------------------------------

#     ### ---------------------- correlation with known degree of heterogeneity ------------------
#     mutations <- data.table::fread('../data/heterogeniety.csv')
#     mutations <- mutations[order(mutations$V1), ]
#     pdata$mutations <- mutations$V2

#     # coral <- cor.test(x=pdata$count, y=pdata$mutations, method = 'spearman')
#     p <- ggplot(pdata, aes(mutations, gl, color=Cancer)) + 
#     geom_point(size=2)+
#     theme(legend.text=element_text(size=12))
#     basesize <- 12
#     p <- p + theme_bw(base_size = basesize * 0.8) +
#     scale_x_continuous(name="Median # of clones", limits=c(0,max(pdata$mutations))) + 
#     scale_y_continuous(name="% of CRPEs in the \ngained or lost category", limits=c(0,max(pdata$gl))) +
#     scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
#     theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
#     axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
#     strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
#     guides(color=guide_legend(title="Cancer type",ncol=2))
#     ggsave(p,filename=paste0("../results/Final_results/Heterogeneity_",net_type[qq],"_gl.png"),width=4.5, height=3, dpi=400)

#      # coral <- cor.test(x=pdata$count, y=pdata$mutations, method = 'spearman')
#     p <- ggplot(pdata, aes(mutations, st, color=Cancer)) + 
#     geom_point(size=2)+
#     theme(legend.text=element_text(size=12))
#     basesize <- 12
#     p <- p + theme_bw(base_size = basesize * 0.8) +
#     scale_x_continuous(name="Median # of clones", limits=c(0,max(pdata$mutations))) + 
#     scale_y_continuous(name="% of CRPEs in the strictly gained + \nstrictly lost categories", limits=c(0,max(pdata$st))) +
#     scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
#     theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
#     axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
#     strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
#     guides(color=guide_legend(title="Cancer type",ncol=2))
#     ggsave(p,filename=paste0("../results/Final_results/Heterogeneity_",net_type[qq],"_strict.png"),width=4.5, height=3, dpi=400)


#     ###----------------------------------------------------------------

# }
