
##--- Generate data details --------

rm(list=ls())
library(ggplot2)
library(Rcpp)
library('biomaRt')
source("eein_cancer_util.r")

outdir <- '../data/Final_results'

if(!dir.exists(outdir)){
    dir.create(outdir, recursive=TRUE)
}

cppFunction("List setAss(CharacterVector ex1, CharacterVector ex2, CharacterVector exon1, CharacterVector exon2, CharacterVector id){

    int loop1 = ex1.size();
    int loop2 = exon1.size();
    CharacterVector allids;

    for(int k=0; k<loop1; k++){
        for(int j=0; j<loop2; j++){
            if((ex1[k] == exon1[j]) & (ex2[k] == exon2[j])){
                allids.push_back(id[j]);
                break;
            }
            if((ex1[k] == exon2[j]) & (ex2[k] == exon1[j])){
                allids.push_back(id[j]);
                break;
            }
        }
    }
    List L = List::create(allids);
    return L;
}")

pval_thres <- 0.05 # threshold to consider an edge to be survival correlated
cancer_type <- c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 
    'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA')

##------------ select networks for each method ----------------------------------------------------------------------------
##--- contact-based -------------------------------------------------------------------------------------------------------
cpm_threshold <- c(0.05, 0.1, 0.2, 0.5, 1, 2)

for(thr in 1:length(cpm_threshold)){

    indir <- paste0('../data/CONTACT_survival/threshold_',cpm_threshold[thr])
    allnet1 <- data.frame(matrix(ncol=2, nrow=0))

    contact_pert <- c()
    back_net <- data.table::fread(paste0('../data/CONTACT_networks/CONTACT_net_6_1.txt'))
    back_net <- igraph::graph_from_data_frame(back_net[,c(3,4)], directed=FALSE)

    for(k in 1:length(cancer_type)){

        templ <- data.table::fread(paste0(indir,'/',cancer_type[k],'_Lost_Surv_.txt'))
        tempg <- data.table::fread(paste0(indir,'/',cancer_type[k],'_Gained_Surv_.txt'))

        templx1 <- templ[(templ$pval <= pval_thres) ,]
        tempgx1 <- tempg[(tempg$pval <= pval_thres) ,]

        tempn1 <- rbind(templx1[,c(1,2)], tempgx1[,c(1,2)])
        allnet1 <- rbind(allnet1, tempn1)

        ##-- percentage kept
        temp_net <- igraph::graph_from_data_frame(tempn1, directed=FALSE)
        contact_pert <- c(contact_pert, igraph::ecount(temp_net)/igraph::ecount(back_net))

    }
        
    allnet_final <- igraph::graph_from_data_frame(allnet1, directed=FALSE)
    allnet_final <- igraph::simplify(allnet_final, remove.multiple=TRUE)
    allnet_final1 <- igraph::as_data_frame(allnet_final)
    data.table::fwrite(allnet_final1, paste0('../data/CONTACT_survival/CONTACT_net_final_',cpm_threshold[thr],'.txt'), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

}


##--- pisa-based -------------------------------------------------------------------------------------------------------------------------

for(thr in 1:length(cpm_threshold)){

    indir <- paste0('../data/PISA_survival/threshold_',cpm_threshold[thr])
    allnet1 <- data.frame(matrix(ncol=2, nrow=0))

    pisa_pert <- c()
    back_net <- data.table::fread(paste0('../data/PISA_networks_filtered/PISA_EEIN_0.5.txt'))
    back_net <- igraph::graph_from_data_frame(back_net[,c(1,2)], directed=FALSE)

    for(k in 1:length(cancer_type)){

        templ <- data.table::fread(paste0(indir,'/',cancer_type[k],'_Lost_Surv_.txt'))
        tempg <- data.table::fread(paste0(indir,'/',cancer_type[k],'_Gained_Surv_.txt'))

        templx1 <- templ[(templ$pval <= pval_thres) ,]
        tempgx1 <- tempg[(tempg$pval <= pval_thres) ,]

        tempn1 <- rbind(templx1[,c(1,2)], tempgx1[,c(1,2)])

        allnet1 <- rbind(allnet1, tempn1)

        ##-- percentage kept
        temp_net <- igraph::graph_from_data_frame(tempn1, directed=FALSE)
        pisa_pert <- c(pisa_pert, igraph::ecount(temp_net)/igraph::ecount(back_net))

    }
    
    allnet_final <- igraph::graph_from_data_frame(allnet1, directed=FALSE)
    allnet_final <- igraph::simplify(allnet_final, remove.multiple=TRUE)
    allnet_final1 <- igraph::as_data_frame(allnet_final)
    data.table::fwrite(allnet_final1, paste0('../data/PISA_survival/PISA_net_final_',cpm_threshold[thr],'.txt'), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

}


##--- eppic-based -------------------------------------------------------------------------------------------------------------------------

for(thr in 1:length(cpm_threshold)){

    indir <- paste0('../data/EPPIC_survival/threshold_',cpm_threshold[thr])

    allnet1 <- data.frame(matrix(ncol=2, nrow=0))

    eppic_pert <- c()
    back_net <- data.table::fread('../data/EPPIC_EEIN_filtered.txt')
    back_net <- igraph::graph_from_data_frame(back_net[,c(1,2)], directed=FALSE)

    for(k in 1:length(cancer_type)){

        templ <- data.table::fread(paste0(indir,'/',cancer_type[k],'_Lost_Surv_.txt'))
        tempg <- data.table::fread(paste0(indir,'/',cancer_type[k],'_Gained_Surv_.txt'))

        templx1 <- templ[(templ$pval <= pval_thres) ,]
        tempgx1 <- tempg[(tempg$pval <= pval_thres) ,]

        tempn1 <- rbind(templx1[,c(1,2)], tempgx1[,c(1,2)])

        allnet1 <- rbind(allnet1, tempn1)

        ##-- percentage kept
        temp_net <- igraph::graph_from_data_frame(tempn1, directed=FALSE)
        eppic_pert <- c(eppic_pert, igraph::ecount(temp_net)/igraph::ecount(back_net))

    }
        
    allnet_final <- igraph::graph_from_data_frame(allnet1, directed=FALSE)
    allnet_final <- igraph::simplify(allnet_final, remove.multiple=TRUE)
    allnet_final1 <- igraph::as_data_frame(allnet_final)
    data.table::fwrite(allnet_final1, paste0('../data/EPPIC_survival/EPPIC_net_final_',cpm_threshold[thr],'.txt'), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

}



###----- SELECT THE CPM THRESHOLD -----------------------------------------------------------------------------------------------------------
###--- NETLOW -----
num_edges_nl <- c()
for(thr in 1:length(cpm_threshold)){

    contact <- data.table::fread(paste0('../data/CONTACT_survival/CONTACT_net_final_',cpm_threshold[thr],'.txt'), header=FALSE)
    pisa <- data.table::fread(paste0('../data/PISA_survival/PISA_net_final_',cpm_threshold[thr],'.txt'), header=FALSE)
    eppic <- data.table::fread(paste0('../data/EPPIC_survival/EPPIC_net_final_',cpm_threshold[thr],'.txt'))

    aa_net <- igraph::simplify(igraph::graph_from_data_frame(contact[,1:2], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
    eppic_net <- igraph::simplify(igraph::graph_from_data_frame(eppic[,1:2], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
    PISA_net <- igraph::simplify(igraph::graph_from_data_frame(pisa[,c(1,2)], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)

    entire_net <- igraph::union(igraph::union(PISA_net, aa_net), eppic_net)

    num_edges_nl <- c(num_edges_nl, igraph::ecount(entire_net))

}

###--- NETMEDIUM -----
num_edges_nm <- c()
for(thr in 1:length(cpm_threshold)){

    contact <- data.table::fread(paste0('../data/CONTACT_survival/CONTACT_net_final_',cpm_threshold[thr],'.txt'), header=FALSE)
    pisa <- data.table::fread(paste0('../data/PISA_survival/PISA_net_final_',cpm_threshold[thr],'.txt'), header=FALSE)
    eppic <- data.table::fread(paste0('../data/EPPIC_survival/EPPIC_net_final_',cpm_threshold[thr],'.txt'))

    aa_net <- igraph::simplify(igraph::graph_from_data_frame(contact[,1:2], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
    eppic_net <- igraph::simplify(igraph::graph_from_data_frame(eppic[,1:2], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
    PISA_net <- igraph::simplify(igraph::graph_from_data_frame(pisa[,c(1,2)], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)

    entire_net <- igraph::union(igraph::union(igraph::intersection(PISA_net,eppic_net,keep.all.vertices=FALSE),
        igraph::intersection(aa_net,eppic_net,keep.all.vertices=FALSE)), 
        igraph::intersection(aa_net,PISA_net,keep.all.vertices=FALSE))


    num_edges_nm <- c(num_edges_nm, igraph::ecount(entire_net))

}

###--- NETHIGH -----
num_edges_nh <- c()
for(thr in 1:length(cpm_threshold)){

    contact <- data.table::fread(paste0('../data/CONTACT_survival/CONTACT_net_final_',cpm_threshold[thr],'.txt'), header=FALSE)
    pisa <- data.table::fread(paste0('../data/PISA_survival/PISA_net_final_',cpm_threshold[thr],'.txt'), header=FALSE)
    eppic <- data.table::fread(paste0('../data/EPPIC_survival/EPPIC_net_final_',cpm_threshold[thr],'.txt'))

    aa_net <- igraph::simplify(igraph::graph_from_data_frame(contact[,1:2], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
    eppic_net <- igraph::simplify(igraph::graph_from_data_frame(eppic[,1:2], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
    PISA_net <- igraph::simplify(igraph::graph_from_data_frame(pisa[,c(1,2)], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)

    entire_net <- igraph::intersection(igraph::intersection(PISA_net,eppic_net,keep.all.vertices=FALSE),
    aa_net,keep.all.vertices=FALSE) 

    num_edges_nh <- c(num_edges_nh, igraph::ecount(entire_net))

}


##--- DRAW FIGURE -----------------
pdata1 <- data.frame(v1=cpm_threshold, v2=num_edges_nl)
pdata1$flag <- rep('NETLOW',length(cpm_threshold))
pdata2 <- data.frame(v1=cpm_threshold, v2=num_edges_nm)
pdata2$flag <- rep('NETMEDIUM',length(cpm_threshold))
pdata3 <- data.frame(v1=cpm_threshold, v2=num_edges_nh)
pdata3$flag <- rep('NETHIGH',length(cpm_threshold))

pdata4 <- rbind(pdata1, pdata2)
pdata <- rbind(pdata4, pdata3)

basesize  <- 10
ppx <- ggplot(data = pdata3, aes(x = v1, y=v2, color=flag)) + geom_point()+geom_line()+
xlab("Exon expression threshold\n(normalized CPM value)")+
scale_color_manual(values=c('#377eb8','#4daf4a','#e41a1c'))+
scale_y_continuous(name="# of cancer relevant edges") +
theme_bw(base_size = basesize * 0.8)+
theme(axis.text.x = element_text(size = 0.8*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 0.8*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"), 
    axis.title = element_text(size=0.8*basesize), 
    legend.text=element_text(size=0.8*basesize))+
guides(color="none")

##--- guides(color=guide_legend(title="Network"))
ggsave(ppx,filename=paste0("../results/Final_results/NETHIGH_cpm_thres_cancer_relevant_edges.png"),width=3.5, height=2.5, dpi=400)


basesize  <- 10
ppx <- ggplot(data = pdata2, aes(x = v1, y=v2, color=flag)) + geom_point()+geom_line()+
xlab("Exon expression threshold\n(normalized CPM value)")+
scale_color_manual(values=c('#377eb8','#4daf4a','#e41a1c'))+
scale_y_continuous(name="# of cancer relevant edges") +
theme_bw(base_size = basesize * 0.8)+
theme(axis.text.x = element_text(size = 0.8*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 0.8*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"), 
    axis.title = element_text(size=0.8*basesize), 
    legend.text=element_text(size=0.8*basesize))+
guides(color="none")

# guides(color=guide_legend(title="Network"))
ggsave(ppx,filename=paste0("../results/Final_results/NETMEDIUM_cpm_thres_cancer_relevant_edges.png"),width=3.5, height=2.5, dpi=400)


basesize  <- 10
ppx <- ggplot(data = pdata1, aes(x = v1, y=v2, color=flag)) + geom_point()+geom_line()+
xlab("Exon expression threshold\n(normalized CPM value)")+
scale_color_manual(values=c('#377eb8','#4daf4a','#e41a1c'))+
scale_y_continuous(name="# of cancer relevant edges") +
theme_bw(base_size = basesize * 0.8)+
theme(axis.text.x = element_text(size = 0.8*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
    axis.text.y = element_text(size = 0.8*basesize, angle = 0, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"), 
    axis.title = element_text(size=0.8*basesize), 
    legend.text=element_text(size=0.8*basesize))+
guides(color="none")

# guides(color=guide_legend(title="Network"))
ggsave(ppx,filename=paste0("../results/Final_results/NETLOW_cpm_thres_cancer_relevant_edges.png"),width=3.5, height=2.5, dpi=400)






best_thr <- pdata[[1]][which(pdata$v2 == max(pdata$v2))]

###----- SAVE FINAL NETWORKS ----------------------------------------------------------------------------------------------------------------

##------ 1. Details of network size ------
contact <- data.table::fread(paste0('../data/CONTACT_survival/CONTACT_net_final_',best_thr,'.txt'), header=FALSE)
net_file <- data.table::fread(paste0('../data/CONTACT_networks/CONTACT_net_6_1.txt'))
net_file[[1]] <- unlist(lapply(strsplit(net_file[[1]], '[_]'), '[[', 1))
net_file[[2]] <- unlist(lapply(strsplit(net_file[[2]], '[_]'), '[[', 1))
contact <- mapProtein(contact[[1]], contact[[2]], net_file[,c(3,4,1,2)])
contact_eein <- igraph::simplify(igraph::graph_from_data_frame(contact[,c(1,2)], directed=FALSE))
contact_ppin <- igraph::simplify(igraph::graph_from_data_frame(contact[,c(3,4)], directed=FALSE))

pisa <- data.table::fread(paste0('../data/PISA_survival/PISA_net_final_',best_thr,'.txt'), header=FALSE)
net_file <- data.table::fread(paste0('../data/PISA_networks_filtered/PISA_EEIN_0.5.txt'))
pisa <- mapProtein(pisa[[1]], pisa[[2]], net_file[,c(1,2,5,6)])
pisa_eein <- igraph::simplify(igraph::graph_from_data_frame(pisa[,c(1,2)], directed=FALSE))
pisa_ppin <- igraph::simplify(igraph::graph_from_data_frame(pisa[,c(3,4)], directed=FALSE))

eppic <- data.table::fread(paste0('../data/EPPIC_survival/EPPIC_net_final_',best_thr,'.txt'))
net_file <- data.table::fread(paste0('../data/EPPIC_EEIN_filtered.txt'))
eppic <- mapProtein(eppic[[1]], eppic[[2]], net_file[,c(1,2,5,6)])
eppic_eein <- igraph::simplify(igraph::graph_from_data_frame(eppic[,c(1,2)], directed=FALSE))
eppic_ppin <- igraph::simplify(igraph::graph_from_data_frame(eppic[,c(3,4)], directed=FALSE))

# ##---- 2. Overlap of the different networks ---------------------------------------------------
aa_net <- igraph::simplify(igraph::graph_from_data_frame(contact[,1:2], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
eppic_net <- igraph::simplify(igraph::graph_from_data_frame(eppic[,1:2], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
PISA_net <- igraph::simplify(igraph::graph_from_data_frame(pisa[,c(1,2)], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)

entire_net <- igraph::as_data_frame(igraph::union(igraph::union(PISA_net, aa_net), eppic_net))
entire_net$ID <- paste0('Id',seq(1,length(entire_net[[1]])))
##---------------------------------------------------------------------------------------------
##-------- 3. save the three different networks -----------------------------------------------
out_dir <- paste0('../data/Final_networks/level1')
dir.create(out_dir, recursive=TRUE)

##-- entire net -- one or more network definitions --
data.table::fwrite(entire_net[,c(1,2)], paste0(out_dir,'/atleast_1.txt'), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
##-- two or more network definitions --
int_net1 <- igraph::intersection(PISA_net, eppic_net)
int_net2 <- igraph::intersection(PISA_net, aa_net)
int_net3 <- igraph::intersection(aa_net, eppic_net)
fnet <- igraph::union(igraph::union(int_net1, int_net2), int_net3)
data.table::fwrite(igraph::as_data_frame(fnet), paste0(out_dir,'/atleast_2.txt'), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
##-- all three network definitions ---
int_net <- igraph::intersection(igraph::intersection(PISA_net, eppic_net), aa_net)
data.table::fwrite(igraph::as_data_frame(int_net), paste0(out_dir,'/atleast_3.txt'), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

##----- onerlap of networks for different confidence ---
net_file <- data.table::fread(paste0('../data/CONTACT_networks/CONTACT_net_6_1.txt'))
net_file$p1 <- unlist(lapply(strsplit(net_file[[1]], '[_]'), '[[', 1))
net_file$p2 <- unlist(lapply(strsplit(net_file[[2]], '[_]'), '[[', 1))
contact <- net_file[,c(3,4,12,13)]
contact_eein <- igraph::simplify(igraph::graph_from_data_frame(contact[,c(1,2)], directed=FALSE))
contact_ppin <- igraph::simplify(igraph::graph_from_data_frame(contact[,c(3,4)], directed=FALSE))

net_file <- data.table::fread(paste0('../data/PISA_networks_filtered/PISA_EEIN_0.5.txt'))
pisa <- net_file[,c(1,2,5,6)]
pisa_eein <- igraph::simplify(igraph::graph_from_data_frame(pisa[,c(1,2)], directed=FALSE))
pisa_ppin <- igraph::simplify(igraph::graph_from_data_frame(pisa[,c(3,4)], directed=FALSE))

net_file <- data.table::fread(paste0('../data/EPPIC_EEIN_filtered.txt'))
eppic <- net_file[,c(1,2,5,6)]
eppic_eein <- igraph::simplify(igraph::graph_from_data_frame(eppic[,c(1,2)], directed=FALSE))
eppic_ppin <- igraph::simplify(igraph::graph_from_data_frame(eppic[,c(3,4)], directed=FALSE))

# ##---- 2. Overlap of the different networks ---------------------------------------------------
aa_net <- igraph::simplify(igraph::graph_from_data_frame(contact[,1:2], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
eppic_net <- igraph::simplify(igraph::graph_from_data_frame(eppic[,1:2], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
PISA_net <- igraph::simplify(igraph::graph_from_data_frame(pisa[,c(1,2)], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)

aa_net1 <- igraph::as_data_frame(aa_net)
eppic_net1 <- igraph::as_data_frame(eppic_net)
PISA_net1 <- igraph::as_data_frame(PISA_net)

#--- Overlap ------
entire_net <- igraph::as_data_frame(igraph::union(igraph::union(PISA_net, aa_net), eppic_net))
entire_net$ID <- paste0('Id',seq(1,length(entire_net[[1]])))
setA_contact <- setAss(aa_net1[[1]], aa_net1[[2]], entire_net[[1]], entire_net[[2]], entire_net[[3]])
setB_eppic <- setAss(eppic_net1[[1]], eppic_net1[[2]], entire_net[[1]], entire_net[[2]], entire_net[[3]])
setC_pisa <- setAss(PISA_net1[[1]], PISA_net1[[2]], entire_net[[1]], entire_net[[2]], entire_net[[3]])

x <- list(`Contact-based` = setA_contact[[1]], `Evolution-based` = setB_eppic[[1]], `Energy-based` = setC_pisa[[1]])
p <- ggvenn::ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),stroke_size = 0.5, set_name_size = 6)

png(filename=paste0(outdir,"/Network_overlaps.png"), width=5, height=5, units="in", res=500)
print(p)
dev.off()



