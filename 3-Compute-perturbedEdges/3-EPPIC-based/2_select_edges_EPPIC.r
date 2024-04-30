##--- Analysis for differences among the cancer types ------
##--- Only run survival analysis to select the edges --------

# rm(list=ls())
args=commandArgs(TRUE)

library(data.table)
library(pvclust)
library('biomaRt')
library(GenomicDataCommons)
library("survival")
library("survminer")
source("eein_cancer_util.r")
library(Rcpp)

c_type <- args[1]#'ESCA'#
outdir3 <- args[2]#'../results/CONTACT_survival'#
pval_thres <- 0.05
cpm_threshold <- c(0.05, 0.1, 0.2, 0.5, 1, 2)

patient <- 0
rnd_run <- 100

in_dir <- paste0('../data/EPPIC_weighted_networks')
net_file <- data.table::fread('../data/EPPIC_EEIN_filtered.txt')
nodes <- union(net_file[[1]], net_file[[2]])

all_expr_data <- data.table::fread(paste0('../data/normalized_exons_sv/', c_type,'/normalized_all.txt'), header=TRUE)
all_survival <- get_survival(c_type, 0)

my.t.test.p.value <- function(...) {
   obj<-try(t.test(...), silent=TRUE)
   if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

for(thr in 1:length(cpm_threshold)){

    indir <- paste0(in_dir,'/threshold_',cpm_threshold[thr])
    gained <- data.table::fread(paste0(indir,'/',c_type,'_gained.txt'))
    ##-- edges with perturbations ----
    wh <- which(gained$patient > patient)
    gained_pert <- gained[wh, ]
    lost <- data.table::fread(paste0(indir,'/',c_type,'_lost.txt'))
    ##-- edges with perturbations ----
    wh <- which(lost$patient > patient)
    lost_pert <- lost[wh, ]
    ##-- only gained / lost / both ---
    gained_pert_g <- igraph::graph_from_data_frame(gained_pert, directed=FALSE)
    lost_pert_g <- igraph::graph_from_data_frame(lost_pert, directed=FALSE)
    only_gained <- gained_pert_g
    only_lost <- lost_pert_g

    ##-- get the condition samples only ---------------------------------------------------------------------------------------
    temp1 <- all_expr_data[all_expr_data$EXON %in% nodes, ]
    exons <- temp1$EXON
    temp2 <- temp1[,EXON:=NULL]
    temp2 <- as.data.frame(temp2)
    temp_brk <- unlist(lapply(strsplit(colnames(temp2), '[_]'), '[[', 2))
    sample_names <- unlist(lapply(strsplit(colnames(temp2), '[_]'), '[[', 1))
    control_pos <- which(temp_brk == '11')
    temp_condition <- temp2[,-control_pos]
    samplenames <- sample_names[-control_pos]

    ##--- check which patients are not present in the survival data but not in the expression data and vice versa
    wh1 <- setdiff(samplenames, all_survival$submitter_id)
    wh2 <- setdiff(all_survival$submitter_id, samplenames)
    if(length(wh1) > 0){
        wh3 <- which(samplenames %in% wh1)
        samplenames <- samplenames[-wh3]
        temp_condition <- temp_condition[,-wh3]
    }
    if(length(wh2) > 0){
        wh3 <- which(all_survival$submitter_id %in% wh2)
        all_survival <- all_survival[-wh3, ]
    }
    
    ##--- call survival analysis -------
    sgpv <- survival_qval_num(only_gained, cpm_threshold[thr], temp_condition, all_survival, exons, rnd_run) ##-- real + label randomization
    slpv <- survival_qval_num(only_lost, cpm_threshold[thr], temp_condition, all_survival, exons, rnd_run)

    save_GL_EEI(only_gained,paste0('Gained_'), outdir3, cpm_threshold[thr], c_type, gained, net_file)
    save_GL_EEI(only_lost,paste0('Lost_'), outdir3, cpm_threshold[thr], c_type, lost, net_file)  
    
    ###--- expr random survival analysis -------------------------------------------------------------------------------------
    sgpv_rnd <- survival_qval_random(only_gained, cpm_threshold[thr], temp_condition, all_survival, exons, rnd_run) ##-- expression randomization
    slpv_rnd <- survival_qval_random(only_lost, cpm_threshold[thr], temp_condition, all_survival, exons, rnd_run)

    ##--- random tests -------------------------------------------------------------------------------------------------------
    fracg1_test <- c()
    fracg2_test <- c()
    for(j in 1:length(sgpv[[1]])){
        wh <- which(sgpv_rnd[[1]] == sgpv[[1]][j])
        rdis <- unlist(lapply(sgpv_rnd[[2]], '[[', wh))
        fracg1_test <- c(fracg1_test, length(which(rdis <= sgpv[[2]][j]))/rnd_run)
        fracg2_test <- c(fracg2_test, length(which(sgpv[[3]][[j]] <= sgpv[[2]][j]))/rnd_run)
    }

    fracl1_test <- c()
    fracl2_test <- c()
    for(j in 1:length(slpv[[1]])){
        wh <- which(slpv_rnd[[1]] == slpv[[1]][j])
        rdis <- unlist(lapply(slpv_rnd[[2]], '[[', wh))
        fracl1_test <- c(fracl1_test, length(which(rdis <= slpv[[2]][j]))/rnd_run)
        fracl2_test <- c(fracl2_test, length(which(slpv[[3]][[j]] <= slpv[[2]][j]))/rnd_run)
    }
    ##-------------------------------------------------------------------------------------------------------------------------

    save_selected_EEI(sgpv[[1]],paste0('Gained_Surv_'),outdir3,cpm_threshold[thr],c_type,sgpv[[2]],gained,net_file, fracg1_test, fracg2_test)
    save_selected_EEI(slpv[[1]],paste0('Lost_Surv_'),outdir3,cpm_threshold[thr],c_type,slpv[[2]],lost,net_file, fracl1_test, fracl2_test)

}

