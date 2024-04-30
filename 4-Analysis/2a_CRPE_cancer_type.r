##############################################################################################
# Purpose: Preprocess existing studies data
## Survival analyze data from --> Genome-wide identification and analysis of prognostic features in human cancers, Cell Reports, 2022
##############################################################################################

rm(list=ls())
library(data.table)
library('biomaRt')

cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))

## read excel file of survival -----
gexpr <- readxl::read_excel('../data/survival_analysis.xlsx',4)
gexpr1 <- as.data.frame(gexpr[1:length(gexpr[[1]]),])
allgenes <- substr(gexpr1[[1]], 2,100)
gexpr1$Gene <- allgenes
wh <- which(colnames(gexpr1) %in% union(cancer_type,'Gene'))
gexpr2 <- gexpr1[,wh]
gexpr3 <- gexpr2[,-1]
signif_can <- sapply(gexpr3, function(x) length(which(x > 1.96 | x < -1.96)))
data.table::fwrite(gexpr2, '../data/survival_analysis_all.txt', sep='\t', quote=FALSE, row.names=FALSE)


##-- map genes to uniprot ids ---
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
attrs <- c('external_gene_name', 'uniprotswissprot')
alld <- getBM(attributes=attrs, filters='external_gene_name', values=gexpr1$Gene, mart=ensembl)
alld1 <- alld[alld$uniprotswissprot != '',]

uco <- plyr::count(alld1[[2]])
whc <- which( uco$freq > 1) ## which uniprotids have been used more than once
ucd <- alld1[alld1$uniprotswissprot %in% uco[[1]][whc], ]

alluni <- rep('', length(gexpr2[[1]]))
for(k in 1:length(allgenes)){
	wh <- which(alld1[[1]] == allgenes[k])
	if(length(wh) > 0){
		alluni[k] <- alld1[[2]][wh[1]]
	}
}

gexpr2$Uniprot <- alluni

##-- entries only with valid uniprots ---------------
gexpr4 <- gexpr2[gexpr2$Uniprot != '', ]
data.table::fwrite(gexpr4, '../data/survival_analysis_processed.txt', sep='\t', quote=FALSE, row.names=FALSE)


###------ protein expression ----
##--- read excel file of survival -----
gexpr <- readxl::read_excel('../data/survival_analysis.xlsx',7)
gexpr1 <- as.data.frame(gexpr[1:length(gexpr[[1]]),])
allgenes <- substr(gexpr1[[1]], 2,100)
gexpr1$Gene <- allgenes
wh <- which(colnames(gexpr1) %in% union(cancer_type,'Gene'))
gexpr2 <- gexpr1[,wh]
gexpr3 <- gexpr2[,-1]
signif_can <- sapply(gexpr3, function(x) length(which(x > 1.96 | x < -1.96)))
data.table::fwrite(gexpr2, '../data/survival_analysis_all_protein.txt', sep='\t', quote=FALSE, row.names=FALSE)


###------ miRNA expression ----
## read excel file of survival -----
gexpr <- readxl::read_excel('../data/survival_analysis.xlsx',5)
gexpr1 <- as.data.frame(gexpr[1:length(gexpr[[1]]),])
allgenes <- substr(gexpr1[[1]], 2,100)
gexpr1$Gene <- allgenes
wh <- which(colnames(gexpr1) %in% union(cancer_type,'Gene'))
gexpr2 <- gexpr1[,wh]
gexpr3 <- gexpr2[,-1]
signif_can <- sapply(gexpr3, function(x) length(which(x > 1.96 | x < -1.96)))
data.table::fwrite(gexpr2, '../data/survival_analysis_all_mirna.txt', sep='\t', quote=FALSE, row.names=FALSE)


###------ CNA ----
## read excel file of survival -----
gexpr <- readxl::read_excel('../data/survival_analysis.xlsx',2)
gexpr1 <- as.data.frame(gexpr[1:length(gexpr[[1]]),])
allgenes <- substr(gexpr1[[1]], 2,100)
gexpr1$Gene <- allgenes
wh <- which(colnames(gexpr1) %in% union(cancer_type,'Gene'))
gexpr2 <- gexpr1[,wh]
gexpr3 <- gexpr2[,-1]
signif_can <- sapply(gexpr3, function(x) length(which(x > 1.96 | x < -1.96)))
data.table::fwrite(gexpr2, '../data/survival_analysis_all_CNA.txt', sep='\t', quote=FALSE, row.names=FALSE)


###------ methylation ----
## read excel file of survival -----
gexpr <- readxl::read_excel('../data/survival_analysis.xlsx',3)
gexpr1 <- as.data.frame(gexpr[1:length(gexpr[[1]]),])
allgenes <- substr(gexpr1[[1]], 2,100)
gexpr1$Gene <- allgenes
wh <- which(colnames(gexpr1) %in% union(cancer_type,'Gene'))
gexpr2 <- gexpr1[,wh]
gexpr3 <- gexpr2[,-1]
signif_can <- sapply(gexpr3, function(x) length(which(x > 1.96 | x < -1.96)))
data.table::fwrite(gexpr2, '../data/survival_analysis_all_Methylation.txt', sep='\t', quote=FALSE, row.names=FALSE)


###------ mutations ----
## read excel file of survival -----
gexpr <- readxl::read_excel('../data/survival_analysis.xlsx',6)
gexpr1 <- as.data.frame(gexpr[1:length(gexpr[[1]]),])
allgenes <- substr(gexpr1[[1]], 2,100)
gexpr1$Gene <- allgenes
wh <- which(colnames(gexpr1) %in% union(cancer_type,'Gene'))
gexpr2 <- gexpr1[,wh]
gexpr3 <- gexpr2[,-1]
signif_can <- sapply(gexpr3, function(x) length(which(x > 1.96 | x < -1.96)))
data.table::fwrite(gexpr2, '../data/survival_analysis_all_Mutation.txt', sep='\t', quote=FALSE, row.names=FALSE)
