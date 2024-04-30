##############################################################################################
# Purpose: RNA-seq normalization for the survival analysis
##############################################################################################
args=commandArgs(TRUE)

library(edgeR)
library(ggplot2)
library(data.table)
library(plyr)
library(Rsubread)
library(Rsamtools) 
# turning warnings into errors
options(warn=2)

## all cancer types
cancer_type <- args[1] #'ESCA'

out_dir <- paste0('../data/featureCounts_exons/', cancer_type)
out_dir_normalization <- paste0('../data/normalized_exons_sv/', cancer_type)
if(!dir.exists(out_dir_normalization)){
	dir.create(out_dir_normalization, recursive=TRUE)
}

temp_final1 <- fread(paste0('../data/',cancer_type,'_manifest_final.txt'))
temp_final2 <- fread(paste0('../data/',cancer_type,'_manifest_final_allCancer.txt'))
temp_final <- rbind(temp_final1, temp_final2)

##--------- Count summary -------------------------------------------------------------------------------
all_nids <- temp_final$nid
all_sample <- temp_final$sample_id

##--- save the exon list
input_file <- paste0(out_dir,'/',all_nids[1],'_',all_sample[1],'.tsv')
if(!file.exists(input_file)){
	input_file <- paste0(out_dir,'/',all_nids[1],'_0',all_sample[1],'.tsv')
}
cmd <- paste(c('cut -d"\t" -f1 ', input_file, ' > ',paste0(out_dir,'/allexons.tmp')), collapse="")
system(cmd)
##-- remove the first and second lines
system(paste0("tail -n+3 ",out_dir,"/allexons.tmp > ",out_dir,"/allexons.lst"))

##--- create empty file
file.create(paste0(out_dir,'/allExpression.exp'))

##--- saving the first and the seventh column for each file
not_present <- c()
counter <- 1
for(k in 1:length(all_nids)){

	input_file <- paste0(out_dir,'/',all_nids[k],'_',all_sample[k],'.tsv')
	if(!file.exists(input_file)){
		input_file <- paste0(out_dir,'/',all_nids[k],'_0',all_sample[k],'.tsv')
		if(!file.exists(input_file)){
			not_present <- c(not_present, paste0(all_nids[k],'_',all_sample[k]))
			next
		}
	}

	cmd <- paste0('cut -f7 ',input_file,' > ',out_dir,'/temp1.tmp')
	system(cmd)

	##-- remove the first and second lines
	system(paste0("tail -n+3 ",out_dir,"/temp1.tmp > ",out_dir,"/temp_",counter,".tmp"))

	counter <- counter+1
	cat('Sample',k, ' of ', length(all_nids), ' done\n')

}

##--- combine all expressions -------
STEP <- floor(counter/4)
STEP1 <- STEP*2
STEP2 <- STEP*3

cmd <- paste0('paste -d"\t" ',out_dir,'/temp_{1..',STEP,'}.tmp > ',out_dir,'/allExpression1.tmp')
system(cmd)

cmd <- paste0('paste -d"\t" ',out_dir,'/temp_{',STEP+1,'..',STEP1,'}.tmp > ',out_dir,'/allExpression2.tmp')
system(cmd)

cmd <- paste0('paste -d"\t" ',out_dir,'/temp_{',STEP1+1,'..',STEP2,'}.tmp > ',out_dir,'/allExpression3.tmp')
system(cmd)

cmd <- paste0('paste -d"\t" ',out_dir,'/temp_{',STEP2+1,'..',counter-1,'}.tmp > ',out_dir,'/allExpression4.tmp')
system(cmd)

cmd <- paste0('paste -d"\t" ',out_dir,'/allExpression{1..4}.tmp > ',out_dir,'/allExpression.exp')
system(cmd)


##--- remove temp files
system(paste0('rm ',out_dir,'/*.tmp'))

##--- read the expression file
temp <- fread(paste0(out_dir,'/allExpression.exp'), sep='\t', header=FALSE)

n_col2 <- length(temp)
## remove gene names columns
group <- rep("Condition",n_col2)

## convert to DGEList format
x.edger <- DGEList(counts=temp,group=group)

## TMM normalization
x.edger.tmm <- calcNormFactors(x.edger, method='TMM')

## compute cpm
x.norm.tmm.cpm <- cpm(x.edger.tmm, log = FALSE, normalized.lib.sizes=TRUE) 

## counts table
counts_tab <- as.data.frame(x.norm.tmm.cpm)
all_IDs <- paste0(all_nids,'_',all_sample)
wh <- which(all_IDs %in% not_present)
if(length(wh) != 0){
	all_IDs <- all_IDs[-wh]
}
colnames(counts_tab) <- all_IDs
allexons <- data.table::fread(paste0(out_dir,'/allexons.lst'), header=FALSE)
counts_tab$EXON <- allexons[[1]]
# counts_tab_filt <- counts_tab[counts_tab$EXON %in% nodes, ]
fwrite(counts_tab, paste0(out_dir_normalization,'/normalized_all.txt'), sep='\t', row.names=FALSE, quote=FALSE)


