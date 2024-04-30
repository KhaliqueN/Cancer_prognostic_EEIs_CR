#!/usr/bin/bash

cancertype=( STAD ESCA KICH BLCA BRCA KIRC HNSC KIRP LIHC LUAD LUSC UCEC THCA COAD PRAD )

##--- for these three cancer types, at least one .bam file is not getting downloaded
## due to "Division by Zero" from the server --
for k in `seq 0 14`
do
	Rscript 1_downloadBAM_alignReads.r ${cancertype[k]}
	echo ${cancertype[k]}

done



##--- data copy to local drive

## rsync -rt newazkha@max-display.desy.de:/beegfs/desy/user/newazkha/project2-exon-expression/data/featureCounts_exons .