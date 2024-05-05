#!/bin/sh

##############################################################################################
# Purpose: Run the survival analysis
##############################################################################################

outDir=../data/EPPIC_survival
cancertype=( BLCA BRCA KIRC HNSC KIRP LIHC LUAD LUSC UCEC THCA COAD PRAD KICH STAD ESCA )

Rscript 1_weighted_network_EPPIC.r

if [ ! -d $outDir ]; then
	mkdir $outDir2
fi


for f in `seq 0 14`; do

  screen -dm Rscript 2_select_edges_EPPIC.r ${cancertype[f]} $outDir

done




