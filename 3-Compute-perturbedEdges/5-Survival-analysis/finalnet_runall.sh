#!/bin/sh

############################################################################################--
# Purpose: Run the survival analysis
##############################################################--#############################--

outDir=../data/Final_survival
cancertype=( BLCA BRCA KIRC HNSC KIRP LIHC LUAD LUSC UCEC THCA COAD PRAD KICH STAD ESCA )
LEVEL=( level1 level2 level3 )

if [ ! -d $outDir ]; then
	mkdir $outDir
fi

for l in `seq 0 0`; do

	Rscript 1_weighted_network_final.r ${LEVEL[l]}

done


for l in `seq 0 0`; do

	outDir1=$outDir/${LEVEL[l]}

	if [ ! -d $outDir1 ]; then
		mkdir $outDir1
	fi

	for f in `seq 0 14`; do

	  screen -dm Rscript 2_analysis_individual_final.r ${cancertype[f]} $outDir1 ${LEVEL[l]}
	  # Rscript 2_analysis_individual_final.r ${cancertype[f]} $outDir1 ${LEVEL[l]}


	done

done







