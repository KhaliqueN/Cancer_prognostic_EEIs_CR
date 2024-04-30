#!/usr/bin/bash

## cancer types
# cancertype=( ESCA STAD )
cancertype=( STAD ESCA KICH BLCA BRCA KIRC HNSC KIRP LIHC LUAD LUSC UCEC THCA COAD PRAD )

for k in `seq 0 14`
do
	Rscript 3_normalizeCounts_all.r ${cancertype[k]}
	echo ${cancertype[k]}

done

