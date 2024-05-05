#!/usr/bin/bash

cancertype=( STAD ESCA KICH BLCA BRCA KIRC HNSC KIRP LIHC LUAD LUSC UCEC THCA COAD PRAD )

Rscript 0_preprocess_data.r

for k in `seq 0 14`
do
	Rscript 1_downloadBAM_alignReads.r ${cancertype[k]}
	echo ${cancertype[k]}

done

