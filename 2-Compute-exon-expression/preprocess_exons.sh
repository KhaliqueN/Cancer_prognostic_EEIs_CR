#!/usr/bin/sh

# select only exons
awk -F '\t' '$3 == "exon"' ../data/Homo_sapiens.GRCh38.105.gtf > ../data/Homo_sapiens.GRCh38.105_exons.gtf
awk -F '\t' '$3 == "gene"' ../data/Homo_sapiens.GRCh38.105.gtf > ../data/Homo_sapiens.GRCh38.105_genes.gtf

