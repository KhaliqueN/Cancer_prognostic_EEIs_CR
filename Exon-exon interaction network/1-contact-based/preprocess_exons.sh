#!/usr/bin/sh

# select only exons
awk -F '\t' '$3 == "exon"' ../data/Homo_sapiens.GRCh38.105.gtf > ../data/processed_exons.tmp
awk -F '\t' '$3 == "CDS"' ../data/Homo_sapiens.GRCh38.105.gtf > ../data/processed_cds.tmp

