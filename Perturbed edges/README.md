**This folder can be used to identify perturbed exon-exon interactions (EEIs) given the exon expression and an EEI network**
- Put the path to the token in line 14 of the "1_downloadBAM_alignReads.r" script
- Run the script "run_all.sh" as sh ./run_all.sh to download the .bam files and otain counts of exon expressions
- After the above script is finished running, run sh ./run_all_nm.sh to normalize the exon counts
- The results will be generated in the "data" folder
- Folder "../data/featureCounts_exons/" will contain exon counts
- Folder "../data/normalized_exons_sv/" will contain normalized exon counts

