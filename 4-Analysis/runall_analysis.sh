#!/bin/sh

############################################################################################
# Purpose: Run all analysis
###########################################################################################--

Rscript 1_cpm_variation.r

Rscript 2a_CRPE_cancer_type.r

Rscript 2b_CRPE_cancer_type.r

Rscript 2c_CRPE_cancer_type.r

Rscript 2d_CRPE_cancer_type.r

Rscript 3_CRPE_category.r

Rscript 4a_CRPE_biomarker.r

Rscript 4c_CRPE_biomarker.r

Rscript 5a_PPI_vs_EEI.r

Rscript 6_Supplementary_Table_S1.r

