# Alternative polyadenylation site (APA) analysis for Ras mutants

This repo contains essential processed data and scripts to reproduce RNA-seq analysis and relevant figures in manuscript: **Alternative polyadenylation is a determinant of oncogenic Ras function** (Subramanian et al. 2021).

**celegans_data_analysis**: Scripts and data for C.elegans 3'UTR-seq data analysis:

Author: huayun.hou@sickkids.ca


  - **PA_sites_identification.html**: code for identifying polyadenylation (PA) sites from 3'UTR-seq data in celegans
  - **differential_PA.html**: code for performing differential PA site usage between conditions of interest and for reproducing all figures (for celegans data) in manuscript. 
  - **derfinder_ce11_utr_functions.r**: customized functions used in the analysis
  - **gdata**: folder contains celegans genomic annotations used 
  - **mapped_bigwig**: contains bigwig files of aligned reads, used for PA site identification 

  RData and RDS files are included for convenience of reproducing figures. 
 
 **human_data_analysis**: Scripts and data for human RNA-seq data analysis. 
 
 Author: marat.mufteev@mail.utoronto.ca 

