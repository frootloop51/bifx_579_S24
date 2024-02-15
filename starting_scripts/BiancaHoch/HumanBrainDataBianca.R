############################################
# Title: Aquisition of normalized RNAseq count data using recount2 - Immunological
# Author: Miranda Darby  (darby@hood.edu)
# Date: updated on 5/12/19
# Recount2 paper: https://www.biorxiv.org/content/10.1101/247346v2
# Recount2 published workflow: https://f1000research.com/articles/6-1558/v1
# Recount2 website: https://jhubiostatistics.shinyapps.io/recount/
# Recount-brain paper for future use: https://www.biorxiv.org/content/10.1101/618025v1
###########################################

### Install BiocManager. Only needs to be done once.
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

### Install recount2. Also only needs to be done once
# BiocManager::install("recount")

library(recount)

unique(brainStudies$tissue_site_1)
# [1] "Cerebral cortex"    NA                   "Mixed"              "Whole brain"       
# [5] "Cerebellum"         "Brainstem"          "Corpus callosum"    "Nucleus accumbens" 
# [9] "Dura mater"         "Putamen"            "Substantia nigra"   "Hippocampus"       
# [13] "Lumbar spinal cord" "Frontal Cortex"     "Caudate" 

# to add brain specific metadata to Recount2 
brainStudies <- add_metadata()
## Subset out studies of interest

#### Find studies that include samples for Nucelus Accumbens, Putamen, Substantia Nigra, Hippocampus, Frontal Cortex, Caudate
NucAcc <- which(brainStudies[, "tissue_site_1"] == "Nucleus accumbens")
Put <- which(brainStudies$tissue_site_1=="Putamen")
SubNig <- which(brainStudies$tissue_site_1=="Substantia nigra")
Hip <- which(brainStudies$tissue_site_1=="Hippocampus")
FC <- which(brainStudies$tissue_site_1=="Frontal Cortex")
Caud <- which(brainStudies$tissue_site_1=="Caudate")

#### Find studies that include samples for dorsolateral prefrontal cortex, prefontal cortex
unique(brainStudies$tissue_site_3)
# [1] "Superior temporal gyrus"       
# [2] "Dorsolateral prefrontal cortex"
# [3] NA                              
# [4] "Motor cortex"                  
# [5] "Superior frontal gyrus"        
# [6] "Anterior cingulate gyrus"      
# [7] "Somatosensory cortex"          
# [8] "Prefrontal cortex"             
# [9] "Anterior prefrontal cortex"    
DPC <-which(brainStudies$tissue_site_3=="Dorsolateral prefrontal cortex")
PC <- which(brainStudies$tissue_site_3=="Prefrontal cortex")

AllRegions <- c(NucAcc, Put, SubNig, Hip, FC, Caud, DPC, PC)
(AllRegions <- unique(AllRegions))
(BestRegions <- c(NucAcc, SubNig, Hip, FC))

#Add step for pulling SRA numbers

#### Print info on workspace environment
sessionInfo()
# R version 3.5.1 (2018-07-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.6
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices
# [6] utils     datasets  methods   base     
# 
# other attached packages:
#   [1] recount_1.8.2              
# [2] SummarizedExperiment_1.12.0
# [3] DelayedArray_0.8.0         
# [4] BiocParallel_1.16.6        
# [5] matrixStats_0.54.0         
# [6] Biobase_2.42.0             
# [7] GenomicRanges_1.34.0       
# [8] GenomeInfoDb_1.18.2        
# [9] IRanges_2.16.0             
# [10] S4Vectors_0.20.1           
# [11] BiocGenerics_0.28.0        
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-6              bit64_0.9-7              
# [3] RColorBrewer_1.1-2        progress_1.2.0           
# [5] httr_1.4.0                GenomicFiles_1.18.0      
# [7] doRNG_1.7.1               tools_3.5.1              
# [9] backports_1.1.4           R6_2.4.0                 
# [11] rpart_4.1-15              Hmisc_4.2-0              
# [13] DBI_1.0.0                 lazyeval_0.2.2           
# [15] colorspace_1.4-1          nnet_7.3-12              
# [17] withr_2.1.2               tidyselect_0.2.5         
# [19] gridExtra_2.3             prettyunits_1.0.2        
# [21] bit_1.1-14                compiler_3.5.1           
# [23] htmlTable_1.13.1          derfinder_1.16.1         
# [25] xml2_1.2.0                pkgmaker_0.27            
# [27] rtracklayer_1.42.2        scales_1.0.0             
# [29] checkmate_1.9.3           readr_1.3.1              
# [31] stringr_1.4.0             digest_0.6.18            
# [33] Rsamtools_1.34.1          foreign_0.8-71           
# [35] rentrez_1.2.2             GEOquery_2.50.5          
# [37] XVector_0.22.0            base64enc_0.1-3          
# [39] pkgconfig_2.0.2           htmltools_0.3.6          
# [41] bibtex_0.4.2              limma_3.38.3             
# [43] BSgenome_1.50.0           htmlwidgets_1.3          
# [45] rlang_0.3.4               rstudioapi_0.10          
# [47] RSQLite_2.1.1             jsonlite_1.6             
# [49] acepack_1.4.1             dplyr_0.8.0.1            
# [51] VariantAnnotation_1.28.13 RCurl_1.95-4.12          
# [53] magrittr_1.5              GenomeInfoDbData_1.2.0   
# [55] Formula_1.2-3             Matrix_1.2-17            
# [57] Rcpp_1.0.1                munsell_0.5.0            
# [59] stringi_1.4.3             yaml_2.2.0               
# [61] zlibbioc_1.28.0           qvalue_2.14.1            
# [63] bumphunter_1.24.5         plyr_1.8.4               
# [65] grid_3.5.1                blob_1.1.1               
# [67] crayon_1.3.4              lattice_0.20-38          
# [69] Biostrings_2.50.2         splines_3.5.1            
# [71] GenomicFeatures_1.34.8    hms_0.4.2                
# [73] derfinderHelper_1.16.1    locfit_1.5-9.1           
# [75] knitr_1.22                pillar_1.3.1             
# [77] rngtools_1.3.1            reshape2_1.4.3           
# [79] codetools_0.2-16          biomaRt_2.38.0           
# [81] XML_3.98-1.19             glue_1.3.1               
# [83] downloader_0.4            latticeExtra_0.6-28      
# [85] data.table_1.12.2         foreach_1.4.4            
# [87] gtable_0.3.0              purrr_0.3.2              
# [89] tidyr_0.8.3               assertthat_0.2.1         
# [91] ggplot2_3.1.1             xfun_0.6                 
# [93] xtable_1.8-4              survival_2.44-1.1        
# [95] tibble_2.1.1              iterators_1.0.10         
# [97] registry_0.5-1            GenomicAlignments_1.18.1 
# [99] AnnotationDbi_1.44.0      memoise_1.1.0            
# [101] cluster_2.0.9 

