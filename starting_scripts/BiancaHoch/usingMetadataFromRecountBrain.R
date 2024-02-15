############################################
# Title: Aquisition of normalized RNAseq count data using recount2 - Immunological
# Author: Miranda Darby  (darby@hood.edu)
# Modified by: Alyssa Klein (amk19@hood.edu)
# Date: updated on 5/12/19
# Date: Updated by Alyssa Klein on 10/11/2019
# Recount2 paper: https://www.biorxiv.org/content/10.1101/247346v2
# Recount2 published workflow: https://f1000research.com/articles/6-1558/v1
# Recount2 website: https://jhubiostatistics.shinyapps.io/recount/
# Recount-brain paper for future use: https://www.biorxiv.org/content/10.1101/618025v1
###########################################

### Install BiocManager. Only needs to be done once.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

### Install recount2. Also only needs to be done once
# BiocManager::install("recount")

library(recount)

# to add brain specific metadata to Recount2 
# Recount-brain paper for future use: https://www.biorxiv.org/content/10.1101/618025v1
brainStudies <- add_metadata()

### Alyssa, you can use a similar searches to the ones below by 
### finding the right categories to access in this dataframe. 
summary(brainStudies)

### For example, you might try:
unique(brainStudies$disease)
# [1] "Schizophrenia"                              
# [2] NA                                           
# [3] "brain tumor unspecified"                    
# [4] "Hutchinson-Gilford progeria syndrome"       
# [5] "Bipolar disorder"                           
# [6] "Cortical ischemic stroke tissue"            
# [7] "Autism spectrum disorder"                   
# [8] "Epilepsy"                                   
# [9] "Huntington's disease"                       
# [10] "Spinal muscular atrophy"                    
# [11] "Alzheimer’s disease"                        
# [12] "Parkinson’s disease"                        
# [13] "Rett syndrome"                              
# [14] "Amyotrophic lateral sclerosis"              
# [15] "Parkinson's disease"                        
# [16] "ZNF804A Knockdown"                          
# [17] "Angelman syndrome"                          
# [18] "Dup15q syndrome"                            
# [19] "Embryonal tumors with multilayered rosettes"
# [20] "Primitive neuroectodermal tumor"            
# [21] "Primary Tumor"                              
# [22] "Recurrent Tumor"                            
# [23] "Solid Tissue Normal"                        
# [24] "ill_expected"                               
# [25] "violent_fast"                               
# [26] "fast_natural"                               
# [27] "ill_unexpected"                             
# [28] "ventilator"  

#### Find studies that include samples for specific diseases to see what data lies in recount
#### Picked the diseases or syndromes we would be most interested in at this time
# generates an ouput of indexes for what studies include samples where the sample is == to the specific disease
schizophrenia <- which(brainStudies[, "disease"] == "Schizophrenia")
# command below does same as command above but in a different way
# schizophrenia_2 <- which(brainStudies$disease == "Schizophrenia")

##### Schizophrenia #####
# generates an output for what study IDs include samples where the sample is == to the specific disease
sczStudies <-brainStudies[schizophrenia,]
# scz_tissue <- brainStudies[,c("tissue_site_1", "tissue_site_2", "tissue_site_3")]
# select for specific columns of interest
sczStudies[,c("Study_full", "sra_sample_s", "sra_study_s", "count_file_identifier")]
sczStudies[,c("Study_full", "tissue_site_1", "tissue_site_2", "tissue_site_3")]
# identify the tissue sites for the samples
scz_site_1 <- unique(sczStudies$tissue_site_1)
scz_site_2 <- unique(sczStudies$tissue_site_2)
scz_site_3 <- unique(sczStudies$tissue_site_3)

# Find unique SRA IDs for each disease
unique(sczStudies$Study_full)

# place the tissue sites into a data frame for easy viewing
table_of_tissue_scz <- list(scz_site_1, scz_site_2, scz_site_3)

brain_tumor <- which(brainStudies[, "disease"] == "brain tumor unspecified")
tumorStudies<-brainStudies[brain_tumor,]

##### Bipolar disorder #####

# generates an output for what study IDs include samples where the sample is == to the specific disease
bipolar <- which(brainStudies[, "disease"] == "Bipolar disorder")
bipolarStudies<-brainStudies[bipolar,]
# select for specific columns of interest
bipolarStudies[,c("Study_full", "sra_sample_s", "sra_study_s", "count_file_identifier")]
bipolarStudies[,c("Study_full", "tissue_site_1", "tissue_site_2", "tissue_site_3")]

bipol_site_1 <- unique(bipolarStudies$tissue_site_1)
bipol_site_2 <- unique(bipolarStudies$tissue_site_2)
bipol_site_3 <- unique(bipolarStudies$tissue_site_3)

# Find unique SRA IDs for each disease
unique(bipolarStudies$Study_full)

# place the tissue sites into a data frame for easy viewing
table_of_tissue_bipol <- list(bipol_site_1, bipol_site_2, bipol_site_3)

##### Autism Spectrum Disorder #####
autism_spec <- which(brainStudies[, "disease"] == "Autism spectrum disorder")
autismStudies <- brainStudies[autism_spec,]
autismStudies[,c("Study_full", "sra_sample_s", "sra_study_s", "count_file_identifier")]
autismStudies[,c("Study_full", "tissue_site_1", "tissue_site_2", "tissue_site_3")]

autism_site_1 <- unique(autismStudies$tissue_site_1)
autism_site_2 <- unique(autismStudies$tissue_site_2)
autism_site_3 <- unique(autismStudies$tissue_site_3)

# Find unique SRA IDs for each disease
unique(autismStudies$Study_full)

# place the tissue sites into a data frame for easy viewing
table_of_tissue_autism <- list(autism_site_1, autism_site_2, autism_site_3)

##### Huntington's Disease #####
huntington <- which(brainStudies[, "disease"] == "Huntington's disease")
huntingtonStudies<-brainStudies[huntington,]
huntingtonStudies[,c("Study_full", "sra_sample_s", "sra_study_s", "count_file_identifier")]
huntingtonStudies[,c("Study_full", "tissue_site_1", "tissue_site_2", "tissue_site_3")]

# Look for unique studies to be downloaded
huntingtonStudies[,"Study_full"]

huntington_site_1 <- unique(huntingtonStudies$tissue_site_1)
huntington_site_2 <- unique(huntingtonStudies$tissue_site_2)
huntington_site_3 <- unique(huntingtonStudies$tissue_site_3)

# Find unique SRA IDs for each disease
unique(huntingtonStudies$Study_full)

# place the tissue sites into a data frame for easy viewing
table_of_tissue_hunt <- list(huntington_site_1, huntington_site_2, huntington_site_3)

##### Alzheimer's disease #####
alz_disease <- which(brainStudies[, "disease"] == "Alzheimer’s disease")
alz_diseaseStudies<-brainStudies[alz_disease,]
alz_diseaseStudies[,c("Study_full", "sra_sample_s", "sra_study_s", "count_file_identifier")]
alz_diseaseStudies[,c("Study_full", "tissue_site_1", "tissue_site_2", "tissue_site_3")]

alz_site_1 <- unique(alz_diseaseStudies$tissue_site_1)
alz_site_2 <- unique(alz_diseaseStudies$tissue_site_2)
alz_site_3 <- unique(alz_diseaseStudies$tissue_site_3)

# Find unique SRA IDs for each disease
unique(alz_diseaseStudies$Study_full)

# place the tissue sites into a data frame for easy viewing
table_of_tissue_alz <- list(alz_site_1, alz_site_2, alz_site_3)

##### Rett syndrome ######
rett_syn <- which(brainStudies[, "disease"] == "Rett syndrome")
rettStudies<-brainStudies[rett_syn,]
rettStudies[,c("Study_full", "sra_sample_s", "sra_study_s", "count_file_identifier")]
rettStudies[,c("Study_full", "tissue_site_1", "tissue_site_2", "tissue_site_3")]

rettStudies_site_1 <- unique(rettStudies$tissue_site_1)
rettStudies_site_2 <- unique(rettStudies$tissue_site_2)
rettStudies_site_3 <- unique(rettStudies$tissue_site_3)

# Find unique SRA IDs for each disease
unique(rettStudies$Study_full)

# place the tissue sites into a data frame for easy viewing
table_of_tissue_rett <- list(rettStudies_site_1, rettStudies_site_2, rettStudies_site_3)

##### Parkinson's disease #####
park_disease <- which(brainStudies[, "disease"] == "Parkinson's disease")
parkStudies<-brainStudies[park_disease,]
parkStudies[,c("Study_full", "sra_sample_s", "sra_study_s", "count_file_identifier")]
parkStudies[,c("Study_full", "tissue_site_1", "tissue_site_2", "tissue_site_3")]

parkStudies_site_1 <- unique(parkStudies$tissue_site_1)
parkStudies_site_2 <- unique(parkStudies$tissue_site_2)
parkStudies_site_3 <- unique(parkStudies$tissue_site_3)

# Find unique SRA IDs for each disease
unique(parkStudies$Study_full)

# place the tissue sites into a data frame for easy viewing
table_of_tissue_park <- list(parkStudies_site_1, parkStudies_site_2, parkStudies_site_3)

znf804a_knockdown <- which(brainStudies[, "disease"] == "ZNF804A Knockdown")
znf804a_Studies<-brainStudies[znf804a_knockdown,]

dup15q_syn <- which(brainStudies[, "disease"] == "Dup15q syndrome")
dup15q_synStudies<-brainStudies[dup15q_syn,]

prim_neur_tumor <- which(brainStudies[, "disease"] == "Primitve neuroectodermal tumor")
prim_neur_tumorStudies<-brainStudies[prim_neur_tumor,]

##### Hutchinson-Gilford progeria syndrome #####

hutchinson_gilf <- which(brainStudies[, "disease"] == "Hutchinson-Gilford progeria syndrome")
hutchinson_gilfStudies<-brainStudies[hutchinson_gilf,]
hutchinson_gilfStudies[,c("Study_full", "sra_sample_s", "sra_study_s", "count_file_identifier")]
hutchinson_gilfStudies[,c("Study_full", "tissue_site_1", "tissue_site_2", "tissue_site_3")]

hutchinson_gilfStudies_site_1 <- unique(rettStudies$tissue_site_1)
hutchinson_gilfStudies_site_2 <- unique(rettStudies$tissue_site_2)
hutchinson_gilfStudies_site_3 <- unique(rettStudies$tissue_site_3)

# Find unique SRA IDs for each disease
unique(hutchinson_gilfStudies$Study_full)

# place the tissue sites into a data frame for easy viewing
table_of_tissue_hutch <- list(hutchinson_gilfStudies_site_1, hutchinson_gilfStudies_site_2, hutchinson_gilfStudies_site_3)

cortical_stroke <- which(brainStudies[, "disease"] == "Cortical ischemic stroke tissue")
cortical_strokeStudies<-brainStudies[cortical_stroke,]

epilepsy <- which(brainStudies[, "disease"] == "Epilepsy")
epilepsy_synStudies<-brainStudies[epilepsy,]

spinal_musc_atr <- which(brainStudies[, "disease"] == "Spinal muscular atrophy")
spinal_musc_atrStudies<-brainStudies[spinal_musc_atr,]

amyo_lat_scler <- which(brainStudies[, "disease"] == "Amyotrophic lateral sclerosis")
amyo_lat_sclerStudies<-brainStudies[amyo_lat_scler,]

brain_tumor <- which(brainStudies[, "disease"] == "Brain tumor")
brain_tumorStudies<-brainStudies[brain_tumor,]

##### Angelmann syndrome #####
angelman_syn <- which(brainStudies[, "disease"] == "Angelman syndrome")
angelman_synStudies<-brainStudies[angelman_syn,]
angelman_synStudies[,c("Study_full", "sra_sample_s", "sra_study_s", "count_file_identifier")]
angelman_synStudies[,c("Study_full", "tissue_site_1", "tissue_site_2", "tissue_site_3")]

angelman_synStudies_site_1 <- unique(angelman_synStudies$tissue_site_1)
angelman_synStudies_site_2 <- unique(angelman_synStudies$tissue_site_2)
angelman_synStudies_site_3 <- unique(angelman_synStudies$tissue_site_3)

# Find unique SRA IDs for each disease
unique(angelman_synStudies$Study_full)

# place the tissue sites into a data frame for easy viewing
table_of_tissue_angelman <- list(angelman_synStudies_site_1, angelman_synStudies_site_2, angelman_synStudies_site_3)

embryonal_rosettes <- which(brainStudies[, "disease"] == "Embryonal tumors with multilayered rosettes")
embryonal_rosettesStudies<-brainStudies[embryonal_rosettes,]

#### Find studies that include samples for dorsolateral prefrontal cortex, prefontal cortex
# looks like the tissue sites become more specific from tissue site 1 being the broadest to tissue site 3 being more specific

unique(brainStudies$tissue_site_1)
# [1] "Cerebral cortex"                 
# [2] NA                                
# [3] "Mixed"                           
# [4] "Whole brain"                     
# [5] "Cerebellum"                      
# [6] "Brainstem"                       
# [7] "Corpus callosum"                 
# [8] "Nucleus accumbens"               
# [9] "Dura mater"                      
# [10] "Putamen"                         
# [11] "Substantia nigra"                
# [12] "Hippocampus"                     
# [13] "Lumbar spinal cord"              
# [14] "Frontal Cortex"                  
# [15] "Caudate"                         
# [16] "LGG"                             
# [17] "GBM"                             
# [18] "Cortex"                          
# [19] "Anterior cingulate cortex (BA24)"
# [20] "Hypothalamus"                    
# [21] "Spinal cord (cervical c-1)"      
# [22] "Cerebellar Hemisphere"           
# [23] "Amygdala"

unique(brainStudies$tissue_site_2)
# [1] "Temporal lobe"   "Frontal lobe"    NA               
# [4] "Cingulate gyrus" "Medulla"         "Insular cortex" 
# [7] "Parietal lobe"   "Pons"            "Occipital lobe" 

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

# AllRegions <- c(NucAcc, Put, SubNig, Hip, FC, Caud, DPC, PC)
# (AllRegions <- unique(AllRegions))
# (BestRegions <- c(NucAcc, SubNig, Hip, FC))

#Add step for pulling SRA numbers
# change the SRP to match those output for studies for each disease

######################################################################
# Downloading the data for each study in four steps:
##### 1) Download the data if it is not there
##### 2) Check that the file was downloaded
##### 3) Load the data
##### 4) Use colData() function to look at pheno data

# Download the data if it is not there
# Input the study_full names here that would need to be downloaded

### IDs for the studies with specific disease types ###
# Schizophrenia: SRP035524, ERP001304
# Bipolar disorder: SRP035524, SRP043684, SRP033725
# Autism spectrum disorder: SRP048683, SRP007483
# Huntington's disease: SRP051844, ERP004592
# Alzheimer's disease: SRP056604, SRP056863, SRP026562
# Rett syndrome: SRP066394
# Parkinson's disease: SRP015668
# Hutchinson-Gilford progeria syndrome:SRP033078
# Angelmann syndrome: SRP044749

if(!file.exists(file.path("SRP044749", "rse_gene.Rdata"))) {
  download_study("SRP044749", type = "rse-gene")
}

# Check that the file was downloaded
file.exists(file.path("SRP044749", "rse_gene.Rdata"))

# Load the data
loaded <- load(file.path("SRP044749", "rse_gene.Rdata"))

# Try and add something to generate table for disease, region of brain, number of samples, and the studies (SRA numbers)
sampleData <- colData(rse_gene)

####################################################################

# downloading the RSE_v2 will download the actual data (R summarized experiments)
# the following command downloads the url for the data in the project specified
# srp035524 <- download_study("SRP035524", type= "rse-gene", download= TRUE, version= 2)

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

