############################################
#
# Title: Exploring available data in recount2 and recount-brain, downloading and combining studies for gene expression analysis.
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 06/28/2019
#
###########################################

rm(list=ls()) #Clear R memory
gc()

###Set path variables that contain the different directories:###

###Install packages:
#install.packages(c("BiocManager","xfun")
#BiocManager::install(c("recount","dplyr","downloader","janitor", "rtracklayer"))


### Load libraries:
library('recount')
library('dplyr')
library('downloader')
library('janitor')
library('xfun')
#install_github('LieberInstitute/jaffelab')
library('jaffelab')
library('rtracklayer')

#Explore available samples.
#Search among GTEx samples, and samples in recount-brain to determine which samples and which
#studies I will use for analysis

#Get recount_brain metadata:
recount_brain_v1 <- recount::add_metadata(source = 'recount_brain_v1')

#Loading objects:
brain <- recount_brain_v1

#Get GTEx metadata:
gtex <- as.data.frame(recount::all_metadata('GTEx'))

#Find GTEx samples that are derived from brains
gtex %>% filter(smts == "Brain") %>% group_by(project, smtsd) %>% summarise(n = n())

  # # A tibble: 13 x 3
  # # Groups:   project [1]
  # project   smtsd                                         n
  # <chr>     <chr>                                     <int>
  #   1 SRP012682 Brain - Amygdala                             81
  # 2 SRP012682 Brain - Anterior cingulate cortex (BA24)     99
  # 3 SRP012682 Brain - Caudate (basal ganglia)             134
  # 4 SRP012682 Brain - Cerebellar Hemisphere               118
  # 5 SRP012682 Brain - Cerebellum                          145
  # 6 SRP012682 Brain - Cortex                              132
  # 7 SRP012682 Brain - Frontal Cortex (BA9)                120
  # 8 SRP012682 Brain - Hippocampus                         103
  # 9 SRP012682 Brain - Hypothalamus                        104
  # 10 SRP012682 Brain - Nucleus accumbens (basal ganglia)   123
  # 11 SRP012682 Brain - Putamen (basal ganglia)             103
  # 12 SRP012682 Brain - Spinal cord (cervical c-1)           76
  # 13 SRP012682 Brain - Substantia nigra                     71

#Find recount-brain samples that are taken from tissue types listed above, and are not cancerous:
unique(brain$tissue_site_1)

# [1] "Cerebral cortex"    NA                   "Mixed"             
# [4] "Whole brain"        "Cerebellum"         "Brainstem"         
# [7] "Corpus callosum"    "Nucleus accumbens"  "Dura mater"        
# [10] "Putamen"            "Substantia nigra"   "Hippocampus"       
# [13] "Lumbar spinal cord" "Frontal Cortex"     "Caudate"  

(cortex <- brain %>% 
    filter(is.na(tumor_type)) %>% 
    filter(tissue_site_1 == "Cerebral cortex") %>%
    group_by(sra_study_s) %>%
    summarize(n=n())
           )
sum(cortex$n)
#cortex has 634 samples

(hippo <- brain %>% 
    filter(is.na(tumor_type)) %>% 
    filter(tissue_site_1 == "Hippocampus") %>%
    group_by(sra_study_s) %>%
    summarize(n=n())
)
sum(hippo$n)
#hippocampus has 25 samples


(cerebellum <- brain %>% 
    filter(is.na(tumor_type)) %>% 
    filter(tissue_site_1 == "Cerebellum") %>%
    group_by(sra_study_s) %>%
    summarize(n=n())
)
sum(cerebellum$n)
#cerebellum has 29 samples

(nucacc <- brain %>% 
    filter(is.na(tumor_type)) %>% 
    filter(tissue_site_1 == "Nucleus accumbens") %>%
    group_by(sra_study_s) %>%
    summarize(n=n())
)
sum(nucacc$n)
#Nucleus accumbens has 1 sample

(putamen <- brain %>% 
    filter(is.na(tumor_type)) %>% 
    filter(tissue_site_1 == "Putamen") %>%
    group_by(sra_study_s) %>%
    summarize(n=n())
)
sum(putamen$n)
#Putamen has 6 samples

(subn <- brain %>% 
    filter(is.na(tumor_type)) %>% 
    filter(tissue_site_1 == "Substantia nigra") %>%
    group_by(sra_study_s) %>%
    summarize(n=n())
)
sum(subn$n)
#Substantia nigra has 1 sample

#Now I will download the expression data (bigwig files) from both GTEx and Recount2 samples:
#First, download expression data from recount2 studies:
#Combine studies that contain our tissues of interest shown above.
tissue_table <- bind_rows(cortex,hippo,cerebellum,nucacc,putamen,subn)
#List out the unique studies
recount_list <- unique(tissue_table$sra_study_s)

#Download recount-brain sample bigwig files using the following code:
for(i in 1:length(recount_list)){
  if(!file.exists(file.path(recount_list[i]))) {
    download_study(recount_list[i], type = "samples", outdir = "~/Documents/splice_project/data/rawdata/bigwigs/recount2_bw", download = TRUE)
  }
}

#Now I will download expression data from the GTEx study:
if(!file.exists(file.path('SRP012682'))) {
       download_study('SRP012682', type = "samples", outdir = "~/Documents/splice_project/data/rawdata/bigwigs/gtex_bw", download = TRUE)
     }








