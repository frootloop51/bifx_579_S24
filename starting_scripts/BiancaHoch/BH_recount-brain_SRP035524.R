############################################
#
# Title: recount-brain metadata applied to the SRA study SRP035524.
# Abstract: The DNA methylome and transcriptome of different brain regions in schizophrenia and bipolar disorder.
# Total number of samples: 35.
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 05/21/2019
#
###########################################

rm(list=ls()) #Clear R memory

###Set path variables that contain the different directories:###

###Install packages:
#install.packages("BiocManager")
#BiocManager::install(c("recount", "Gviz"))

### Load libraries:
library(recount)
library(Gviz)

#Retrieve recount-brain metadata:
meta <- add_metadata()

#Find unique studies:
unique(meta$disease)

#  [1] "Schizophrenia"                              
#  [2] NA                                           
#  [3] "brain tumor unspecified"                    
#  [4] "Hutchinson-Gilford progeria syndrome"       
#  [5] "Bipolar disorder"                           
#  [6] "Cortical ischemic stroke tissue"            
#  [7] "Autism spectrum disorder"                   
#  [8] "Epilepsy"                                   
#  [9] "Huntington's disease"                       
# [10] "Spinal muscular atrophy"                    
# [11] "Alzheimerâ???Ts disease"                      
# [12] "Parkinsonâ???Ts disease"                      
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

#Find biploar disorder and schizophrenia projects:
BPD <-which(meta$disease=="Bipolar disorder")
SCZ <- which(meta$disease=="Schizophrenia")

#Combine projects of interest:
BPDSCZ <- c(BPD,SCZ)
BPDSCZ

# [1] 125 126 127 128 129 130 131 142
#  [9] 143 144 145 146 157 158 200 201
# [17] 202 203 204 205 206 207 208 209
# [25] 210 211 604 605 606 607 608 609
# [33] 610 611 612 613 614 615 616 617
# [41] 618 619 620 621 622 623 624 625
# [49] 626 627 628 654 655 656 657 658
# [57] 659   1   2   3   4   5   6   7
# [65]   8   9 137 138 139 140 141 152
# [73] 153 154 155 159

#Retrieve metadata of projects of interest:
BPDSCZmeta <- (meta[(BPDSCZ),])
head(BPDSCZmeta)

#See which projects are present using SRA study ID:
SRAs <- unique(BPDSCZmeta$sra_study_s)
SRAs 
# [1] "SRP035524" "SRP043684" "SRP033725"
# [4] "ERP001304"


#Download the first project ("SRP035524"):
download_study(SRAs[1])

?download_study
# 2019-05-29 13:46:28 downloading file rse_gene.Rdata to SRP035524
# trying URL 'http://duffel.rail.bio/recount/v2/SRP035524/rse_gene.Rdata'
# Content type 'application/octet-stream' length 5275872 bytes (5.0 MB)
# downloaded 5.0 MB

#Load the data:
load(file.path((SRAs[1]), 'rse_gene.Rdata'))

#View the RangedSummarizedExperiment object:
rse_gene
str(rse_gene)

#Total number of samples used:
length(unique(rse_gene$sample))
#35

#Add sample metadata from recount-brain to "SRP035524":
RSEmeta <- add_metadata(rse_gene)
RSEmeta

#View sample phenotype data provided by the recount project:
colData(RSEmeta)

#View distribution of samples among conditions:
summary(as.factor(colData(RSEmeta)$disease))

#View genes and base pair lengths:
rowData(RSEmeta)

#Once we have identified the study of interest, we can use the browse_study() function to browse the study at the SRA website:
browse_study(SRAs[1])

#Preparing for differential expression analysis:
#Scale counts by taking into account the total coverage per sample:
rse <- scale_counts(rse_gene)
rse

##### Details about counts #####

#Scale counts to 40 million mapped reads. Not all mapped reads are in exonic
#sequences, so the total is not necessarily 40 million.
colSums(assays(rse)$counts) / 1e6

# SRR1130143 SRR1130243 SRR1130448 SRR1130508 
#   32.35457   35.27974   34.10197   36.77475 
# SRR1130605 SRR1130768 SRR1130770 SRR1130772 
#   29.93161   35.67637   35.24074   34.46396 
# SRR1130774 SRR1130827 SRR1130927 SRR1131110 
#   34.69101   34.66357   34.66702   32.18319 
# SRR1131112 SRR1131114 SRR1131116 SRR1131119 
#   33.70177   30.51355   33.08400   33.54209 
# SRR1131625 SRR1131659 SRR1131661 SRR1131663 
#   34.10075   33.82674   34.92719   37.06550 
# SRR1131665 SRR1131667 SRR1131745 SRR1131920 
#   31.52163   33.53092   33.21323   35.16038 
# SRR1131929 SRR1131931 SRR1131935 SRR1132102 
#   37.91320   30.65899   34.11697   32.68606 
# SRR1132101 SRR1132104 SRR1132106 SRR1132108 
#   32.56337   34.77098   33.09669   28.78073 
# SRR1132110 SRR1132112 SRR1132122 
#   33.92096   32.64020   33.12940 

#Compute read counts:
rse_read_counts <- read_counts(rse_gene)

#Difference between read counts and number of reads downloaded by Rail-RNA:
colSums(assays(rse_read_counts)$counts) / 1e6 -
    colData(rse_read_counts)$reads_downloaded / 1e6


#STOPPED HERE 05/29/2019
#Define expressed regions for study SRP035524:
regions <- expressed_regions('SRP035524', cutoff = 5L, 
    maxClusterGap = 3000L)

?expressed_regions

