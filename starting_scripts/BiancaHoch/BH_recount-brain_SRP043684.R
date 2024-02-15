############################################
#
#DO NOT USE
# Title: recount-brain metadata applied to the SRA study SRP043684.
# Abstract: Transcriptome Analysis.
# The development of induced pluripotent stem cell (iPSC) technology has provided such a new approach. 
# Here, we developed a human BD iPSC model and investigated the cellular phenotypes of hippocampal dentate gyrus neurons derived from the patient iPSCs. Using patch clamp recording, somatic Ca2+ imaging and RNA-seq techniques, 
# we found that the neurons derived from BD patients exhibited hyperactive action potential (AP) firing, 
# up-regulated expression of PKA/PKC/AP and mitochondria-related genes. 
# Moreover, lithium selectively reversed these alterations in the neurons of patients who responded to lithium treatment. 
# Therefore, hyper-excitability is one endophenotype of BD that is probably achieved through enhancement in the PKA/PKC 
# and Na+ channel signaling systems, and our BD iPSC model can be used to develop new therapies and drugs aimed at clinical 
# treatment of this disease. Overall design: total RNAseq from neurons generated from BD patient-specific iPS cells.
# Total number of samples: 18
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 05/30/2019
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


#Download the first project ("SRP043684"):
download_study(SRAs[2])

?download_study
# 2019-05-29 13:46:28 downloading file rse_gene.Rdata to SRP043684
# trying URL 'http://duffel.rail.bio/recount/v2/SRP035524/rse_gene.Rdata'
# Content type 'application/octet-stream' length 5275872 bytes (5.0 MB)
# downloaded 5.0 MB

#Load the data:
load(file.path((SRAs[2]), 'rse_gene.Rdata'))

#View the RangedSummarizedExperiment object:
rse_gene
str(rse_gene)

#Total number of samples used:
length(unique(rse_gene$sample))
#18

#Add sample metadata from recount-brain to "SRP043684":
RSEmeta <- add_metadata(rse_gene)
RSEmeta

#View sample phenotype data provided by the recount project:
colData(RSEmeta)

#View genes and base pair lengths:
rowData(RSEmeta)

#Once we have identified the study of interest, we can use the browse_study() function to browse the study at the SRA website:
browse_study(SRAs[2])

#Preparing for differential expression analysis:
#Scale counts by taking into account the total coverage per sample:
rse <- scale_counts(rse_gene)
rse

##### Details about counts #####

#Scale counts to 40 million mapped reads. Not all mapped reads are in exonic
#sequences, so the total is not necessarily 40 million.
colSums(assays(rse)$counts) / 1e6

# 32.05871   31.28918   33.00268   29.08565   27.32943   30.57827   15.18811   23.13391 
# SRR1501373 SRR1501374 SRR1501375 SRR1501376 SRR1501377 SRR1501378 SRR1501379 SRR1501380 
# 29.62647   32.96646   32.30376   35.47805   39.79535   31.02983   29.57096   32.18764 
# SRR1501381 SRR1501382 
# 31.28560   41.31374

#Compute read counts:
rse_read_counts <- read_counts(rse_gene)

#Difference between read counts and number of reads downloaded by Rail-RNA:
colSums(assays(rse_read_counts)$counts) / 1e6 -
  colData(rse_read_counts)$reads_downloaded / 1e6


#STOPPED HERE 05/29/2019
#Define expressed regions for study SRP043684:
regions <- expressed_regions('SRP043684', cutoff = 5L, 
                             maxClusterGap = 3000L)

?expressed_regions