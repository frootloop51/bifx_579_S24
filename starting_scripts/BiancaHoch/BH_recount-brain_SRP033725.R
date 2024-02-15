############################################
#
# Title: recount-brain metadata applied to the SRA study SRP033725.
# Abstract: We used RNA-seq to survey the brain transcriptome in high-quality post-mortem dorsolateral prefrontal 
# cortex from 11 individuals diagnosed with bipolar disorder (BD) and from 11 age- and gender-matched controls. 
# Deep sequencing was performed, with over 350 million reads per specimen. 
# At a false-discovery rate of <5%, we detected 5 differentially-expressed (DE) genes and 12 DE transcripts, 
# most of which have not been previously implicated in BD. 
# Among these, PROM1/CD133 and ABCG2 play important roles in neuroplasticity. 
# We also show for the first time differential expression of long non-coding RNAs (lncRNAs) in BD. 
# DE transcripts include those of SRSF5 and RFX4, which along with lncRNAs play a role in mammalian circadian rhythms. 
# The DE genes were significantly enriched for several Gene Ontology (GO) categories. 
# Of these, genes involved with GTPase binding were also enriched for BD-associated SNPs from previous genome-wide 
# association studies, suggesting that differential expression of these genes is not simply a consequence of BD or 
# its treatment. Many of these findings were replicated by microarray in an independent sample of 60 cases and controls. 
# These results highlight common pathways for inherited and non-inherited influences on disease risk that may constitute 
# good targets for novel therapies. Overall design: Brain transcriptome in bipolar disorder
# Total number of samples: 22
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


#Download the first project ("SRP033725"):
download_study(SRAs[3])

?download_study
# 2019-05-29 13:46:28 downloading file rse_gene.Rdata to SRP033725:
# trying URL 'http://duffel.rail.bio/recount/v2/SRP035524/rse_gene.Rdata'
# Content type 'application/octet-stream' length 5275872 bytes (5.0 MB)
# downloaded 5.0 MB

#Load the data:
load(file.path((SRAs[3]), 'rse_gene.Rdata'))

#View the RangedSummarizedExperiment object:
rse_gene
str(rse_gene)

#Total number of samples used:
length(unique(rse_gene$sample))
#22

#Add sample metadata from recount-brain to "SRP033725":
RSEmeta <- add_metadata(rse_gene)
RSEmeta

#View sample phenotype data provided by the recount project:
colData(RSEmeta)

#View genes and base pair lengths:
rowData(RSEmeta)

#Once we have identified the study of interest, we can use the browse_study() function to browse the study at the SRA website:
browse_study(SRAs[3])

#Preparing for differential expression analysis:
#Scale counts by taking into account the total coverage per sample:
rse <- scale_counts(rse_gene)
rse

##### Details about counts #####

#Scale counts to 40 million mapped reads. Not all mapped reads are in exonic
#sequences, so the total is not necessarily 40 million.
colSums(assays(rse)$counts) / 1e6

# SRR1047813 SRR1047814 SRR1047815 SRR1047816 SRR1047817 SRR1047818 SRR1047819 SRR1047820 
# 33.53576   33.70302   33.69540   33.68758   33.05451   32.91335   33.16262   33.16866 
# SRR1047821 SRR1047822 SRR1047823 SRR1047824 SRR1047825 SRR1047826 SRR1047827 SRR1047828 
# 33.03853   32.75003   33.51883   33.50838   33.72702   33.70987   33.55104   31.61073 
# SRR1047829 SRR1047830 SRR1047831 SRR1047832 SRR1047833 SRR1047834 SRR1047835 SRR1047836 
# 31.42060   31.58233   31.56462   31.50677   29.11526   29.42787   29.15107   29.15551 
# SRR1047837 SRR1047838 SRR1047839 SRR1047840 SRR1047841 SRR1047842 SRR1047843 SRR1047844 
# 29.12833   34.66215   34.46921   34.77220   34.66775   34.66298   33.43065   33.15786 
# SRR1047845 SRR1047846 SRR1047847 SRR1047848 SRR1047849 SRR1047850 SRR1047851 SRR1047852 
# 33.43665   33.41124   32.89196   32.85440   32.72315   32.67064   32.75450   32.78701 
# SRR1047853 SRR1047854 SRR1047855 SRR1047856 SRR1047857 SRR1047858 SRR1047859 SRR1047860 
# 31.94776   32.14998   32.14926   31.94375   31.77387   34.68209   34.54113   34.69250 
# SRR1047861 SRR1047863 SRR1047864 SRR1047865 SRR1047866 SRR1047867 SRR1047868 SRR1047869 
# 34.67021   28.74245   26.99074   29.61936   26.30670   24.40949   28.24896   31.14360 
# SRR1047870 SRR1047871 SRR1047872 SRR1047873 SRR1047874 
# 27.78519   29.60858   30.24933   31.12525   32.03458 

#Compute read counts:
rse_read_counts <- read_counts(rse_gene)

#Difference between read counts and number of reads downloaded by Rail-RNA:
colSums(assays(rse_read_counts)$counts) / 1e6 -
  colData(rse_read_counts)$reads_downloaded / 1e6


#STOPPED HERE 05/29/2019
#Define expressed regions for study SRP033725:
regions <- expressed_regions('SRP033725', cutoff = 5L, 
                             maxClusterGap = 3000L)

?expressed_regions