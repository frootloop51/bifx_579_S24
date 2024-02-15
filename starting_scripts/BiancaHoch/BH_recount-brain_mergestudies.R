############################################
#
# Title: Attempt to merge 2 SRA studies, combinining data from only controls and schizophrenia patients.
#
# Workflows and references used:
# BigWigs to count matrix- This protocol explains how to create a feature count matrix from coverage data stored in BigWig files.
# This feature count matrix can then be used for differential expression analyses. http://lcolladotor.github.io/protocols/bigwig_DEanalysis/
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 06/18/2019
#
###########################################

rm(list=ls()) #Clear R memory

###Set path variables that contain the different directories:###

###Install packages:
#install.packages("BiocManager")
#BiocManager::install(c("ggplot2", "pheatmap", "RColorBrewer", "DESeq2", "rtracklayer","recount", "Gviz", "derfinder", "derfinderData", "GenomicRanges"))


### Load libraries:
library('recount')
library('Gviz')
library('derfinder')
library('derfinderData')
library('GenomicRanges')
library('rtracklayer')
library('dplyr')
library('DESeq2')
library('pheatmap')
library('RColorBrewer')
library('ggplot2')


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
# [11] "Alzheimer????Ts disease"                      
# [12] "Parkinson????Ts disease"                      
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

#The following doesn't work for some reason.
# #Download the fourth project ("ERP001304"):
# ERP <- download_study(SRAs[4], download = TRUE)
# ERP <- load(file.path(SRAs[4], 'rse_gene.Rdata'))
# 
# #Download the first project ("SRP035524"):
# SRP <- download_study(SRAs[1], download = TRUE)
# SRP <- load(file.path(SRAs[1], 'rse_gene.Rdata'))

#Download the sample BigWig files for the projects:
#download_study('ERP001304', type = "samples", outdir = "~/Documents/splice_project/data/rawdata/bigwigs", download = TRUE)
#download_study('SRP035524', type = "samples", outdir = "~/Documents/splice_project/data/rawdata/bigwigs/SRP035524", download = TRUE)

#Construct full paths to the BigWig files and fix the names:
bwfilesERP <- rawFiles(datadir = "~/Documents/splice_project/data/rawdata/bigwigs", samplepatt = '*.bw' , fileterm = NULL)
names(bwfilesERP) <- gsub('.bw', '', names(bwfilesERP))

bwfilesSRP <- rawFiles(datadir = "~/Documents/splice_project/data/rawdata/bigwigs/SRP035524/bw", samplepatt = '*.bw' , fileterm = NULL)
names(bwfilesSRP) <- gsub('.bw', '', names(bwfilesSRP))

#Load the GRanges object that lists the new exons that we are interested in:
newGR <- ("~/Documents/splice_project/data/analysis/newGRanges_06112019.Rdata")
load(newGR)
newGR

#In the following code, I will lift the features from the GRanges list (hg19) to hg38.
chain <- import.chain("~/Documents/splice_project/data/hg19ToHg38.over.chain")
newGRhg38 <- unlist(liftOver(newGR, chain))

#Get a list of chromosomes from the GRanges list containing the new exons:
chromosomes <- levels(seqnames(newGRhg38))

#Import data and create a count matrix for all chromosomes present in the GRanges list (WITH CLOOJ):
bw <- BigWigFileList(bwfilesSRP)
exons <- newGRhg38
counts_exons <- matrix(NA, nrow = length(exons), ncol = length(bw))
colnames(counts_exons) <- names(bw)
for(i in seq_len(length(bw))) {
  for(chr in chromosomes) {
    eval(parse(text=paste0("coverage <- import(bw[[i]], as = 'RleList')$",chr)))
    counts_exons[, i] <- sum(Views(coverage, ranges(exons)))
    filename <- paste("~/Documents/splice_project/data/analysis/counts/SRP035524", chr, ".csv", sep = "")
    write.csv(counts_exons, file = filename)
  }
}





