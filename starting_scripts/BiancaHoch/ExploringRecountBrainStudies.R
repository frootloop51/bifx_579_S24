############################################
#
# Title: Exploring recount2 brain studies (available patient ages and sample size).
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 06/19/2019
#
###########################################

rm(list=ls()) #Clear R memory

###Set path variables that contain the different directories:###

###Install packages:
#install.packages("BiocManager")
BiocManager::install(c("ggplot2", "pheatmap", "RColorBrewer", "DESeq2", "rtracklayer","recount", "Gviz", "derfinder", "derfinderData", "GenomicRanges"))


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

#SAMPLE SIZES-----------

#Retrieve recount-brain metadata:
meta <- add_metadata()

#Find out which studies contain the most samples:
sort(table(meta$sra_study_s), decreasing = TRUE)[1:10]

# SRP025982 SRP012682 SRP027383 SRP044668 SRP058181 SRP045638 SRP051844 
# 2898      1409       274        94        73        72        69 
# SRP033725 SRP056477 SRP043364 
# 62        53        52 

#Find unique studies:

#SUBJECT AGE INFORMATION-----------

#Available ages:
sort(unique(meta$age))

#Age units (years, months, etc.)
unique(meta$age_units)
# [1] "Years"                 "Months"                NA                     
# [4] "Post Conception Weeks" "Weeks"                 "Days"                 

#See which studies have age_units for "Post Conception Weeks"
Fetal <-which(meta$age_units == "Post Conception Weeks")
Fetalmeta <- (meta[(Fetal),])
head(Fetalmeta)
#See which projects are present using SRA study ID:
SRAs <- unique(Fetalmeta$sra_study_s)
SRAs 
# [1] "SRP045638"
#Description for SRP045638: RNAseq data of 36 samples across human brain development by age group from LIBD.




