###########################################
#
# Title: Lifting Putative Novel Exon regions from hg19 to hg38 for first pass analysis.
#
# In this script I am lifting the information from the new supplementary table - NovelExonsWithStrandAndProteinTranslation.csv.
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 10/23/2019
#
###########################################

rm(list=ls()) #Clear R memory

#Download packages:
BiocManager::install(c("GenomicRanges","rtracklayer","Repitools"))

#Load Libraries:
library(GenomicRanges)
library(rtracklayer)
library(Repitools)

#Load dataframe containing putative novel exons:
PEx <- read.csv("~/Documents/splice_project/data/rawdata/NovelExonsWithStrandAndProteinTranslation.csv", header = TRUE)

#Variables:
novelEx <- PEx$putativeExon
charEx <- as.character(novelEx)

#Get the chromosome information for each REL exon position:
chr <- sapply(strsplit(charEx,":"),"[[",1)

#Get the starting point of the REL exon:
theStart <- sapply(strsplit(charEx,":"),"[[",2)
theStart <- sapply(strsplit(theStart,"-"),"[[",1)
theStart <- as.integer(theStart)

#Get the end of the REL exon:
theEnd <- sapply(strsplit(charEx,":"),"[[",2)
theEnd <- sapply(strsplit(theEnd,"-"),"[[",2)
theEnd <- as.integer(theEnd)

#Make the hg19 GRanges lists:
PutativeExhg19 <- GRanges(seqnames= chr, IRanges(start = theStart, end = theEnd), strand = PEx$strand)

#In the following code, I will lift the features from the GRanges list (hg19) to hg38.
chain <- import.chain("~/Documents/splice_project/data/hg19ToHg38.over.chain")
PutativeExhg38 <- unlist(liftOver(PutativeExhg19, chain))

#Now the GRanges list will be converted into a dataframe, 
#so that the information can be saved as a .csv file to load into the makeSnapScripts.R script.
DFhg38 <- annoGR2DF(PutativeExhg38)
DFhg38_2 <- paste(DFhg38$chr,DFhg38$start, sep = ":")
DFhg38_3 <- paste(DFhg38_2,DFhg38$end, sep = "-")
DFhg38_3 <- as.data.frame(DFhg38_3)

#The region and strand information will be loaded in the following dataframe:
PutativeNovelExons_hg38 <- bind_cols(DFhg38_3,DFhg38[5])
write.csv(PutativeNovelExons_hg38, file="~/Documents/splice_project/data/rawdata/PutativeNovelExons_hg38.csv")

#There are atleast 6 more regions in the hg38 GRanges list than the hg19 GRanges list.