###########################################
#
# Title: Picking Novel Exons for First Pass Analysis with New Supplementary Table 20, and
#        keeping all novel exons.
#
# In this script I have added a step that inserts the column containing 
# strand specific information (from the new supplementary table - NovelExonsWithStrandAndProteinTranslation.csv) into
# the previous supplementary table that did not have strand information (SupplementaryTable_19.csv).
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 07/29/2019
#
###########################################

rm(list=ls()) #Clear R memory

#Download packages:
BiocManager::install(c("GenomicRanges","rtracklayer"))

#Load Libraries:
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Establish Paths:

nE <- read.csv("~/Documents/splice_project/data/rawdata/SupplementaryTable 19.csv", header = TRUE)
nE <- as.data.frame(nE)

nE2 <- read.csv("~/Documents/splice_project/data/rawdata/NovelExonsWithStrandAndProteinTranslation.csv", header = TRUE)
nE2 <- as.data.frame(nE2)

#The columns in both dataframes containing the novel exon position do not match:
identical(nE[,3],nE2[,1])

#Create a new dataframe that only contains rows in which first table has matching keys with the second table:
newDF <- merge(x=nE,y=nE2,by="putativeExon")

#The new dataframe now contains information about the strand.

#Variables:
newEx <- newDF$putativeExon
charEx <- as.character(newEx)

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

# # #Create regions for the start and ends of the REL, respectively:
# startjunction <- paste(theStart,theStart,sep = "-", collapse = NULL)
# endjunction <- paste(theEnd,theEnd,sep = "-", collapse = NULL)

#Make 2 GRanges lists:
regions1 <- GRanges(seqnames= chr, IRanges(start = theStart, end = theStart), strand = newDF$strand)
regions2 <- GRanges(seqnames= chr, IRanges(start = theEnd, end = theEnd), strand = newDF$strand)

#Add metadata:
values(regions1) <- DataFrame(Gene_ID <- newDF$GENE_ID.x)
values(regions2) <- DataFrame(Gene_ID <- newDF$GENE_ID.x)

#Save the GRanges objects:
save(regions1, file="~/Documents/splice_project/data/analysis/GRanges_regions1.Rdata")
save(regions2, file="~/Documents/splice_project/data/analysis/GRanges_regions2.Rdata")