###########################################
#
# Title: Picking Novel Exons for First Pass Analysis with New Supplementary Table 20
#
# In this script I have added a stept that inserts the column containing 
# strand specific information (from the new supplementary table - NovelExonsWithStrandAndProteinTranslation.csv) into
# the previous supplementary table that did not have strand information (SupplementaryTable_19.csv).
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 05/30/2019
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

nE <- read.csv("~/Documents/splice_project/data/rawdata/SupplementaryTable_19.csv", header = TRUE)
nE <- as.data.frame(nE)

nE2 <- read.csv("~/Documents/splice_project/data/rawdata/NovelExonsWithStrandAndProteinTranslation.csv", header = TRUE)
nE2 <- as.data.frame(nE2)

#The columns in both dataframes containing the novel exon position do not match:
identical(nE[,3],nE2[,1])

#Create a new dataframe that only contains rows in which first table has matching keys with the second table:
newDF <- merge(x=nE,y=nE2,by="REL.Exon_Position")

#The new dataframe now contains information about the strand.

#Variables:
newEx <- newDF$REL.Exon_Position
ex1 <- newDF$Exon.1_Position
ex2 <- newDF$Exon.2_Position

#make the user defined function to make the Genomic Ranges:
makeGR <- function(regions)
{
  a=as.character(regions)
  eNames=sapply(strsplit(a,":"),"[[",1)
  theCoords=sapply(strsplit(a,":"),"[[",2)
  #theStrand=substr(theCoords,nchar(theCoords),nchar(theCoords))
  #theCoords=substr(theCoords,1,nchar(theCoords)-2)
  theStart=sapply(strsplit(theCoords,"-"),"[[",1)
  theEnd=sapply(strsplit(theCoords,"-"),"[[",2)
  theStart=as.integer(theStart)
  theEnd=as.integer(theEnd)
  #regions=GRanges(seqnames=eNames,ranges=IRanges(start=theStart,end=theEnd),strand=theStrand)
  regions=GRanges(seqnames=eNames,ranges=IRanges(start=theStart,end=theEnd), strand = newDF$strand)
}

#Make genomic ranges for new exons using the makeGR function:
newExons <- makeGR(newEx)

#Find width of new exons:
newWidth <- width(newExons)

#Divide coverage by width:
normCovNewEx <- newDF$REL.Exon_Coverage/newWidth

#Make genomic ranges for exon 1s using the makeGR function:
ex1gr <- makeGR(ex1)

#Find width of exon 1s:
ex1width <- width(ex1gr)

#Divide coverage by width:
normEx1 <- newDF$Exon.1_Coverage/ex1width

#Make genomic ranges for exon 2s using the makeGR function:
ex2gr <- makeGR(ex2)

#Find width of exon 2s:
ex2width <- width(ex2gr)

#Divide coverage by width:
normEx2 <- newDF$Exon.2_Coverage/ex2width

#Determine which novel exons are at least 10% the size of normalized exon 1
goodEx1 <- which(normCovNewEx >= 0.1*normEx1)

#Determine which novel exons are at least 10% the size of normalized exon 2
goodEx2 <- which(normCovNewEx >= 0.1*normEx2)

#Merge the exons that pass the criteria:
goodEx <- unique(goodEx1,goodEx2)

#Create a subset of the newExons GRanges object that only contains the novel exons that meet our desired criteria:
newGR <- newExons[goodEx]

#Save the GRanges object:
save(newGR, file="~/Documents/splice_project/data/analysis/newGRanges_06112019.Rdata")
