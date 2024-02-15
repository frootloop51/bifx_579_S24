###########################################
#
# Title: Picking Novel Exons for First Pass Analysis
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 05/30/2019
#
###########################################

rm(list=ls()) #Clear R memory

#Install packages
install.packages("readxl")

#Load Libraries:
library(GenomicRanges)
library(readxl)

#Establish Paths:

nE <- read.csv("~/Documents/splice_project/data/rawdata/SupplementaryTable_19.csv", header = TRUE)
nE <- as.data.frame(nE)

#Set the title of the first column to avoid errors on reading in the data:
colnames(nE)[1] <- "putativeExon"

#Variables:
newEx <- nE$REL.Exon_Position
ex1 <- nE$Exon.1_Position
ex2 <- nE$Exon.2_Position

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
  regions=GRanges(seqnames=eNames,ranges=IRanges(start=theStart,end=theEnd),strand="*")
}

#Make genomic ranges for new exons using the makeGR function:
newExons <- makeGR(newEx)

#Find width of new exons:
newWidth <- width(newExons)

#Divide coverage by width:
normCovNewEx <- nE$REL.Exon_Coverage/newWidth

#Make genomic ranges for exon 1s using the makeGR function:
ex1gr <- makeGR(ex1)

#Find width of exon 1s:
ex1width <- width(ex1gr)

#Divide coverage by width:
normEx1 <- nE$Exon.1_Coverage/ex1width

#Make genomic ranges for exon 2s using the makeGR function:
ex2gr <- makeGR(ex2)

#Find width of exon 2s:
ex2width <- width(ex2gr)

#Divide coverage by width:
normEx2 <- nE$Exon.2_Coverage/ex2width

#Determine which novel exons are at least 10% the size of normalized exon 1
goodEx1 <- which(normCovNewEx >= 0.1*normEx1)

#Determine which novel exons are at least 10% the size of normalized exon 2
goodEx2 <- which(normCovNewEx >= 0.1*normEx2)

#Merge the exons that pass the criteria:
goodEx <- unique(goodEx1,goodEx2)

#Create a subset of the newExons GRanges object that only contains the novel exons that meet our desired criteria:
newGR <- newExons[goodEx]

#Save the GRanges object:
save(newGR, file="~/splice_project/data/analysis/newGRanges.Rdata")

