###########################################
#
# Title: Creating GRanges lists of Exon-Exon junctions for First Pass Analysis with New Supplementary Table 11, and
#        keeping all novel junctions.
#
# Description: I am using data from supplementary table 11. Table 11 contains acceptor and donor splice sites located in intronic and intergenic 
# repeat elements in the human brain. The rownames indicate the location of the splice junction that forms between an RE and an 
# annotated exon.  Columns indicate the number of reads spanning the splice junction in a sample.
#
# Make 2 new GRanges list:
# One that contains exon-exon junctions upstream from the novel exons
# Another that contains exon-exon junctions downstream from the novel exons
#
# Load in Supplementary table 11, add strand information as well. Replace table 19.
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 07/22/2019
#
###########################################

rm(list=ls()) #Clear R memory

#Download packages:
BiocManager::install(c("GenomicRanges","rtracklayer","TxDb.Hsapiens.UCSC.hg19.knownGene"))

#Load Libraries:
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Establish paths:
#Intron acceptor site
asiteintron <- read.csv("~/Documents/splice_project/data/rawdata/ST11acceptorsiteintron.csv", header = TRUE)
asiteintron <- as.data.frame(asiteintron)

#Intron donor site
dsiteintron <- read.csv("~/Documents/splice_project/data/rawdata/ST11donorsiteintron.csv", header = TRUE)
dsiteintron <- as.data.frame(dsiteintron)

#Intergenic acceptor site
asiteintergenic <- read.csv("~/Documents/splice_project/data/rawdata/ST11acceptorsiteintergenic.csv", header = TRUE)
asiteintergenic <- as.data.frame(asiteintergenic)

#Intergenic donor site
dsiteintergenic <- read.csv("~/Documents/splice_project/data/rawdata/ST11donorsiteintergenic.csv", header = TRUE)
dsiteintergenic <- as.data.frame(dsiteintergenic)

#List of chromosome positions for the first GRanges (intron acceptor site):
acceptor1 <- asiteintron$chromosomeposition[2:19979]
#List of chromosome positions for the second GRanges (intron donor site):
donor1 <- dsiteintron$chromosomeposition[2:18962]
#List of chromosome positions for the third GRanges (intergenic acceptor site):
acceptor2 <- asiteintergenic$chromosomeposition[2:4708]
#List of chromosome positions for the fourth GRanges (intergenic donor site):
donor2 <- dsiteintergenic$chromosomeposition[2:2626]

#make the user defined function to make the Genomic Ranges:
makeGR <- function(regions)
{
  a=as.character(regions)
  eNames=sapply(strsplit(a,":"),"[[",1)
  theCoords=sapply(strsplit(a,":"),"[[",2)
  theStrand=sapply(strsplit(theCoords,"_"),"[[",2)
  #theCoords=substr(theCoords,1,nchar(theCoords)-2)
  theStart=sapply(strsplit(theCoords,"-"),"[[",1)
  theEnd=sapply(strsplit(theCoords,"-"),"[[",2)
  theEnd=sapply(strsplit(theEnd,"_"),"[[",1)
  theStart=as.integer(theStart)
  theEnd=as.integer(theEnd)
  #regions=GRanges(seqnames=eNames,ranges=IRanges(start=theStart,end=theEnd),strand=theStrand)
  regions=GRanges(seqnames=eNames,ranges=IRanges(start=theStart,end=theEnd), strand = theStrand)
}

#Make genomic ranges for exon-exon junctions using the makeGR function:
asiteintronjunc <- makeGR(acceptor1)
dsiteintronjunc <- makeGR(donor1)
asiteintergenicjunc <- makeGR(acceptor2)
dsiteintergenicjunc <- makeGR(donor2)

#Save the GRanges objects:
save(asiteintronjunc, file="~/Documents/splice_project/data/analysis/GRanges_asiteintronjunc.Rdata")
save(dsiteintronjunc, file="~/Documents/splice_project/data/analysis/GRanges_dsiteintronjunc.Rdata")
save(asiteintergenicjunc, file="~/Documents/splice_project/data/analysis/GRanges_asiteintergenicjunc.Rdata")
save(dsiteintergenicjunc, file="~/Documents/splice_project/data/analysis/GRanges_dsiteintergenicjunc.Rdata")
