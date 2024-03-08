#Install Bioconductor for the first time
#if (!require("BiocManager", quietly = TRUE))
#  +   install.packages("BiocManager")
#BiocManager::install(version = "3.18")
  
#Install Bioconductor Packages for the first time
#BiocManager::install("GenomicRanges", version = "3.18")
#Load Packages after initial installation
library(GenomicRanges)
library(readr)

NovelCommonMajorIsoforms <- read_csv("Dropbox/Hood/BiancaThesis/NovelCommonMajorIsoforms.csv")

#Variables:
novelExons <- NovelCommonMajorIsoforms$PEx
(charEx <- as.character(novelExons))

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

#Make the Common Major Isoform Novel Exon GRanges lists:
novelExons_CM <- GRanges(seqnames= chr, IRanges(start = theStart, end = theEnd), strand = NovelCommonMajorIsoforms$strand.downstream)
novelExons_CM

#A Genomic Ranges Object is an array of vectors. Treat it like a vector most of the time as in unique() below.
(novelExons_CM <- unique(novelExons_CM))

#But you can also use functions to call the vectors
(chromsome <- as.character(seqnames(novelExons_CM)))
(start <- as.numeric(start(novelExons_CM)))
(end <- as.numeric(end(novelExons_CM)))
(strand <- as.character(strand(novelExons_CM)))

#Make a new csv with just the exons.
(NovelCommonMajorExons <-cbind(chromsome,start,end,strand))
write.csv(NovelCommonMajorExons,"Dropbox/Hood/BiancaThesis/NovelCommonMajorExons.csv",col.names = TRUE, row.names = FALSE)
