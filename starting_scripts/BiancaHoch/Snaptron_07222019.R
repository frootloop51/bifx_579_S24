############################################
#
# Title: Querying a list of novel exons against publicly available databases using the Snaptron search engine.
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 07/22/2019
#
###########################################

rm(list=ls()) #Clear R memory

###Set path variables that contain the different directories:###

###Install packages:
#BiocManager::install("recount")

### Load libraries:
library(recount)
library(rtracklayer)

#----CODE GOES BELOW----

#Load the GRanges object that lists the new exons that we are interested in (based on hg19 coordinates):
asiteintronjunc <- get(load("~/Documents/splice_project/data/analysis/GRanges_asiteintronjunc.Rdata"))
dsiteintronjunc <- get(load("~/Documents/splice_project/data/analysis/GRanges_dsiteintronjunc.Rdata"))
asiteintergenicjunc <- get(load("~/Documents/splice_project/data/analysis/GRanges_asiteintergenicjunc.Rdata"))
dsiteintergenicjunc <- get(load("~/Documents/splice_project/data/analysis/GRanges_dsiteintergenicjunc.Rdata"))
regions1 <- get(load("~/Documents/splice_project/data/analysis/GRanges_regions1.Rdata"))
regions2 <- get(load("~/Documents/splice_project/data/analysis/GRanges_regions2.Rdata"))

#In the following code, I will lift the features from the GRanges list (hg19) to hg38.
chain <- import.chain("~/Documents/splice_project/data/hg19ToHg38.over.chain")
asiteintron38 <- unlist(liftOver(asiteintronjunc, chain))
dsiteintron38 <- unlist(liftOver(dsiteintronjunc, chain))
asiteintergenic38 <- unlist(liftOver(asiteintergenicjunc, chain))
dsiteintergenic38 <- unlist(liftOver(dsiteintergenicjunc, chain))
regions138 <- unlist(liftOver(regions1, chain))
regions238 <- unlist(liftOver(regions2, chain))

#Execute snaptron query with the hg38 GRanges:
snaptron_query(asiteintron38, version = 'gtex')
# 2019-07-24 12:08:26 querying Snaptron
# 2019-07-24 12:19:52 processing results
# 2019-07-24 12:19:52 found no exon-exon junctions in Intropolis version gtex matching your query: this version uses hg38 coordinates.

snaptron_query(dsiteintron38, version = 'gtex')
# 2019-07-24 12:22:44 querying Snaptron
# 2019-07-24 12:33:12 processing results
# 2019-07-24 12:33:12 found no exon-exon junctions in Intropolis version gtex matching your query: this version uses hg38 coordinates.
# NULL

snaptron_query(asiteintergenic38, version = 'gtex')
# 2019-07-24 12:37:41 querying Snaptron
# 2019-07-24 12:44:26 processing results
# 2019-07-24 12:44:26 found no exon-exon junctions in Intropolis version gtex matching your query: this version uses hg38 coordinates.
# NULL

snaptron_query(dsiteintergenic38, version = 'gtex')
# 2019-07-24 12:44:53 querying Snaptron
# 2019-07-24 12:50:46 processing results
# 2019-07-24 12:50:46 found no exon-exon junctions in Intropolis version gtex matching your query: this version uses hg38 coordinates.
# NULL

snaptron_query(asiteintron38, version = 'srav2')
# 2019-07-24 12:44:53 querying Snaptron
# 2019-07-24 12:50:46 processing results
# 2019-07-24 12:50:46 found no exon-exon junctions in Intropolis version gtex matching your query: this version uses hg38 coordinates.
# NULL

snaptron_query(dsiteintron38, version = 'srav2')
# 2019-07-24 13:20:40 querying Snaptron
# 2019-07-24 13:37:48 processing results
# 2019-07-24 13:37:48 found no exon-exon junctions in Intropolis version srav2 matching your query: this version uses hg38 coordinates.
# NULL

snaptron_query(asiteintergenic38, version = 'srav2')
# 2019-07-24 13:38:25 querying Snaptron
# 2019-07-24 13:56:11 processing results
# 2019-07-24 13:56:11 found no exon-exon junctions in Intropolis version srav2 matching your query: this version uses hg38 coordinates.
# NULL

snaptron_query(dsiteintergenic38, version = 'srav2')
# 2019-07-24 13:56:38 querying Snaptron
# 2019-07-24 14:11:34 processing results
# 2019-07-24 14:11:34 found no exon-exon junctions in Intropolis version srav2 matching your query: this version uses hg38 coordinates.
# NULL

snaptron_query(asiteintronjunc, version = 'srav1')
# 2019-07-24 14:13:11 querying Snaptron
# 2019-07-24 14:19:31 processing results
# 2019-07-24 14:19:31 found no exon-exon junctions in Intropolis version srav1 matching your query: this version uses hg19 coordinates.
# NULL

#Above queries are all returning NULL. Trying a new snaptron_query with the new GRanges (regions1, 07/29/2019):
snaptron_query(regions138, version = 'srav2')
# 2019-07-29 16:33:42 querying Snaptron
# 2019-07-29 16:38:53 processing results
# 2019-07-29 16:38:53 found no exon-exon junctions in Intropolis version srav2 matching your query: this version uses hg38 coordinates.
# NULL

#Check to see if GLRB regions match those in the Snaptron example dataset.
# subset(regions138, Gene_ID....newDF.GENE_ID.x == "GLRB")
# GRanges object with 2 ranges and 1 metadata column:
#   seqnames    ranges strand | Gene_ID....newDF.GENE_ID.x
# <Rle> <IRanges>  <Rle> |                   <factor>
#   [1]     chr4 157146414      + |                       GLRB
# [2]     chr4 157146419      + |                       GLRB
# -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

#Yes, they do:
# chr4:157146414-157146414	2	coverage_sum>=1&strand=+	GLRB_2
# chr4:157146419-157146419	2	coverage_sum>=1&strand=+	GLRB_1
# chr4:157146481-157146481	1	coverage_sum>=1&strand=+	GLRB_1
# chr4:157146481-157146481	1	coverage_sum>=1&strand=+	GLRB_2

snaptron_query(regions138, version = 'gtex')
# 2019-07-29 17:05:18 querying Snaptron
# 2019-07-29 17:10:27 processing results
# 2019-07-29 17:10:27 found no exon-exon junctions in Intropolis version gtex matching your query: this version uses hg38 coordinates.
# NULL

snaptron_query(regions238, version = 'gtex')




