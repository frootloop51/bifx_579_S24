############################################
#
# Title: recount-brain metadata applied to the SRA study ERP001304.
# Abstract: The RNA was extracted from post-mortem cortical grey matter of STG from the left hemisphere of 9 pairs 
# of male subjects with schizophrenia and matched non-psychiatric controls. 
# Each sample was subjected to 76 cycles of sequencing from a single end in one lane of Illumina Genome Analyzer II.
# Total number of samples: 18
#
# Workflows and references used:
# This code taken from the derfinder quick start vignette gets the region matrix from a set of BigWig
# files included in derfinderData.http://lcolladotor.github.io/rfigshare_test/create_data.html
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 05/21/2019
#
###########################################

rm(list=ls()) #Clear R memory

###Set path variables that contain the different directories:###

###Install packages:
#install.packages("BiocManager")
BiocManager::install(c("recount", "Gviz", "derfinder", "derfinderData", "GenomicRanges"))


### Load libraries:
library('recount')
library('Gviz')
library('derfinder')
library('derfinderData')
library('GenomicRanges')

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


#Download the first project ("ERP001304"):
download_study(SRAs[4])

?download_study
# 2019-05-29 13:46:28 downloading file rse_gene.Rdata to ERP001304
# trying URL 'http://duffel.rail.bio/recount/v2/SRP035524/rse_gene.Rdata'
# Content type 'application/octet-stream' length 5275872 bytes (5.0 MB)
# downloaded 5.0 MB

#Load the data:
load(file.path((SRAs[4]), 'rse_gene.Rdata'))

#View the RangedSummarizedExperiment object:
rse_gene
str(rse_gene)

#Total number of samples used:
length(unique(rse_gene$sample))
#18

#Add sample metadata from recount-brain to "ERP001304":
RSEmeta <- add_metadata(rse_gene)
RSEmeta

#View sample phenotype data provided by the recount project:
colData(RSEmeta)

#View genes and base pair lengths:
rowData(RSEmeta)

#Once we have identified the study of interest, we can use the browse_study() function to browse the study at the SRA website:
browse_study(SRAs[4])

#Preparing for differential expression analysis:
#Scale counts by taking into account the total coverage per sample:
rse <- scale_counts(rse_gene)
rse

##### Details about counts #####

#Scale counts to 40 million mapped reads. Not all mapped reads are in exonic
#sequences, so the total is not necessarily 40 million.
colSums(assays(rse)$counts) / 1e6

# ERR103423 ERR103435 ERR103437 ERR103434 ERR103424 ERR103422 ERR103438 ERR103436 ERR103432 
# 36.69383  36.27202  33.81209  38.89148  35.47921  38.30259  32.61524  38.55911  36.97086 
# ERR103433 ERR103428 ERR103431 ERR103426 ERR103430 ERR103427 ERR103429 ERR103425 ERR103421 
# 38.21954  38.09043  38.74372  38.00316  39.04242  34.75546  38.26299  38.86450  38.74521 

#Compute read counts:
rse_read_counts <- read_counts(rse_gene)

#Difference between read counts and number of reads downloaded by Rail-RNA:
colSums(assays(rse_read_counts)$counts) / 1e6 -
  colData(rse_read_counts)$reads_downloaded / 1e6

#06/03/2019
#Load the GRanges object that lists the new exons that we are interested in:
newGR <- ("~/Documents/splice_project/data/analysis/newGRanges.Rdata")
load(newGR)
newGR

#Download the sample BigWig files for the project ("ERP001304"):
#download_study('ERP001304', type = "samples", outdir = "~Documents/splice_project/data/rawdata/bigwigs", download = TRUE)

#Construct full paths to the BigWig files and fix the names:
bwfiles <- rawFiles(datadir = "~/Documents/splice_project/data/rawdata/bigwigs/", samplepatt = '*.bw' , fileterm = NULL)
names(bwfiles) <- gsub('.bw', '', names(bwfiles))

#Get a list of chromosomes from the GRanges list containing the new exons:
chromosomes <- levels(seqnames(newGR))

#Load the unfiltered coverage information from a group of bigwig files and a list of chromosomes.
#For a group of samples this function reads the coverage information for several chromosomes directly from the BAM files. 
#Per chromosome, it merges the unfiltered coverage by sample into a DataFrame. 
#The end result is a list with one such DataFrame objects per chromosome.

#fullCov <- fullCoverage(files = bwfiles, chrs = chromosomes, verbose = FALSE)
#save(fullCov, file="~/Documents/splice_project/data/analysis/fullCov.Rdata")
fullCov <- ("~/Documents/splice_project/data/analysis/fullCov.Rdata")
load(fullCov)

#Identify regions data by a coverage filter and get a count matrix.
#Given a set of un-filtered coverage data (see fullCoverage), 
#create candidate regions by applying a cutoff on the coverage values, 
#and obtain a count matrix where the number of rows corresponds to the number of candidate regions 
#and the number of columns corresponds to the number of samples. 
#The values are the mean coverage for a given sample for a given region.

#regionMat <- regionMatrix(fullCov, cutoff = 30, L = 76, verbose = FALSE)
#save(regionMat, file="~/Documents/splice_project/data/analysis/regionMat.Rdata")
regionMat <- ("~/Documents/splice_project/data/analysis/regionMat.Rdata")
load(regionMat)

#The following code creates a count table, and a RangedSummarixedExperiment object:
#Round matrix
counts <- round(regionMat$chr1$coverageMatrix)

#Reorganize RSEmeta so that it can be used as colData for the new RSE (rseNew)
orMeta <- colData(RSEmeta)
orMeta <- orMeta[order(rownames(orMeta)),]
orMeta

#Check to see if orMeta is identical to the colnames of counts
identical(rownames(orMeta),colnames(counts))

#Create a RSE
rseNew <- SummarizedExperiment(assays = list('counts' = counts), colData = orMeta, rowRanges = regionMat$chr1$regions)




