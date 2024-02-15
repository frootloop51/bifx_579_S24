######################################################

# Differential Expression of the putative exonic regions using recount2 database - SRP035524
# Author: Alyssa Klein
# Date: December 30, 2019
# Modified: March 25, 2020
# Created based off of the following workflows:
#   derfinder:https://bioconductor.org/packages/release/bioc/vignettes/derfinder/inst/doc/derfinder-users-guide.html#611_find_ders_with_deseq2
#   recount: http://master.bioconductor.org/packages/release/workflows/vignettes/recountWorkflow/inst/doc/recount-workflow.html#23_coverage_counts_provided_by_recount2
#   BigWigs to count matrix how-to: http://lcolladotor.github.io/protocols/bigwig_DEanalysis/#BigWigs_to_count_matrix

#####################################################

# The purpose of this script is to conduct a Differential Expression Analysis of a set of putative exons.  In this case,
# the putative exons are being observed for differential expression when regarding psychiatric disease, including 
# schizophrenia and bipolar disorder. There are three steps to complete the analysis.  The script begins with generating 
# a GenomicRanges object of the putative exons that the user wishes to check for differential expression. The second step
# includes calculating the coverage of the desired regions, scaling these counts, adding necessary colData, and combining
# all of the coverage matrices into one overall coverage matrix.  The thirs step includes Exploratory Analysis and 
# the Differential Expression Analysis.

#####################################################
# Step One: Inputting the Regions to be tested for DE Through Creation of GenomicRanges Object

# load the necessary packages needed to complete the differential expression analysis
library(GenomicRanges)
library(recount)
library(DESeq2)

# It is important to note that the coordinates in recount are from the hg38 assembly. If yours are assembly hg19,
# be sure to convert them to hg38. The UCSC Genome Browser has created a great utility for converting coordinates 
# between assemblies.  The link to the converter can be found here: https://genome.ucsc.edu/cgi-bin/hgLiftOver. 

# Create a GenomicRanges object of the putative regions
# This is equivalent to the recount package 'expressed_regions()' function, which is used to define a set of regions based on the mean bigWig file for a given project
# read in the file that contains the information necessary to create the GRanges object
region_info <- read.csv("putative_region_info_for_GRanges.csv", header= TRUE, sep = ",")

# variables for each of the GRanges components- creating the GenomicRanges object is the alternative to 
# using the 'expressed_regions' function in recount

# start and stop positions
# start positions
start_pos <- region_info$start

# end positions
end_pos <- region_info$end

# seqnames
seqnames <- rep(as.character(region_info$seqnames,length.out= length(start_pos)))

# strand
strand <- rep(as.character(region_info$strand,length.out= length(start_pos)))

# Use the variables created above to create the GRanges object
granges_putative_exons <- GRanges(seqnames= seqnames,ranges=IRanges(start= start_pos,end=end_pos),strand= strand)

#####################################################

# Step Two: Import individual bigWig data for each sample, and use it to generate a count matrix
# for the regions of interest

# Load in the recount project rse_gene in order to generate the correct bigWig URLs

if(!file.exists(file.path("SRP035524", "rse_gene_SRP035524.Rdata"))) {
  download_study("SRP035524", type = "rse-gene")
}

# Check that the file was downloaded using the file.exists() function
file.exists(file.path("SRP035524", "rse_gene_SRP035524.Rdata"))

load(file.path("SRP035524", "rse_gene_SRP035524.Rdata"))

# rename the rse_gene object to be named specifically for the study
rse_gene_SRP035524 <- rse_gene

# Construct the list of bigWig URLs
# They have the following form:
# http://duffel.rail.bio/recount/
# project id
# /bw/
# sample run id
# .bw

bws_SRP035524 <- paste0("http://duffel.rail.bio/recount/SRP035524/bw/",
                        colData(rse_gene_SRP035524)$bigwig_file)

## note that they are also present in the recount_url data.frame; alternative to preceding line of code
# bws_SRP035524 <- recount_url$url[match(colData(rse_gene)$bigwig_file,
#                             recount_url$file_name)]

# use the sample run IDs as the sample names
names(bws_SRP035524) <- colData(rse_gene_SRP035524)$run

# load the derfinder package to complete the determine region coverage
library(derfinder)
library(rtracklayer)

# Obtain coverage data for each chromosome in our regions of interest
# loop for each chromosome that contains region(s) of interest to import the bigWig file information
# needed for those specific regions
# Each loop creates an empty matric, and then it is filled in with the data from each bigWig from each sample
# Loop modified from following: http://lcolladotor.github.io/protocols/bigwig_DEanalysis/#BigWigs_to_count_matrix

# chr1
counts_chr1_035524 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr1")]), ncol = length(bws_SRP035524))
colnames(counts_chr1_035524) <- names(bws_SRP035524)
i <- 1
for(bw_file in bws_SRP035524) {
  coverage <- import(bws_SRP035524[[i]], as = 'RleList')$chr1
  counts_chr1_035524[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr1")])))
  i <- i + 1
}
# chr10
counts_chr10_035524 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr10")]), ncol = length(bws_SRP035524))
colnames(counts_chr10_035524) <- names(bws_SRP035524)
i <- 1
for(bw_file in bws_SRP035524) {
  coverage <- import(bws_SRP035524[[i]], as = 'RleList')$chr10
  counts_chr10_035524[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr10")])))
  i <- i + 1
}
# chr11
counts_chr11_035524 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr11")]), ncol = length(bws_SRP035524))
colnames(counts_chr11_035524) <- names(bws_SRP035524)
i <- 1
for(bw_file in bws_SRP035524) {
  coverage <- import(bws_SRP035524[[i]], as = 'RleList')$chr11
  counts_chr11_035524[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr11")])))
  i <- i + 1
}

# chr12
counts_chr12_035524 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr12")]), ncol = length(bws_SRP035524))
colnames(counts_chr12_035524) <- names(bws_SRP035524)
i <- 1
for(bw_file in bws_SRP035524) {
  coverage <- import(bws_SRP035524[[i]], as = 'RleList')$chr12
  counts_chr12_035524[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr12")])))
  i <- i + 1
}

# chr13
counts_chr13_035524 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr13")]), ncol = length(bws_SRP035524))
colnames(counts_chr13_035524) <- names(bws_SRP035524)
i <- 1
for(bw_file in bws_SRP035524) {
  coverage <- import(bws_SRP035524[[i]], as = 'RleList')$chr13
  counts_chr13_035524[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr13")])))
  i <- i + 1
}

# chr14
counts_chr14_035524 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr14")]), ncol = length(bws_SRP035524))
colnames(counts_chr14_035524) <- names(bws_SRP035524)
i <- 1
for(bw_file in bws_SRP035524) {
  coverage <- import(bws_SRP035524[[i]], as = 'RleList')$chr14
  counts_chr14_035524[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr14")])))
  i <- i + 1
}

# chr15
counts_chr15_035524 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr15")]), ncol = length(bws_SRP035524))
colnames(counts_chr15_035524) <- names(bws_SRP035524)
i <- 1
for(bw_file in bws_SRP035524) {
  coverage <- import(bws_SRP035524[[i]], as = 'RleList')$chr15
  counts_chr15_035524[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr15")])))
  i <- i + 1
}

# chr16
counts_chr16_035524 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr16")]), ncol = length(bws_SRP035524))
colnames(counts_chr16_035524) <- names(bws_SRP035524)
i <- 1
for(bw_file in bws_SRP035524) {
  coverage <- import(bws_SRP035524[[i]], as = 'RleList')$chr16
  counts_chr16_035524[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr16")])))
  i <- i + 1
}

# Write each of the initial count matrices to csv files JUST IN CASE!
write.csv(counts_chr1_035524, file = "counts_chr1_035524.csv", row.names= TRUE)
write.csv(counts_chr10_035524, file = "counts_chr10_035524.csv", row.names= TRUE)
write.csv(counts_chr11_035524, file = "counts_chr11_035524.csv", row.names= TRUE)
write.csv(counts_chr12_035524, file = "counts_chr12_035524.csv", row.names= TRUE)
write.csv(counts_chr13_035524, file = "counts_chr13_035524.csv", row.names= TRUE)
write.csv(counts_chr14_035524, file = "counts_chr14_035524.csv", row.names= TRUE)
write.csv(counts_chr15_035524, file = "counts_chr15_035524.csv", row.names= TRUE)
write.csv(counts_chr16_035524, file = "counts_chr16_035524.csv", row.names= TRUE)

# Need to divide by read length and round to integer numbers
# Determine the read length for each individual sample and divide by that specific read length
# Load in the phenotype data and use the "avg_read_length" column

if(!file.exists(file.path("SRP035524", "phenotype.tsv"))) {
  download_study("SRP035524", type = "phenotype")
}

# Check that the file was downloaded using the file.exists() function
file.exists(file.path("SRP035524", "SRP035524.tsv"))

# Read in the phenotype file as a csv-type (convert tsv to csv manually)
phenotype_SRP035524 <- read.csv("~/Dropbox/Hood/Alyssa Thesis/SRP035524.csv", header = TRUE, sep = ",")

# Now the counts can be divided by each sample's read length and rounded
# Completed for each count matrix for each chromosome that had a region of interest

# chr1
counts_chr1_rounded_035524 <- counts_chr1_035524
i <- 1
for (each_col in counts_chr1_035524){
  read_length <- phenotype_SRP035524$avg_read_length[which(colnames(counts_chr1_035524) == colnames(counts_chr1_035524)[i])]
  counts_chr1_rounded_035524[,i] <- round(counts_chr1_035524[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr1_rounded_035524) + 1){
    break
  }
}

# chr10
counts_chr10_rounded_035524 <- counts_chr10_035524
i <- 1
for (each_col in counts_chr10_035524){
  read_length <- phenotype_SRP035524$avg_read_length[which(colnames(counts_chr10_035524) == colnames(counts_chr10_035524)[i])]
  counts_chr10_rounded_035524[,i] <- round(counts_chr10_035524[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr10_rounded_035524_035524) + 1){
    break
  }
}

# chr11
counts_chr11_rounded_035524 <- counts_chr11_035524
i <- 1
for (each_col in counts_chr11_035524){
  read_length <- phenotype_SRP035524$avg_read_length[which(colnames(counts_chr11_035524) == colnames(counts_chr11_035524)[i])]
  counts_chr11_rounded_035524[,i] <- round(counts_chr11_035524[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr11_rounded_035524) + 1){
    break
  }
}

# chr12
counts_chr12_rounded_035524 <- counts_chr12_035524
i <- 1
for (each_col in counts_chr12_035524){
  read_length <- phenotype_SRP035524$avg_read_length[which(colnames(counts_chr12_035524) == colnames(counts_chr12_035524)[i])]
  counts_chr12_rounded_035524[,i] <- round(counts_chr12_035524[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr12_rounded_035524) + 1){
    break
  }
}

# chr13
counts_chr13_rounded_035524 <- counts_chr13_035524
i <- 1
for (each_col in counts_chr13_035524){
  read_length <- phenotype_SRP035524$avg_read_length[which(colnames(counts_chr13_035524) == colnames(counts_chr13_035524)[i])]
  counts_chr13_rounded_035524[,i] <- round(counts_chr13_035524[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr13_rounded_035524) + 1){
    break
  }
}

# chr14
counts_chr14_rounded_035524 <- counts_chr14_035524
i <- 1
for (each_col in counts_chr14_035524){
  read_length <- phenotype_SRP035524$avg_read_length[which(colnames(counts_chr14_035524) == colnames(counts_chr14_035524)[i])]
  counts_chr14_rounded_035524[,i] <- round(counts_chr14_035524[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr14_rounded_035524) + 1){
    break
  }
}

# chr15
counts_chr15_rounded_035524 <- counts_chr15_035524
i <- 1
for (each_col in counts_chr15_035524){
  read_length <- phenotype_SRP035524$avg_read_length[which(colnames(counts_chr15_035524) == colnames(counts_chr15_035524)[i])]
  counts_chr15_rounded_035524[,i] <- round(counts_chr15_035524[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr15_rounded_035524) + 1){
    break
  }
}

# chr16
counts_chr16_rounded_035524 <- counts_chr16_035524
i <- 1
for (each_col in counts_chr16_035524){
  read_length <- phenotype_SRP035524$avg_read_length[which(colnames(counts_chr16_035524) == colnames(counts_chr16_035524)[i])]
  counts_chr16_rounded_035524[,i] <- round(counts_chr16_035524[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr16_rounded_035524) + 1){
    break
  }
}

# Update the row names for each count matrix to be the region of interest
rownames(counts_chr1_rounded_035524) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr1'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))])), length.out= nrow(counts_chr1_rounded_035524)))
rownames(counts_chr10_rounded_035524) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr10'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))])), length.out= nrow(counts_chr10_rounded_035524)))
rownames(counts_chr11_rounded_035524) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr11'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))])), length.out= nrow(counts_chr11_rounded_035524)))
rownames(counts_chr12_rounded_035524) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr12'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))])), length.out= nrow(counts_chr12_rounded_035524)))
rownames(counts_chr13_rounded_035524) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr13'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))])), length.out= nrow(counts_chr13_rounded_035524)))
rownames(counts_chr14_rounded_035524) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr14'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))])), length.out= nrow(counts_chr14_rounded_035524)))
rownames(counts_chr15_rounded_035524) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr15'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))])), length.out= nrow(counts_chr15_rounded_035524)))
rownames(counts_chr16_rounded_035524) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr16'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))])), length.out= nrow(counts_chr16_rounded_035524)))

# Write each of the rounded count matrices to csv files JUST IN CASE!
write.csv(counts_chr1_rounded_035524, file = "counts_chr1_rounded_035524.csv", row.names= TRUE)
write.csv(counts_chr10_rounded_035524, file = "counts_chr10_rounded_035524.csv", row.names= TRUE)
write.csv(counts_chr11_rounded_035524, file = "counts_chr11_rounded_035524.csv", row.names= TRUE)
write.csv(counts_chr12_rounded_035524, file = "counts_chr12_rounded_035524.csv", row.names= TRUE)
write.csv(counts_chr13_rounded_035524, file = "counts_chr13_rounded_035524.csv", row.names= TRUE)
write.csv(counts_chr14_rounded_035524, file = "counts_chr14_rounded_035524.csv", row.names= TRUE)
write.csv(counts_chr15_rounded_035524, file = "counts_chr15_rounded_035524.csv", row.names= TRUE)
write.csv(counts_chr16_rounded_035524, file = "counts_chr16_rounded_035524.csv", row.names= TRUE)

# Combine all of the rounded matrices into one
counts_all_rounded_035524 <- rbind(counts_chr1_rounded_035524, counts_chr10_rounded_035524, counts_chr11_rounded_035524, counts_chr12_rounded_035524,
                                   counts_chr13_rounded_035524, counts_chr14_rounded_035524, counts_chr15_rounded_035524, counts_chr16_rounded_035524)

# Write the entire matrix for all regions to a file JUST IN CASE!
write.csv(counts_all_rounded_035524, file = "counts_all_rounded_035524.csv", row.names= TRUE)


# Explore the dimensions of the count matrix to make sure they are what is anticipated
dim(counts_chr1_rounded_035524)
# [1]  9 35
dim(counts_chr10_rounded_035524)
# [1]  3 35
dim(counts_chr11_rounded_035524)
# [1]  7 35
dim(counts_chr12_rounded_035524)
# [1]  6 35
dim(counts_chr13_rounded_035524)
# [1]  2 35
dim(counts_chr14_rounded_035524)
# [1]  4 35
dim(counts_chr15_rounded_035524)
# [1]  8 35
dim(counts_chr16_rounded_035524)
# [1]  1 35

# Save the information from https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP035524
file.exists(file.path("SRP035524", "SraRunTable_SRP035524.txt"))

## Read the table
sra_SRP035524 <- read.table(file.path("~/Dropbox/Hood/Alyssa Thesis/", "SraRunTable_SRP035524.txt"),
                            header = TRUE, sep = ",")

# want to use "disease" and "body_site"
colnames(sra_SRP035524)
# [1] "Run"                "Age"               
# [3] "Assay.Type"         "AvgSpotLen"        
# [5] "Bases"              "BioProject"        
# [7] "BioSample"          "body_site"         
# [9] "Bytes"              "Center.Name"       
# [11] "Consent"            "DATASTORE.filetype"
# [13] "DATASTORE.provider" "DATASTORE.region"  
# [15] "disease"            "Experiment"        
# [17] "Instrument"         "Library.Name"      
# [19] "LibraryLayout"      "LibrarySelection"  
# [21] "LibrarySource"      "MBases"            
# [23] "MBytes"             "Organism"          
# [25] "Platform"           "ReleaseDate"       
# [27] "sample_acc"         "Sample.Name"       
# [29] "sex"                "SRA.Study"  

# save the variables of interest to add to one variable
sra_vars_SRP035524 <- c("disease", "body_site")
sra_vars_SRP035524

# Re-organize the SRA table based on the SRA Run IDs we have
sra_SRP035524 <- sra_SRP035524[match(phenotype_SRP035524$run, sra_SRP035524$Run), ]

# Append the variables of interest to the phenotype data
phenotype_SRP035524 <- cbind(phenotype_SRP035524, sra_SRP035524[, sra_vars_SRP035524])

# check levels and set reference level (control) using relevel() function
levels(phenotype_SRP035524$disease)
# [1] "bipolar disorder (BP)" "Normal Control"       
# [3] "schizophrenia" 

levels(phenotype_SRP035524$disease) <- c("bp", "control", "schizophrenia")
phenotype_SRP035524$disease <- relevel(phenotype_SRP035524$disease, "control")

levels(phenotype_SRP035524$body_site)
# [1] "anterior cingulate" "frontal cortex"

levels(phenotype_SRP035524$body_site) <- c("anterior_cingulate", "frontal_cortex")

nrow(phenotype_SRP035524)
# [1] 35

#####################################################

# Step Three: Exploratory Analysis and Differential Expression Analysis

# Need to add a pseudocount of 1 to all counts in order to DESeq to be able to work
# Removed regions that had 0 for all samples 
# chr12:128980477-128980532_+ : only region that had 0 count for all samples- removed before adding pseudocount of 1

# Read in the adjusted file that removed regions with 0 count for all samples
counts_all_rounded_035524_for_pseudocounts <- read.csv("counts_all_rounded_035524_for_pseudocounts.csv", header = TRUE, sep = ",", row.names = 1)

# Add pseudocount of 1 to all counts
i <- 1
for (each_count in counts_all_rounded_035524_for_pseudocounts){
  counts_all_rounded_035524_for_pseudocounts[,i] <- (counts_all_rounded_035524_for_pseudocounts[,i] + 1)
  i <- i + 1
}

# Print the first six lines to make sure that there are no longer any zeroes
head(counts_all_rounded_035524_for_pseudocounts)

write.csv(counts_all_rounded_035524_for_pseudocounts, file = "counts_all_035524_pseudocount_added.csv", row.names= TRUE)

#### MD added this to test####
counts_all_rounded_035524_for_pseudocounts <- read.csv("~/Dropbox/Hood/Alyssa Thesis/counts_all_035524_pseudocount_added.csv", row.names = 1)
#################################

# set up DESeqDataSetfromMatrix 
# countdata object to store the coverage matrix that includes all of the chromosomes
countData <- counts_all_rounded_035524_for_pseudocounts
# DESeq2 needs integers in order to analyze the data
# coldata object to store the colData for rse_gene_SRP035524
colData <- phenotype_SRP035524

# countdata: a table with the fragment counts
# coldata: a table with information about the samples

# design for SRP035524 based on disease_state and tissue
# for this particular study, the DE will be based only on disease_state since there is only one tissue site 
ddsMat_SRP035524 <- DESeqDataSetFromMatrix(countData,
                                           colData,
                                           design = ~ body_site + disease)

nrow(ddsMat_SRP035524)
# [1] 39

# keep only the rows where the counts are greater than 1
keep <- rowSums(counts(ddsMat_SRP035524)) > 1
table(keep)
# TRUE
# 39      

ddsMat_SRP035524 <- ddsMat_SRP035524[keep,]
nrow(ddsMat_SRP035524)
# [1] 39

# vst transformation to do exploratory analysis
rld_SRP035524 <- rlog(ddsMat_SRP035524, blind = FALSE)
head(assay(rld_SRP035524), 6)

# observe the distances between samples visualized through a heatmap
sampleDists <- dist(t(assay(rld_SRP035524)))
sampleDists

# generate the heatmap
library(pheatmap)
library(RColorBrewer)

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste0(rld_SRP035524$disease, rld_SRP035524$body_site)
colnames(sampleDistMatrix) <- paste0(rld_SRP035524$disease, rld_SRP035524$body_site)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# principal component analysis for additional exploratory analysis
pcaData_SRP035524 <- DESeq2::plotPCA(rld_SRP035524, intgroup = c("disease","body_site"), returnData= TRUE)
pcaData_SRP035524

# generate the PCA Plot
percentVar_SRP035524 <- round(100 * attr(pcaData_SRP035524, "percentVar"))
percentVar_SRP035524
# [1] 33 15

library(ggplot2)
ggplot(pcaData_SRP035524, aes(x = PC1, y = PC2, color = disease, shape = body_site)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar_SRP035524[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_SRP035524[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with Rlog data")

# set-up the DE analysis using the 'DESeq()' function
ddsMat_SRP035524_deseq <- DESeq(ddsMat_SRP035524)

res_SRP035524 <- results(ddsMat_SRP035524_deseq)
res_SRP035524

summary(res_SRP035524)
# out of 39 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1, 2.6%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)

results(ddsMat_SRP035524_deseq, contrast = c("disease", "control", "bp"))
results(ddsMat_SRP035524_deseq, contrast = c("disease", "control", "schizophrenia"))
results(ddsMat_SRP035524_deseq, contrast = c("disease", "bp", "schizophrenia"))

# observe the number of regions whose alpha was 0.05 or less (and therefore significant in terms of differential expression)
res.05 <- results(ddsMat_SRP035524_deseq, alpha = 0.05)
table(res.05$padj < 0.05)
# FALSE  TRUE 
# 38     1
sum(res_SRP035524$pvalue < 0.05, na.rm=TRUE)
# [1] 1
sum(res_SRP035524$padj < 0.05, na.rm=TRUE)
# [1] 1
sum(!is.na(res_SRP035524$pvalue))
# [1] 39

# observe the number of regions whose alpha was 0.10 of less
res.10 <- results(ddsMat_SRP035524_deseq, alpha = 0.10)
table(res.10$padj < 0.10)
# FALSE  TRUE 
# 38     1
sum(res_SRP035524$padj < 0.1, na.rm=TRUE)
# [1] 1

# regions with up-regulation
resSig <- subset(res_SRP035524, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
# log2 fold change (MLE): disease schizophrenia vs control 
# Wald test p-value: disease schizophrenia vs control 
# DataFrame with 1 row and 6 columns
# baseMean
# <numeric>
#   chr10:119757894-119757977_+ 18.5327242588834
# log2FoldChange
# <numeric>
#   chr10:119757894-119757977_+ 4.39299793331817
# lfcSE
# <numeric>
#  chr10:119757894-119757977_+ 0.830525198195625
# stat
# <numeric>
#   chr10:119757894-119757977_+ 5.28942161280871
# pvalue
# <numeric>
#   chr10:119757894-119757977_+ 1.22703754191094e-07
# padj
# <numeric>
#   chr10:119757894-119757977_+ 4.78544641345267e-06

# regions with down-regulation
# head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

# Counts Plots for each of the regions 
rownames_de_regions_035524 <-rownames(res_SRP035524)[which(res_SRP035524$padj < 0.05)]
rownames_de_regions_035524
# [1] "chr10:119757894-119757977_+"

# for loop to generate a counts plot for each individual differentially expressed region identified from the study
i <- 1
for (each in rownames_de_regions_035524){
  # Load the ggbeeswarm package in order to use the plotCounts function to plot
  # each of the differentially expressed regions
  # "gene" argument needs to be changed each time when changing the region being plotted
  library("ggbeeswarm")
  geneCounts <- plotCounts(ddsMat_SRP035524_deseq, gene = rownames_de_regions_035524[[i]], intgroup = c("disease"),
                           returnData = TRUE)
  pdf(paste0(rownames_de_regions_035524[[i]], ".pdf"))
  print(ggplot(geneCounts, aes(x = disease, y = count, color = disease)) +
          scale_y_log10() +  geom_beeswarm(cex = 3) + ggtitle(paste0("Counts Plot:", rownames_de_regions_035524[[i]])))
  dev.off()
  i <- i + 1
}


