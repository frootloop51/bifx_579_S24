######################################################

# Differential Expression of the putative exonic regions using recount2 database - SRP033725
# Author: Alyssa Klein
# Date: December 30, 2019
# Modified: March 22, 2020
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

# Construct the list of bigWig URLs
# They have the following form:
# http://duffel.rail.bio/recount/
# project id
# /bw/
# sample run id
# .bw

bws_SRP033725 <- paste0("http://duffel.rail.bio/recount/SRP033725/bw/",
              colData(rse_gene)$bigwig_file)

## note that they are also present in the recount_url data.frame; alternative to preceding line of code
# bws_SRP033725 <- recount_url$url[match(colData(rse_gene)$bigwig_file,
#                             recount_url$file_name)]

# use the sample run IDs as the sample names
names(bws_SRP033725) <- colData(rse_gene_SRP033725)$run

# load the derfinder package to complete the determine region coverage
library(derfinder)
library(rtracklayer)

# Obtain coverage data for each chromosome in our regions of interest
# loop for each chromosome that contains region(s) of interest to import the bigWig file information
# needed for those specific regions
# Each loop creates an empty matric, and then it is filled in with the data from each bigWig from each sample
# Loop modified from following: http://lcolladotor.github.io/protocols/bigwig_DEanalysis/#BigWigs_to_count_matrix

# chr1
counts_chr1 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr1")]), ncol = length(bws_SRP033725))
colnames(counts_chr1) <- names(bws_SRP033725)
i <- 1
for(bw_file in bws_SRP033725) {
  coverage <- import(bws_SRP033725[[i]], as = 'RleList')$chr1
  counts_chr1[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr1")])))
  i <- i + 1
}
# chr10
counts_chr10 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr10")]), ncol = length(bws_SRP033725))
colnames(counts_chr10) <- names(bws_SRP033725)
i <- 1
for(bw_file in bws_SRP033725) {
  coverage <- import(bws_SRP033725[[i]], as = 'RleList')$chr10
  counts_chr10[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr10")])))
  i <- i + 1
}
# chr11
counts_chr11 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr11")]), ncol = length(bws_SRP033725))
colnames(counts_chr11) <- names(bws_SRP033725)
i <- 1
for(bw_file in bws_SRP033725) {
  coverage <- import(bws_SRP033725[[i]], as = 'RleList')$chr11
  counts_chr11[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr11")])))
  i <- i + 1
}

# chr12
counts_chr12 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr12")]), ncol = length(bws_SRP033725))
colnames(counts_chr12) <- names(bws_SRP033725)
i <- 1
for(bw_file in bws_SRP033725) {
  coverage <- import(bws_SRP033725[[i]], as = 'RleList')$chr12
  counts_chr12[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr12")])))
  i <- i + 1
}

# chr13
counts_chr13 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr13")]), ncol = length(bws_SRP033725))
colnames(counts_chr13) <- names(bws_SRP033725)
i <- 1
for(bw_file in bws_SRP033725) {
  coverage <- import(bws_SRP033725[[i]], as = 'RleList')$chr13
  counts_chr13[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr13")])))
  i <- i + 1
}

# chr14
counts_chr14 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr14")]), ncol = length(bws_SRP033725))
colnames(counts_chr14) <- names(bws_SRP033725)
i <- 1
for(bw_file in bws_SRP033725) {
  coverage <- import(bws_SRP033725[[i]], as = 'RleList')$chr14
  counts_chr14[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr14")])))
  i <- i + 1
}

# chr15
counts_chr15 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr15")]), ncol = length(bws_SRP033725))
colnames(counts_chr15) <- names(bws_SRP033725)
i <- 1
for(bw_file in bws_SRP033725) {
  coverage <- import(bws_SRP033725[[i]], as = 'RleList')$chr15
  counts_chr15[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr15")])))
  i <- i + 1
}

# chr16
counts_chr16 <- matrix(NA, nrow = length(granges_putative_exons[which(seqnames == "chr16")]), ncol = length(bws_SRP033725))
colnames(counts_chr16) <- names(bws_SRP033725)
i <- 1
for(bw_file in bws_SRP033725) {
  coverage <- import(bws_SRP033725[[i]], as = 'RleList')$chr16
  counts_chr16[, i] <- sum(Views(coverage, ranges(granges_putative_exons[which(seqnames == "chr16")])))
  i <- i + 1
}

# Write each of the initial count matrices to csv files JUST IN CASE!
write.csv(counts_chr1, file = "counts_chr1.csv", row.names= TRUE)
write.csv(counts_chr10, file = "counts_chr10.csv", row.names= TRUE)
write.csv(counts_chr11, file = "counts_chr11.csv", row.names= TRUE)
write.csv(counts_chr12, file = "counts_chr12.csv", row.names= TRUE)
write.csv(counts_chr13, file = "counts_chr13.csv", row.names= TRUE)
write.csv(counts_chr14, file = "counts_chr14.csv", row.names= TRUE)
write.csv(counts_chr15, file = "counts_chr15.csv", row.names= TRUE)
write.csv(counts_chr15, file = "counts_chr16.csv", row.names= TRUE)

# Need to divide by read length and round to integer numbers
# Determine the read length for each individual sample and divide by that specific read length
# Load in the phenotype data and use the "avg_read_length" column

if(!file.exists(file.path("SRP033725", "phenotype.tsv"))) {
  download_study("SRP033725", type = "phenotype")
}

# Check that the file was downloaded using the file.exists() function
file.exists(file.path("SRP033725", "SRP033725.tsv"))

phenotype_SRP033725 <- read.csv("SRP033725/SRP033725.csv", header = TRUE, sep = ",")

# Now the counts can be divided by each sample's read length and rounded
# Completed for each count matrix for each chromosome that had a region of interest

# chr1
counts_chr1_rounded <- counts_chr1
i <- 1
for (each_col in counts_chr1){
  read_length <- phenotype_SRP033725$avg_read_length[which(colnames(counts_chr1) == colnames(counts_chr1)[i])]
  counts_chr1_rounded[,i] <- round(counts_chr1[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr1_rounded) + 1){
    break
  }
}

# chr10
counts_chr10_rounded <- counts_chr10
i <- 1
for (each_col in counts_chr10){
  read_length <- phenotype_SRP033725$avg_read_length[which(colnames(counts_chr10) == colnames(counts_chr10)[i])]
  counts_chr10_rounded[,i] <- round(counts_chr10[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr10_rounded) + 1){
    break
  }
}

# chr11
counts_chr11_rounded <- counts_chr11
i <- 1
for (each_col in counts_chr11){
  read_length <- phenotype_SRP033725$avg_read_length[which(colnames(counts_chr11) == colnames(counts_chr11)[i])]
  counts_chr11_rounded[,i] <- round(counts_chr11[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr11_rounded) + 1){
    break
  }
}

# chr12
counts_chr12_rounded <- counts_chr12
i <- 1
for (each_col in counts_chr12){
  read_length <- phenotype_SRP033725$avg_read_length[which(colnames(counts_chr12) == colnames(counts_chr12)[i])]
  counts_chr12_rounded[,i] <- round(counts_chr12[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr12_rounded) + 1){
    break
  }
}

# chr13
counts_chr13_rounded <- counts_chr13
i <- 1
for (each_col in counts_chr13){
  read_length <- phenotype_SRP033725$avg_read_length[which(colnames(counts_chr13) == colnames(counts_chr13)[i])]
  counts_chr13_rounded[,i] <- round(counts_chr13[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr13_rounded) + 1){
    break
  }
}

# chr14
counts_chr14_rounded <- counts_chr14
i <- 1
for (each_col in counts_chr14){
  read_length <- phenotype_SRP033725$avg_read_length[which(colnames(counts_chr14) == colnames(counts_chr14)[i])]
  counts_chr14_rounded[,i] <- round(counts_chr14[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr14_rounded) + 1){
    break
  }
}

# chr15
counts_chr15_rounded <- counts_chr15
i <- 1
for (each_col in counts_chr15){
  read_length <- phenotype_SRP033725$avg_read_length[which(colnames(counts_chr15) == colnames(counts_chr15)[i])]
  counts_chr15_rounded[,i] <- round(counts_chr15[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr15_rounded) + 1){
    break
  }
}

# chr16
counts_chr16_rounded <- counts_chr16
i <- 1
for (each_col in counts_chr16){
  read_length <- phenotype_SRP033725$avg_read_length[which(colnames(counts_chr16) == colnames(counts_chr16)[i])]
  counts_chr16_rounded[,i] <- round(counts_chr16[,i] / read_length, 0)
  i <- i + 1
  if(i == ncol(counts_chr16_rounded) + 1){
    break
  }
}

# Update the row names for each count matrix to be the region of interest
rownames(counts_chr1_rounded) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr1'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))])), length.out= nrow(counts_chr1_rounded)))
rownames(counts_chr10_rounded) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr10'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))])), length.out= nrow(counts_chr10_rounded)))
rownames(counts_chr11_rounded) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr11'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))])), length.out= nrow(counts_chr11_rounded)))
rownames(counts_chr12_rounded) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr12'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))])), length.out= nrow(counts_chr12_rounded)))
rownames(counts_chr13_rounded) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr13'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))])), length.out= nrow(counts_chr13_rounded)))
rownames(counts_chr14_rounded) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr14'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))])), length.out= nrow(counts_chr14_rounded)))
rownames(counts_chr15_rounded) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr15'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))])), length.out= nrow(counts_chr15_rounded)))
rownames(counts_chr16_rounded) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr16'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))])), length.out= nrow(counts_chr16_rounded)))

# Write each of the rounded count matrices to csv files JUST IN CASE!
write.csv(counts_chr1_rounded, file = "counts_chr1_rounded.csv", row.names= TRUE)
write.csv(counts_chr10_rounded, file = "counts_chr10_rounded.csv", row.names= TRUE)
write.csv(counts_chr11_rounded, file = "counts_chr11_rounded.csv", row.names= TRUE)
write.csv(counts_chr12_rounded, file = "counts_chr12_rounded.csv", row.names= TRUE)
write.csv(counts_chr13_rounded, file = "counts_chr13_rounded.csv", row.names= TRUE)
write.csv(counts_chr14_rounded, file = "counts_chr14_rounded.csv", row.names= TRUE)
write.csv(counts_chr15_rounded, file = "counts_chr15_rounded.csv", row.names= TRUE)
write.csv(counts_chr16_rounded, file = "counts_chr16_rounded.csv", row.names= TRUE)

# Combine all of the rounded matrices into one
counts_all_rounded <- rbind(counts_chr1_rounded, counts_chr10_rounded, counts_chr11_rounded, counts_chr12_rounded,
                            counts_chr13_rounded, counts_chr14_rounded, counts_chr15_rounded, counts_chr16_rounded)

# Write the entire matrix for all regions to a file JUST IN CASE!
write.csv(counts_all_rounded, file = "counts_all_rounded.csv", row.names= TRUE)


# Explore the dimensions of the count matrix to make sure they are what is anticipated
dim(counts_chr1_rounded)
# [1]  9 61
dim(counts_chr10_rounded)
# [1]  3 61
dim(counts_chr11_rounded)
# [1]  7 61
dim(counts_chr12_rounded)
# [1]  6 61
dim(counts_chr13_rounded)
# [1]  2 61
dim(counts_chr14_rounded)
# [1]  4 61
dim(counts_chr15_rounded)
# [1]  8 61
dim(counts_chr16_rounded)
# [1]  1 61

# Save the information from https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP033725
file.exists(file.path("SRP033725", "SraRunTable_SRP033725.txt"))

# Read the table
sra_SRP033725 <- read.table(file.path("SRP033725", "SraRunTable_SRP033725.txt"),
                            header = TRUE, sep = ",")


# want to use "disease_state" and "tissue" since we are interested in differential expression in disease state and tissue site
colnames(sra_SRP033725)
# [1] "Run"                "Assay.Type"        
# [3] "AvgSpotLen"         "Bases"             
# [5] "BioProject"         "BioSample"         
# [7] "Bytes"              "Center.Name"       
# [9] "Consent"            "DATASTORE.filetype"
# [11] "DATASTORE.provider" "DATASTORE.region"  
# [13] "disease_state"      "Experiment"        
# [15] "GEO_Accession"      "Instrument"        
# [17] "LibraryLayout"      "LibrarySelection"  
# [19] "LibrarySource"      "MBases"            
# [21] "MBytes"             "Organism"          
# [23] "Platform"           "ReleaseDate"       
# [25] "sample_acc"         "Sample.Name"       
# [27] "source_name"        "SRA.Study"         
# [29] "tissue"

# save the variables of interest to add to one variable
sra_vars_SRP033725 <- c("disease_state", "tissue")
sra_vars_SRP033725

# Re-organize the SRA table based on the SRA Run IDs we have
sra_SRP033725 <- sra_SRP033725[match(phenotype_SRP033725$run, sra_SRP033725$Run), ]

# Double check the order to make sure that the phenotype run run column is equivalent to the sra Run column
identical(phenotype_SRP033725$run, sra_SRP033725$Run)

# Append the variables of interest to the phenotype data
phenotype_SRP033725 <- cbind(phenotype_SRP033725, sra_SRP033725[, sra_vars_SRP033725])

# check levels and set reference level (control) using relevel() function
levels(phenotype_SRP033725$disease_state)
# [1] "BD"      "Control" 

levels(phenotype_SRP033725$disease_state) <- c("bp", "control")
phenotype_SRP033725$disease_state <- relevel(phenotype_SRP033725$disease_state, "control")
levels(phenotype_SRP033725$disease_state)
# [1] "control" "bp"     

levels(phenotype_SRP033725$tissue)
# [1] "dorsolateral prefrontal cortex"

levels(phenotype_SRP033725$tissue) <- ("dorsolateral_prefrontal_cortex")

levels(phenotype_SRP033725$tissue)
# [1] "dorsolateral_prefrontal_cortex"

nrow(phenotype_SRP033725)
# [1] 62

#####################################################

# Step Three: Exploratory Analysis and Differential Expression Analysis

# set up DESeqDataSetfromMatrix 
# countdata object to store the coverage matrix that includes all of the chromosomes
countData <- counts_all_rounded
# DESeq2 needs integers in order to analyze the data
# coldata object to store the colData for rse_gene_SRP037725
colData <- phenotype_SRP033725

# need to remove SRR1047862 because it does not have counts/bigWig file
which(colData$run == "SRR1047862")
# [1] 62

# remove row from colData
colData <- colData[-c(62),]
nrow(colData)
# [1] 61
# countdata: a table with the fragment counts
# coldata: a table with information about the samples

# design for SRP033725 based on disease_state and tissue
# for this particular study, the DE will be based only on disease_state since there is only one tissue site 
ddsMat_SRP033725 <- DESeqDataSetFromMatrix(countData,
                                           colData,
                                           design = ~ disease_state)

nrow(ddsMat_SRP033725)
# [1] 40

# keep only the rows where the counts are greater than 1
keep <- rowSums(counts(ddsMat_SRP033725)) > 1
table(keep)
# TRUE 
# 40

ddsMat_SRP033725 <- ddsMat_SRP033725[keep,]
nrow(ddsMat_SRP033725)
# [1] 40

# vst transformation to do exploratory analysis
rld_SRP033725 <- rlog(ddsMat_SRP033725, blind = FALSE)
head(assay(rld_SRP033725), 6)

# observe the distances between samples visualized through a heatmap
sampleDists <- dist(t(assay(rld_SRP033725)))
sampleDists

# generate the heatmap
library(pheatmap)
library(RColorBrewer)

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste0(rld_SRP033725$disease_state, rld_SRP033725$run)
colnames(sampleDistMatrix) <- paste0(rld_SRP033725$disease_state, rld_SRP033725$run)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# principal component analysis for additional exploratory analysis
pcaData_SRP033725 <- DESeq2::plotPCA(rld_SRP033725, intgroup = "disease_state", returnData= TRUE)
pcaData_SRP033725

# generate the PCA Plot
percentVar_SRP033725 <- round(100 * attr(pcaData_SRP033725, "percentVar"))
percentVar_SRP033725
# [1] 27  10

library(ggplot2)
ggplot(pcaData_SRP033725, aes(x = PC1, y = PC2, color = disease_state, shape = disease_state)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar_SRP033725[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_SRP033725[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with Rlog data")

# set-up the DE analysis using the 'DESeq()' function
ddsMat_SRP033725_deseq <- DESeq(ddsMat_SRP033725)

res_SRP033725 <- results(ddsMat_SRP033725_deseq)
res_SRP033725

summary(res_SRP033725)
# out of 40 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1, 2.5%
# LFC < 0 (down)     : 3, 7.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)

# observe the number of regions whose alpha was 0.05 or less (and therefore significant in terms of differential expression)
res.05 <- results(ddsMat_SRP033725_deseq, alpha = 0.05)
table(res.05$padj < 0.05)
# FALSE  TRUE 
# 36     4
sum(res_SRP033725$pvalue < 0.05, na.rm=TRUE)
# [1] 7
sum(res_SRP033725$padj < 0.05, na.rm=TRUE)
# [1] 4
sum(!is.na(res_SRP033725$pvalue))
# [1] 40

# observe the number of regions whose alpha was 0.10 of less
res.10 <- results(ddsMat_SRP033725_deseq, alpha = 0.10)
table(res.10$padj < 0.10)
# FALSE  TRUE 
# 36     4
sum(res_SRP033725$padj < 0.1, na.rm=TRUE)
# [1] 4

# regions with up-regulation
resSig <- subset(res_SRP033725, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
# log2 fold change (MLE): disease state bp vs control 
# Wald test p-value: disease state bp vs control 
# DataFrame with 4 rows and 6 columns
# baseMean
# <numeric>
# chr11:63746946-6807088_+   6.34615667249411
# chr15:75007660-155866737_+ 20.9699268812706
# chr15:75007587-155310453_+ 30.8024523543537
# chr15:75351300-180052316_+ 5.95426189025667
# log2FoldChange
# <numeric>
# chr11:63746946-6807088_+   -0.609296498800728
# chr15:75007660-155866737_+ -0.414103371682752
# chr15:75007587-155310453_+ -0.378107453337308
# chr15:75351300-180052316_+  0.519656959266494
# lfcSE
# <numeric>
#  chr11:63746946-6807088_+   0.166409945015499
# chr15:75007660-155866737_+ 0.123808856698154
# chr15:75007587-155310453_+ 0.105811292869366
# chr15:75351300-180052316_+  0.18485354717531
# stat
# <numeric>
# chr11:63746946-6807088_+   -3.66141878566199
# chr15:75007660-155866737_+ -3.34469910090792
# chr15:75007587-155310453_+ -3.57341303639601
# chr15:75351300-180052316_+  2.81118197192974
# pvalue
# <numeric>
# chr11:63746946-6807088_+   0.000250822373668652
# chr15:75007660-155866737_+ 0.000823718831696048
# chr15:75007587-155310453_+ 0.000352358299948455
# chr15:75351300-180052316_+  0.00493598606949876
# padj
# <numeric>
# chr11:63746946-6807088_+   0.00704716599896909
# chr15:75007660-155866737_+  0.0109829177559473
# chr15:75007587-155310453_+ 0.00704716599896909
# chr15:75351300-180052316_+  0.0493598606949876

# regions with down-regulation
# head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

# Counts Plots for each of the regions 
rownames_de_regions <-rownames(res_SRP033725)[which(res_SRP033725$padj < 0.05)]
# [1] "chr11:63746946-6807088_+"   "chr15:75007587-155310453_+"
# [3] "chr15:75007660-155866737_+" "chr15:75351300-180052316_+"

# for loop to generate a counts plot for each individual differentially expressed region identified from the study
i <- 1
for (each in rownames_de_regions){
  # Load the ggbeeswarm package in order to use the plotCounts function to plot
  # each of the differentially expressed regions
  # "gene" argument needs to be changed each time when changing the region being plotted
  library("ggbeeswarm")
  geneCounts <- plotCounts(ddsMat_SRP033725_deseq, gene = rownames_de_regions[[i]], intgroup = c("disease_state"),
                           returnData = TRUE)
  pdf(paste0(rownames_de_regions[[i]], ".pdf"))
  print(ggplot(geneCounts, aes(x = disease_state, y = count, color = disease_state)) +
    scale_y_log10() +  geom_beeswarm(cex = 3) + ggtitle(paste0("Counts Plot:", rownames_de_regions[[i]])))
  dev.off()
  i <- i + 1
}


