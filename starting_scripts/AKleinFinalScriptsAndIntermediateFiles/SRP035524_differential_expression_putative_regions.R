######################################################

# Differential Expression of the putative exonic regions using recount2 database - SRP035524
# Author: Alyssa Klein
# Date: December 30, 2019
# Modified: February 14, 2020
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
# This is equivalent to expressed_regions(), which is used to define a set of regions based on the mean bigWig file for a given project
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

# Convert the GenomicRanges object to a dataframe
# Make the rownames of the dataframe reflect the region (chr#: coordinates_strand) to be used when creating coverage matrix

#####################################################

# Step Two: Calculating coverage of the desired regions, scaling the counts, adding necessary colData, and combining
# all of the coverage matrices

# Load in the recount project to be used in the differential expression analysis

if(!file.exists(file.path("SRP035524", "rse_gene_SRP035524.Rdata"))) {
  download_study("SRP035524", type = "rse-gene")
}

# Check that the file was downloaded using the file.exists() function
file.exists(file.path("SRP035524", "rse_gene_SRP035524.Rdata"))

load(file.path("SRP035524", "rse_gene_SRP035524.Rdata"))

rse_gene_SRP035524 <- rse_gene

# Calculate the coverage over the putative regions (using coverage_matrix() function in recount)
rse_er_SRP035524_chr1 <- coverage_matrix("SRP035524", "chr1", granges_putative_exons,
                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP035524_chr10 <- coverage_matrix("SRP035524", "chr10", granges_putative_exons,
                                         verboseLoad = FALSE, scale = FALSE)
rse_er_SRP035524_chr11 <- coverage_matrix("SRP035524", "chr11", granges_putative_exons,
                                         verboseLoad = FALSE, scale = FALSE)
rse_er_SRP035524_chr12 <- coverage_matrix("SRP035524", "chr12", granges_putative_exons,
                                         verboseLoad = FALSE, scale = FALSE)
rse_er_SRP035524_chr13 <- coverage_matrix("SRP035524", "chr13", granges_putative_exons,
                                         verboseLoad = FALSE, scale = FALSE)
rse_er_SRP035524_chr14 <- coverage_matrix("SRP035524", "chr14", granges_putative_exons,
                                         verboseLoad = FALSE, scale = FALSE)
rse_er_SRP035524_chr15 <- coverage_matrix("SRP035524", "chr15", granges_putative_exons,
                                         verboseLoad = FALSE, scale = FALSE)
rse_er_SRP035524_chr16 <- coverage_matrix("SRP035524", "chr16", granges_putative_exons,
                                         verboseLoad = FALSE, scale = FALSE)

# rename the columns names for chromosome 16- for some reason they are not named with the same schema as the others
colnames(rse_er_SRP035524_chr16) <- rse_er_SRP035524_chr16$run


# scale the counts for each of the coverage matrices for each chromosome
rse_er_scaled_SRP035524_chr1 <- scale_counts(rse_er_SRP035524_chr1)
rse_er_scaled_SRP035524_chr10 <- scale_counts(rse_er_SRP035524_chr10)
rse_er_scaled_SRP035524_chr11 <- scale_counts(rse_er_SRP035524_chr11)
rse_er_scaled_SRP035524_chr12 <- scale_counts(rse_er_SRP035524_chr12)
rse_er_scaled_SRP035524_chr13 <- scale_counts(rse_er_SRP035524_chr13)
rse_er_scaled_SRP035524_chr14 <- scale_counts(rse_er_SRP035524_chr14)
rse_er_scaled_SRP035524_chr15 <- scale_counts(rse_er_SRP035524_chr15)
rse_er_scaled_SRP035524_chr16 <- scale_counts(rse_er_SRP035524_chr16)

# Use the expanded metadata we built for the gene model- here we are simply updating the colData for the scaled
# counts with the colData for the original rse_gene(s)
colData(rse_er_scaled_SRP035524_chr1) <- colData(rse_gene_SRP035524)
colData(rse_er_scaled_SRP035524_chr10) <- colData(rse_gene_SRP035524)
colData(rse_er_scaled_SRP035524_chr11) <- colData(rse_gene_SRP035524)
colData(rse_er_scaled_SRP035524_chr12) <- colData(rse_gene_SRP035524)
colData(rse_er_scaled_SRP035524_chr13) <- colData(rse_gene_SRP035524)
colData(rse_er_scaled_SRP035524_chr14) <- colData(rse_gene_SRP035524)
colData(rse_er_scaled_SRP035524_chr15) <- colData(rse_gene_SRP035524)
colData(rse_er_scaled_SRP035524_chr16) <- colData(rse_gene_SRP035524)

# Save the information from https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP035524
file.exists(file.path("SRP035524", "SraRunTable_SRP035524.txt"))

## Read the table
sra_SRP035524 <- read.table(file.path("SRP035524", "SraRunTable_SRP035524.txt"),
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
sra_vars <- c("disease", "body_site")
sra_vars

# Re-organize the SRA table based on the SRA Run IDs we have
sra_SRP035524 <- sra_SRP035524[match(colData(rse_er_scaled_SRP035524_chr1)$run, sra_SRP035524$Run), ]

# Double check the order to make sure that the rse_er run column is equivalent to the sra acc column
identical(colData(rse_er_scaled_SRP035524_chr1)$run, as.character(sra_SRP035524$Run))
identical(colData(rse_er_scaled_SRP035524_chr10)$run, as.character(sra_SRP035524$Run))
identical(colData(rse_er_scaled_SRP035524_chr11)$run, as.character(sra_SRP035524$Run))
identical(colData(rse_er_scaled_SRP035524_chr12)$run, as.character(sra_SRP035524$Run))
identical(colData(rse_er_scaled_SRP035524_chr13)$run, as.character(sra_SRP035524$Run))
identical(colData(rse_er_scaled_SRP035524_chr14)$run, as.character(sra_SRP035524$Run))
identical(colData(rse_er_scaled_SRP035524_chr15)$run, as.character(sra_SRP035524$Run))
identical(colData(rse_er_scaled_SRP035524_chr16)$run, as.character(sra_SRP035524$Run))

# Append the variables of interest to the colData for the rse_er (specifically rse_er_scaled)
colData(rse_gene_SRP035524) <- cbind(colData(rse_gene_SRP035524), sra_SRP035524[, sra_vars])
colData(rse_er_scaled_SRP035524_chr1) <- cbind(colData(rse_er_SRP035524_chr1), sra_SRP035524[, sra_vars])
colData(rse_er_scaled_SRP035524_chr10) <- cbind(colData(rse_er_SRP035524_chr10), sra_SRP035524[, sra_vars])
colData(rse_er_scaled_SRP035524_chr11) <- cbind(colData(rse_er_SRP035524_chr11), sra_SRP035524[, sra_vars])
colData(rse_er_scaled_SRP035524_chr12) <- cbind(colData(rse_er_SRP035524_chr12), sra_SRP035524[, sra_vars])
colData(rse_er_scaled_SRP035524_chr13) <- cbind(colData(rse_er_SRP035524_chr13), sra_SRP035524[, sra_vars])
colData(rse_er_scaled_SRP035524_chr14) <- cbind(colData(rse_er_SRP035524_chr14), sra_SRP035524[, sra_vars])
colData(rse_er_scaled_SRP035524_chr15) <- cbind(colData(rse_er_SRP035524_chr15), sra_SRP035524[, sra_vars])
colData(rse_er_scaled_SRP035524_chr16) <- cbind(colData(rse_er_SRP035524_chr16), sra_SRP035524[, sra_vars])

# check levels and set reference level (control) using relevel() function
levels(rse_gene_SRP035524$disease)
# [1] "bipolar disorder (BP)" "Normal Control"       
# [3] "schizophrenia" 

levels(rse_gene_SRP035524$disease) <- c("bp", "control", "schizophrenia")
rse_gene_SRP035524$disease <- relevel(rse_gene_SRP035524$disease, "control")

levels(rse_gene_SRP035524$body_site)
# [1] "anterior cingulate" "frontal cortex"

levels(rse_gene_SRP035524$body_site) <- c("anterior_cingulate", "frontal_cortex")
#head(colData(rse_er_scaled))

assayNames(rse_er_scaled_SRP035524_chr1)

# Set the rownames for each rse_scaled to be equivalent to the chromosome #, genomic coordinates of the novel region, and the strand
rownames(rse_er_scaled_SRP035524_chr1) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr1'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))])), length.out= nrow(rse_er_scaled_SRP035524_chr1)))
rownames(rse_er_scaled_SRP035524_chr10) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr10'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))])), length.out= nrow(rse_er_scaled_SRP035524_chr10)))
rownames(rse_er_scaled_SRP035524_chr11) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr11'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))])), length.out= nrow(rse_er_scaled_SRP035524_chr11)))
rownames(rse_er_scaled_SRP035524_chr12) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr12'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))])), length.out= nrow(rse_er_scaled_SRP035524_chr12)))
rownames(rse_er_scaled_SRP035524_chr13) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr13'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))])), length.out= nrow(rse_er_scaled_SRP035524_chr13)))
rownames(rse_er_scaled_SRP035524_chr14) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr14'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))])), length.out= nrow(rse_er_scaled_SRP035524_chr14)))
rownames(rse_er_scaled_SRP035524_chr15) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr15'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))])), length.out= nrow(rse_er_scaled_SRP035524_chr15)))
rownames(rse_er_scaled_SRP035524_chr16) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr16'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))])), length.out= nrow(rse_er_scaled_SRP035524_chr16)))

# Coverage matrices for each chromosome
coverage_SRP035524_chr1 <- assay(rse_er_scaled_SRP035524_chr1)
coverage_SRP035524_chr10 <- assay(rse_er_scaled_SRP035524_chr10)
coverage_SRP035524_chr11 <- assay(rse_er_scaled_SRP035524_chr11)
coverage_SRP035524_chr12 <- assay(rse_er_scaled_SRP035524_chr12)
coverage_SRP035524_chr13 <- assay(rse_er_scaled_SRP035524_chr13)
coverage_SRP035524_chr14 <- assay(rse_er_scaled_SRP035524_chr14)
coverage_SRP035524_chr15 <- assay(rse_er_scaled_SRP035524_chr15)
coverage_SRP035524_chr16 <- assay(rse_er_scaled_SRP035524_chr16)

# Put all of the coverage matrices together using rbind()
coverage_SRP035524 <- rbind(coverage_SRP035524_chr1, coverage_SRP035524_chr10, coverage_SRP035524_chr11, coverage_SRP035524_chr12, coverage_SRP035524_chr13,
      coverage_SRP035524_chr14, coverage_SRP035524_chr15, coverage_SRP035524_chr16)

#####################################################

# Step Three: Exploratory Analysis and Differential Expression Analysis

# set up DESeqDataSetfromMatrix
# countdata object to store the coverage matrix that includes all of the chromosomes
countdata <- coverage_SRP035524
# coldata object to store the colData for rse_gene_SRP037725
coldata <- colData(rse_gene_SRP035524)
# countdata: a table with the fragment counts
# coldata: a table with information about the samples

# design for SRP035524 based on disease and body_site
# every other study will be just based on disease since there is only one tissue site
ddsMat_SRP035524 <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ body_site + disease)
nrow(ddsMat_SRP035524)
# [1] 40
# keep only the rows where the counts are greater than 1
keep <- rowSums(counts(ddsMat_SRP035524)) > 1
table(keep)
#FALSE  TRUE 
#  2    38 

# observe which rows contain no counts for any of the samples
which(keep == FALSE)
# chr10:124606586-124606647_+ chr12:128980477-128980532_+ 
ddsMat_SRP035524 <- ddsMat_SRP035524[keep,]
nrow(ddsMat_SRP035524)
# [1] 38

# rlog transformation to do exploratory analysis
rld_SRP035524 <- rlog(ddsMat_SRP035524, blind = FALSE)
head(assay(rld_SRP035524), 6)

# observe the distances between samples visualized through a heatmap
sampleDists <- dist(t(assay(rld_SRP035524)))
sampleDists

# generate the heatmap
library(pheatmap)
library(RColorBrewer)

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld_SRP035524$disease, rld_SRP035524$body_site, sep = " - " )
colnames(sampleDistMatrix) <- paste( rld_SRP035524$disease, rld_SRP035524$body_site, sep = " - " )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# principal component analysis for additional exploratory analysis
pcaData <- DESeq2::plotPCA(rld_SRP035524, intgroup = c( "disease", "body_site"), returnData= TRUE)
pcaData

# generate the PCA Plot
percentVar <- round(100 * attr(pcaData, "percentVar"))
# [1] 31 11

library(ggplot2)
ggplot(pcaData, aes(x = PC1, y = PC2, color = disease, shape = body_site)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with Rlog data")

# set-up the DE analysis using the 'DESeq()' function
ddsMat_SRP035524 <- DESeq(ddsMat_SRP035524)

res <- results(ddsMat_SRP035524)
res

summary(res)

results(ddsMat_SRP035524, contrast = c("disease", "control", "bp"))
results(ddsMat_SRP035524, contrast = c("disease", "control", "schizophrenia"))
results(ddsMat_SRP035524, contrast = c("disease", "bp", "schizophrenia"))

# observe the number of regions whose alpha was 0.05 or less (and therefore significant in terms of differential expresssion)
res.05 <- results(ddsMat_SRP035524, alpha = 0.05)
table(res.05$padj < 0.05)
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))

# observe the number of regions whose alpha was 0.10 of less
res.10 <- results(ddsMat_SRP035524, alpha = 0.10)
table(res.10$padj < 0.10)
sum(res$padj < 0.1, na.rm=TRUE)

# regions with up-regulation
resSig <- subset(res, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])

# regions with down-regulation
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

# If there were significant differentially expressed regions, then a counts plot could be generated.  
# But for these regions and this particular study, there were not any significant regions identified.