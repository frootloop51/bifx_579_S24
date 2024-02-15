######################################################

# Differential Expression of the putative exonic regions using recount2 database - SRP043684
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

if(!file.exists(file.path("SRP043684", "rse_gene_SRP043684.Rdata"))) {
  download_study("SRP043684", type = "rse-gene")
}

# Check that the file was downloaded using the file.exists() function
file.exists(file.path("SRP043684", "rse_gene_SRP043684.Rdata"))

load(file.path("SRP043684", "rse_gene_SRP043684.Rdata"))

rse_gene_SRP043684 <- rse_gene

# Calculate the coverage over the putative regions (using coverage_matrix() function in recount)
rse_er_SRP043684_chr1 <- coverage_matrix("SRP043684", "chr1", granges_putative_exons,
                                         verboseLoad = FALSE, scale = FALSE)
rse_er_SRP043684_chr10 <- coverage_matrix("SRP043684", "chr10", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP043684_chr11 <- coverage_matrix("SRP043684", "chr11", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP043684_chr12 <- coverage_matrix("SRP043684", "chr12", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP043684_chr13 <- coverage_matrix("SRP043684", "chr13", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP043684_chr14 <- coverage_matrix("SRP043684", "chr14", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP043684_chr15 <- coverage_matrix("SRP043684", "chr15", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP043684_chr16 <- coverage_matrix("SRP043684", "chr16", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)

# rename the columns names for chromosome 16- for some reason they are not named with the same schema as the others
colnames(rse_er_SRP043684_chr16) <- rse_er_SRP043684_chr16$run


# scale the counts for each of the coverage matrices for each chromosome
rse_er_scaled_SRP043684_chr1 <- scale_counts(rse_er_SRP043684_chr1)
rse_er_scaled_SRP043684_chr10 <- scale_counts(rse_er_SRP043684_chr10)
rse_er_scaled_SRP043684_chr11 <- scale_counts(rse_er_SRP043684_chr11)
rse_er_scaled_SRP043684_chr12 <- scale_counts(rse_er_SRP043684_chr12)
rse_er_scaled_SRP043684_chr13 <- scale_counts(rse_er_SRP043684_chr13)
rse_er_scaled_SRP043684_chr14 <- scale_counts(rse_er_SRP043684_chr14)
rse_er_scaled_SRP043684_chr15 <- scale_counts(rse_er_SRP043684_chr15)
rse_er_scaled_SRP043684_chr16 <- scale_counts(rse_er_SRP043684_chr16)

# Use the expanded metadata we built for the gene model- here we are simply updating the colData for the scaled
# counts with the colData for the original rse_gene(s)
colData(rse_er_scaled_SRP043684_chr1) <- colData(rse_gene_SRP043684)
colData(rse_er_scaled_SRP043684_chr10) <- colData(rse_gene_SRP043684)
colData(rse_er_scaled_SRP043684_chr11) <- colData(rse_gene_SRP043684)
colData(rse_er_scaled_SRP043684_chr12) <- colData(rse_gene_SRP043684)
colData(rse_er_scaled_SRP043684_chr13) <- colData(rse_gene_SRP043684)
colData(rse_er_scaled_SRP043684_chr14) <- colData(rse_gene_SRP043684)
colData(rse_er_scaled_SRP043684_chr15) <- colData(rse_gene_SRP043684)
colData(rse_er_scaled_SRP043684_chr16) <- colData(rse_gene_SRP043684)

# Save the information from https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP043684
file.exists(file.path("SRP043684", "SraRunTable_SRP043684.txt"))

## Read the table
sra_SRP043684 <- read.table(file.path("SRP043684", "SraRunTable_SRP043684.txt"),
                            header = TRUE, sep = ",")

# want to use "bd_patient", "Cell_type", and "lithium_treated"
colnames(sra_SRP043684)
# [1] "Run"                "Assay.Type"        
# [3] "AvgSpotLen"         "Bases"             
# [5] "BioProject"         "BioSample"         
# [7] "Bytes"              "Cell_type"         
# [9] "Center.Name"        "Consent"           
# [11] "DATASTORE.filetype" "DATASTORE.provider"
# [13] "DATASTORE.region"   "Experiment"        
# [15] "GEO_Accession"      "Instrument"        
# [17] "LibraryLayout"      "LibrarySelection"  
# [19] "LibrarySource"      "MBases"            
# [21] "MBytes"             "Organism"          
# [23] "Platform"           "ReleaseDate"       
# [25] "sample_acc"         "Sample.Name"       
# [27] "source_name"        "SRA.Study"         
# [29] "bd_patient"         "lithium_treated"   
# [31] "prox1gfp"   

# save the variables of interest to add to one variable
sra_vars_SRP043684 <- c("bd_patient", "Cell_type", "lithium_treated")
sra_vars_SRP043684

# Re-organize the SRA table based on the SRA Run IDs we have
sra_SRP043684 <- sra_SRP043684[match(colData(rse_er_scaled_SRP043684_chr1)$run, sra_SRP043684$Run), ]

# Double check the order to make sure that the rse_er run column is equivalent to the sra acc column
identical(colData(rse_er_scaled_SRP043684_chr1)$run, as.character(sra_SRP043684$Run))
identical(colData(rse_er_scaled_SRP043684_chr10)$run, as.character(sra_SRP043684$Run))
identical(colData(rse_er_scaled_SRP043684_chr11)$run, as.character(sra_SRP043684$Run))
identical(colData(rse_er_scaled_SRP043684_chr12)$run, as.character(sra_SRP043684$Run))
identical(colData(rse_er_scaled_SRP043684_chr13)$run, as.character(sra_SRP043684$Run))
identical(colData(rse_er_scaled_SRP043684_chr14)$run, as.character(sra_SRP043684$Run))
identical(colData(rse_er_scaled_SRP043684_chr15)$run, as.character(sra_SRP043684$Run))
identical(colData(rse_er_scaled_SRP043684_chr16)$run, as.character(sra_SRP043684$Run))

# Append the variables of interest to the colData for the rse_er (specifically rse_er_scaled)
colData(rse_gene_SRP043684) <- cbind(colData(rse_gene_SRP043684), sra_SRP043684[, sra_vars_SRP043684])
colData(rse_er_scaled_SRP043684_chr1) <- cbind(colData(rse_er_SRP043684_chr1), sra_SRP043684[, sra_vars_SRP043684])
colData(rse_er_scaled_SRP043684_chr10) <- cbind(colData(rse_er_SRP043684_chr10), sra_SRP043684[, sra_vars_SRP043684])
colData(rse_er_scaled_SRP043684_chr11) <- cbind(colData(rse_er_SRP043684_chr11), sra_SRP043684[, sra_vars_SRP043684])
colData(rse_er_scaled_SRP043684_chr12) <- cbind(colData(rse_er_SRP043684_chr12), sra_SRP043684[, sra_vars_SRP043684])
colData(rse_er_scaled_SRP043684_chr13) <- cbind(colData(rse_er_SRP043684_chr13), sra_SRP043684[, sra_vars_SRP043684])
colData(rse_er_scaled_SRP043684_chr14) <- cbind(colData(rse_er_SRP043684_chr14), sra_SRP043684[, sra_vars_SRP043684])
colData(rse_er_scaled_SRP043684_chr15) <- cbind(colData(rse_er_SRP043684_chr15), sra_SRP043684[, sra_vars_SRP043684])
colData(rse_er_scaled_SRP043684_chr16) <- cbind(colData(rse_er_SRP043684_chr16), sra_SRP043684[, sra_vars_SRP043684])

# check levels and set reference level (control) using relevel() function
levels(rse_gene_SRP043684$bd_patient)
# [1] ""    "no"  "yes" 

# START HERE
levels(rse_gene_SRP043684$bd_patient) <- c("", "control", "bp")
rse_gene_SRP043684$bd_patient <- relevel(rse_gene_SRP043684$bd_patient, "control")

levels(rse_gene_SRP043684$Cell_type)
# [1] "Hippocampal dentate gyrus (DG) granule neurons"
# [2] "iPSC-derived neurons" 

levels(rse_gene_SRP043684$source_name)
# [1] "iPSC-derived hippocampal dentate gyrus granule neurons"
# [2] "iPSC-derived neurons" 

assayNames(rse_er_scaled_SRP043684_chr1)

# Set the rownames for each rse_scaled to be equivalent to the chromosome #, genomic coordinates of the novel region, and the strand
rownames(rse_er_scaled_SRP043684_chr1) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr1'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))])), length.out= nrow(rse_er_scaled_SRP043684_chr1)))
rownames(rse_er_scaled_SRP043684_chr10) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr10'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))])), length.out= nrow(rse_er_scaled_SRP043684_chr10)))
rownames(rse_er_scaled_SRP043684_chr11) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr11'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))])), length.out= nrow(rse_er_scaled_SRP043684_chr11)))
rownames(rse_er_scaled_SRP043684_chr12) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr12'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))])), length.out= nrow(rse_er_scaled_SRP043684_chr12)))
rownames(rse_er_scaled_SRP043684_chr13) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr13'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))])), length.out= nrow(rse_er_scaled_SRP043684_chr13)))
rownames(rse_er_scaled_SRP043684_chr14) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr14'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))])), length.out= nrow(rse_er_scaled_SRP043684_chr14)))
rownames(rse_er_scaled_SRP043684_chr15) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr15'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))])), length.out= nrow(rse_er_scaled_SRP043684_chr15)))
rownames(rse_er_scaled_SRP043684_chr16) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr16'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))])), length.out= nrow(rse_er_scaled_SRP043684_chr16)))

# Coverage matrices for each chromosome
coverage_SRP043684_chr1 <- assay(rse_er_scaled_SRP043684_chr1)
coverage_SRP043684_chr10 <- assay(rse_er_scaled_SRP043684_chr10)
coverage_SRP043684_chr11 <- assay(rse_er_scaled_SRP043684_chr11)
coverage_SRP043684_chr12 <- assay(rse_er_scaled_SRP043684_chr12)
coverage_SRP043684_chr13 <- assay(rse_er_scaled_SRP043684_chr13)
coverage_SRP043684_chr14 <- assay(rse_er_scaled_SRP043684_chr14)
coverage_SRP043684_chr15 <- assay(rse_er_scaled_SRP043684_chr15)
coverage_SRP043684_chr16 <- assay(rse_er_scaled_SRP043684_chr16)

# Put all of the coverage matrices together using rbind()
coverage_SRP043684 <- rbind(coverage_SRP043684_chr1, coverage_SRP043684_chr10, coverage_SRP043684_chr11, coverage_SRP043684_chr12, coverage_SRP043684_chr13,
                            coverage_SRP043684_chr14, coverage_SRP043684_chr15, coverage_SRP043684_chr16)

#####################################################

# Step Three: Exploratory Analysis and Differential Expression Analysis

# set up DESeqDataSetfromMatrix 
# countdata object to store the coverage matrix that includes all of the chromosomes
countdata <- coverage_SRP043684
# coldata object to store the colData for rse_gene_SRP037725
coldata <- colData(rse_gene_SRP043684)
# countdata: a table with the fragment counts
# coldata: a table with information about the samples

# does not allow me to add lithium_treated variable, so just used bd_patient in the design
# could go back and test for batch effects to somehow include the lithium_treated variable
# Error in checkFullRank(modelMatrix) : 
#  the model matrix is not full rank, so the model cannot be fit as specified.
# One or more variables or interaction terms in the design formula are linear
# combinations of the others and must be removed.
ddsMat_SRP043684 <- DESeqDataSetFromMatrix(countData = countdata,
                                           colData = coldata,
                                           design = ~ bd_patient)
nrow(ddsMat_SRP043684)
# [1] 40
# keep only the rows where the counts are greater than 1
keep <- rowSums(counts(ddsMat_SRP043684)) > 1
table(keep)
#FALSE  TRUE 
#  2    38 

# observe which rows contain no counts for any of the samples
which(keep == FALSE)
# chr10:124606586-124606647_+ chr12:128980477-128980532_+
ddsMat_SRP043684 <- ddsMat_SRP043684[keep,]
nrow(ddsMat_SRP043684)
# [1] 38

# rlog transformation to do exploratory analysis
rld_SRP043684 <- rlog(ddsMat_SRP043684, blind = FALSE)
head(assay(rld_SRP043684), 6)

# observe the distances between samples visualized through a heatmap
sampleDists <- dist(t(assay(rld_SRP043684)))
sampleDists

# generate the heatmap
library(pheatmap)
library(RColorBrewer)

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld_SRP043684$bd_patient
colnames(sampleDistMatrix) <- rld_SRP043684$bd_patient
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# principal component analysis for additional exploratory analysis
pcaData <- DESeq2::plotPCA(rld_SRP043684, intgroup = c( "bd_patient"), returnData= TRUE)
pcaData

# generate the PCA Plot
percentVar <- round(100 * attr(pcaData, "percentVar"))
# [1] 30 14

library(ggplot2)
ggplot(pcaData, aes(x = PC1, y = PC2, color = bd_patient, shape = bd_patient)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with Rlog data")

# set-up the DE analysis using the 'DESeq()' function
ddsMat_SRP043684 <- DESeq(ddsMat_SRP043684)

res_SRP043684 <- results(ddsMat_SRP043684)
res_SRP043684

summary(res_SRP043684)

results(ddsMat_SRP043684, contrast = c("bd_patient", "control", "bp"))

# observe the number of regions whose alpha was 0.05 or less (and therefore significant in terms of differential expresssion)
res.05 <- results(ddsMat_SRP043684, alpha = 0.05)
table(res.05$padj < 0.05)
sum(res_SRP043684$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res_SRP043684$pvalue))

# observe the number of regions whose alpha was 0.10 of less
res.10 <- results(ddsMat_SRP043684, alpha = 0.10)
table(res.10$padj < 0.10)
sum(res_SRP043684$padj < 0.1, na.rm=TRUE)

# regions with up-regulation
resSig <- subset(res_SRP043684, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])

# regions with down-regulation
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

# If there were significant differentially expressed regions, then a counts plot could be generated.  
# But for these regions and this particular study, there were not any significant regions identified.