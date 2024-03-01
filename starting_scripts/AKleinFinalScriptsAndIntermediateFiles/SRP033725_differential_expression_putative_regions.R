######################################################

# Differential Expression of the putative exonic regions using recount2 database - SRP033725
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

# Convert the GenomicRanges object to a dataframe
# Make the rownames of the dataframe reflect the region (chr#: coordinates_strand) to be used when creating coverage matrix

#####################################################

# Step Two: Calculating coverage of the desired regions, scaling the counts, adding necessary colData, and combining
# all of the coverage matrices

# Load in the recount project to be used in the differential expression analysis

if(!file.exists(file.path("SRP033725", "rse_gene_SRP033725.Rdata"))) {
  download_study("SRP033725", type = "rse-gene")
}

# Check that the file was downloaded using the file.exists() function
file.exists(file.path("SRP033725", "rse_gene_SRP033725.Rdata"))

load(file.path("SRP033725", "rse_gene_SRP033725.Rdata"))

# rename the rse_gene object to be named specifically for the study
rse_gene_SRP033725 <- rse_gene

# Calculate the coverage over the putative regions (using coverage_matrix() function in recount)
rse_er_SRP033725_chr1 <- coverage_matrix("SRP033725", "chr1", granges_putative_exons,
                                         verboseLoad = FALSE, scale = FALSE)
rse_er_SRP033725_chr10 <- coverage_matrix("SRP033725", "chr10", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP033725_chr11 <- coverage_matrix("SRP033725", "chr11", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP033725_chr12 <- coverage_matrix("SRP033725", "chr12", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP033725_chr13 <- coverage_matrix("SRP033725", "chr13", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP033725_chr14 <- coverage_matrix("SRP033725", "chr14", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP033725_chr15 <- coverage_matrix("SRP033725", "chr15", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)
rse_er_SRP033725_chr16 <- coverage_matrix("SRP033725", "chr16", granges_putative_exons,
                                          verboseLoad = FALSE, scale = FALSE)

# rename the columns names for chromosome 16- for some reason they are not named with the same schema as the others
colnames(rse_er_SRP033725_chr16) <- rse_er_SRP033725_chr16$run

# scale the counts for each of the coverage matrices for each chromosome
rse_er_scaled_SRP033725_chr1 <- scale_counts(rse_er_SRP033725_chr1)
rse_er_scaled_SRP033725_chr10 <- scale_counts(rse_er_SRP033725_chr10)
rse_er_scaled_SRP033725_chr11 <- scale_counts(rse_er_SRP033725_chr11)
rse_er_scaled_SRP033725_chr12 <- scale_counts(rse_er_SRP033725_chr12)
rse_er_scaled_SRP033725_chr13 <- scale_counts(rse_er_SRP033725_chr13)
rse_er_scaled_SRP033725_chr14 <- scale_counts(rse_er_SRP033725_chr14)
rse_er_scaled_SRP033725_chr15 <- scale_counts(rse_er_SRP033725_chr15)
rse_er_scaled_SRP033725_chr16 <- scale_counts(rse_er_SRP033725_chr16)

# Use the expanded metadata we built for the gene model- here we are simply updating the colData for the scaled
# counts with the colData for the original rse_gene(s)
colData(rse_er_scaled_SRP033725_chr1) <- colData(rse_gene_SRP033725)
colData(rse_er_scaled_SRP033725_chr10) <- colData(rse_gene_SRP033725)
colData(rse_er_scaled_SRP033725_chr11) <- colData(rse_gene_SRP033725)
colData(rse_er_scaled_SRP033725_chr12) <- colData(rse_gene_SRP033725)
colData(rse_er_scaled_SRP033725_chr13) <- colData(rse_gene_SRP033725)
colData(rse_er_scaled_SRP033725_chr14) <- colData(rse_gene_SRP033725)
colData(rse_er_scaled_SRP033725_chr15) <- colData(rse_gene_SRP033725)
colData(rse_er_scaled_SRP033725_chr16) <- colData(rse_gene_SRP033725)

# Save the information from https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP033725
file.exists(file.path("SRP033725", "SraRunTable_SRP033725.txt"))

# Read the table
sra_SRP033725 <- read.table(file.path("SRP033725", "SraRunTable_SRP033725.txt"),
                            header = TRUE, sep = ",")

# want to use "disease" and "body_site" since we are interested in differential expression in disease state and tissue site
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
sra_SRP033725 <- sra_SRP033725[match(colData(rse_er_scaled_SRP033725_chr1)$run, sra_SRP033725$Run), ]

# Double check the order to make sure that the rse_er run column is equivalent to the sra acc column
identical(colData(rse_er_scaled_SRP033725_chr1)$run, as.character(sra_SRP033725$Run))
identical(colData(rse_er_scaled_SRP033725_chr10)$run, as.character(sra_SRP033725$Run))
identical(colData(rse_er_scaled_SRP033725_chr11)$run, as.character(sra_SRP033725$Run))
identical(colData(rse_er_scaled_SRP033725_chr12)$run, as.character(sra_SRP033725$Run))
identical(colData(rse_er_scaled_SRP033725_chr13)$run, as.character(sra_SRP033725$Run))
identical(colData(rse_er_scaled_SRP033725_chr14)$run, as.character(sra_SRP033725$Run))
identical(colData(rse_er_scaled_SRP033725_chr15)$run, as.character(sra_SRP033725$Run))
identical(colData(rse_er_scaled_SRP033725_chr16)$run, as.character(sra_SRP033725$Run))

# Append the variables of interest to the colData for the rse_er (specifically rse_er_scaled)
colData(rse_gene_SRP033725) <- cbind(colData(rse_gene_SRP033725), sra_SRP033725[, sra_vars_SRP033725])
colData(rse_er_scaled_SRP033725_chr1) <- cbind(colData(rse_er_SRP033725_chr1), sra_SRP033725[, sra_vars_SRP033725])
colData(rse_er_scaled_SRP033725_chr10) <- cbind(colData(rse_er_SRP033725_chr10), sra_SRP033725[, sra_vars_SRP033725])
colData(rse_er_scaled_SRP033725_chr11) <- cbind(colData(rse_er_SRP033725_chr11), sra_SRP033725[, sra_vars_SRP033725])
colData(rse_er_scaled_SRP033725_chr12) <- cbind(colData(rse_er_SRP033725_chr12), sra_SRP033725[, sra_vars_SRP033725])
colData(rse_er_scaled_SRP033725_chr13) <- cbind(colData(rse_er_SRP033725_chr13), sra_SRP033725[, sra_vars_SRP033725])
colData(rse_er_scaled_SRP033725_chr14) <- cbind(colData(rse_er_SRP033725_chr14), sra_SRP033725[, sra_vars_SRP033725])
colData(rse_er_scaled_SRP033725_chr15) <- cbind(colData(rse_er_SRP033725_chr15), sra_SRP033725[, sra_vars_SRP033725])
colData(rse_er_scaled_SRP033725_chr16) <- cbind(colData(rse_er_SRP033725_chr16), sra_SRP033725[, sra_vars_SRP033725])

# check levels and set reference level (control) using relevel() function
levels(rse_gene_SRP033725$disease_state)
# [1] "BD"      "Control" 

levels(rse_gene_SRP033725$disease_state) <- c("bp", "control")
rse_gene_SRP033725$disease_state <- relevel(rse_gene_SRP033725$disease_state, "control")
levels(rse_gene_SRP033725$disease_state)

levels(rse_gene_SRP033725$tissue)
# [1] "dorsolateral prefrontal cortex"

levels(rse_gene_SRP033725$tissue) <- ("dorsolateral_prefrontal_cortex")

assayNames(rse_er_scaled_SRP033725_chr1)

# Set the rownames for each rse_scaled to be equivalent to the chromosome #, genomic coordinates of the novel region, and the strand
rownames(rse_er_scaled_SRP033725_chr1) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr1'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr1'))])), length.out= nrow(rse_er_scaled_SRP033725_chr1)))
rownames(rse_er_scaled_SRP033725_chr10) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr10'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr10'))])), length.out= nrow(rse_er_scaled_SRP033725_chr10)))
rownames(rse_er_scaled_SRP033725_chr11) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr11'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr11'))])), length.out= nrow(rse_er_scaled_SRP033725_chr11)))
rownames(rse_er_scaled_SRP033725_chr12) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr12'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr12'))])), length.out= nrow(rse_er_scaled_SRP033725_chr12)))
rownames(rse_er_scaled_SRP033725_chr13) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr13'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr13'))])), length.out= nrow(rse_er_scaled_SRP033725_chr13)))
rownames(rse_er_scaled_SRP033725_chr14) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr14'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr14'))])), length.out= nrow(rse_er_scaled_SRP033725_chr14)))
rownames(rse_er_scaled_SRP033725_chr15) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr15'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr15'))])), length.out= nrow(rse_er_scaled_SRP033725_chr15)))
rownames(rse_er_scaled_SRP033725_chr16) <- c(rep(paste0(seqnames[(which(seqnames(granges_putative_exons) =='chr16'))],":",  start(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))]), "-", end(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))]), "_", strand(granges_putative_exons[(which(seqnames(granges_putative_exons) =='chr16'))])), length.out= nrow(rse_er_scaled_SRP033725_chr16)))

# Coverage matrices for each chromosome
coverage_SRP033725_chr1 <- assay(rse_er_scaled_SRP033725_chr1)
coverage_SRP033725_chr10 <- assay(rse_er_scaled_SRP033725_chr10)
coverage_SRP033725_chr11 <- assay(rse_er_scaled_SRP033725_chr11)
coverage_SRP033725_chr12 <- assay(rse_er_scaled_SRP033725_chr12)
coverage_SRP033725_chr13 <- assay(rse_er_scaled_SRP033725_chr13)
coverage_SRP033725_chr14 <- assay(rse_er_scaled_SRP033725_chr14)
coverage_SRP033725_chr15 <- assay(rse_er_scaled_SRP033725_chr15)
coverage_SRP033725_chr16 <- assay(rse_er_scaled_SRP033725_chr16)

# Put all of the coverage matrices together using rbind()
coverage_SRP033725 <- rbind(coverage_SRP033725_chr1, coverage_SRP033725_chr10, coverage_SRP033725_chr11, coverage_SRP033725_chr12, coverage_SRP033725_chr13,
                            coverage_SRP033725_chr14, coverage_SRP033725_chr15, coverage_SRP033725_chr16)

#####################################################

# Step Three: Exploratory Analysis and Differential Expression Analysis

# set up DESeqDataSetfromMatrix 
# countdata object to store the coverage matrix that includes all of the chromosomes
countdata <- coverage_SRP033725
# coldata object to store the colData for rse_gene_SRP037725
coldata <- colData(rse_gene_SRP033725)
# countdata: a table with the fragment counts
# coldata: a table with information about the samples

# design for SRP033725 based on disease and body_site
# for this particular study, the DE will be based only on disease_state since there is only one tissue site 
ddsMat_SRP033725 <- DESeqDataSetFromMatrix(countData = countdata,
                                           colData = coldata,
                                           design = ~ disease_state)
nrow(ddsMat_SRP033725)
# [1] 40
# keep only the rows where the counts are greater than 1
keep <- rowSums(counts(ddsMat_SRP033725)) > 1
table(keep)
# TRUE 
#  40 

# observe which rows contain no counts for any of the samples
# which(keep == FALSE)
ddsMat_SRP033725 <- ddsMat_SRP033725[keep,]
nrow(ddsMat_SRP033725)
# [1] 40

# rlog transformation to do exploratory analysis
rld_SRP033725 <- rlog(ddsMat_SRP033725, blind = FALSE)
head(assay(rld_SRP033725), 6)

# observe the distances between samples visualized through a heatmap
sampleDists <- dist(t(assay(rld_SRP033725)))
sampleDists

# generate the heatmap
library(pheatmap)
library(RColorBrewer)

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld_SRP033725$disease_state
colnames(sampleDistMatrix) <- rld_SRP033725$disease_state
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
# [1] 25  9

library(ggplot2)
ggplot(pcaData_SRP033725, aes(x = PC1, y = PC2, color = disease_state, shape = disease_state)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar_SRP033725[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_SRP033725[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with Rlog data")

# set-up the DE analysis using the 'DESeq()' function
ddsMat_SRP033725 <- DESeq(ddsMat_SRP033725)

res_SRP033725 <- results(ddsMat_SRP033725)
res_SRP033725

summary(res_SRP033725)
# help with LFC explanation > 0 etc
# out of 40 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1, 2.5%
# LFC < 0 (down)     : 3, 7.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)

results(ddsMat_SRP033725, contrast = c("disease_state", "control", "bp"))

# observe the number of regions whose alpha was 0.05 or less (and therefore significant in terms of differential expresssion)
res.05 <- results(ddsMat_SRP033725, alpha = 0.05)
table(res.05$padj < 0.05)
# FALSE  TRUE 
# 36     4
sum(res_SRP033725$pvalue < 0.05, na.rm=TRUE)
sum(res_SRP033725$padj < 0.05, na.rm=TRUE)
sum(!is.na(res_SRP033725$pvalue))

# observe the number of regions whose alpha was 0.10 of less
res.10 <- results(ddsMat_SRP033725, alpha = 0.10)
table(res.10$padj < 0.10)
sum(res_SRP033725$padj < 0.1, na.rm=TRUE)

# regions with up-regulation
resSig <- subset(res_SRP033725, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])

# regions with down-regulation
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
# log2 fold change (MLE): disease state bp vs control 
# Wald test p-value: disease state bp vs control 
# DataFrame with 4 rows and 6 columns
#                                 baseMean
#                               <numeric>
# chr15:75351300-75351376_+ 5.94100193473498
# chr15:75007587-75007728_+ 31.7388463003704
# chr15:75007660-75007728_+ 21.6305109421663
# chr11:63746946-63747024_+ 6.59300730494249
#                               log2FoldChange
#                                 <numeric>
# chr15:75351300-75351376_+  0.743854728015551
# chr15:75007587-75007728_+ -0.356402532110837
# chr15:75007660-75007728_+ -0.386254486197221
# chr11:63746946-63747024_+ -0.597331937653378
#                                   lfcSE
#                                 <numeric>
# chr15:75351300-75351376_+  0.172502316064305
# chr15:75007587-75007728_+ 0.0991368368300148
# chr15:75007660-75007728_+  0.113801780712846
# chr11:63746946-63747024_+  0.157785845235945
#                                   stat
#                                 <numeric>
# chr15:75351300-75351376_+   4.3121434250092
# chr15:75007587-75007728_+ -3.59505652497207
# chr15:75007660-75007728_+ -3.39409878982341
# chr11:63746946-63747024_+ -3.78571307686161
                                  # pvalue
#                                 <numeric>
# chr15:75351300-75351376_+ 1.61679494899326e-05
# chr15:75007587-75007728_+ 0.000324321152991125
# chr15:75007660-75007728_+ 0.000688548492843674
# chr11:63746946-63747024_+ 0.000153268439148034
#                                      padj
#                                    <numeric>
# chr15:75351300-75351376_+ 0.000646717979597306
# chr15:75007587-75007728_+  0.00432428203988167
# chr15:75007660-75007728_+  0.00688548492843674
# chr11:63746946-63747024_+  0.00306536878296067

# Counts Plot
# indexed for each of the regions according to which row it was indexed as in rownames(res_SRP037725)
region_chr15_75351300 <- rownames(res_SRP033725)[which.min(res_SRP033725$padj)]

region_chr15_75007587 <- rownames(res_SRP033725[37,])

region_chr15_75007660 <- rownames(res_SRP033725[38,])

region_chr11_63746946 <- rownames(res_SRP033725[14,])

# Load the ggbeeswarm package in order to use the plotCounts function to plot
# each of the differentially expressed regions
# "gene" argument needs to be changed each time when changing the region being plotted
library("ggbeeswarm")
geneCounts <- plotCounts(ddsMat_SRP033725, gene = region_chr15_75351300, intgroup = c("disease_state"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = disease_state, y = count, color = disease_state)) +
  scale_y_log10() +  geom_beeswarm(cex = 3) + ggtitle("Counts Plot for chr15:75351300-75351376_+")
