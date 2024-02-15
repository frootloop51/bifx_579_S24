########################################################################
#
# RNA-seq workflow example
# BIFX 572: Computational Genomics and Proteomics
# Miranda M. Darby, Ph.D.
# Date: 11/5/2019
#
# Based on workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# Title of workflow: RNA-seq workflow: gene-level exploratory analysis and differential expression
#
#######################################################################
#
#   #This workflow is an example of an end-to-end gene-level RNA-seq differential expression.
#   #The data used for this workflow can be found in the airway package.  This contains data
#   #in which airway smooth muscle cells were treated or untreated with the synthetic
#   #glucocorticoid steroid dexamthasone- which is supposed to reduce inflammation.
#   #Four cell lines of airway smooth muscle cells were treated, and there is a treated
#   #and an untreated sample for each of the four cell lines.
#
#####################################################################
#
# Preparing quantification input to DESeq2
#
# input data is expected to be in the form of a matrix of un-normalized counts
# gene is denoted as i and sample as j
# each value in the i-th row and j-th column is the number of reads for a gene in a particular sample
# never provide counts that were pre-normalized as model works best when applied to unnormalized counts
#
# 2.1 Transcript quantification and tximport/tximeta
#
# workflow shows importing transcript-level quantification data that is then
# aggregated to the gene-level with tximport/tximeta
# Salmon or RSEM- mapping to reference transcripts and outputting estimated counts per transcript
# previous workflow aligned reads to the genome and then counted them

# tximeta is an extension of tximport by allowiing for addition of annotation metadata for
# commonly used transcriptomes
# tximeta produces SummarizedExperiment object that can be loaded into DESeq2

# 2.2 Quantifying with Salmon
# Salmon used too quantify files in the workflow
# Can be run on a cluster using Snakemake workflow management system
# The following Snakemake file was used to quantify each of the sample (8)
# I created a separate Snakemake file with the name "Snakemake_RNA-seq.R"
#
# Samples downloaded from the SRA (SRR is the run identifier)
# DATASETS = ["SRR1039508",
#             "SRR1039509",
#             "SRR1039512",
#             "SRR1039513",
#             "SRR1039516",
#             "SRR1039517",
#             "SRR1039520",
#             "SRR1039521"]
# 
# provides the path to Salmon
# SALMON = "/path/to/salmon_0.14.1/bin/salmon"
# 
# rule all:
#   input: expand("quants/{dataset}/quant.sf", dataset=DATASETS)
# 
# rule salmon_quant:
#   input:
#       r1 = "fastq/{sample}_1.fastq.gz", # uses fasta.gz files for each sample as input
#       r2 = "fastq/{sample}_2.fastq.gz",
#       index = "/path/to/gencode.v29_salmon_0.14.1" # where the index was created for the Snakemake file to use
#   output:
#     "quants/{sample}/quant.sf" # output files to be generated with quantification data
#   params:
#       dir = "quants/{sample}"
# following command runs Salmon 
# quantifies using specific index, using 6 threads, validateMappings setting
# GC bias correction, and writing out 20 Gibbs samples
# - o denotes output directory, and -1 and -2 denote read files (input files)
#   shell:
#       "{SALMON} quant -i {input.index} -l A -p 6 --validateMappings \
#       --gcBias --numGibbsSamples 20 -o {params.dir} \
#       -1 {input.r1} -2 {input.r2}"
#
# Creation of the Salmon index
salmon index -t trnascripts.fa.gz -i name_of_index
#
# Snakemake file created executed on cluster 
snakemake -j 4 --latency-wait 30 --cluster "sbatch -N -n 6"
#
# 2.3 Reading in data with tximeta
# 
# airway package has two quantification directories output by Salmon
# quant.sf files have been gzipped to make the data package smaller
# says to use quant.sf on out own machine?

# install the airway package if not previously installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("airway")

# load the airway data package that contains the airway smooth muscle cell data
library(airway)

# Use of system.file function to find out where computer files from a package have been installed
# asks for the path to extdata directory of the airway package
dir <- system.file("extdata", packaage="airway", mustWorkk=TRUE)

# list the files that are now stored in the dir object after running the system.file function
# number of files, including BAM files and files in the quants directory (which is our focus)
# list the files that are found only in the quants directory
list.files(file.path(dir, "quants"))

# load csv file that contains a table with info about each of the samples that links them to FASTQ
# and Salmon directories
# find the sample_table.csv file path
# read in the table csv file
csvfile <- file.path(dir, "sample_table.csv")
# read in the table csv file using the file path identified
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
# output what is now stored in coldata
coldata

# Loading Salmon quantification data into R
# use only first and second row samples (found in the quants directory)
coldata <- coldata[1:2,]
# create names column from coldata$Run
coldat$names <- coldata$Run
# create files column from the file identified for the "quant.sf.gz" file
# check to see if the coldata$files object exists (and therefore the file.path function worked)
file.exists(coldata$files)

# install the tximeta package if not previously installed 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("tximeta")


# Load the tximeta package
library(tximeta)

# Run the main function of the tximeta package
se <- tximeta(coldata)
# tximeta imports data at the transcript level

# Check the dimension of the se object created
dim(se)
# Output the first six rownames of the se object
head(rownames(se))

# Summarize the transcript-level quantifications to the gene-level
gse <- summarizeToGene(se)

# Check that the summarizeToGene function worked by checking dimensions of the gse object
# Checking that row IDs are now gene IDs (and so dimensions should be reduced)
dim(gse)

# Output the first six rownames of gse object
head(rownames(gse))

# 2.4 DESeq2 import functions
# Following tools can be used to generate or compile count data for use with DESeq:
# tximport, tximeta, htseq-count, featureCounts, and summarizeOverlaps
#
# 2.5 SummarizedExperiment
# assay contains matrix of counts
# rowRanges contains info about genomic ranges (GRanges of the genes)
# colData contains info about the samples
# 
# three matrices can be found in gse object created above:
# counts- fragment counts for each gene and sample
# abundance- estimated transcript abundances
# length- effective gene lengths

# load the full count matrix corresponding to all samples and all data
data(gse)
# Output the gse object
gse

# Output names of the assays
assayNames(gse)
# Output the first three rows of assay(gse)
# First matrix contains the counts
head(assay(gse), 3)

# Output the column sums for each of the counts per SRR identifier
colSums(assay(gse))

# Output rowRanges to show the ranges for the first five and last five genes
rowRanges(gse)

# Use the seqinfo function to print metadata about the sequences in the gse object
seqinfo(rowRanges(gse))

# colData reflection the datframe that was imported to tximeta function
colData(gse)


# 2.6 Branching point
# Could now use a variety of Bioconductor packages in analysis and DE
# This workflow will continue with the use of DESeq2

# The DESeqDataSet object, sample information and the design formula
# custom class in DESeq2 is DESeqDataSet; built on top of SummarizedExperiment class

# examine the columns of the colData of gse by useing $ directly
# on SummarizedExperiment or DESEqDataSet
gse$donor

gse$condition

# rename the variables to make easier for analysis
# use cell to denote donor cell line and dex to denote treatment condition
gse$cell <- gse$donor
gse$dex <- gse$condition

# change the names of the levels (be sure not to change the order!)
# rename "Untreated" as "untrt" and "Dexamethasone" as "trt"
levels(gse$dex)
levels(gse$dex) <- c("untrt", "trt")
#
# first level in R prefered to be the reference level (control or untreated)
# to relevel if necessary (though not in this case)
# load the magrittr library
# library(magrittr)
# relevel the first level of the factor to be the untreated samples
# gse$dex %<>% relevel ("untrt)
# gse$dex

############################
### Workflow modification to use prepared GSE
### Alternative workflow uses summarized experiment from airway package
### Commands pulled from F1000Rsearch 2016, 4:1070 (Love et al 2016)
#############################
library(airway)
data("airway")
gse <- airway
# 
# # 3.1 Starting from SummarizedExperiment
# # Check the millions of fragments that could be mapped by Salmon to the genes
# # use of round function to specify how many decimal points to keep (1 in this case)
# round(colSums(assay(gse))/ 1e6 , 1 )
# 
# change the names of the levels (be sure not to change the order!)
# rename "Untreated" as "untrt" and "Dexamethasone" as "trt"
levels(gse$dex)
levels(gse$dex) <- c("untrt", "trt")


# construct a DESeqDataSet object 
# add an appropriate design for the analysis 
# design in this case is to test for the effect of dex controlling for the effect of the different cell lines
# load the DESeq2 library 
library(DESeq2)
dds <- DESeqDataSet(gse, design = ~ cell + dex)
#
# 3.2 Starting from count matrices
# Showing how to build a DESeqDataSet if there is only a count matrix and a table of sample info
# Skip is SummarizedExperiment object has been prepared
# To see the fragment counts, use assay accessor function
# countdata <- round(assays(gse)[["counts"]])
# 
# # Output the first three lines of countdata
# # head(countdata, 3)
# 
# 
# # Important to check that the columns of the count matrix
# # correspond to the rows of the sample information table
# coldata <_ colData(gse)
# 
# # construct the DESeqDataSet object from the matrix of counts and the sample info table
# # ddsMat <- DESeqDataSetFromMatrix(countData = countData,
# #                                  colData = coldata,
# #                                  design = ~ cell + dex)
# 
# # Exploratory analysis and visualization
# 
# # 4.1 Pre-filtering the dataset
# # removing rows of the DESeqDataSet that have no counts, or only a single count across all samples
# # check for number of rows in dds object
# nrow(dds)
# 
# tabulate row sums for the counts in dds where the number of counts is greater than 1
keep <- rowSums(counts(dds)) > 1
# modify the dds object to keep only rows where the number of counts is greater than 1 (and all columns)
dds <- dds[keep, ]
# check the number of rows in the dds object to make sure rows were removed

# # Alternative subset
# # keep <- rowSums(counts(dds) >= 10) >= 3
# 
# # 4.2 The variance stabilizing transformation and the rlog
# # In RNA-seq, expected variance grows with the mean
# # to avoid plots revolving around genes with highest counts, take the log
# # of the normalized count values plus a pseudocount of 1 (though lowest counts
# # could then cause noise depending on pseudocount value chosen)
# 
# # Simulated data to show the abobve scenario in lines 271-2743
# # Plot the Sd of each row (genes) against the mean
# lambda <- 10^seq(from = -1, to = 2, length = 1000)
# cts <- matrix(rpois(1000*100, lambda), ncol = 100)
# library(vsn)
# meanSdPlot(cts, ranks = FALSE)
# 
# # Log-transformed counts
# lot.cts.one <- log2(cts + 1)
# meanSdPlot(log.cts.one, ranks = FALSE)
# 
# Two transformations in DESeq2 to stabilize variance across the mean
# Variance stabilizing transformation and regularized-log transformation
# Both return a DESeqTransform object based on SummarizedExperiment class

# blind = FALSE means that differences between cell lines and treatment will not
# # contribute to the expected variance-mean trend of the experiment
# # use of VST transformation (recommended for datasets where n > 30)
# vsd <- vst(dds, blind = FALSE)
# # # output first three rows of vsd
# # head(assay(vsd), 3)
# # 
# # # access colData attached to dds
# # colData(vsd)
# # 
# # use of the rlog transformation (recommended for datasets where n < 30)
# rld <- rlog(dds, blind = FALSE)
# # output the first three rows of rld
# # head(assay(rld), 3)
# 
# # Show the effect of the transformation just performed
# # Plotting of the first sample against the second using log2, VST, and rlog
# 
# load the dplyr package
library(dplyr)
# load the ggplot2 package 
library(ggplot2)

# need to estimate size factors to account for sequencing depth in log2 approach
# dds <- estimateSizeFactors(dds)
# 
# df <- bind_rows(
#   as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
#         mutate(transformation = "log2(x + 1)"),
#   as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
#   as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
# 
# # set columns names of df
# colnames(df)[1:2] <- c("x", "y")
# 
# # plot df using gglot 
# # this will plot all three transformation examples side by side
# ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + 
#   coord_fixed() + facet_grid( .~ transformation)
# 
# # 4.3 Sample distances
# # Calculate Euclidean distance between samples using dist function
# # Transpose matrix values so samples are in the rows
# sampleDists <- dist(t(assay(vsd)))
# sampleDists
# 
# # Visualization of the distances in heatmap
# # load the pheatmap package
# library(pheatmap)
# # load the RColorBrewer package
# library(RColorBrewer)
# 
# # Manually provide sampleDists created above to clustering_distance arg of pheatmap function
# # create object to store sampleDists as a matrix
# sampleDistMatrix <- as.matrix( sampleDists )
# # Sets rownames as treated or untreated, and cell line
# rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
# # Sets column names to null
# colnames(sampleDistMatrix) <- NULL
# # Sets coloration of the heatmap to be created
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# # open pdf device to save heatmap
# pdf(file="sampleDistanceHeatmap.pdf")
# # use of the pheatmap function to generate the heatmap
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows = sampleDists,
#          clustering_distance_cols = sampleDists,
#          col = colors)
# # turn off pdf device
# dev.off()
# 
# # Calculating sample distanes using Poisson Distance
# # load the PoiClaClu package
# library(PoiClaClu)
# # Calculate distance using PoissonDistance function 
# poisd <- PoissonDistance(t(counts(dds)))
# 
# # Plot the heatmap of the distances calculated
# # create object to store sampleDists as a matrix
# samplePoisDistMatrix <- as.matrix( poisd$dd )
# # Sets rownames as treated or untreated, and cell line
# rownames(samplePoisDistMatrix) <- paste ( dds$dex, dds$cell, sep= " - ")
# # Sets column names to null
# colnames(samplePoisDistMatrix) <- NULL
# # open pdf device to save heatmap
# pdf(file="poissonHeatmap.pdf")
# 
# # use of the pheatmap function to generate the heatmap
# pheatmap(samplePoisDistMatrix,
#          clustering_distance_rows = poisd$dd,
#          clustering_distance_cols = poisd$dd,
#          col = colors)
# # turn off pdf device
# dev.off()
# 
# 4.4 PCA Plot
# PCA is yet another way to visualize sampe-to-sample distance
# open pdf device to save heatmap
# pdf(file="vsdPCA.pdf")
# 
# plotPCA(vsd, intgroup = c("dex", "cell"))
# # In the plot generated, sample are projected onto 2D plane
# # x axis is the drxn that separates the data points the most
# # y axis is drxn that separates the data points the second most
# 
# # Can also request that plotPCA return the data rather than building the plot
# pcaData <- plotPCA(vsd, intgroup = c("dex", "cel"), returnData = TRUE)
# pcaData
# 
# # Can now use the data returned to make PCA plot from scratch using ggplot2 package
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
#   geom_point(size =3) +
#   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
#   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
#   coord_fixed() +
#   ggtitle("PCA with VST data")

# Can see from that plot generated that the differences between cells is considerable
# but not stronger than the difference in treated vs. untreated

# PCA Plot using Generalized PCA
# Technique for performing dimension reduction on data not normally distributed
# load the glmpca package
# library(glmpca)
# # use of glmpca function
# gpca <- glmpca(counts(dds), L=2)
# gpca.dat <- gpca$factors
# gpca.dat$dex <- dds$dex
# gpca.dat$cell <- dds$cell
# 
# # use of ggplot to generate plot for glmpca
# ggplot(gpca.dat, aes(x = dim1, y = dim2, color = dex, shape = cell)) +
#   geom_point(size = 3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
# 
# # 4.6 MDS plot
# # Multidimensional scaling function in base R
# # Good to use when you only have matrix of distances rather than data
# 
# # Compute MDS for distances in VST data
# mds <- as.data.frame(colData(vsd)) %>%
#   cbind(cmdscale(sampleDistMatrix))
# # Plot mds using ggplot
# ggplot(mds, aes(x = '1', y = '2', color = dex, shape = cell)) +
#   geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
# 
# # Generation of same plot but for PoissonDistance
# mdsPois <- as.data.frame(colData(dds)) %>%
#   cbind(cmdscale(samplePoisDistMatrix))
# ggplot(mdsPois, aes(x = '1', y = '2', color = dex, shape = cell)) +
#   geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")
# 
# Differential Expression Analysis

# 5.1 Running the differential expression pipeline
# Can run DE pipeline on the raw counts with a single call to the DESeq function
dds <- DESeq(dds)

# 5.2 Building the results table
# Calling results without any args extracts estimated log2 fold changes and p values
# for the last variable in the design formula
res <- results(dds)
res

# Produce the same results table when calling result with a more specific command
res <- reults(dds, contrast=c("dex", "trt", "untrt"))

# Since res is a DataFrame object, it carries metadata with info on the meaning of the columns
mcols(res, use.names = TRUE)

# Summarize the results using summary() function
summary(res)

# Ways to be more strict in decision as to which genes are significant
# Lower the false discovery rate threshold
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

# Raise the log2 fold change threshold from 0 
# to show more substantial changes due to treatment 
# use the lfcThreshold argument
# setting of lfcThreshold equal to 1 tests for genes with significant effects 
# of treatment on gene counts more than doubling or less than halving
resLFC1 <- resuts(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

# 5.3 Other Comparisons
# Use contrast arg in results function to extract results for comparions between any
# two variable levels
# contrast needs to specify: 1) name of variable, 2) name of level  for numerator, 
# and 3) name of level for denominator
# extract results for log2 of fold change of one cell line over another
results(dds, contrast = c("cell", "N061011", "N61311"))

# 5.4 Multiple Testing
# In biology, use the p values to correct for multiple testing
# Generates number of genes with p-value <0.05
sum(res$pvalue < 0.05, na.rm = TRUE)

# Generates number of genes where a p-value was reported
sum(!is.na(res$pvalue))

# Consider here a fraction of 10% false positives acceptable
# Therefore consider all genes with adj p-value below 10% as significant
sum(res$padj < 0.1, na.rm = TRUE)

# Subset results table for genes with adj p-value less than 10%
# and then sort by log2 fold change estimates 
# This will output significant genes with the strongest down-regulation
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

# Output significant genes with the strongest up-regulation
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

# Plotting Results

# 6.1 Counts plot
# Quick way to visualize counts for a particular gene is plotCounts function
# args for plotCounts are: 1) DeSeqDataSet, 2) gene name, and 3) group over
# which to plot the counts
# Select for specific genes to be considered in plotCounts function
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dex"))

# Can also make custom plots using ggplot function from ggplot2
# load the ggbeeswarm library
library(ggbeeswarm)
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex", "cell"),
                         returnData = TRUE)
# Plot using geom_beeswarm- offsets points within categories to reduce overplotting
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() + geom_beeswarm(cex = 3)

# Plot that will display lines connecting cell lines
ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()

# 6.2 MA-plot
# also known as a mean-difference plot
# MA-plot provides overview for distribution of estimated coefficients in the model
# y-axis, M stands for minus (subtraction of log vals equivalent to the log of the ratio)
# x-axis, A stands for average

# load the apeglm library
library(apeglm)
# Need to shrink log2 fold changes for comparison of dex vs. untreated samples
# shrink using the lfcShrink function
resultsNames(dds)

res <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = apeglm)
# Generation of the MA plot
# each dot represents a gene
# genes with adj p-value below threshold (0.1) are in red
plotMA(res, ylim = c(-5, 5))

# Generation of plot example without using lfcShrink first
res.noshr <- results(dds, name= "dex_trt_vs_untrt")
plotMA(res.noshr, ylim = c(-5,5))

# Label individual points on the MA-plot using the with function
plotMA(res, ylim = c(-5, 5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

# Generate histogram of p-values (excludes genes with small counts that could generate spikes)
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

# 6.3 Gene Clustering
# Previous heatmap showed clustering of the samples

# load the genefilter library
library(genefilter)

# selects the 20 genes with the highest variance across samples
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

# More interesting to look at amount by which each gene deviates in a sample from the 
# gene's average across all samples
mat <- assay(vsd) [topVarGenes, ]
# center each genes' values acoss samples
mat <- mat - rowMeans(mat)
# as.data.frame function used to tell pheatmap how to label columns
anno <- as.data.frame(colData(vsd)[, c("cell", "dex")])
# call of pheatmap function to plot heatmap
pheatmap(mat, annotation_col = anno)

# Is this necessary since DESeq2 automatically peforms independent filtering?
# 6.4 Independent filtering
# Removing low count genes from the input to FDR procedure allows for more genes to be significant
# among those kept- improves the power of the test

# Generation of a bar plot by examining ratio of small p-values for genes binned by mean 
# normalized count
# create bins by using the quantile function
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
# bin genes by base mean using cut
bins <- cut(resLFC1$baseMean, qs)
# rename the levels of bins using middle point
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
# calculate the ratio of p values less than 0.05 for each bin
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
                          mean(p < .05, na.rm = TRUE))
# plot the ratios
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")
# plot generated shows that genes with very low mean count have little or now power
# and should therefore be excluded from testing

# 6.5 Independent Hypothesis Weighting
# weight hypotheses to optimize power
# IHW- Independent Hypothesis Weighting used in lieu of independent filtering in 6.4

# load the IHW package
library(IHW)
res.ihw <- results(dds, filterFun= ihw)

# Annotating and exporting results
# result table so far is composed of Ensembl gene IDs (alternative gene names may be more helpful)

# load the AnnotationDbi package
library(AnnotationDbi)

# load the org.Hs.ef.db package
# this is the annotation package for Homo sapiens that uses Entrez Gene IDs as primary key
library(org.Hs.eg.db)

# access the available key types
columns(org.Hs.eg.db)

# Add individual columns to the results table using mapIds function
# provide row names of our results table as a kay
# Add the gene symbol and Entrez ID to our results table (call mapIds function twice)
# column argument states what info is desired
# multiVals tells mapIds what value to use if there are multiple values for a single input value
# substrings of rownames for character 1-15?
ens.str <- substr(rownames(res), 1, 15)
# Adding the gene symbol to the table
res$symbol <- mapIds(org.Hs.eg.db, 
                     keys = ens.str,
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# Adding the entrez id to the table
res$entrez <- mapIds(org.Hs.eg.db, 
                     keys = ens.str,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# Results now has the following external gene IDs
resOrdered <- res[order(res$pvalue), ]
# Output first six rows of resOrdered
head(resOrdered)

# 7.1 Exporting Results
# convert the DataFrame object (IRanges package) to a data.frame object that can be
# processed by write.csv
# use of the top 100 genes for demonstration
resOrderedDF <- as.data.frame(resOrdered[1:100, ])
# save the results in a csv file named "results.csv"
write.csv(resOrderedDF, file = "results.csv")

# another way to export results is through use of the ReportingTools package
# ReportingTools automatically generated HTML documents with links to external databases using gene identifiers
# and boxplots summarizing normalized counts across groups

# load the ReportingTools library
library(ReportingTools)

htmlRep <- HTMLReport(shortName = "report", title = "My report",
                      reportDirectory = "./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)

# 7.2 Plotting Fold Changes in Genomic Space
# use of lfcShrink function to return a GRanges object
# arguments in lfcShrink include the dds dataframe, the coefficient as the treatment type
# use of the apeglm method for shrinking, and format to return a GRanges object
resGR <- lfcShrink(dds, coef= "dex_trt_vs_untrt", type = "apeglm", format = "GRanges")
# output of the GRanges object created using the lfcShrink function
resGR

# Add in the symol for labeling the genes on the plot
# substrings of rownames for character 1-15?
ens.str <- substr(names(resGR), 1, 15)
# use of mapIds function to add "SYMBOL" column to resGR object
resGR$symbol <- mapIds(org.Hs.eg.db, ens.str, "SYMBOL", "ENSEMBL")

# Use of the Gviz package to plot GRanges and associated metadata
# load the Gviz package
library(Gviz)

# Specification of a window of 1 million bp upstream and downstream 
# from the gene with the smallest p-value
window <- resGR[topGene] + 1e6
# specify the strand (either strand in this case)
strand(window) <- "*"
# creation of a subset of the full results
resGRsub <- resGR[resGR %over% window]
# check if the gene symbol exists or is duplicated
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
# add the gene symbol as a name if the symbol exists and is not duplicated in the subset
# done through use of naOrDup object created in previous line
resGRsub$group <- ifelse(naOrCup, names(resGRsub), resGRsub$symbol)

# creation of a vector specifying if the genes in the subset had a low value of padj
status <- factor(ifelse(resGRsub$padj < 0.5 & !is.na(resGRsub$padj), "sig", "notsig"))

# Plot the results using Gviz functions
# specify to not use use UCSC chromosome names
options(ucscChromosomeNames = FALSE)
# creation of axis track to specify location in the genome
g <- GenomeAxisTrack()
# creation of axis track to show the genes and their names
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
# creation of a data track to display vertical bars showing the moderated log fold change
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
# genes and their names colored by significance in the plot generated
plotTracks(list(g, d, a), groupAnnotation = "group", notsig = "grey", sig = "hotpink")

# Removing Hidden Batch Effects

# can use statistical methods designed for RNA-seq from the sva or RUVSeq package 
# to detect groupings that could be causing hidden or unwanted variation 
# example in our case would be if we were unaware that there were different cell lines used

# 8.1 Using SVA with DESeq2

# load the sva packages
library(sva)

# normalized counts from dds put into dat object
dat <- counts(dds, normalized = TRUE)
# rowMeans in dat with average count across samples being larger than 1
idx <- rowMeans(dat) > 1
# recreate dat object to include only those samples in the idx object (count greater than 1)
dat <- dat[idx, ]
# use of a full model matrix with th dex variable
mod <- model.matrix(~ dex, colData(dds))
# use of a reduced/null model matric with only an intercept form
mod0 <- model.matrix(~ 1, colData(dds))
# specify when using the svaseq function that we want 2 surrogate variables estimated
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

# output the sv values
svseq$sv

# Since we know the cell lines, we can observe how well the SVA method did in recovering
# the variables using the figure below
# set 2 rows and 1 column using mfrow
# set the margins of the plot using mar
par(mfrow = c(2,1), mar = c(3, 5, 3, 1))
# for loop to generate chart for SV1 and SV2
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$cell, vertical = TRUE, main = paste0("sv", i))
  abline(h = 0)
}

# remove any effects on the counts from the surrogate variables
# add the variables as columns to DESeqDataSet and then add them to the design
ddsva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV2 + SV2 + dex

# 8.2 Using RUV with DESeq2

# Another method to detect hidden batch effects

# load the RUVSeq package
library(RUVSeq)

# creation of neqSeqExpressionSet object of the dds counts
set <- newSeqExpressionSet(counts(dds))
# rowSums where counts are greater than 5 and then greater than or equal to 2
idx <- rowSums(counts(set) > 5) >= 2
# recreate set object to include only those samples in the idx object
set <- set[idx, ]
# recreate set object and use the betweenLaneNormalization function
set <- betweenLaneNormalization(set, which = "upper")
# pulling out rownames where the p-value is greater than 0.1 (and therefore not significant)
not.sig <- rownames(res)[which(res$pvalue > .1)]
# creation of empirical object to store set of empirical control genes (those that do not have a small p-value)
empirical <- rownames(set)[ rownames,(set) %in% not.sig ]
# use of the RUVg function to estimate factors or unwanted variation
# same as SVA's surrogate variables
set <- RUVg(set, empirical, k=2)
# display pData of set object after performing RUV method
pData(set)

# Plot the factors estimated by RUV method
# set 2 rows and 1 column using mfrow
# set the margins of the plot using mar
par(mfrow = c(2,1), ma = c(3,5,3,1))
# for loop to generate chart for W1 and W2
for (i in 1:2) {
  stripchart(pData(set)[, i] ~ dds$cell, vertical = TRUE, main = paste0("w", i))
  abline(h = 0)
}

# to control for these factors, add to DEseqDataSet and then th design
ddsruv <- dds
ddsruv$W1 <- set$W_1
ddsruuv$W2 <- set$W_2
design(ddsruv) <- ~ W1 + W2 + dex

# Time Course Experiments

# DESeq2 can be used to analyse time course experiments to find genes that react
# in a condition-specific manner overtime compared to baseline samples
#################################################################
#
#   #Example today is use of the fission package that contains gene counts for RNA-seq
#   #time course of fission yeast. The yeast were exposed to oxidatibe stress, and half
#   #of the samples have a delection of atf21 gene
#
#################################################################
# load the fission package
library(fission)
# load in the data
data(fission)
# Design formula using fission data, strain difference at 0 (strain), difference over time (minute), 
# and any strain-specific differences over time (strain:minute)
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)

# likelihood ratio test
# removal of the strain -specific differencs over time
# small p-values are those which at one or more time points after time 0 showed a strain-specific effect
ddsTC <- DESeq(ddxTC, test = "LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
# addition of column to resTC with symbol from ddsTC
resTC$symbol <- mcols(ddsTC)$symbol
# output first for rows in order of increasing adjusted p-value
head(resTC[order(resTC$padj), ], 4)

# plot the counts for the two groups using ggplot2 package
# use of plotCounts function to aid in visualize counts
# specify group over which to plot the counts to be minute and strain
fiss <- plotCounts(ddsTC, which.min(resTC$padj),
                   intgroup = c("minute", "strain"), returnData = TRUE)
# creation of minute column in fiss to be numeric
fiss$minute <- as.numeric(as.character(fiss$minute))
# plot generated using fiss counts, colored and grouped by strain
ggplot(fiss, aes(x= minute, y= count, color= strain, group= strain)) +
  geom_point() + stat_summary(fun.y= mean, geom= "line") + scale+y_log10()

# check the names of what can be accessed in ddsTC results function output
resultsNames(dds)

# Wald test using the test argument in the results function to test for
# log2 fold changes at individual time points
res30 <- results(ddsTC, name = "strainmut.minute30", test= "Wald")
res30[which.min(resTC$padj), ]

# cluster significant genes by their profiles
# extract matrix of shrunken log2 fold changes using coef function
betas <- coef(ddsTC)
# observe column names in betas
colnames(betas)

# plot the log2 fold changes in a heatmap
# extract genes with smallest adjusted p-values
topGenes <- head(order(resTC$padj), 20)
mat <- betas[topGenes, -c(1,2)]
# establish threshold of 3
thr <- 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks= seq(from= -thr, to= thr, length= 101), cluster_col= FALSE)

# Session Information

# call sessionInfo function to report versions of R and al packakges used in the session
sessionInfo()