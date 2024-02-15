############################################
#
# Title: recount-brain metadata applied to the SRA study ERP001304. WITH CLOOJ.
# Abstract: The RNA was extracted from post-mortem cortical grey matter of STG from the left hemisphere of 9 pairs 
# of male subjects with schizophrenia and matched non-psychiatric controls. 
# Each sample was subjected to 76 cycles of sequencing from a single end in one lane of Illumina Genome Analyzer II.
# Total number of samples: 18
#
# Workflows and references used:
# BigWigs to count matrix- This protocol explains how to create a feature count matrix from coverage data stored in BigWig files.
# This feature count matrix can then be used for differential expression analyses. http://lcolladotor.github.io/protocols/bigwig_DEanalysis/
#
# Author: Bianca Hoch  (bih1@hood.edu)
# Date: 05/21/2019
#
###########################################

rm(list=ls()) #Clear R memory

###Set path variables that contain the different directories:###

###Install packages:
#install.packages("BiocManager")
BiocManager::install(c("ggplot2", "pheatmap", "RColorBrewer", "DESeq2", "rtracklayer","recount", "Gviz", "derfinder", "derfinderData", "GenomicRanges"))


### Load libraries:
library('recount')
library('Gviz')
library('derfinder')
library('derfinderData')
library('GenomicRanges')
library('rtracklayer')
library('dplyr')
library('DESeq2')
library('pheatmap')
library('RColorBrewer')
library('ggplot2')


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

#Details about counts

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

#Preparing for Differential Expression Analysis:

#Load the GRanges object that lists the new exons that we are interested in:
newGR <- ("~/Documents/splice_project/data/analysis/newGRanges_06112019.Rdata")
load(newGR)
newGR

#In the following code, I will lift the features from the GRanges list (hg19) to hg38.
chain <- import.chain("~/Documents/splice_project/data/hg19ToHg38.over.chain")
newGRhg38 <- unlist(liftOver(newGR, chain))

#Download the sample BigWig files for the project ("ERP001304"):
#download_study('ERP001304', type = "samples", outdir = "~/Documents/splice_project/data/rawdata/bigwigs", download = TRUE)

#Construct full paths to the BigWig files and fix the names:
bwfiles <- rawFiles(datadir = "~/Documents/splice_project/data/rawdata/bigwigs", samplepatt = '*.bw' , fileterm = NULL)
names(bwfiles) <- gsub('.bw', '', names(bwfiles))

#Get a list of chromosomes from the GRanges list containing the new exons:
chromosomes <- levels(seqnames(newGRhg38))

# This only needs to be run once.
# Import data and create a count matrix for all chromosomes present in the GRanges list (WITH CLOOJ):
# bw <- BigWigFileList(bwfiles)
# exons <- newGRhg38
# counts_exons <- matrix(NA, nrow = length(exons), ncol = length(bw))
# colnames(counts_exons) <- names(bw)
# for(i in seq_len(length(bw))) {
#   for(chr in chromosomes) {
#     eval(parse(text=paste0("coverage <- import(bw[[i]], as = 'RleList')$",chr)))
#     counts_exons[, i] <- sum(Views(coverage, ranges(exons)))
#     filename <- paste("~/Documents/splice_project/data/analysis/counts", chr, ".csv", sep = "")
#     write.csv(counts_exons, file = filename)
#   }
# }

# Take the count .csv files that were output from the function above and rbind them:
# Read in the .csv files
t1	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr1.csv")
t2	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr10.csv")
t3	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr11.csv")
t4	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countsChr12.csv")
t5	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr13.csv")
t6	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr14.csv")
t7	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr15.csv")
t8	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr16.csv")
t9	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr17.csv")
t10	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr18.csv")
t11	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr19.csv")
t12	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr2.csv")
t13	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr20.csv")
t14	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr21.csv")
t15	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr22.csv")
t16	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr3.csv")
t17	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr4.csv")
t18	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr5.csv")
t19	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr6.csv")
t20	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr7.csv")
t21	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr8.csv")
t22	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschr9.csv")
t23	<- read.csv(file="/Users/biancahoch/Documents/splice_project/data/analysis/counts/countschrX.csv")

#Rename the rownames column to "Novel Exon" for each counts dataframe
names(t1)[1] <- "NovelExon"
names(t2)[1] <- "NovelExon"
names(t3)[1] <- "NovelExon"
names(t4)[1] <- "NovelExon"
names(t5)[1] <- "NovelExon"
names(t6)[1] <- "NovelExon"
names(t7)[1] <- "NovelExon"
names(t8)[1] <- "NovelExon"
names(t9)[1] <- "NovelExon"
names(t10)[1] <- "NovelExon"
names(t11)[1] <- "NovelExon"
names(t12)[1] <- "NovelExon"
names(t13)[1] <- "NovelExon"
names(t14)[1] <- "NovelExon"
names(t15)[1] <- "NovelExon"
names(t16)[1] <- "NovelExon"
names(t17)[1] <- "NovelExon"
names(t18)[1] <- "NovelExon"
names(t19)[1] <- "NovelExon"
names(t20)[1] <- "NovelExon"
names(t21)[1] <- "NovelExon"
names(t22)[1] <- "NovelExon"
names(t23)[1] <- "NovelExon"

#bind them
megaTable <- bind_rows(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23) %>%
  group_by(NovelExon) %>%
  summarise_each(funs(sum))

#Turn the megaTable into a dataframe
mt <- as.data.frame(megaTable)

#Divide by read length and round to integer numbers
mtr <- round(mt[,2:19] / 76, 0)

#Add novel exon position as rownames to the table
rownames(mtr) <- make.names(exNames, unique = TRUE)

#Explore a little portion of the count matrix
dim(mtr)
mtr[1425:1500, 1:8]

#Filter out rows that have no or nearly no information about the amount of gene expression.
mtr <- mtr[ rowSums(mtr) > 1, ]
  
#Now that the count matrix is prepared, we can proceed with differential expression.
#We will use DESeq2 to find differentially expressed exons between schizophrenic and nonschizophrenic individuals.

#Reorganize RSEmeta so that it can be used as colData for the new RSE (rseNew)
orMeta <- colData(RSEmeta)
orMeta <- orMeta[order(rownames(orMeta)),]
orMeta

#Check to see if orMeta is identical to the colnames of counts
identical(rownames(orMeta),colnames(mtr))

#Round matrix and specify design
cd <- orMeta
cd <- DataFrame(cd)

#Replace NAs in "disease" variable with "control"
cd$disease[1:9] <- "control"
dse <- DESeqDataSetFromMatrix(mtr,cd, ~ disease)
dse <- dse[rowSums(counts(dse))>1,]

#Perform DE analysis
dse <- DESeq(dse, test = 'Wald')
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 692 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing

#Visualizing sample distances:

#------------06/17/2019------------#

#Conducting an rlog transformation for exploratory analysis of multidimensional data.
#rld <-rlog(dse, blind = FALSE)
# Warning message:
#   In sparseTest(counts(object, normalized = TRUE), 0.9, 100, 0.1) :
#   the rlog assumes that data is close to a negative binomial distribution, an assumption
# which is sometimes not compatible with datasets where many genes have many zero counts
# despite a few very large counts.
# In this data, for 71.7% of genes with a sum of normalized counts above 100, it was the case 
# that a single sample's normalized count made up more than 90% of the sum over all samples.
# the threshold for this warning is 10% of genes. See plotSparsity(dds) for a visualization of this.
# We recommend instead using the varianceStabilizingTransformation or shifted log (see vignette).

#Using varianceStabilizingTransformation
vst1 <- varianceStabilizingTransformation(dse, blind = FALSE, fitType = "parametric")
plot(assay(vst1)[,1:2], pch=16, cex=0.3)

#Plotting sample distances
sampleDists <- dist(t(assay(vst1)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(orMeta$run, orMeta$disease, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9,"Blues"))) (255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors)

#Another way to visualize sample-to-sample distances is a principle componenets analysis (PCA).
plotPCA(vst1, intgroup = c("run", "disease"))

#Multidimensional scaling (MDS) plot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vst1)))
ggplot(mds, aes(X1,X2,color=run,shape=disease)) + geom_point(size=3)

#According to the PCA, MDS, and Sample Distance plots, it seems that the sample ERR103438 (schizophrenia) is an outlier.
#Extract results
results <- results(dse)

#Explore results:
results

#How many have a significant p-value?
sum(results$padj < 0.1, na.rm = TRUE)
#A total of 332 potential novel exons.

#Subset the results which have a padj higher than 0.1.
sig <- results[which(results$padj < 0.1),] 

#Look at metdata with information on the meaning of the columns.
mcols(sig, use.names = TRUE)

#We can also summarize the results with the following line of code, which reports some addition information.
summary(sig)

#Here we look at the significant novel exons with the strongest down-regulation:
head(sig[order(sig$log2FoldChange), ])

#These are the significant novel exons with the strongest up-regulation:
head(sig[order(sig$log2FoldChange, decreasing = TRUE), ])

#Exploratory analysis and visualization.

#Plotting the results
topExon <- rownames(sig)[which.min(sig$padj)]
exonCounts <- plotCounts(dse, gene = topExon, intgroup = c("disease","run"), returnData = TRUE)
ggplot(exonCounts, aes(x=disease, y=count, color=run)) + scale_y_log10() +
  geom_point(position = position_jitter(width=.1,height=0), size=3)

#According to the PCA, MDS, and Sample Distance plots, it seems that the sample ERR103438 (schizophrenia) is an outlier.
#I will remove this sample and reconstruct the plots.

#Remove the ERR103438:Schizophrenia sample
mtr2 <- mtr[,1:17]

#Remove all rows that have mostly zeros or little to no information on expression
mtr2 <- mtr2[ rowSums(mtr2) > 1, ]

#Now that the count matrix is prepared, we can proceed with differential expression.
#We will use DESeq2 to find differentially expressed exons between schizophrenic and nonschizophrenic individuals.

#Remove the ERR103438:Schizophrenia sample from RSEmeta so that it can be used as 
#colData for the new RSE.
orMeta2 <- orMeta[1:17,]

#Check to see if orMeta2 is identical to the colnames of counts
identical(rownames(orMeta2),colnames(mtr2))

#Round matrix and specify design
cd2 <- orMeta2
cd2 <- DataFrame(cd2)

#Replace NAs in "disease" variable with "control"
cd2$disease[1:9] <- "control"
dse2 <- DESeqDataSetFromMatrix(mtr2,cd2, ~ disease)
dse2 <- dse2[rowSums(counts(dse2))>1,]

#Perform DE analysis
dse2 <- DESeq(dse2, test = 'Wald')
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 215 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing

#Visualizing sample distances:

#------------06/18/2019------------#

#Conducting an rlog transformation for exploratory analysis of multidimensional data.
rld2 <-rlog(dse2, blind = FALSE)

#Plotting the transformation.
plot(assay(rld2)[,1:2], pch=16, cex=0.3)

#Plotting sample distances.
sampleDists <- dist(t(assay(rld2)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(orMeta2$run, orMeta2$disease, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9,"Blues"))) (255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors)

#Another way to visualize sample-to-sample distances is a principle componenets analysis (PCA).
plotPCA(rld2, intgroup = c("run", "disease"))

#Multidimensional scaling (MDS) plot.
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld2)))
ggplot(mds, aes(X1,X2,color=run,shape=disease)) + geom_point(size=3)

#Extract results
results <- results(dse2)

#Explore results:
results

#How many have a significant p-value?
sum(results$padj < 0.1, na.rm = TRUE)
#A total of 0 potential novel exons.

#Subset the results which have a padj higher than 0.1.
sig2 <- results[which(results$padj < 0.1),] 

#Look at metdata with information on the meaning of the columns.
mcols(sig2, use.names = TRUE)

#We can also summarize the results with the following line of code, which reports some addition information.
summary(sig2)

#Here we look at the significant novel exons with the strongest down-regulation:
head(sig[order(sig$log2FoldChange), ])

#These are the significant novel exons with the strongest up-regulation:
head(sig[order(sig$log2FoldChange, decreasing = TRUE), ])

#Exploratory analysis and visualization.

#Plotting the results
topExon <- rownames(sig)[which.min(sig$padj)]
exonCounts <- plotCounts(dse, gene = topExon, intgroup = c("disease","run"), returnData = TRUE)
ggplot(exonCounts, aes(x=disease, y=count, color=run)) + scale_y_log10() +
  geom_point(position = position_jitter(width=.1,height=0), size=3)
