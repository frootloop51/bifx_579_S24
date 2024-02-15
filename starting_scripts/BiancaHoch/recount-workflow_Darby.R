##############################################################
#### Author: Miranda Darby modification of script by Leonardo Collado Torres
#### Link to original script: http://master.bioconductor.org/packages/release/workflows/vignettes/recountWorkflow/inst/doc/recount-workflow.R
#### Recount2 demonstration workflow
#### http://master.bioconductor.org/packages/release/workflows/vignettes/recountWorkflow/inst/doc/recount-workflow.html
#### Last modified: November 26, 2019
##############################################################

###### ONLY NEEDED ONCE
#  ## Install packages from Bioconductor
#  install.packages("BiocManager")
  # BiocManager::install(c("recount", "GenomicRanges", "limma", "edgeR", "DESeq2",
  #     "regionReport", "clusterProfiler", "org.Hs.eg.db", "gplots",
  #     "derfinder", "rtracklayer", "GenomicFeatures", "bumphunter",
  #     "derfinderPlot", "sessioninfo"))

## ----"load libraries", message = FALSE, warning = FALSE--------------------
library("recount")
library("GenomicRanges")
library("limma")
library("edgeR")
library("DESeq2")
library("regionReport")
library("clusterProfiler")
library("org.Hs.eg.db")
library("gplots")
library("derfinder")
library("rtracklayer")
library("GenomicFeatures")
library("bumphunter")
library("derfinderPlot")
library("sessioninfo")

################################################################
### Section 2.3 Demonstrate how exons are counted by recount2
################################################################
#QUESTION 1:What does the disjoin function do?

#QUESTION 2: How is calculating the AUC (area under the curve) different from simply counting the number of reads that overlap a region (which is how most gene count matrices are generated)?

#QUESTION 3: How can can the AUC values in the recount RSEs be converted to scaled read counts (in other words, what math is being done)?
################################################################

library("GenomicRanges")
exons <- GRanges("seq", IRanges(start = c(1, 1, 13), end = c(5, 8, 15)))
exons
disjoin(exons)

#### Demonstrate how coverage of each base is summed to give overall count 
library("GenomicRanges")
reads <- GRanges("seq", IRanges(
    start = rep(
        c(1, 2, 3, 4, 5, 7, 8, 9, 10, 13, 14), 
        c(3, 1, 2, 1, 2, 1, 2, 1, 2, 4, 1)
    ), width = rep(
        c(1, 3, 2, 3, 1, 2, 1, 3, 2, 3, 2, 1, 3),
        c(1, 4, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1)
    )
))
## Get the base-level genome coverage curve
cov <- as.integer(coverage(reads)$seq)

## AUC = AREA UNDER THE CURVE
sum(cov)

#### Demonstrate how to scale the reads so that analyses can be combined across studies in recount2 
#### Default scaling is for 40 million reads per sample with read length = 100bp
?scale_counts()

## Check that the number of reads is less than or equal to 40 million
## after scaling.
library("recount")
rse_scaled <- scale_counts(rse_gene_SRP009615, round = FALSE)
summary(colSums(assays(rse_scaled)$counts)) / 1e6

#### Nice discussion of annotated exons versus intron retention versus novel exon

################################################################
### Section 2.4 Analyzing differential expression of pre-counted genes
################################################################
#QUESTION 4: How can you find the total number of reads and the read length to adjust the scale_counts() function?

#QUESTION 5: What are the different functions that you can use to automatically find and add metadata to your dataset?

#QUESTION 6: How can you add a new metadata column that you created to the colData for an RSE?
################################################################

library("recount")

## Find the project ID by searching abstracts of studies
abstract_search("human brain development by age")

## Download the data if it is not there
if(!file.exists(file.path("SRP045638", "rse_gene.Rdata"))) {
    download_study("SRP045638", type = "rse-gene")
}

## Check that the file was downloaded
file.exists(file.path("SRP045638", "rse_gene.Rdata"))

## Load the data
load(file.path("SRP045638", "rse_gene.Rdata"))

## ----colData---------------------------------------------------------------
## One row per sample, one column per phenotype variable
dim(colData(rse_gene))

## Mostly technical variables are included
colnames(colData(rse_gene))

## ----"explore colData"-----------------------------------------------------
## Input reads: number reported by SRA might be larger than number
## of reads Rail-RNA downloaded
colData(rse_gene)[,
    c("read_count_as_reported_by_sra", "reads_downloaded")]
summary(
    colData(rse_gene)$proportion_of_reads_reported_by_sra_downloaded
)

## AUC information used by scale_counts() by default
head(colData(rse_gene)$auc)

## Alternatively, scale_scounts() can use the number of mapped reads
## and other information
colData(rse_gene)[, c("mapped_read_count", "paired_end",
    "avg_read_length")]

## ----sharq-----------------------------------------------------------------
## SHARQ tissue predictions: not present for all studies
head(colData(rse_gene)$sharq_beta_tissue)
head(colData(rse_gene_SRP009615)$sharq_beta_tissue)

## ----characteristics-------------------------------------------------------
## GEO information was absent for the SRP045638 data set
colData(rse_gene)[, c("geo_accession", "title", "characteristics")]

## GEO information for the SRP009615 data set
head(colData(rse_gene_SRP009615)$geo_accession)
head(colData(rse_gene_SRP009615)$title, 2)
head(colData(rse_gene_SRP009615)$characteristics, 2)

## Similar but not exactly the same wording used for two different samples
colData(rse_gene_SRP009615)$characteristics[[1]]
colData(rse_gene_SRP009615)$characteristics[[11]]

## Extract the target information
target <- sapply(colData(rse_gene_SRP009615)$characteristics, "[", 2)
target

## Build a useful factor vector, set the reference level and append the result 
## to the colData() slot
target_factor <- sapply(strsplit(target, "targeting "), "[", 2)
target_factor[is.na(target_factor)] <- "none"
target_factor <- factor(target_factor)
target_factor <- relevel(target_factor, "none")
target_factor
colData(rse_gene_SRP009615)$target_factor <- target_factor

## ----"add_predictions"-----------------------------------------------------
## Before adding predictions
dim(colData(rse_gene))

## Add the predictions
rse_gene <- add_predictions(rse_gene)

## After adding the predictions
dim(colData(rse_gene))

## Explore the variables
colData(rse_gene)[, 22:ncol(colData(rse_gene))]

## ----"sra_run_table"-------------------------------------------------------
## Save the information from 
## https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP045638
## to a table. We saved the file as SRP045638/SraRunTable.txt.
file.exists(file.path("SRP045638", "SraRunTable.txt"))

## Read the table
sra <- read.table(file.path("SRP045638", "SraRunTable.txt"),
    header = TRUE, sep = "\t")

## Explore it
head(sra)

## We will remove some trailing '_s' from the variable names
colnames(sra) <- gsub("_s$", "", colnames(sra))

## Choose some variables we want to add
sra_vars <- c("sex", "race", "RIN", "age", "isolate", "disease",
    "tissue")

## Re-organize the SRA table based on the SRA Run IDs we have
sra <- sra[match(colData(rse_gene)$run, sra$Run), ]

## Double check the order
identical(colData(rse_gene)$run, as.character(sra$Run))

## Append the variables of interest
colData(rse_gene) <- cbind(colData(rse_gene), sra[, sra_vars])

## Final dimensions
dim(colData(rse_gene))

## Explore result
colData(rse_gene)[, 34:ncol(colData(rse_gene))]

## ----"sex preds"-----------------------------------------------------------
table("Predicted" = colData(rse_gene)$predicted_sex,
    "Observed" = colData(rse_gene)$sex)

## ----"age_groups"----------------------------------------------------------
## Create the original 6 age groups
age_bins <- cut( colData(rse_gene)$age, c(-1, 0, 1, 10, 20, 50, Inf),
    include.lowest=TRUE )
levels( age_bins ) <- c("prenatal", "infant", "child", "teen", "adult",
    "late life")
colData(rse_gene)$age_group <- age_bins

## ----"prenatal_factor"-----------------------------------------------------
## Create prenatal factor
colData(rse_gene)$prenatal <- factor(
    ifelse(colData(rse_gene)$age_group == "prenatal", "prenatal",
        "postnatal"),
    levels = c("prenatal", "postnatal"))

## ----"scale_counts"--------------------------------------------------------
## Scale counts
rse_gene_scaled <- scale_counts(rse_gene)

## To highlight that we scaled the counts
rm(rse_gene)

## ----"filter_low"----------------------------------------------------------
## Extract counts and filter out lowly expressed geens
counts <- assays(rse_gene_scaled)$counts
filter <- rowMeans(counts) > 0.5

########## THIS IS ANOTHER WAY TO DO DIFFERENTIAL EXPRESSION OR YOU CAN CAN ADD DESIGN FORMULA TO RSE AND CREATE DESEQ DATASET
######## skip if desired
# ## ----"limmade1", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Multi-dimensional scaling plot of the gene-level data by age group."----
# library("limma")
# library("edgeR")
# 
# ## Build DGEList object
# dge <- DGEList(counts = counts[filter, ])
# 
# ## Calculate normalization factors
# dge <- calcNormFactors(dge)
# 
# ## Explore the data
# plotMDS(dge, labels = substr(colData(rse_gene_scaled)$prenatal, 1, 2) )
# 
# ## ----"limmade2", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Multi-dimensional scaling plot of the gene-level data by sex."----
# plotMDS(dge, labels = substr(colData(rse_gene_scaled)$sex, 1, 1) )
# tapply(colData(rse_gene_scaled)$RIN, colData(rse_gene_scaled)$prenatal,
#     summary)
# 
# ## Specify our design matrix
# design <- with(colData(rse_gene_scaled),
#     model.matrix(~ sex + RIN + prenatal))
# 
# ## ----"limmade3", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "voom mean-variance plot of the gene-level data."----
# ## Run voom
# v <- voom(dge, design, plot = TRUE)
# 
# ## Run remaining parts of the DE analysis
# fit <- lmFit(v, design)
# fit <- eBayes(fit)
# 
# ## ----"limmaplots1", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "MA plot of the gene-level data. Testing for prenatal and postnatal DE adjusting for sex and RIN."----
# ## Visually explore DE results
# limma::plotMA(fit, coef = 4)
# 
# ## ----"limmaplots2", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Volcano plot of the gene-level data. Testing for prenatal and postnatal DE adjusting for sex and RIN."----
# limma::volcanoplot(fit, coef = 4)
# 
# ## ----"report_setup"--------------------------------------------------------
# ## Extract data from limma-voom results
# top <- topTable(fit, number = Inf, sort.by = "none",
#     coef = "prenatalpostnatal")
################################


## Build a DESeqDataSet with the count data and model we used
library("DESeq2")
dds <- DESeqDataSet(rse_gene_scaled[filter, ], ~ sex + RIN + prenatal)

## Add gene names keeping only the Ensembl part of the Gencode IDs
rownames(dds) <- gsub("\\..*", "", rownames(dds))
####### At this point you can enter the RNAseq workflow that we already completed.

######## skipped code


######### This is a way to do GO term enrichment. You can also use DESeq results for this. 
#########  You do not need to run this code, but know that it exists.

library("clusterProfiler")
library("org.Hs.eg.db")

## Remember that limma_res had ENSEMBL IDs for the genes, you could also use DESeq2 results after converting the gene names to ensembl
head(rownames(limma_res))

## Perform enrichment analysis for Biological Process (BP)
## Note that the argument is keytype instead of keyType in Bioconductor 3.5
enrich_go <- enrichGO(
    gene = rownames(limma_res)[limma_res$padj < 0.001],
    OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
    universe = rownames(limma_res))

## Visualize enrichment results
dotplot(enrich_go, font.size = 7)

#####################################
### This shows how you can dowload and analyze a count matrix that is based on disjoint exons rather than gene counts.
#########  You do not need to run this code, but know that it exists.

######## skipped code

################################################################
### Section 2.5 Use bigWig Files and derfinder to analyze differential expression of "expressed regions"
################################################################
#QUESTION 7: The expressed_regions() function divides the genome into chucks of bases that are expressed above a certain threshold. How many expressed regions are on chr21 in this study?
################################################################


## Define expressed regions for study SRP045638, only for chromosome 21
regions <- expressed_regions("SRP045638", "chr21", cutoff = 5L,
    maxClusterGap = 3000L)
    
## Explore the resulting expressed regions
regions
summary(width(regions))
table(width(regions) >= 100)

## Keep only the ones that are at least 100 bp long
regions <- regions[width(regions) >= 100]
length(regions)

## Compute coverage matrix for study SRP045638, only for chromosome 21
## Takes about 4 minutes
rse_er <- coverage_matrix("SRP045638", "chr21", regions,
    chunksize = 2000, verboseLoad = FALSE, scale = FALSE)

# ## Use the expanded metadata we built for the gene model
# colData(rse_er) <- colData(rse_gene_scaled)

################# Alternatively, add colData from RSE for the same study
load(file.path("SRP045638", "rse_gene.Rdata"))
colData(rse_er) <- colData(rse_gene)
################

## Scale the coverage matrix
rse_er_scaled <- scale_counts(rse_er)

## To highlight that we scaled the counts
rm(rse_er)

########## You can skip the differential expression analysis of the regions -- you can use DESeq2 instead or try this code.
# ## ----"erdeanalysis1", eval = .Platform$OS.type != "windows", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Multi-dimensional scaling plot of the expressed regions level data by age group."----
# ## Build DGEList object
# dge_er <- DGEList(counts = assays(rse_er_scaled)$counts)
# 
# ## Calculate normalization factors
# dge_er <- calcNormFactors(dge_er)
# 
# ## Explore the data
# plotMDS(dge_er, labels = substr(colData(rse_er_scaled)$prenatal, 1, 2) )
# 
# ## ----"erdeanalysis2", eval = .Platform$OS.type != "windows", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Multi-dimensional scaling plot of the expressed regions level data by sex."----
# plotMDS(dge_er, labels = substr(colData(rse_er_scaled)$sex, 1, 1) )
# 
# ## ----"erdeanalysis3", eval = .Platform$OS.type != "windows", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "voom mean-variance plot of the expressed regions level data."----
# ## Run voom
# v_er <- voom(dge_er, design, plot = TRUE)
# 
# ## Run remaining parts of the DE analysis
# fit_er <- lmFit(v_er, design)
# fit_er <- eBayes(fit_er)
# 
# ## ----"erdeanalysis4", eval = .Platform$OS.type != "windows", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Volcano plot of the expressed regions level data. Testing for prenatal and postnatal DE adjusting for sex and RIN."----
# ## Visually explore the results
# limma::volcanoplot(fit_er, coef = 4)
# 
# ## Number of DERs
# top_er <- topTable(fit_er, number = Inf, sort.by = "none",
#     coef = "prenatalpostnatal")
# table(top_er$adj.P.Val < 0.001)
# 
# ## ----"sort_qvalue", eval = .Platform$OS.type != "windows"------------------
# ## Sort regions by q-value
# regions_by_padj <- regions[order(top_er$adj.P.Val, decreasing = FALSE)]
# 
# ## Look at the top 10
# regions_by_padj[1:10]
# width(regions_by_padj[1:10])
# 


################################################################
### Section 2.5.3 Plotting genomic regions of bigWigs
################################################################
#QUESTION 8: Which is the only bigWig in a project that does not have to be scaled before visualization?

################################################################
## Construct the list of bigWig URLs
## They have the following form:
## http://duffel.rail.bio/recount/
## project id
## /bw/
## sample run id
## .bw
################# Alternatively, if rse_er_scaled is not in memory, add colData from RSE for the same study
load(file.path("SRP045638", "rse_gene.Rdata"))
colData(rse_gene) ### replaced colData(rse_er_scaled) with colData(rse_gene) in all subsequent code to make this work
# add fake "prenatal" category
prenatal <- c(rep("prenatal", 36),rep("postnatal",36)) ## aribtrarily set first half as prenatal and second half as postnatal
prenatal <- as.factor(prenatal)
colData(rse_gene)$prenatal <- prenatal
################

bws <- paste0("http://duffel.rail.bio/recount/SRP045638/bw/",
    colData(rse_gene)$bigwig_file)

## Note that they are also present in the recount_url data.frame
bws <- recount_url$url[match(colData(rse_gene)$bigwig_file,
    recount_url$file_name)]

## Use the sample run IDs as the sample names
names(bws) <- colData(rse_gene)$run

## Skip this if you did not run DE on regions
# # ## Add 100 bp padding on each side
# regions_resized <- resize(regions_by_padj[1:10],
#     width(regions_by_padj[1:10]) + 200, fix = "center")

############### Alternative code if you did not run DE on regions
regions <- expressed_regions("SRP045638", "chr21", cutoff = 5L,
                             maxClusterGap = 3000L)
regions_resized <- resize(regions[1:10],
    width(regions[1:10]) + 200, fix = "center")
####################

## ----"regionCov", eval = .Platform$OS.type != "windows"--------------------
## Get the bp coverage data for the plots
library("derfinder")
regionCov <- getRegionCoverage(regions = regions_resized, files = bws,
    targetSize = 40 * 1e6 * 100,
    totalMapped = colData(rse_gene)$auc,
    verbose = FALSE)

## ----"gencode_txdb", eval = .Platform$OS.type != "windows"-----------------
## Import the Gencode v25 hg38 gene annotation
library("rtracklayer")
gencode_v25_hg38 <- import(paste0(
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human",
    "/release_25/gencode.v25.annotation.gff3.gz"))
            
## Keep only the chr21 info
gencode_v25_hg38 <- keepSeqlevels(gencode_v25_hg38, "chr21",
    pruning.mode = "coarse")

## Get the chromosome information for hg38
library("GenomicFeatures")
chrInfo <- getChromInfoFromUCSC("hg38")
chrInfo$chrom <- as.character(chrInfo$chrom)
chrInfo <- chrInfo[chrInfo$chrom %in% seqlevels(regions), ]
chrInfo$isCircular <- FALSE

## Assign the chromosome information to the object we will use to
## create the txdb object
si <- with(chrInfo, Seqinfo(as.character(chrom), length, isCircular,
    genome = "hg38"))
seqinfo(gencode_v25_hg38) <- si

## Switch from Gencode gene IDs to Ensembl gene IDs
gencode_v25_hg38$gene_id <- gsub("\\..*", "", gencode_v25_hg38$gene_id)

## Create the TxDb object
gencode_v25_hg38_txdb <- makeTxDbFromGRanges(gencode_v25_hg38)

## Explore the TxDb object
gencode_v25_hg38_txdb

## ----"bump_ann", eval = .Platform$OS.type != "windows"---------------------
library("bumphunter")
## Annotate all transcripts for gencode v25 based on the TxDb object
## we built previously.
ann_gencode_v25_hg38 <- annotateTranscripts(gencode_v25_hg38_txdb,
    annotationPackage = "org.Hs.eg.db",
    mappingInfo = list("column" = "ENTREZID", "keytype" = "ENSEMBL",
    "multiVals" = "first"), simplifyGeneID = TRUE)
    
## Annotate the regions of interest
## Note that we are using the original regions, not the resized ones
# skip this line and use next if regions_by_padj is not in memory
# nearest_ann <- matchGenes(regions_by_padj[1:10], ann_gencode_v25_hg38)
# alternative for if regions_by_padj is not in memory
nearest_ann <- matchGenes(regions[1:10], ann_gencode_v25_hg38)

## ----"make_gs", eval = .Platform$OS.type != "windows"----------------------
## Create the genomic state object using the gencode TxDb object
gs_gencode_v25_hg38 <- makeGenomicState(gencode_v25_hg38_txdb,
    chrs = seqlevels(regions))
    
## Annotate the original regions
regions_ann <- annotateRegions(regions_resized,
    gs_gencode_v25_hg38$fullGenome)

## ----"regionplots", eval = .Platform$OS.type != "windows", out.width=ifelse(on.bioc, "100%", "75%"), fig.align="center", fig.cap = 'Base-pair resolution plot of differentially expressed region 2.'----
library("derfinderPlot")
pdf('region_plots.pdf')
plotRegionCoverage(regions = regions_resized, regionCoverage = regionCov, 
   groupInfo = colData(rse_gene)$prenatal,
   nearestAnnotation = nearest_ann, 
   annotatedRegions = regions_ann,
   txdb = gencode_v25_hg38_txdb,
   scalefac = 1, ylab = "Coverage (RP40M, 100bp)",
   ask = FALSE, verbose = FALSE)
dev.off()

# ## Visualize DER #2
# plotRegionCoverage(regions = regions_resized, regionCoverage = regionCov, 
#    groupInfo = colData(rse_er_scaled)$prenatal,
#    nearestAnnotation = nearest_ann, 
#    annotatedRegions = regions_ann,
#    txdb = gencode_v25_hg38_txdb,
#    scalefac = 1, ylab = "Coverage (RP40M, 100bp)",
#    ask = FALSE, verbose = FALSE, whichRegions = 2)

# ## ----sessionInfo----------------------------------------------------------------------------------
# ## Final list of files created
# dir("SRP045638")
# 
# ## Pandoc information
# library("rmarkdown")
# pandoc_version()
# 
# ## Time for reproducing this workflow, in minutes
# round(proc.time()[3] / 60, 1)
# 
# options(width = 100)
# library("sessioninfo")
# session_info()
# 
