# DESEq2 vignette and paper explanations of dispersions and independent filtering


# From the DESeq2 paper:
# 
# DESeq2 uses the average expression strength of each gene, across all samples,
# as its filter criterion, and it omits all genes with mean
# normalized counts below a filtering threshold from multiple testing adjustment. DESeq2 by default will choose a
# threshold that maximizes the number of genes found at
# a user-specified target FDR.

# "alpha" is the significance cutoff used for optimizing the independent filtering (by default 0.1). 
# If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

# From the DESeq2 vignette:

# "alpha" for the results function:
#   By default, independent filtering is performed to select a set of genes for multiple test correction which maximizes the number of adjusted p-values less than a given critical value alpha (by default 0.1).
# 
# In DESeq2 version >= 1.10, the threshold that is chosen is the lowest quantile of the filter for which the number
# of rejections is close to the peak of a curve fit to the number of rejections over the filter quantiles.
# ’Close to’ is defined as within 1 residual standard deviation. The adjusted p-values for the genes
# which do not pass the filter threshold are set to NA.
# 
# By default, results assigns a p-value of NA to genes containing count outliers, as identified using
# Cook’s distance. See the cooksCutoff argument for control of this behavior. Cook’s distances for
# 42 results
# each sample are accessible as a matrix "cooks" stored in the assays() list. This measure is useful
# for identifying rows where the observed counts might not fit to a Negative Binomial distribution.

## See section 6.4 for an explanation of independent filtering. 

# To examine how this works:
library(DESeq2)
library(dplyr)
library(ggplot2)
library(airway)
data("airway")
gse <- airway
levels(gse$dex)
levels(gse$dex) <- c("untrt", "trt")
dds <- DESeqDataSet(gse, design = ~ cell + dex)
dds <- DESeq(dds)
res <- results(dds)
(table(res$padj < 0.1))
(table(res$padj < 0.05))
# Lower the false discovery rate threshold
res.05 <- results(dds, alpha = 0.05)
(table(res.05$padj < 0.1))
(table(res.05$padj < 0.05))

# Increase lfcThreshold
resLFC1 <- results(dds, lfcThreshold=1)
(table(resLFC1$padj < 0.1))
(table(resLFC1$padj < 0.05))

# Lower FDR and increase LFC
# Increase lfcThreshold
resLFC1a05 <- results(dds, alpha= 0.05, lfcThreshold=1)
(table(resLFC1a05$padj < 0.1))
(table(resLFC1a05$padj < 0.05))



(nhits <- sum(res$pvalue < 0.05, na.rm=TRUE))

(ngenes <-sum(!is.na(res$pvalue)))

(fivePercentOfGenes <- 0.05*ngenes)

(percentFalsePos <- fivePercentOfGenes/nhits)

# DESeq2 uses the Benjamini-Hochberg (BH) adjustment (Benjamini and Hochberg 1995) as implemented in the base R p.adjust function; in brief, this method calculates for each gene an adjusted p value that answers the following question: if one called significant all genes with an adjusted p value less than or equal to this gene’s adjusted p value threshold, what would be the fraction of false positives (the false discovery rate, FDR) among them, in the sense of the calculation outlined above?

(adhits <- sum(res$padj < 0.1, na.rm=TRUE))
(falsePos <- 0.1*adhits)
(adjustedPercentFalsePos <- falsePos/adhits)

resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

results(dds, contrast = c("cell", "N061011", "N61311"))

#plotting counts
topGene <- rownames(res)[which.min(res$padj)]
pdf(file="topGeneScatter.pdf")
plotCounts(dds, gene = topGene, intgroup=c("dex"))
dev.off()

library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                         returnData = TRUE)
pdf(file="colorScatter.pdf")
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

pdf(file="colorAndLineScatter.pdf")
ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()
dev.off()

#command to transfer plots from server
#scp darby@144.175.88.21:/home/darby/*Scatter.pdf ~/Dropbox/Hood/BIFX_572_2019/

#MA plot
library("apeglm")
resultsNames(dds)
pdf(file="shrinkMA.pdf")
res <- lfcShrink(dds, coef="dex_trt_vs_untrt", type="apeglm")
plotMA(res, ylim = c(-5, 5))
dev.off()

pdf(file="noShrinkMA.pdf")
res.noshr <- results(dds, name="dex_trt_vs_untrt")
plotMA(res.noshr, ylim = c(-5, 5))
dev.off()

pdf(file="geneLabelMA.pdf")
plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()

#histogram of p values
pdf(file="phist.pdf")
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
dev.off()

#Gene Clustering
vsd <- vst(dds, blind = FALSE)
library("genefilter")
library("pheatmap")
library("RColorBrewer")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])

pdf(file="mostVariableGenesHM.pdf")
pheatmap(mat, annotation_col = anno)
dev.off()

## Brief discussion of independent filtering
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))

pdf(file="pvalVsCountsHist.pdf")
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")
dev.off()

#skip IHW for today

# annotating and exporting
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

ens.str <- substr(rownames(res), 1, 15)
#substr(rownames(res), 1, 15) <- c(1:15)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resOrdered <- res[order(res$pvalue),]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results.csv")

library("ReportingTools")
library("htmltools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
# url is the location/name of the file, transfer folder and file to desktop to open with webbrowser

resGR <- lfcShrink(dds, coef="dex_trt_vs_untrt", type="apeglm", format="GRanges")
resGR
ens.str <- substr(names(resGR), 1, 15)
resGR$symbol <- mapIds(org.Hs.eg.db, ens.str, "SYMBOL", "ENSEMBL")

library("Gviz")

window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

status <- factor(ifelse(resGRsub$padj < 0.05 & !is.na(resGRsub$padj),
                        "sig", "notsig"))

options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
pdf(file="gvisPlot.pdf")
plotTracks(list(g, d, a), groupAnnotation = "group",
           notsig = "grey", sig = "hotpink")
dev.off()

library("sva")
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv

pdf(file="surrogateVariableInterceptsVsCellLine.pdf")
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$cell, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}
dev.off()

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex
 #how would we get the results?

library("RUVSeq")
set <- newSeqExpressionSet(counts(dds))
idx  <- rowSums(counts(set) > 5) >= 2
set  <- set[idx, ]
set <- betweenLaneNormalization(set, which="upper")
not.sig <- rownames(res)[which(res$pvalue > .1)]
empirical <- rownames(set)[ rownames(set) %in% not.sig ]
set <- RUVg(set, empirical, k=2)
pData(set)
pdf(file="RUVplot.pdf")
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(pData(set)[, i] ~ dds$cell, vertical = TRUE, main = paste0("W", i))
  abline(h = 0)
}
dev.off()

ddsruv <- dds
ddsruv$W1 <- set$W_1
ddsruv$W2 <- set$W_2
design(ddsruv) <- ~ W1 + W2 + dex

library("fission")
data("fission")
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),], 4)

fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("minute","strain"), returnData = TRUE)
fiss$minute <- as.numeric(as.character(fiss$minute))
pdf(file="timecourseCounts.pdf")
ggplot(fiss,
       aes(x = minute, y = count, color = strain, group = strain)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()
dev.off()

resultsNames(ddsTC)
res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")
res30[which.min(resTC$padj),]

betas <- coef(ddsTC)
colnames(betas)

topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr

pdf("timecourseHeatmap.pdf")
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)
dev.off()




