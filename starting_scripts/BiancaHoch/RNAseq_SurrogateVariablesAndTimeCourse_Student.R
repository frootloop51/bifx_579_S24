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
library("pheatmap")
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




