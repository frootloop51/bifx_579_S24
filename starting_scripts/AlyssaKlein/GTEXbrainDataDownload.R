############################################
# Title: Aquisition of normalized RNAseq count data using recount2
# Author: Miranda Darby  (darby@hood.edu)
# Date: created on 7/22/19
# Recount2 paper: https://www.biorxiv.org/content/10.1101/247346v2
# Recount2 published workflow: https://f1000research.com/articles/6-1558/v1
# Recount2 website: https://jhubiostatistics.shinyapps.io/recount/
# Recount-brain paper: https://www.biorxiv.org/content/10.1101/618025v1
# Working directory 
###########################################
setwd("~/Dropbox/WomenHackathon/")
## Install BiocManager. Only needs to be done once.
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

### Install recount2. Also only needs to be done once
# BiocManager::install("recount")

library(recount)
### actually found online download easier. Clicked "brain.Rdata" under "Download RSE by Tissue" 
### on "GTEx" page from https://jhubiostatistics.shinyapps.io/recount/
### file location is: "/Users/darby/Downloads/rse_fc_SRP012682_brain.Rdata"
### download was on 7/22/19

load("/Users/darby/Downloads/rse_fc_SRP012682_brain.Rdata")

sampleData <- colData(rse_gene_SRP009615)

### access sample data and column data from full study
### used this to see what the data was like and what accessors were there
# sampleData <- colData(rse_fc)
# # dim(sampleData)
# # [1] 1409   82
# 
# counts <- assays(rse_fc)$counts
# # dim(counts)
# # [1] 109873   1409

#subset out samples with RIN greater than 7.5
# ranged summarized experiments are relational databases that can be subset based 
# on the information in the colData just as if they were a dataframe

gtex_brain <- rse_fc[,which(rse_fc$smrin >= 7.5)]

#find different brain regions
regions <- as.factor(gtex_brain$smtsd)
summary(regions)
# Brain - Amygdala 
# 8 
# Brain - Anterior cingulate cortex (BA24) 
# 22 
# Brain - Caudate (basal ganglia) 
# 82 
# Brain - Cerebellar Hemisphere 
# 76 
# Brain - Cerebellum 
# 22 
# Brain - Cortex 
# 14 
# Brain - Frontal Cortex (BA9) 
# 53 
# Brain - Hippocampus 
# 23 
# Brain - Hypothalamus 
# 33 
# Brain - Nucleus accumbens (basal ganglia) 
# 52 
# Brain - Putamen (basal ganglia) 
# 47 
# Brain - Spinal cord (cervical c-1) 
# 34 
# Brain - Substantia nigra 
# 14 

amygdala <- gtex_brain[,which(gtex_brain$smtsd == "Brain - Amygdala")]
amygdala_counts <- assays(amygdala)$counts
amygdala_sampleData <- colData(amygdala)
write.csv(amygdala_counts, file="amygdala_counts.csv")
write.csv(amygdala_sampleData, file="amygdala_sampleData.csv")

BA24 <- gtex_brain[,which(gtex_brain$smtsd == "Brain - Anterior cingulate cortex (BA24)")]
BA24_counts <- assays(BA24)$counts
BA24_sampleData <- colData(BA24)
write.csv(BA24_counts, file="BA24_counts.csv")
write.csv(BA24_sampleData, file="BA24_sampleData.csv")

caudate <- gtex_brain[,which(gtex_brain$smtsd == "Brain - Caudate (basal ganglia)")]
caudate_counts <- assays(caudate)$counts
caudate_sampleData <- colData(caudate)
write.csv(caudate_counts, file="caudate_counts.csv")
write.csv(caudate_sampleData, file="caudate_sampleData.csv")

cerebellum <- gtex_brain[,which(gtex_brain$smtsd == "Brain - Cerebellum")]
cerebellum_counts <- assays(cerebellum)$counts
cerebellum_sampleData <- colData(cerebellum)
write.csv(cerebellum_counts, file="cerebellum_counts.csv")
write.csv(cerebellum_sampleData, file="cerebellum_sampleData.csv")

BA9 <- gtex_brain[,which(gtex_brain$smtsd == "Brain - Frontal Cortex (BA9)")]
BA9_counts <- assays(BA9)$counts
BA9_sampleData <- colData(BA9)
write.csv(BA9_counts, file="BA9_counts.csv")
write.csv(BA9_sampleData, file="BA9_sampleData.csv")

hippocampus <- gtex_brain[,which(gtex_brain$smtsd == "Brain - Hippocampus")]
hippocampus_counts <- assays(hippocampus)$counts
hippocampus_sampleData <- colData(hippocampus)
write.csv(hippocampus_counts, file="hippocampus_counts.csv")
write.csv(hippocampus_sampleData, file="hippocampus_sampleData.csv")

hypothalamus <- gtex_brain[,which(gtex_brain$smtsd == "Brain - Hypothalamus")]
hypothalamus_counts <- assays(hypothalamus)$counts
hypothalamus_sampleData <- colData(hypothalamus)
write.csv(hypothalamus_counts, file="hypothalamus_counts.csv")
write.csv(hypothalamus_sampleData, file="hypothalamus_sampleData.csv")

nucleus_accumbens <- gtex_brain[,which(gtex_brain$smtsd == "Brain - Nucleus accumbens (basal ganglia)")]
nucleus_accumbens_counts <- assays(nucleus_accumbens)$counts
nucleus_accumbens_sampleData <- colData(nucleus_accumbens)
write.csv(nucleus_accumbens_counts, file="nucleus_accumbens_counts.csv")
write.csv(nucleus_accumbens_sampleData, file="nucleus_accumbens_sampleData.csv")

putamen <- gtex_brain[,which(gtex_brain$smtsd == "Brain - Putamen (basal ganglia)")]
putamen_counts <- assays(putamen)$counts
putamen_sampleData <- colData(putamen)
write.csv(putamen_counts, file="putamen_counts.csv")
write.csv(putamen_sampleData, file="putamen_sampleData.csv")

substantia_nigra <- gtex_brain[,which(gtex_brain$smtsd == "Brain - Substantia nigra")]
substantia_nigra_counts <- assays(substantia_nigra)$counts
substantia_nigra_sampleData <- colData(substantia_nigra)
write.csv(substantia_nigra_counts, file="substantia_nigra_counts.csv")
write.csv(substantia_nigra_sampleData, file="substantia_nigra_sampleData.csv")

