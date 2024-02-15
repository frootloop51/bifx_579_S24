
##############################################################
#### Author: Miranda Darby modification of script by Leonardo Collado Torres
#### Link to original script: http://www.bioconductor.org/packages/release/bioc/vignettes/recount/inst/doc/recount-quickstart.R
#### Recount2 quick start workflow
#### http://www.bioconductor.org/packages/release/bioc/vignettes/recount/inst/doc/recount-quickstart.html
#### Last modified: November 26, 2019
##############################################################

## Load recount R package
library('recount')

## ----'search_abstract'-----------------------------------------------------
## Find a project of interest
project_info <- abstract_search('GSE32465')

## Explore info
project_info

## ----'download'------------------------------------------------------------
## Download the gene-level RangedSummarizedExperiment data
download_study(project_info$project)

## Load the data
load(file.path(project_info$project, 'rse_gene.Rdata'))

## Delete it if you don't need it anymore
unlink(project_info$project, recursive = TRUE)

## ----'explore_rse'---------------------------------------------------------
rse_gene

## This is the sample phenotype data provided by the recount project
colData(rse_gene)

## At the gene level, the row data includes the gene Gencode ids, the gene
## symbols and the sum of the disjoint exons widths, which can be used for 
## taking into account the gene length.
rowData(rse_gene)

## At the exon level, you can get the gene Gencode ids from the names of:
# rowRanges(rse_exon)

## ----'browse'--------------------------------------------------------------
## Browse the project at SRA
browse_study(project_info$project)

## ----'sample_info', warning = FALSE----------------------------------------
## View GEO ids
colData(rse_gene)$geo_accession

## Extract the sample characteristics
geochar <- lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)

## Note that the information for this study is a little inconsistent, so we
## have to fix it.
geochar <- do.call(rbind, lapply(geochar, function(x) {
    if('cells' %in% colnames(x)) {
        colnames(x)[colnames(x) == 'cells'] <- 'cell.line'
        return(x)
    } else {
        return(x)
    }
}))

## We can now define some sample information to use
sample_info <- data.frame(
    run = colData(rse_gene)$run,
    group = ifelse(grepl('uninduced', colData(rse_gene)$title), 'uninduced', 'induced'),
    gene_target = sapply(colData(rse_gene)$title, function(x) { strsplit(strsplit(x,
        'targeting ')[[1]][2], ' gene')[[1]][1] }),
    cell.line = geochar$cell.line
)

## ----'scale_counts'--------------------------------------------------------
## Scale counts by taking into account the total coverage per sample
rse <- scale_counts(rse_gene)

##### Details about counts #####

## Scale counts to 40 million mapped reads. Not all mapped reads are in exonic
## sequence, so the total is not necessarily 40 million.
colSums(assays(rse)$counts) / 1e6

## Compute read counts
rse_read_counts <- read_counts(rse_gene)
## Difference between read counts and number of reads downloaded by Rail-RNA
colSums(assays(rse_read_counts)$counts) / 1e6 -
    colData(rse_read_counts)$reads_downloaded / 1e6
## Check the help page for read_counts() for more details

## ----'add_sample_info'-----------------------------------------------------
## Add sample information for DE analysis
colData(rse)$group <- sample_info$group
colData(rse)$gene_target <- sample_info$gene_target

## ----'de_analysis'---------------------------------------------------------
## Perform differential expression analysis with DESeq2
library('DESeq2')

## Specify design and switch to DESeq2 format
dds <- DESeqDataSet(rse, ~ gene_target + group)

## Perform DE analysis
dds <- DESeq(dds, test = 'LRT', reduced = ~ gene_target, fitType = 'local')
res <- results(dds)

## ----'maplot', fig.cap = "MA plot for group differences adjusted by gene target for SRP009615 using gene-level data."----
## Explore results
plotMA(res, main="DESeq2 results for SRP009615")

## ----'make_report', eval = FALSE-------------------------------------------
#  ## Make a report with the results
#  library('regionReport')
#  report <- DESeq2Report(dds, res = res, project = 'SRP009615',
#      intgroup = c('group', 'gene_target'), outdir = '.',
#      output = 'SRP009615-results', nBest = 10, nBestFeatures = 2)

## ----'make_report_real', echo = FALSE, results = 'hide'--------------------
library('regionReport')

## Make it so that the report will be available as a vignette
original <- readLines(system.file('DESeq2Exploration', 'DESeq2Exploration.Rmd',
    package = 'regionReport'))
vignetteInfo <- c(
    'vignette: >',
    '  %\\VignetteEngine{knitr::rmarkdown}',
    '  %\\VignetteIndexEntry{Basic DESeq2 results exploration}',
    '  %\\VignetteEncoding{UTF-8}'
)
new <- c(original[1:16], vignetteInfo, original[17:length(original)])
writeLines(new, 'SRP009615-results-template.Rmd')

## Now create the report
report <- DESeq2Report(dds, res = res, project = 'SRP009615',
    intgroup = c('group', 'gene_target'), outdir = '.',
    output = 'SRP009615-results', device = 'png', template = 'SRP009615-results-template.Rmd', nBest = 10, nBestFeatures = 2)
    
## Clean up
file.remove('SRP009615-results-template.Rmd')

## ----'geneSymbols'---------------------------------------------------------
## Load required library
library('org.Hs.eg.db')

## Extract Gencode gene ids
gencode <- gsub('\\..*', '', names(recount_genes))

## Find the gene information we are interested in
gene_info <- select(org.Hs.eg.db, gencode, c('ENTREZID', 'GENENAME', 'SYMBOL',
    'ENSEMBL'), 'ENSEMBL')

## Explore part of the results
dim(gene_info)
head(gene_info)

## ----'define_ers', eval = .Platform$OS.type != 'windows'-------------------
## Define expressed regions for study SRP009615, only for chromosome Y
regions <- expressed_regions('SRP009615', 'chrY', cutoff = 5L, 
    maxClusterGap = 3000L)

## Briefly explore the resulting regions
regions

## ----'compute_covMat', eval = .Platform$OS.type != 'windows'---------------
## Compute coverage matrix for study SRP009615, only for chromosome Y
system.time( rse_ER <- coverage_matrix('SRP009615', 'chrY', regions) )

## Explore the RSE a bit
dim(rse_ER)
rse_ER

## ----'to_integer', eval = .Platform$OS.type != 'windows'-------------------
## Round the coverage matrix to integers
covMat <- round(assays(rse_ER)$counts, 0)

## ----'phenoData', eval = .Platform$OS.type != 'windows'--------------------
## Get phenotype data for study SRP009615
pheno <- colData(rse_ER)

## Complete the phenotype table with the data we got from GEO
m <- match(pheno$run, sample_info$run)
pheno <- cbind(pheno, sample_info[m, 2:3])

## Explore the phenotype data a little bit
head(pheno)

## ----'ers_dds', eval = .Platform$OS.type != 'windows'----------------------
## Build a DESeqDataSet
dds_ers <- DESeqDataSetFromMatrix(countData = covMat, colData = pheno,
    design =  ~ gene_target + group)

## ----'de_analysis_ers', eval = .Platform$OS.type != 'windows'--------------
## Perform differential expression analysis with DESeq2 at the ER-level
dds_ers <- DESeq(dds_ers, test = 'LRT', reduced = ~ gene_target,
    fitType = 'local')
res_ers <- results(dds_ers)

## ----'maploters', eval = .Platform$OS.type != 'windows', fig.cap = "MA plot for group differences adjusted by gene target for SRP009615 using ER-level data (just chrY)."----
## Explore results
plotMA(res_ers, main="DESeq2 results for SRP009615 (ER-level, chrY)")

## ----'report2', eval = FALSE-----------------------------------------------
#  ## Create the report
#  report2 <- DESeq2Report(dds_ers, res = res_ers,
#      project = 'SRP009615 (ER-level, chrY)',
#      intgroup = c('group', 'gene_target'), outdir = '.',
#      output = 'SRP009615-results-ER-level-chrY')

## ---- 'gene_annotation'----------------------------------------------------
## Gene annotation information
rowRanges(rse_gene_SRP009615)

## Also accessible via
rowData(rse_gene_SRP009615)

## ---- 'exon_annotation'----------------------------------------------------
## Get the rse_exon object for study SRP009615
download_study('SRP009615', type = 'rse-exon')
load(file.path('SRP009615', 'rse_exon.Rdata'))

## Delete it if you don't need it anymore
unlink('SRP009615', recursive = TRUE)

## Annotation information
rowRanges(rse_exon)

## Add a gene_id column
rowRanges(rse_exon)$gene_id <- rownames(rse_exon)

## Example subset
rse_exon_subset <- subset(rse_exon, subset = gene_id == 'ENSG00000000003.14')
rowRanges(rse_exon_subset)

## ---- eval = .Platform$OS.type != 'windows'--------------------------------
## Get the disjoint exons based on EnsDb.Hsapiens.v79 which matches hg38
exons <- reproduce_ranges('exon', db = 'EnsDb.Hsapiens.v79')

## Change the chromosome names to match those used in the BigWig files
library('GenomeInfoDb')
seqlevelsStyle(exons) <- 'UCSC'

## Get the count matrix at the exon level for chrY
exons_chrY <- keepSeqlevels(exons, 'chrY', pruning.mode = 'coarse')
exonRSE <- coverage_matrix('SRP009615', 'chrY', unlist(exons_chrY),
    chunksize = 3000)
exonMatrix <- assays(exonRSE)$counts
dim(exonMatrix)
head(exonMatrix)

## Summary the information at the gene level for chrY
exon_gene <- rep(names(exons_chrY), elementNROWS(exons_chrY))
geneMatrix <- do.call(rbind, lapply(split(as.data.frame(exonMatrix),
    exon_gene), colSums))
dim(geneMatrix)
head(geneMatrix)

## ----'gene_fusions'--------------------------------------------------------
library('recount')

## Download and load RSE object for SRP009615 study
download_study('SRP009615', type = 'rse-jx')
load(file.path('SRP009615', 'rse_jx.Rdata'))

## Delete it if you don't need it anymore
unlink('SRP009615', recursive = TRUE)

## Exon-exon junctions by class
table(rowRanges(rse_jx)$class)

## Potential gene fusions by chromosome
fusions <- rowRanges(rse_jx)[rowRanges(rse_jx)$class == 'fusion']
fusions_by_chr <- sort(table(seqnames(fusions)), decreasing = TRUE)
fusions_by_chr[fusions_by_chr > 0]

## Genes with the most fusions
head(sort(table(unlist(fusions$symbol_proposed)), decreasing = TRUE))

## ----snaptron--------------------------------------------------------------
library('GenomicRanges')
junctions <- GRanges(seqnames = 'chr2', IRanges(
    start = c(28971711, 29555082, 29754983),
    end = c(29462418, 29923339, 29917715)))

snaptron_query(junctions)

## ----'rse-fc example'------------------------------------------------------
## Example use of download_study()
## for downloading the FANTOM-CAT `rse_fc`
## file for a given study
download_study('DRP000366', type = 'rse-fc', download = FALSE)

## ----'recount-brain'-------------------------------------------------------
## recount-brain citation details
citation('recount')[5]

## ----'downloadAll'---------------------------------------------------------
## Get data.frame with all the URLs
library('recount')
dim(recount_url)

## Explore URLs
head(recount_url$url)

## ----'actuallyDownload', eval = FALSE--------------------------------------
#  ## Download all files (will take a while! It's over 8 TB of data)
#  sapply(unique(recount_url$project), download_study, type = 'all')

## ----Figure1, out.width="100%", fig.align="center", fig.cap = "SIgn up to SciServer", echo = FALSE----
knitr::include_graphics("https://raw.githubusercontent.com/leekgroup/recount-website/master/sciserver_demo/register.png")

## ----Figure2, out.width="100%", fig.align="center", fig.cap = "Click on SciServer Compute,", echo = FALSE----
knitr::include_graphics("https://raw.githubusercontent.com/leekgroup/recount-website/master/sciserver_demo/sciserver-compute.png")

## ----Figure3, out.width="100%", fig.align="center", fig.cap = "Select the appropriate image and load the recount public volume.", echo = FALSE----
knitr::include_graphics("https://raw.githubusercontent.com/leekgroup/recount-website/master/sciserver_demo/new-container.png")

## ----Figure4, out.width="100%", fig.align="center", fig.cap = "Create a R notebook.", echo = FALSE----
knitr::include_graphics("https://raw.githubusercontent.com/leekgroup/recount-website/master/sciserver_demo/new-R.png")

## ----Figure5, out.width="100%", fig.align="center", fig.cap = "Install recount and DESeq2 on SciServer.", echo = FALSE----
knitr::include_graphics("https://raw.githubusercontent.com/leekgroup/recount-website/master/sciserver_demo/demo1.png")

## ----Figure6, out.width="100%", fig.align="center", fig.cap = "recount files are available locally so remmeber to use the local data when working on SciServer.", echo = FALSE----
knitr::include_graphics("https://raw.githubusercontent.com/leekgroup/recount-website/master/sciserver_demo/demo2.png")

## ----'ER-sciserver', eval = FALSE------------------------------------------
#  ## Expressed-regions SciServer example
#  regions <- expressed_regions('SRP009615', 'chrY', cutoff = 5L,
#      maxClusterGap = 3000L, outdir = file.path(scipath, 'SRP009615'))
#  regions

## ----createVignette, eval=FALSE--------------------------------------------
#  ## Create the vignette
#  library('rmarkdown')
#  system.time(render('recount-quickstart.Rmd', 'BiocStyle::html_document'))
#  
#  ## Extract the R code
#  library('knitr')
#  knit('recount-quickstart.Rmd', tangle = TRUE)

## ----createVignette2-------------------------------------------------------
## Clean up
file.remove('quickstartRef.bib')

## ----reproduce1, echo=FALSE------------------------------------------------
## Date the vignette was generated
Sys.time()

## ----reproduce2, echo=FALSE------------------------------------------------
## Processing time in seconds
totalTime <- diff(c(startTime, Sys.time()))
round(totalTime, digits=3)

## ----reproduce3, echo=FALSE-------------------------------------------------------------------------------------------
## Session info
library('sessioninfo')
options(width = 120)
session_info()

## ----'datasetup', echo = FALSE, eval = FALSE--------------------------------------------------------------------------
#  ## Code for re-creating the data distributed in this package
#  
#  ## Genes/exons
#  library('recount')
#  recount_both <- reproduce_ranges('both', db = 'Gencode.v25')
#  recount_exons <- recount_both$exon
#  save(recount_exons, file = '../data/recount_exons.RData', compress = 'xz')
#  recount_genes <- recount_both$gene
#  save(recount_genes, file = '../data/recount_genes.RData', compress = 'xz')
#  
#  ## URL table
#  load('../../recount-website/fileinfo/upload_table.Rdata')
#  recount_url <- upload_table
#  ## Create URLs
#  is.bw <- grepl('[.]bw$', recount_url$file_name)
#  recount_url$url <- NA
#  recount_url$url[!is.bw] <-paste0('http://duffel.rail.bio/recount/',
#      recount_url$project[!is.bw], '/', recount_url$file_name[!is.bw])
#  recount_url$url[!is.bw & recount_url$version2] <-
#      paste0('http://duffel.rail.bio/recount/v2/',
#      recount_url$project[!is.bw & recount_url$version2], '/',
#      recount_url$file_name[!is.bw & recount_url$version2])
#  recount_url$url[is.bw] <- paste0('http://duffel.rail.bio/recount/',
#      recount_url$project[is.bw], '/bw/', recount_url$file_name[is.bw])
#  save(recount_url, file = '../data/recount_url.RData', compress = 'xz')
#  
#  ## Add FC-RC2 data
#  load('../data/recount_url.RData')
#  load('../../recount-website/fileinfo/fc_rc_url_table.Rdata')
#  recount_url <- rbind(recount_url, fc_rc_url_table)
#  rownames(recount_url) <- NULL
#  save(recount_url, file = '../data/recount_url.RData', compress = 'xz')
#  
#  
#  ## Abstract info
#  load('../../recount-website/website/meta_web.Rdata')
#  recount_abstract <- meta_web[, c('number_samples', 'species', 'abstract')]
#  recount_abstract$project <- gsub('.*">|</a>', '', meta_web$accession)
#  Encoding(recount_abstract$abstract) <- 'latin1'
#  recount_abstract$abstract <- iconv(recount_abstract$abstract, 'latin1', 'UTF-8')
#  save(recount_abstract, file = '../data/recount_abstract.RData',
#      compress = 'bzip2')
#  
#  ## Example rse_gene file
#  system('scp e:/dcl01/leek/data/recount-website/rse/rse_sra/SRP009615/rse_gene.Rdata .')
#  load('rse_gene.Rdata')
#  rse_gene_SRP009615 <- rse_gene
#  save(rse_gene_SRP009615, file = '../data/rse_gene_SRP009615.RData',
#      compress = 'xz')
#  unlink('rse_gene.Rdata')

## ----'cleanup', echo = FALSE, results = 'hide'------------------------------------------------------------------------
## Clean up
unlink('SRP009615', recursive = TRUE)

## ----vignetteBiblio, results = 'asis', echo = FALSE, warning = FALSE, message = FALSE---------------------------------
## Print bibliography
bibliography()

