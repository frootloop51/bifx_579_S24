#PLOT ALL CONDITIONS

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ggbio)
data(genesymbol, package = "biovizBase")

#load("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/flippedStrandControlRanges.rda")
#load("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/flippedStrandBipolarRanges.rda")
#load("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/flippedStrandDepressionRanges.rda")
#load("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/flippedStrandSchizophreniaRanges.rda")

setwd("/dcs01/stanley/work/miranda/human/orbFrontal_Repeats/orbFront_wGT1_SarvenRep/CoveragePlots")
load("combinedGRwAnnotation.rda")
#Plot genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
transcripts <- transcripts(txdb)

#Find all reads in GRanges that overlap region
findReads <- function(GR,range){
  ovrlp <- queryHits(findOverlaps(GR, range, ignore.strand=TRUE))
  ovrlp <- GR[ovrlp]
}

ggbiosave <- function (filename = default_name(plot), plot = last_plot(),
                       device = default_device(filename), path = NULL,
                       scale = 1,
                       width = par("din")[1], height = par("din")[2],
                       units = c("in",
                                 
                                 "cm", "mm"), dpi = 300, limitsize = TRUE, ...)
{
  ## simply comment out the check part
  ## if (!inherits(plot, "ggplot"))
  ##     stop("plot should be a ggplot2 plot")
  eps <- ps <- function(..., width, height) grDevices::postscript(...,
                                                                  width = width, height = height, onefile = FALSE, horizontal =
                                                                    FALSE,
                                                                  paper = "special")
  tex <- function(..., width, height) grDevices::pictex(...,
                                                        width = width, height = height)
  pdf <- function(..., version = "1.4") grDevices::pdf(...,
                                                       version = version)
  svg <- function(...) grDevices::svg(...)
  wmf <- function(..., width, height) grDevices::win.metafile(...,
                                                              width = width, height = height)
  emf <- function(..., width, height) grDevices::win.metafile(...,
                                                              width = width, height = height)
  png <- function(..., width, height) grDevices::png(..., width = width,
                                                     height = height, res = dpi, units = "in")
  jpg <- jpeg <- function(..., width, height) grDevices::jpeg(...,
                                                              width = width, height = height, res = dpi, units = "in")
  bmp <- function(..., width, height) grDevices::bmp(..., width = width,
                                                     height = height, res = dpi, units = "in")
  tiff <- function(..., width, height) grDevices::tiff(...,
                                                       width = width, height = height, res = dpi, units = "in")
  default_name <- function(plot) {
    paste(digest.ggplot(plot), ".pdf", sep = "")
  }
  default_device <- function(filename) {
    pieces <- strsplit(filename, "\\.")[[1]]
    ext <- tolower(pieces[length(pieces)])
    match.fun(ext)
  }
  units <- match.arg(units)
  convert_to_inches <- function(x, units) {
    x <- switch(units, `in` = x, cm = x/2.54, mm = x/2.54/10)
  }
  convert_from_inches <- function(x, units) {
    x <- switch(units, `in` = x, cm = x * 2.54, mm = x *
                  2.54 * 10)
  }
  if (!missing(width)) {
    width <- convert_to_inches(width, units)
  }
  if (!missing(height)) {
    height <- convert_to_inches(height, units)
  }
  if (missing(width) || missing(height)) {
    message("Saving ", prettyNum(convert_from_inches(width *
                                                       scale, units), digits = 3), " x ",
            prettyNum(convert_from_inches(height *
                                            scale, units), digits = 3), " ", units, " image")
  }
  width <- width * scale
  height <- height * scale
  if (limitsize && (width >= 50 || height >= 50)) {
    stop("Dimensions exceed 50 inches (height and width are specified
         in inches/cm/mm, not pixels).",
         " If you are sure you want these dimensions, use
         'limitsize=FALSE'.")
  }
  if (!is.null(path)) {
    filename <- file.path(path, filename)
  }
  device(file = filename, width = width, height = height, ...)
  on.exit(capture.output(dev.off()))
  print(plot)
  invisible()
}

plotRegion <- function(GR,chr,start,end){
  wh <- GRanges(seqnames = chr, ranges = IRanges(start = start, end =end), strand = "*")
  region <-paste(chr,":",start,"-",end,sep="")
  reads <- findReads(GR=GR,range=wh)
  readsPlot <- autoplot(reads, stat = "coverage", geom = "area", facets = strand~seqnames, aes(fill = strand))
  fileName <- paste(region,".pdf",sep="")
  ggbiosave(filename=fileName, plot=readsPlot)
}
findReadsSS <- function(GR,range){
  ovrlp <- queryHits(findOverlaps(GR, range, ignore.strand=FALSE))
  ovrlp <- GR[ovrlp]
}

plotRegionSS4GR <- function(chr,start,end,strand){
  wh <- GRanges(seqnames = chr, ranges = IRanges(start = start, end =end), strand = strand)
  region <-paste(chr,":",start,"-",end,sep="")
  cntrl <- findReadsSS(GR=flpCntrlGr,range=wh)
  bp <- findReadsSS(GR=flpBpGr,range=wh)
  scz <- findReadsSS(GR=flpSczGr,range=wh)
  dep <- findReadsSS(GR=flpDepGr,range=wh)
  all <- c(cntrl,bp,scz,dep)
  p1 <- autoplot(all, stat = "coverage", geom = "area",facets = value~seqnames)
  fileName <- paste(region,"_AllConditionsSS_facet.pdf",sep="")
  ggbiosave(filename=fileName, tracks(reads=p1))
}
# plotRegion4GR(chr="chr2",start=157183057, end=157183223,strand="-")

findReadsDS <- function(GR,range){
  ovrlp <- queryHits(findOverlaps(GR, range, ignore.strand=TRUE))
  ovrlp <- GR[ovrlp]
}

plotRegionDS4GR <- function(chr,start,end,strand){
  wh <- GRanges(seqnames = chr, ranges = IRanges(start = start, end =end), strand = strand)
  region <-paste(chr,":",start,"-",end,sep="")
  cntrl <- findReadsDS(GR=flpCntrlGr,range=wh)
  bp <- findReadsDS(GR=flpBpGr,range=wh)
  scz <- findReadsDS(GR=flpSczGr,range=wh)
  dep <- findReadsDS(GR=flpDepGr,range=wh)
  all <- c(cntrl,bp,scz,dep)
  p1 <- autoplot(all, stat = "coverage", geom = "area",facets = value~seqnames)
  fileName <- paste(region,"_AllConditionsDS_facet.pdf",sep="")
  ggbiosave(filename=fileName, tracks(reads=p1))
}

plotGene4GR <- function(gene){
  wh <- genesymbol[gene]
  trans <- ggplot() + geom_alignment(txdb, which = wh)
  ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
  cntrl <- findReadsSS(GR=flpCntrlGr,range=wh)
  bp <- findReadsSS(GR=flpBpGr,range=wh)
  scz <- findReadsSS(GR=flpSczGr,range=wh)
  dep <- findReadsSS(GR=flpDepGr,range=wh)
  all <- c(cntrl,bp,scz,dep)
  p1 <- autoplot(all, stat = "coverage", geom = "area",facets = value~seqnames)
  fileName <- paste(gene,"_orbFront_facet.pdf",sep="")
  ggbiosave(filename=fileName, tracks(reads=p1, transcripts=trans, heights=c(7,1.5)))
}

plotRegion4GRwTrans <- function(chr,start,end,strand){
  wh <- GRanges(seqnames = chr, ranges = IRanges(start = start, end =end), strand = strand)
  trans <- ggplot() + geom_alignment(txdb, which = wh)
  ggplot(txdb) + geom_alignment(which = wh, names.expr = "tx_id:::gene_id")
  cntrl <- findReadsSS(GR=flpCntrlGr,range=wh)
  bp <- findReadsSS(GR=flpBpGr,range=wh)
  scz <- findReadsSS(GR=flpSczGr,range=wh)
  dep <- findReadsSS(GR=flpDepGr,range=wh)
  all <- c(cntrl,bp,scz,dep)
  p1 <- autoplot(all, stat = "coverage", geom = "area",facets = value~seqnames)
  region <-paste(chr,":",start,"-",end,sep="")
  fileName <- paste(region,"_AllConditions_facet.pdf",sep="")
  ggbiosave(filename=fileName, tracks(reads=p1, transcripts=trans, heights=c(7,1.5)))
}
