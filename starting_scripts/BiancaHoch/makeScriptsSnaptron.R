###########################################
#
# Creating snap queries that will search for the putative exon region within brain tissue sample lists that were generated using 
# the script in RailtoRun_braintissuegroups.R.
#
###########################################
rm(list=ls()) #Clear R memory
gc() #empty garbage collector

### Load libraries:

#Source putative novel exon regions:
regionFile <- "~/Documents/splice_project/data/rawdata/PutativeNovelExons_hg38.csv"
regionDF <- as.data.frame(read.csv(regionFile))

# The query format contains the folllowing: a path to the snaptron ./qs script, the region of interest, the database being searched against, samples of interest in a list format.
# Example query:
#./qs --region "chr1:6145799-6145853" --datasrc gtex --samples '53969,54029,54069' --filters "strand=+" 1> amygdalatest.tsv

#Setting up the components of the query:
#Call qs script:
#snapDir <-"~/snaptron-experiments/qs"

#Region call:
regionCall <- "--region"

#List of putative exon regions:
regions <- regionDF[,2]

#List of putative exons for output filenames:
regionsNames <- gsub(":",".",regions)

#Database being searched against:
database <- "--datasrc gtex"

#Sample call:
sampleCall <- "--samples"

#Lists of samples which are specific to each tissue:
load("~/Documents/splice_project/data/analysis/acc_rail.Rdata")
acc_rail
load("~/Documents/splice_project/data/analysis/amygdala_rail.Rdata")
amygdala_rail
load("~/Documents/splice_project/data/analysis/caudate_rail.Rdata")
caudate_rail
load("~/Documents/splice_project/data/analysis/cerebellarh_rail.Rdata")
cerebellarh_rail
load("~/Documents/splice_project/data/analysis/cerebellum_rail.Rdata")
cerebellum_rail
load("~/Documents/splice_project/data/analysis/cortex_rail.Rdata")
cortex_rail
load("~/Documents/splice_project/data/analysis/frontalcortex_rail.Rdata")
frontalcortex_rail
load("~/Documents/splice_project/data/analysis/hippoc_rail.Rdata")
hippoc_rail
load("~/Documents/splice_project/data/analysis/hypo_rail.Rdata")
hypo_rail
load("~/Documents/splice_project/data/analysis/nucacc_rail.Rdata")
nucacc_rail
load("~/Documents/splice_project/data/analysis/putamen_rail.Rdata")
putamen_rail
load("~/Documents/splice_project/data/analysis/spinalc_rail.Rdata")
spinalc_rail
load("~/Documents/splice_project/data/analysis/subnigra_rail.Rdata")
subnigra_rail

#Filter call:
filterCall <- "--filters"

#Filtering results by strand specificity:
strandCall <- "strand="
strand <- regionDF[,3]
strands <- paste(strandCall,strand, sep = "")

#Save output to a file:
redirect <- "1>"

#Create a bash script header:
headLines <- paste("#!/bin/bash","set -e","set -o pipefail","set -u","#10/23/2019", sep = "\n")
cat(headLines)

###########TEST EXAMPLE#############
regions5 <- regions[1:5]
strands5 <- strands[1:5]
fileName <- paste("brainregion")
query <- paste("~/snaptron-experiments/qs", regionCall, regions5, database, sampleCall, acc_rail, filterCall, strands5, redirect, sep = " ")
brainregion_fileName <- paste("brainregion_", regions5, ".tsv", sep = "" )
queries2 <- paste(query,brainregion_fileName)
queries3 <- c(headLines,queries2)
#Save the list of queries as a .sh file:
write(queries3, file="~/Documents/splice_project/data/analysis/queries3.sh")
#####################################

#Create a list of queries similar to the example above for each tissue:
##acc_rail
acc_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, acc_rail, filterCall, strands, redirect, sep = " ")
acc_fileName <- paste("~/Documents/snaptronwork/snapout/acc_gtex/acc_", regionsNames, ".tsv", sep = "" )
acc_queries <- paste(acc_query,acc_fileName)
acc_SQ <- c(headLines,acc_queries)
#Save the list of queries as a .sh file:
write(acc_SQ, file="~/Documents/snaptronwork/snapgtex/acc_SQ.sh")

##amygdala_rail
amygdala_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, amygdala_rail, filterCall, strands, redirect, sep = " ")
amygdala_fileName <- paste("~/Documents/snaptronwork/snapout/amygdala_gtex/amygdala_", regionsNames, ".tsv", sep = "" )
amygdala_queries <- paste(amygdala_query,amygdala_fileName)
amygdala_SQ <- c(headLines,amygdala_queries)
#Save the list of queries as a .sh file:
write(amygdala_SQ, file="~/Documents/snaptronwork/snapgtex/amygdala_SQ.sh")

##caudate_rail
caudate_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, caudate_rail, filterCall, strands, redirect, sep = " ")
caudate_fileName <- paste("~/Documents/snaptronwork/snapout/caudate_gtex/caudate_", regionsNames, ".tsv", sep = "" )
caudate_queries <- paste(caudate_query,caudate_fileName)
caudate_SQ <- c(headLines,caudate_queries)
#Save the list of queries as a .sh file:
write(caudate_SQ, file="~/Documents/snaptronwork/snapgtex/caudate_SQ.sh")

##cerebellarh_rail
cerebellarh_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, cerebellarh_rail, filterCall, strands, redirect, sep = " ")
cerebellarh_fileName <- paste("~/Documents/snaptronwork/snapout/cerebellarh_gtex/cerebellarh_", regionsNames, ".tsv", sep = "" )
cerebellarh_queries <- paste(cerebellarh_query,cerebellarh_fileName)
cerebellarh_SQ <- c(headLines,cerebellarh_queries)
#Save the list of queries as a .sh file:
write(cerebellarh_SQ, file="~/Documents/snaptronwork/snapgtex/cerebellarh_SQ.sh")

##cerebellum_rail
cerebellum_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, cerebellum_rail, filterCall, strands, redirect, sep = " ")
cerebellum_fileName <- paste("~/Documents/snaptronwork/snapout/cerebellum_gtex/cerebellum_", regionsNames, ".tsv", sep = "" )
cerebellum_queries <- paste(cerebellum_query,cerebellum_fileName)
cerebellum_SQ <- c(headLines,cerebellum_queries)
#Save the list of queries as a .sh file:
write(cerebellum_SQ, file="~/Documents/snaptronwork/snapgtex/cerebellum_SQ.sh")

##cortex_rail
cortex_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, cortex_rail, filterCall, strands, redirect, sep = " ")
cortex_fileName <- paste("~/Documents/snaptronwork/snapout/cortex_gtex/cortex_", regionsNames, ".tsv", sep = "" )
cortex_queries <- paste(cortex_query,cortex_fileName)
cortex_SQ <- c(headLines,cortex_queries)
#Save the list of queries as a .sh file:
write(cortex_SQ, file="~/Documents/snaptronwork/snapgtex/cortex_SQ.sh")

##frontalcortex_rail
frontalcortex_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, frontalcortex_rail, filterCall, strands, redirect, sep = " ")
frontalcortex_fileName <- paste("~/Documents/snaptronwork/snapout/frontalcortex_gtex/frontalcortex_", regionsNames, ".tsv", sep = "" )
frontalcortex_queries <- paste(frontalcortex_query,frontalcortex_fileName)
frontalcortex_SQ <- c(headLines,frontalcortex_queries)
#Save the list of queries as a .sh file:
write(frontalcortex_SQ, file="~/Documents/snaptronwork/snapgtex/frontalcortex_SQ.sh")

##hippoc_rail
hippoc_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, hippoc_rail, filterCall, strands, redirect, sep = " ")
hippoc_fileName <- paste("~/Documents/snaptronwork/snapout/hippoc_gtex/hippoc_", regionsNames, ".tsv", sep = "" )
hippoc_queries <- paste(hippoc_query,hippoc_fileName)
hippoc_SQ <- c(headLines,hippoc_queries)
#Save the list of queries as a .sh file:
write(hippoc_SQ, file="~/Documents/snaptronwork/snapgtex/hippoc_SQ.sh")

##hypo_rail
hypo_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, hypo_rail, filterCall, strands, redirect, sep = " ")
hypo_fileName <- paste("~/Documents/snaptronwork/snapout/hypo_gtex/hypo_", regionsNames, ".tsv", sep = "" )
hypo_queries <- paste(hypo_query,hypo_fileName)
hypo_SQ <- c(headLines,hypo_queries)
#Save the list of queries as a .sh file:
write(hypo_SQ, file="~/Documents/snaptronwork/snapgtex/hypo_SQ.sh")

##nucacc_rail
nucacc_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, nucacc_rail, filterCall, strands, redirect, sep = " ")
nucacc_fileName <- paste("~/Documents/snaptronwork/snapout/nucacc_gtex/nucacc_", regionsNames, ".tsv", sep = "" )
nucacc_queries <- paste(nucacc_query,nucacc_fileName)
nucacc_SQ <- c(headLines,nucacc_queries)
#Save the list of queries as a .sh file:
write(nucacc_SQ, file="~/Documents/snaptronwork/snapgtex/nucacc_SQ.sh")

##putamen_rail
putamen_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, putamen_rail, filterCall, strands, redirect, sep = " ")
putamen_fileName <- paste("~/Documents/snaptronwork/snapout/putamen_gtex/putamen_", regionsNames, ".tsv", sep = "" )
putamen_queries <- paste(putamen_query,putamen_fileName)
putamen_SQ <- c(headLines,putamen_queries)
#Save the list of queries as a .sh file:
write(putamen_SQ, file="~/Documents/snaptronwork/snapgtex/putamen_SQ.sh")

##spinalc_rail
spinalc_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, spinalc_rail, filterCall, strands, redirect, sep = " ")
spinalc_fileName <- paste("~/Documents/snaptronwork/snapout/spinalc_gtex/spinalc_", regionsNames, ".tsv", sep = "" )
spinalc_queries <- paste(spinalc_query,spinalc_fileName)
spinalc_SQ <- c(headLines,spinalc_queries)
#Save the list of queries as a .sh file:
write(spinalc_SQ, file="~/Documents/snaptronwork/snapgtex/spinalc_SQ.sh")

##subnigra_rail
subnigra_query <- paste("~/snaptron-experiments/qs", regionCall, regions, database, sampleCall, subnigra_rail, filterCall, strands, redirect, sep = " ")
subnigra_fileName <- paste("~/Documents/snaptronwork/snapout/subnigra_gtex/subnigra_", regionsNames, ".tsv", sep = "" )
subnigra_queries <- paste(subnigra_query,subnigra_fileName)
subnigra_SQ <- c(headLines,subnigra_queries)
#Save the list of queries as a .sh file:
write(subnigra_SQ, file="~/Documents/snaptronwork/snapgtex/subnigra_SQ.sh")


