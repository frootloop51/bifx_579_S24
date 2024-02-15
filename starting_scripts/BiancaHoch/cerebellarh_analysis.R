###########################################
#
# Read and Analyze Snaptron output files from the cerebellar hemisphere queries made by makeScriptsSnaptron.R.
#
###########################################
rm(list=ls()) #Clear R memory
gc() #empty garbage collector

###Set working directory:
setwd("/Users/biancahoch/Documents/snaptronwork/snapout/cerebellarh_gtex/")

#Update file path to point toward appropriate folder on the computer
cerebellarh_list <- list.files(pattern=".tsv") 

#Read files:
cerebellarh <- sapply(cerebellarh_list, read.delim, simplify = FALSE)

#Print only the rows in which the number of samples that have one or more reads covering the intron constitute atleast 5% of the total samples searched:

#First find the number of anterior cingulate cortex samples:
load("~/Documents/splice_project/data/analysis/cerebellarh_rail.Rdata")
cerebellarh_rail
splitcerebellarh <- strsplit(cerebellarh_rail, ",")
unlistsplitcerebellarh <- unlist(splitcerebellarh)
print(paste(length(unlistsplitcerebellarh)))

#See which junctions are present in atleast 5% of samples:
f <- function(df){
  (df[which(df$samples_count >= length(unlistsplitcerebellarh) *.05), ])
}

#Apply function to list of dataframes:
dfList <- lapply( X=cerebellarh, FUN=f ) 
#Bind the lists into one dataframe:
allDF <-do.call(rbind,dfList) 
#Write the dataframe to a .csv:
write.csv(allDF, "~/Documents/snaptronwork/snapout_analysis/cerebellarh_fivepercent.csv")
#Check output:
#cerebellarh_fivepercent <- read.csv("~/Documents/snaptronwork/snapout_analysis/cerebellarh_fivepercent.csv", header = TRUE, sep = ",")

#See which junctions are present in atleast 60% of samples:
f2 <- function(df){
  (df[which(df$samples_count >= length(unlistsplitcerebellarh) *.60), ])
}
#Apply function to list of dataframes:
dfList2 <- lapply( X=cerebellarh, FUN=f2 ) 
#Bind the lists into one dataframe:
allDF2 <-do.call(rbind,dfList2) 
#Write the dataframe to a .csv: (NOTE: ADD ROWNAMES)
write.csv(allDF2, "~/Documents/snaptronwork/snapout_analysis/cerebellarh_sixtypercent.csv")
#Check output:
#cerebellarh_sixtypercent <- read.csv("~/Documents/snaptronwork/snapout_analysis/cerebellarh_sixtypercent.csv", header = TRUE, sep = ",")

#############################################################
