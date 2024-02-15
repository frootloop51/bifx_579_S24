###########################################
#
# A Test for Reading and Analyzing Snaptron output files from the queries made by makeScriptsSnaptron.R.
# Finalized methods will be used in all scripts that end with *_analysis.R.
#
# 11/17/2019
###########################################
rm(list=ls()) #Clear R memory
gc() #empty garbage collector

###Set working directory:
setwd("/Users/biancahoch/Documents/snaptronwork/snapout/acc_gtex/")

#Update file path to point toward appropriate folder on the computer
acc_list <- list.files(pattern=".tsv") 

########Starting my analysis on just 5 files as a test########

##Create list of acc filenames, this will be used to verify if the start and end splice sites are present:
acc_list5 <- as.character(acc_list[1:5])

#Read in only 5 files:
acc5 <- sapply(acc_list5, read.delim, simplify = FALSE)

#############################################################
#Checking to see which junctions match the start and end positions of the original 3,911 regions that were queried.

#Remove the first portion of the filenames that indicate tissue:
strsplit1 <- sapply(strsplit(acc_list5,"_"),"[[",2)
strsplit2 <- sapply(strsplit(strsplit1,"[.]"),"[[",2)
RELstart <- sapply(strsplit(strsplit2,"-"),"[[",1)
RELend <- sapply(strsplit(strsplit2,"-"),"[[",2)

#Add a start and end which I know will be found within the first row of the first database in acc5.
#This will confirm whether or not the function can pull a row in which the start and end positions match those of our putative exons.
known1 <- append("9643277", RELstart[2:5])
known2 <- append("10094933", RELend[2:5])

#Begin creating the loop:
#Ensure that the indices of the dataframes and the putative exon start/end region coordinates are the same:
index <- 1

#Search too see if each of the 5 putative exons in acc_list5 can be found within their respective dataframes (acc5):
for (df in acc5){
  if (known1[index] %in% df$start & known2[index] %in% df$end) {
    print(df[which(df$start == known1[index] & df$end == known2[index]), ]) 
  }
  index <- index + 1
}

#############################################################
#Print only the rows in which the number of samples that have one or more reads covering the intron constitute atleast 5% of the total samples searched:

#First find the number of anterior cingulate cortex samples:
load("~/Documents/splice_project/data/analysis/acc_rail.Rdata")
acc_rail
splitacc <- strsplit(acc_rail, ",")
unlistsplitacc <- unlist(splitacc)
print(paste(length(unlistsplitacc)))
#[1] 99

#See which junctions are present in atleast 5% of samples:
f <- function(df){
  (df[which(df$samples_count >= length(unlistsplitacc) *.05), ])
}

#Apply function to list of dataframes:
dfList <- lapply( X=acc5, FUN=f ) 
#Bind the lists into one dataframe:
allDF <-do.call(rbind,dfList) 
#Write the dataframe to a .csv: (NOTE: ADD ROWNAMES)
write.csv(allDF, "~/Documents/snaptronwork/acc5_fivepercent.csv")
#Check output:
acc5_fivepercent <- read.csv("~/Documents/snaptronwork/acc5_fivepercent.csv", header = TRUE, sep = ",")

#############################################################

#See which junctions are present in atleast 60% of samples or are within the 3rd quantile of coverage average:
# f2 <- function(df){
#   (df[which(df$samples_count >= length(unlistsplitacc) *.60) | which(df$coverage_avg >= quantile(df$coverage_avg, 0.75)), ])
# }
## Think about how to feed in variable for naming the csvs to match the tissue type and what other variables might be needed such as being able to feed in thresholds
## Instead of using .60, use variable name and then set the variable as an argument in the function.

f2 <- function(df, var=threshold, tissue etc.){
  sc <- which(df$samples_count >= length(unlistsplitacc) *threshold)
  cov <- which(df$coverage_avg >= quantile(df$coverage_avg, 0.75))
  scDF <- df[sc,]
  write.csv(scDF, "")
  #do same for cov
  combined <- c(sc, cov)
  all <- unique(combined)
  all <- df[all,]
  write.csv(all, "~/putName")
  both <- duplicated(combined)
  both <- df[both,] # and save
}

#Apply function to list of dataframes:
dfList2 <- lapply( X=acc5, FUN=f2 ) 
#Bind the lists into one dataframe:
allDF2 <-do.call(rbind,dfList2) 
#Write the dataframe to a .csv: (NOTE: ADD ROWNAMES)
write.csv(allDF2, "~/Documents/snaptronwork/acc5_sixtypercent.csv")
#Check output:
acc5_sixtypercent <- read.csv("~/Documents/snaptronwork/acc5_sixtypercent.csv", header = TRUE, sep = ",")
%in%


