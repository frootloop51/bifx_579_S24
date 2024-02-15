###########################################
#
# Title: Loading and interpreting Snaptron Client 
#        output after querying junctions that overlap a genomic reagion of interest. 
#
###########################################

rm(list=ls()) #Clear R memory
gc() #empty garbage collector

#Load the file using red.delim
setwd("~/Users/biancahoch/Documents")
data <- ("./snaptest.csv")
data <- read.delim(file = "./snaptest.csv", encoding = "UTF-8")

#Defining each varilabe:
# [1] "DataSource.Type"   "snaptron_id"       "chromosome"        "start"            
# [5] "end"               "length"            "strand"            "annotated"        
# [9] "left_motif"        "right_motif"       "left_annotated"    "right_annotated"  
# [13] "samples"           "samples_count"     "coverage_sum"      "coverage_avg"     
# [17] "coverage_median"   "source_dataset_id"

#DataSource.Type = Differentiates between a return line of type Intron (I), Sample (S), or Gene (G). Also includes the database searched against.
#snaptron_id = unique identifier for the junction.
#chromosome = chromosome where the junction is located.
#start = start position of the intron
#end = end position of the intron
#length = length of the intron
#strand = (+/-)
#annotated = If both ends of the intron are annotated as a splice site in some annotation
#left_motif = Splice site sequence bases at the left end of the intron
#right_motif = Splice site sequence bases at the right end of the intron
#left_annotated = If the left end splice site is annotated or not and which annotations it appears in (maybe more than once)
#right_annotated = If the right end splice site is in an annotated or not, same as left_annotated
#samples = The list of samples which had one or more reads covering the intron and their coverages. IDs are from the IntropolisDB.
#samples_count = Total number of samples that have one or more reads covering this junction
#coverage_sum = Sum of all samples coverage for this junction
#coverage_avg = Average coverage across all samples which had at least 1 read covering the intron in the first pass alignment
#coverage_median = Median coverage across all samples which had at least 1 read covering the intron in the first pass alignment
#source_dataset_id = Snaptron ID for the compilation. GTEx=1, SRAv2=2, TCGA=4)


















