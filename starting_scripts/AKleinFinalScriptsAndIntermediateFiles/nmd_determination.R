######################################################

# Determination of NMD Candidates
# Author: Alyssa Klein
# Date: February 2, 2020
# Modified: February 12, 2020

#####################################################

# This script takes in the amino acid sequences for both the canonical protein sequence and the novel 
# isoform protein sequence, and uses this information to determine if the novel isoform is a potential
# NMD candidate or not. An output text file is generated with the information for each protein, named: 
# "nmd_information_output.txt". There is also a small dataframe generated that includes the information 
# in the initial file read into the script, along with a new column for NMD candidacy.  
# This csv file is named: "intron_positions_for_nmd_with_nmd_candidacy.csv".

#####################################################

# Read in the file that contains the necessary information to determine if a novel isoform is a NMD candidate or not.
# This file contains only the proteins whose isoforms introduced a premature termination codon (PTC)
# See the README file for more information as to what this file needs to contain and where the information in this file comes from.
file_with_intron_position <- read.csv("intron_positions_for_nmd.csv", header=TRUE, sep= ",")

# the files being read in contain two sequences- 1) the canonical sequence for the protein and 2) the novel isoform sequence
# each of the files are named with the following connotation: "Protein_name_Seq_ID_can_iso_seqs.fasta"
# use the paste0() function to store all of the filenames for the files that need to be iterated through in the for loop
name_file_paste <- paste0(file_with_intron_position$Protein_Name, "_",file_with_intron_position$Seq_ID, "_can_iso_seqs.fasta")

# set i equal to 1 to begin the for loop in the following lines
i <- 1
for (files in name_file_paste){

  # Extract the necessary information for determining potential NMD candidates
  # extract the correct strand for each file imported using the base R 'which' function
  strand <- file_with_intron_position$Strand[which(file_with_intron_position$Protein_Name == file_with_intron_position$Protein_Name[[i]])]
  # extract the correct intron position for each file imported using the base R 'which' function
  intron_position <- file_with_intron_position$Final_intron_position[which(file_with_intron_position$Protein_Name == file_with_intron_position$Protein_Name[[i]])]
  # paste0 function to print the protein name and its intron position so it can be displayed to the user
  print(paste0(file_with_intron_position$Protein_Name[[i]], ": ", intron_position))
  # read in each fasta file using the 'read.fasta' function from the seqinr library
  library(seqinr)
  fasta_file <- read.fasta((paste0(file_with_intron_position$Protein_Name[[i]], "_",file_with_intron_position$Seq_ID[[i]], "_can_iso_seqs.fasta")), as.string  =TRUE, forceDNAtolower = FALSE, seqonly = TRUE)
  # print(fasta_file)
  # the following line is an example of a filename for a file being read in
  # fasta_file <- read.fasta("AAMDC_Q9H7C9_can_iso_seqs.fasta", as.string  =TRUE, forceDNAtolower = FALSE, seqonly = TRUE)
  
  # cat function to append the filename to the output text file
  cat(name_file_paste[[i]],file= "nmd_information_output.txt",append=TRUE, sep = "\n")
  # cat function to append the Protein name, intron position, and strand to the output text file
  cat(paste0(file_with_intron_position$Protein_Name[[i]], " Intron position: ", intron_position, " Strand: ", strand),file= "nmd_information_output.txt",append=TRUE, sep = "\n")
  
  # store the canonical sequence from the fasta file as canonical_seq object
  canonical_seq <- fasta_file[[1]]
  
  # store the isoform sequence from the fasta file as isoform_seq object
  isoform_seq <- fasta_file[[2]]
  
  # calculate the number of nucleotides in the canonical sequence (just to observe) by multiplying the number of amino acids by three
  nuc_canonical <- nchar(canonical_seq)*3
  # print the number of nucleotides in the canonical sequence 
  print(nuc_canonical)
  # append the number of nucleotides in the canonical sequence to the output file
  cat(paste0("Nucleotides in canonical sequence: ", nuc_canonical),file= "nmd_information_output.txt",append=TRUE, sep = "\n")
  
  # calculate the number of nucleotides in the isoform sequence and subtract 2 to obtain the base pair number of the start of the stop codon
  nuc_isoform <- nchar(isoform_seq)*3 - 2 
  # print the number of nucleotides in the isoform sequence 
  print(nuc_isoform)
  # append the number of nucleotides in the isoform sequence to the output file
  cat(paste0("Nucleotides in isoform sequence: ", nuc_canonical),file= "nmd_information_output.txt",append=TRUE, sep = "\n")
  
  # make new loop that only does all of the calculations when the PTC is upstream of the final exon-exon junction (intron position) 2/5/2020
  # subtract the final intron position (i.e. final exon-exon junction) from the number of nucleotides in the isoform
  # check to tsee if the intron_position is equal to "NA"
  if (identical(toString(intron_position), "NA") == TRUE){
    # set nmd_candidate variable to "NA"
    nmd_candidate <- "NA"
    # create a new column in the dataframe to store the result of a novel isoform being a nmd candidate or not (or "NA")
    file_with_intron_position$NMD_candidate[[i]] <- nmd_candidate
    # print this information to the console
    print(paste0(file_with_intron_position$Protein_Name[[i]], " did not have junction information available at this time. The potential for this isoform to be an NMD candidate is TBD."))
    # append this information to the output text file
    cat(paste0(file_with_intron_position$Protein_Name[[i]], " did not have junction information available at this time. The potential for this isoform to be an NMD candidate is TBD."), file= "nmd_information_output.txt",append=TRUE, sep = "\n")
  }else{
    # calculation carried out if the intron position is downstream of the end of the isoform nucleotide sequence
    if(intron_position >= nuc_isoform){
      # take the absolute value of the intron position to account for both the + and - strand possibilities
      intron_pos_minus_iso_nt <- abs(intron_position) - nuc_isoform
      if(intron_pos_minus_iso_nt > 55){
        # set nmd_candidate variable to "yes"
        nmd_candidate <- "yes"
        # create a new column in the dataframe to store the result of a novel isoform being a nmd candidate or not (or "NA")
        file_with_intron_position$NMD_candidate[[i]] <- nmd_candidate
        # print this information to the console
        print(paste0("The premature stop codon (PTC) is greater than 55 nucleotides away from the final exon-exon junction. The PTC is ", 
                     intron_pos_minus_iso_nt, " nucleotides away, and therefore is a potential NMD candidate."))
        # append this information to the output text file
        cat(paste0("The premature stop codon (PTC) is greater than 55 nucleotides away from the final exon-exon junction. The PTC is ", 
                    intron_pos_minus_iso_nt, " nucleotides away, and therefore is a potential NMD candidate."), file= "nmd_information_output.txt",append=TRUE, sep = "\n")
      }else{
        # set nmd_candidate variable to "no"
        nmd_candidate <- "no"
        # create a new column in the dataframe to store the result of a novel isoform being a nmd candidate or not (or "NA")
        file_with_intron_position$NMD_candidate[[i]] <- nmd_candidate
        # print this information to the console
        print(paste0("The premature stop codon (PTC) is less than 55 nucleotides away from the final exon-exon junction. The PTC is ", 
                     intron_pos_minus_iso_nt, " nucleotides away, and therefore is not a potential NMD candidate."))
        # append this information to the output text file
        cat(paste0("The premature stop codon (PTC) is less than 55 nucleotides away from the final exon-exon junction. The PTC is ", 
                  intron_pos_minus_iso_nt, " nucleotides away, and therefore is not a potential NMD candidate."), file= "nmd_information_output.txt",append=TRUE, sep = "\n")
    }
    }else{ # the following is if the intron position is upstream of the of the end of the isoform nucleotide sequence
      # set nmd_candidate variable to "no"
      nmd_candidate <- "no"
      # create a new column in the dataframe to store the result of a novel isoform being a nmd candidate or not (or "NA")
      file_with_intron_position$NMD_candidate[[i]] <- nmd_candidate
      # print this information to the console
      print("The position of the final exon-exon junction is before the PTC position. Therefore, this isoform is not a potential NMD candidate.")
      # append this information to the output text file
      cat(paste0("The position of the final exon-exon junction is before the PTC position. Therefore, this isoform is not a potential NMD candidate."), file= "nmd_information_output.txt",append=TRUE, sep = "\n")
    }
  }
  
  # iterate i by one each time through the loop in order to look at the next protein for potential NMD candidacy
  i <- i + 1
  
}

# write the csv file that contains a table with a new column with the name: "NMD_candidate"
# yes- novel protein isoform is a NMD candidate
# no - novel protein isofrom is not an NMD candidate
# NA - not able to be determined as of this time due to lack of information in JuncDB (see README for more information)
write.csv(file_with_intron_position, "intron_positions_for_nmd_with_nmd_candidacy.csv")

