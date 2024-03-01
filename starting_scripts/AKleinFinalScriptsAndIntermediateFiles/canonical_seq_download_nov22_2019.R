# Title: Converting Symbol gene_id names to Uniprot amino acid sequence names
# Author: Alyssa Klein
# Date: November 16, 2019
# Modified: February 12, 2020

##################################################################################

# The purpose of this script is to take the gene ids and convert them to the appropriate 
# UniProt IDs.  There are two steps in this script. The first is to select for the UniProt IDs
# and entry names for each of the gene ids.  This step contains an intermediary step where
# alternative gene names are identified for genes that did not have a UniProt ID
# The IDs identified are then used to download the FASTA canonical protein sequences. 

# The UniProt IDs are the information necessary for the FASTA sequence download.  You will see 
# that I have also identified the UniProt Entry Name in case it was needed later on in research.

##################################################################################

# Step One Part A: Load in the gene ids and determine the UniProt IDs and UniProt Entry Names for each.

# load the UniProt.ws package (this is how the UniProt IDs will be obtained)
library(UniProt.ws)

# create object to call constructor for UniProt.ws and store the taxon as Homo sapiens
uniprot_homo_sapiens <- UniProt.ws(taxId=9606)

# check that the species for the uniprot_homo_sapiens object is Homo sapiens
species(uniprot_homo_sapiens)

# output all possible keytypes for the uniprot_homo_sapiens object so it can be determined which keytype will be used in 
# the select() function later on 
keytypes(uniprot_homo_sapiens)
# [1] "AARHUS/GHENT-2DPAGE"       
# [2] "AGD"                       
# [3] "ALLERGOME"                 
# [4] "ARACHNOSERVER"             
# [5] "BIOCYC"                    
# [6] "CGD"                       
# [7] "CLEANEX"                   
# [8] "CONOSERVER"                
# [9] "CYGD"                      
# [10] "DICTYBASE"                 
# [11] "DIP"                       
# [12] "DISPROT"                   
# [13] "DMDM"                      
# [14] "DNASU"                     
# [15] "DRUGBANK"                  
# [16] "ECHOBASE"                  
# [17] "ECO2DBASE"                 
# [18] "ECOGENE"                   
# [19] "EGGNOG"                    
# [20] "EMBL/GENBANK/DDBJ"         
# [21] "EMBL/GENBANK/DDBJ_CDS"     
# [22] "ENSEMBL"                   
# [23] "ENSEMBL_GENOMES"           
# [24] "ENSEMBL_GENOMES PROTEIN"   
# [25] "ENSEMBL_GENOMES TRANSCRIPT"
# [26] "ENSEMBL_PROTEIN"           
# [27] "ENSEMBL_TRANSCRIPT"        
# [28] "ENTREZ_GENE"               
# [29] "EUHCVDB"                   
# [30] "EUPATHDB"                  
# [31] "FLYBASE"                   
# [32] "GENECARDS"                 
# [33] "GENEFARM"                  
# [34] "GENEID"                    
# [35] "GENETREE"                  
# [36] "GENOLIST"                  
# [37] "GENOMEREVIEWS"             
# [38] "GENOMERNAI"                
# [39] "GERMONLINE"                
# [40] "GI_NUMBER*"                
# [41] "H-INVDB"                   
# [42] "HGNC"                      
# [43] "HOGENOM"                   
# [44] "HPA"                       
# [45] "HSSP"                      
# [46] "KEGG"                      
# [47] "KO"                        
# [48] "LEGIOLIST"                 
# [49] "LEPROMA"                   
# [50] "MAIZEGDB"                  
# [51] "MEROPS"                    
# [52] "MGI"                       
# [53] "MIM"                       
# [54] "MINT"                      
# [55] "NEXTBIO"                   
# [56] "NEXTPROT"                  
# [57] "OMA"                       
# [58] "ORPHANET"                  
# [59] "ORTHODB"                   
# [60] "PATRIC"                    
# [61] "PDB"                       
# [62] "PEROXIBASE"                
# [63] "PHARMGKB"                  
# [64] "PHOSSITE"                  
# [65] "PIR"                       
# [66] "POMBASE"                   
# [67] "PPTASEDB"                  
# [68] "PROTCLUSTDB"               
# [69] "PSEUDOCAP"                 
# [70] "REACTOME"                  
# [71] "REBASE"                    
# [72] "REFSEQ_NUCLEOTIDE"         
# [73] "REFSEQ_PROTEIN"            
# [74] "RGD"                       
# [75] "SGD"                       
# [76] "TAIR"                      
# [77] "TCDB"                      
# [78] "TIGR"                      
# [79] "TUBERCULIST"               
# [80] "UCSC"                      
# [81] "UNIPARC"                   
# [82] "UNIPATHWAY"                
# [83] "UNIPROTKB"                 
# [84] "UNIREF100"                 
# [85] "UNIREF50"                  
# [86] "UNIREF90"                  
# [87] "VECTORBASE"                
# [88] "WORLD-2DPAGE"              
# [89] "WORMBASE"                  
# [90] "WORMBASE_PROTEIN"          
# [91] "WORMBASE_TRANSCRIPT"       
# [92] "XENBASE"                   
# [93] "ZFIN" 

# output the possible columns within the uniprot_homo_sapiens object (to be used later in the select() function)
columns(uniprot_homo_sapiens)
# [1] "3D"                        
# [2] "AARHUS/GHENT-2DPAGE"       
# [3] "AGD"                       
# [4] "ALLERGOME"                 
# [5] "ARACHNOSERVER"             
# [6] "BIOCYC"                    
# [7] "CGD"                       
# [8] "CITATION"                  
# [9] "CLEANEX"                   
# [10] "CLUSTERS"                  
# [11] "COMMENTS"                  
# [12] "CONOSERVER"                
# [13] "CYGD"                      
# [14] "DATABASE(PDB)"             
# [15] "DATABASE(PFAM)"            
# [16] "DICTYBASE"                 
# [17] "DIP"                       
# [18] "DISPROT"                   
# [19] "DMDM"                      
# [20] "DNASU"                     
# [21] "DOMAIN"                    
# [22] "DOMAINS"                   
# [23] "DRUGBANK"                  
# [24] "EC"                        
# [25] "ECHOBASE"                  
# [26] "ECO2DBASE"                 
# [27] "ECOGENE"                   
# [28] "EGGNOG"                    
# [29] "EMBL/GENBANK/DDBJ"         
# [30] "EMBL/GENBANK/DDBJ_CDS"     
# [31] "ENSEMBL"                   
# [32] "ENSEMBL_GENOMES"           
# [33] "ENSEMBL_GENOMES PROTEIN"   
# [34] "ENSEMBL_GENOMES TRANSCRIPT"
# [35] "ENSEMBL_PROTEIN"           
# [36] "ENSEMBL_TRANSCRIPT"        
# [37] "ENTREZ_GENE"               
# [38] "ENTRY-NAME"                
# [39] "EUHCVDB"                   
# [40] "EUPATHDB"                  
# [41] "EXISTENCE"                 
# [42] "FAMILIES"                  
# [43] "FEATURES"                  
# [44] "FLYBASE"                   
# [45] "FUNCTION"                  
# [46] "GENECARDS"                 
# [47] "GENEFARM"                  
# [48] "GENEID"                    
# [49] "GENES"                     
# [50] "GENETREE"                  
# [51] "GENOLIST"                  
# [52] "GENOMEREVIEWS"             
# [53] "GENOMERNAI"                
# [54] "GERMONLINE"                
# [55] "GI_NUMBER*"                
# [56] "GO"                        
# [57] "GO-ID"                     
# [58] "H-INVDB"                   
# [59] "HGNC"                      
# [60] "HOGENOM"                   
# [61] "HPA"                       
# [62] "HSSP"                      
# [63] "ID"                        
# [64] "INTERACTOR"                
# [65] "INTERPRO"                  
# [66] "KEGG"                      
# [67] "KEYWORD-ID"                
# [68] "KEYWORDS"                  
# [69] "KO"                        
# [70] "LAST-MODIFIED"             
# [71] "LEGIOLIST"                 
# [72] "LENGTH"                    
# [73] "LEPROMA"                   
# [74] "MAIZEGDB"                  
# [75] "MEROPS"                    
# [76] "MGI"                       
# [77] "MIM"                       
# [78] "MINT"                      
# [79] "NEXTBIO"                   
# [80] "NEXTPROT"                  
# [81] "OMA"                       
# [82] "ORGANISM"                  
# [83] "ORGANISM-ID"               
# [84] "ORPHANET"                  
# [85] "ORTHODB"                   
# [86] "PATHWAY"                   
# [87] "PATRIC"                    
# [88] "PDB"                       
# [89] "PEROXIBASE"                
# [90] "PHARMGKB"                  
# [91] "PHOSSITE"                  
# [92] "PIR"                       
# [93] "POMBASE"                   
# [94] "PPTASEDB"                  
# [95] "PROTCLUSTDB"               
# [96] "PROTEIN-NAMES"             
# [97] "PSEUDOCAP"                 
# [98] "REACTOME"                  
# [99] "REBASE"                    
# [100] "REFSEQ_NUCLEOTIDE"         
# [101] "REFSEQ_PROTEIN"            
# [102] "REVIEWED"                  
# [103] "RGD"                       
# [104] "SCORE"                     
# [105] "SEQUENCE"                  
# [106] "SGD"                       
# [107] "SUBCELLULAR-LOCATIONS"     
# [108] "TAIR"                      
# [109] "TAXONOMIC-LINEAGE"         
# [110] "TCDB"                      
# [111] "TIGR"                      
# [112] "TOOLS"                     
# [113] "TUBERCULIST"               
# [114] "UCSC"                      
# [115] "UNIPARC"                   
# [116] "UNIPATHWAY"                
# [117] "UNIPROTKB"                 
# [118] "UNIREF100"                 
# [119] "UNIREF50"                  
# [120] "UNIREF90"                  
# [121] "VECTORBASE"                
# [122] "VERSION"                   
# [123] "VIRUS-HOSTS"               
# [124] "WORLD-2DPAGE"              
# [125] "WORMBASE"                  
# [126] "WORMBASE_PROTEIN"          
# [127] "WORMBASE_TRANSCRIPT"       
# [128] "XENBASE"                   
# [129] "ZFIN"

# read in each of the tables that contain the gene ids
gene_ids <- read.table(file = "gene_ids_for_conversion_script.txt", header=FALSE, sep="\n", stringsAsFactors = FALSE)

# see what the rownames are in the gene_ids object
# here it can be seen that the numbering of the rows (not the gene_ids themselves) are the rownames
rownames(gene_ids)

# create a gene_id_names object to store only the column of the gene_ids
gene_id_names <- gene_ids[,1]

# use the gene_id_names object in the select function to select the UniProt ID for each gene name
gene_names_to_uniprot_ids <- select(up,
                  keys = gene_id_names,
                  columns = "ID",
                  keytype = "GENECARDS")

# write a csv file that contains the GENECARDS gene name and the UniProt ID (denoted as "ID")
write.csv(gene_names_to_uniprot_ids, "gene_names_to_ids.csv")

# write another csv file that contains only the unique GENECARDS gene names and the UniProt ID (denoted as "ID")
write.csv(unique(gene_names_to_uniprot_ids), "unique_gene_names_to_ids.csv")

# add the entry names for the proteins as seen on UniProt
gene_names_to_uniprot_entryname<- select(uniprot_homo_sapiens,
                                    keys = gene_id_names,
                                    columns = "ENTRY-NAME",
                                    keytype = "GENECARDS")

# write a csv file that contains the GENECARDS gene name and the entry name (denoted as "ENTRY_NAME")
write.csv(gene_names_to_uniprot_entryname, "gene_names_to_entrynames.csv")

# write a csv file that contains only the unique GENECARDS gene name and the entry name (denoted as "ENTRY_NAME")
write.csv(unique(gene_names_to_uniprot_entryname), "unique_gene_names_to_entrynames.csv")

# Output what genes do not have a UniProt ID 
# These will be genes that are checked to see that they actually aren't protein coding, or if they have an alternative gene name
# Use the %in% function to identify the genes that did not have a UniProt ID assigned to them
genes_no_uniprot_ids <- gene_id_names %in% gene_names_to_uniprot_ids$GENECARDS

# Index for the ids that are not found in the gene_names_to_uniprot object
which(genes_no_uniprot_ids == FALSE)
# [1]    13   23   29   31   38   39   40   49   68   83
# [11]   87   93  103  110  123  128  131  145  184  200
# [21]  201  212  221  222  236  254  261  266  270  276
# [31]  279  286  299  303  315  340  348  358  359  360
# [41]  370  374  378  400  403  412  418  420  423  436
# [51]  472  491  492  495  501  506  526  537  557  575
# [61]  585  588  589  590  593  596  602  616  626  629
# [71]  630  635  642  648  675  677  679  695  707  708
# [81]  725  754  755  758  763  764  778  788  791  796
# [91]  798  804  814  818  819  821  825  831  834  847
# [101]  857  858  867  877  886  897  902  922  924  943
# [111]  946  956  957  968  992 1006 1007 1019 1036 1046
# [121] 1047 1056 1072 1090 1101 1107 1109 1121 1124 1125
# [131] 1129 1135 1139 1153 1159 1173 1183 1190 1191 1192
# [141] 1196 1198 1199 1228 1236 1247 1252 1272 1286 1288
# [151] 1297 1300 1304 1335 1351 1378 1392 1397 1429 1435
# [161] 1462 1478 1483 1487 1499 1500 1521 1535 1538 1541
# [171] 1544 1580 1584 1588 1590 1598 1623 1628 1629 1631
# [181] 1654 1670 1678 1694 1710 1716 1721 1722 1731 1748
# [191] 1749 1759 1760 1766 1791 1803 1810 1825 1830 1866
# [201] 1877 1888 1906 1910 1911 1927 1928 1939 1949 1958
# [211] 1962 1966 1980 1996 2002 2004 2014 2022 2031 2032
# [221] 2049 2050 2061 2062 2072 2093 2106 2119 2125 2128
# [231] 2143 2146 2147 2153 2154 2174 2175 2181 2189 2234
# [241] 2244 2248 2257 2258 2262 2265 2268 2286 2303 2305
# [251] 2314 2321 2338 2339 2346 2350 2362 2363 2376 2385
# [261] 2389 2395 2396 2397 2418 2422 2425 2432 2433 2460
# [271] 2463 2500 2504 2505 2509 2512 2520 2521 2536 2540
# [281] 2549 2563 2570 2589 2598 2599 2605 2608 2609 2680
# [291] 2697 2719 2721 2727 2730 2740 2751 2757 2770 2787
# [301] 2808 2809 2831 2885 2887 2911 2928 2930 2931 2933
# [311] 2961 2966 2972 2979 2983 2984 2995 2997 3020 3030
# [321] 3032 3045 3069 3089 3092 3102 3134 3138 3139 3162
# [331] 3163 3164 3167 3171 3174 3180 3196 3198 3203 3216
# [341] 3235 3236 3251 3252 3265 3276 3280 3282 3283 3288
# [351] 3289 3290 3321 3323 3337 3340 3349 3353 3380 3384
# [361] 3386 3404 3417 3419 3434 3437 3438 3461 3470 3485
# [371] 3500 3509 3511 3525 3528 3534 3539 3560 3568 3572
# [381] 3586 3594 3599 3635 3644 3652 3655 3656 3676 3680
# [391] 3683 3702 3717 3720 3739 3743 3745 3749 3756 3757
# [401] 3819 3824 3825 3826 3840 3865 3866 3878 3885 3911
# [411] 3912 3917 3920 3921 3929 3958 3964 3970 3975 3981
# [421] 3986 3999 4007 4026 4040 4059 4071 4093 4105 4108
# [431] 4109 4110 4120 4124 4146 4164 4187 4190 4191 4193
# [441] 4202 4215 4220 4223 4224 4225 4235 4242 4249 4266
# [451] 4281 4282 4294 4300 4318 4327 4330 4335 4342 4347
# [461] 4372 4379 4380 4409 4413 4440 4441 4442 4481 4493
# [471] 4497 4520 4522 4550 4563 4566 4571

# write all of the gene names that did not identify with a UniProt ID into a file named: "gene_ids_no_uniprot_ids.csv".
write.csv(gene_id_names[c(13, 23, 29, 31, 38, 39, 40, 49, 68, 83, 87, 93, 103,  110,  
                123,  128,  131,  145,  184,  200, 201,  212,  221,  222,  
                236,  254,  261,  266,  270,  276, 279,  286,  299,  303,  
                315,  340,  348,  358,  359,  360, 370,  374,  378,  400,  
                403,  412,  418,  420,  423,  436,472,  491,  492,  495,  
                501,  506,  526,  537,  557,  575, 585,  588,  589,  590,  
                593,  596,  602,  616,  626,  629, 630,  635, 642 , 648 , 
                675,  677,  679,  695,  707,  708,
                725,  754,  755,  758,  763,  764,  778,  788,  791,  796,
                798,  804,  814,  818,  819,  821,  825,  831,  834,  847,
                857,  858,  867,  877,  886,  897,  902,  922, 924,  943,
                946,  956,  957,  968,  992, 1006, 1007, 1019, 1036, 1046,
                1047, 1056, 1072, 1090, 1101, 1107, 1109, 1121, 1124, 1125,
                1129, 1135, 1139, 1153, 1159, 1173, 1183, 1190, 1191, 1192,
                1196, 1198, 1199, 1228, 1236, 1247, 1252, 1272, 1286, 1288,
                1297, 1300, 1304, 1335, 1351, 1378, 1392, 1397, 1429, 1435,
                1462, 1478, 1483, 1487, 1499, 1500, 1521, 1535, 1538, 1541,
                1544, 1580, 1584, 1588, 1590, 1598, 1623, 1628, 1629, 1631,
                1654, 1670, 1678, 1694, 1710, 1716, 1721, 1722, 1731, 1748,
                1749, 1759, 1760, 1766, 1791, 1803, 1810, 1825, 1830, 1866,
                1877, 1888, 1906, 1910, 1911, 1927, 1928, 1939, 1949, 1958,
                1962, 1966, 1980, 1996, 2002, 2004, 2014, 2022, 2031, 2032,
                2049, 2050, 2061, 2062, 2072, 2093, 2106, 2119, 2125, 2128,
                2143, 2146, 2147, 2153, 2154, 2174, 2175, 2181, 2189, 2234,
                2244, 2248, 2257, 2258, 2262, 2265, 2268, 2286, 2303, 2305, 
                2314, 2321, 2338, 2339, 2346, 2350, 2362, 2363,2376, 2385,
                2389, 2395, 2396, 2397, 2418, 2422, 2425, 2432,2433, 2460,
                2463, 2500, 2504, 2505, 2509, 2512, 2520, 2521, 2536, 2540,
                2549, 2563, 2570, 2589, 2598, 2599, 2605, 2608, 2609, 2680,
                2697, 2719, 2721, 2727, 2730, 2740, 2751, 2757, 2770, 2787,
                2808, 2809, 2831, 2885, 2887, 2911, 2928, 2930, 2931, 2933,
                2961, 2966, 2972, 2979, 2983, 2984, 2995, 2997, 3020, 3030,
                3032, 3045, 3069, 3089, 3092, 3102, 3134, 3138, 3139, 3162,
                3163, 3164, 3167, 3171, 3174, 3180, 3196, 3198, 3203, 3216,
                3235, 3236, 3251, 3252, 3265, 3276, 3280, 3282, 3283, 3288,
                3289, 3290, 3321, 3323, 3337, 3340, 3349, 3353, 3380, 3384,
                3386, 3404, 3417, 3419, 3434, 3437, 3438, 3461, 3470, 3485,
                3500, 3509, 3511, 3525, 3528, 3534, 3539, 3560, 3568, 3572,
                3586, 3594, 3599, 3635, 3644, 3652, 3655, 3656, 3676, 3680,
                3683, 3702, 3717, 3720, 3739, 3743, 3745, 3749, 3756, 3757,
                3819, 3824, 3825, 3826, 3840, 3865, 3866, 3878, 3885, 3911,
                3912, 3917, 3920, 3921, 3929, 3958, 3964, 3970, 3975, 3981,
                3986, 3999, 4007, 4026, 4040, 4059, 4071, 4093, 4105, 4108,
                4109, 4110, 4120, 4124, 4146, 4164, 4187, 4190, 4191, 4193,
                4202, 4215, 4220, 4223, 4224, 4225, 4235, 4242, 4249, 4266,
                4281, 4282, 4294, 4300, 4318, 4327, 4330, 4335, 4342, 4347,
                4372, 4379, 4380, 4409, 4413, 4440, 4441, 4442, 4481, 4493,
                4497, 4520,  4522, 4550, 4563, 4566, 4571)], "gene_ids_no_uniprot_ids.csv")

# The following are the gene names that did not have a UniProt ID match:

# [1] "TMEM180"        "FAM66C"         "LINC00507"     
# [4] "TPT1-AS1"       "HERC2P9"        "LOC100996255"  
# [7] "DPH6-AS1"       "SLX1B-SULT1A4"  "MIR548N"       
# [10] "SELO"           "SELT"           "HRASLS"        
# [13] "TMEM161B-AS1"   "LINC00265"      "LINC01000"     
# [16] "PVT1"           "FRG1HP"         "MINOS1-NBL1"   
# [19] "LOC728730"      "LINC00634"      "RPL23AP82"     
# [22] "GTF2H2B"        "TRAF3IP2-AS1"   "TRAF3IP2-AS1"  
# [25] "ARMCX5-GPRASP2" "TMEM5"          "SNORD109B"     
# [28] "SLX1B-SULT1A4"  "DNAJC27-AS1"    "ZFAS1"         
# [31] "FAM19A5"        "LOC645513"      "JPX"           
# [34] "GAREML"         "ATP5F1"         "KIAA1033"      
# [37] "CCDC176"        "C16orf45"       "SLX1B-SULT1A4" 
# [40] "SLX1B-SULT1A4"  "KIAA0195"       "LIPE-AS1"      
# [43] "PQLC3"          "SOX2-OT"        "HRASLS"        
# [46] "GUCY1A3"        "PAPD4"          "KIAA0141"      
# [49] "ADCY10P1"       "LINC00996"      "RQCD1"         
# [52] "SETD8"          "ZNF664-FAM101A" "GAREML"        
# [55] "TMEM161B-AS1"   "C16orf45"       "MINOS1-NBL1"   
# [58] "LINC01140"      "MRVI1-AS1"      "FLJ41278"      
# [61] "KIAA1033"       "SETD8"          "LINC00507"     
# [64] "LINC00507"      "TPT1-AS1"       "RNF219-AS1"    
# [67] "FAM179B"        "WHAMMP2"        "CIRH1A"        
# [70] "USP32P1"        "USP32P1"        "ZNF271P"       
# [73] "C19orf60"       "CENPBD1P1"      "DGCR5"         
# [76] "LINC00634"      "LINC00693"      "LINC01094"     
# [79] "GTF2H2B"        "ZBED3-AS1"      "TRAF3IP2-AS1"  
# [82] "LINC01189"      "LINC01189"      "FAM73B"        
# [85] "JPX"            "ARMCX5-GPRASP2" "BCDIN3D-AS1"   
# [88] "LINC00907"      "MIR548N"        "DGCR5"         
# [91] "WHSC1"          "INMT-FAM188B"   "LOC729970"     
# [94] "DGCR5"          "DGCR5"          "ZBED3-AS1"     
# [97] "MFSD4"          "ZNF664-FAM101A" "PCNX"          
# [100] "DGCR5"          "ZBED3-AS1"      "TMEM161B-AS1"  
# [103] "C1orf228"       "KIAA1033"       "SNORD109B"     
# [106] "DNAJC27-AS1"    "MIR548N"        "ADCY10P1"      
# [109] "TRAF3IP2-AS1"   "C10orf35"       "TPT1-AS1"      
# [112] "LOC100507387"   "ADCY10P1"       "LOC645166"     
# [115] "MIR548N"        "CASC15"         "TRAF3IP2-AS1"  
# [118] "SNAI3-AS1"      "PTCHD3P1"       "KIAA1033"      
# [121] "ZNF664-FAM101A" "SLX1B-SULT1A4"  "AHSA2"         
# [124] "ADAMTS9-AS2"    "PAPD4"          "DOPEY1"        
# [127] "JAZF1-AS1"      "MINOS1-NBL1"    "FAM66C"        
# [130] "LINC00507"      "RBFADN"         "PAPD4"         
# [133] "ARMCX5-GPRASP2" "LOC399815"      "FAM66C"        
# [136] "SNORD109B"      "LOC100287036"   "KIAA1468"      
# [139] "PARD6G-AS1"     "C19orf60"       "LIPE-AS1"      
# [142] "CENPBD1P1"      "LOC100506274"   "C4orf32"       
# [145] "ADCY10P1"       "LINC01189"      "FAM69B"        
# [148] "TRAM2-AS1"      "HIATL1"         "GLUD1P7"       
# [151] "LINC00907"      "LOC100506274"   "QTRTD1"        
# [154] "ZNF664-FAM101A" "KIAA0922"       "ERICH1-AS1"    
# [157] "LINC00493"      "TRAF3IP2-AS1"   "EBLN3"         
# [160] "LOC645166"      "LINC00327"      "MIR548AE2"     
# [163] "C19orf60"       "LINC01122"      "LINC00634"     
# [166] "PVRL3"          "JPX"            "SRP14-AS1"     
# [169] "ABCA17P"        "BRE"            "C7orf55-LUC7L2"
# [172] "DGCR5"          "HRASLS"         "KIAA0922"      
# [175] "ERBB2IP"        "LINC01420"      "ZNF664-FAM101A"
# [178] "SLX1B-SULT1A4"  "SLX1B-SULT1A4"  "C19orf60"      
# [181] "FLJ31104"       "FAM21C"         "RPL23AP82"     
# [184] "LOC389765"      "LOC100506393"   "PCNX"          
# [187] "SLX1B-SULT1A4"  "SLX1B-SULT1A4"  "PQLC3"         
# [190] "DGCR5"          "DGCR5"          "ADCY10P1"      
# [193] "CCDC162P"       "C7orf55-LUC7L2" "ZBED3-AS1"     
# [196] "LINC00467"      "PCNX"           "C4orf22"       
# [199] "FRG1HP"         "NDUFA6-AS1"     "ADCY10P1"      
# [202] "ARMCX5-GPRASP2" "USP32P1"        "DNAJC27-AS1"   
# [205] "LOC728730"      "JPX"            "ADCK3"         
# [208] "NUPL2"          "FLJ41278"       "LINC00189"     
# [211] "JAZF1-AS1"      "FRG1HP"         "DPH6-AS1"      
# [214] "RPL23AP82"      "SKIV2L2"        "ZBED3-AS1"     
# [217] "PTCHD3P1"       "TPTE2P3"        "BRE"           
# [220] "LOC728730"      "LOC100506207"   "ADCY10P1"      
# [223] "CCDC94"         "BCRP2"          "CCDC64"        
# [226] "LOC100507387"   "HERC2P9"        "ERBB2IP"       
# [229] "GLUD1P7"        "FLJ41278"       "CENPBD1P1"     
# [232] "SEPT2"          "SEPT2"          "LOC340017"     
# [235] "LOC340017"      "ZNF271P"        "KIAA1468"      
# [238] "DGCR5"          "GDNF-AS1"       "SKIV2L2"       
# [241] "KIAA1033"       "ZFAS1"          "LOC440040"     
# [244] "USP32P1"        "ADAMTS9-AS2"    "LOC100506207"  
# [247] "ARMCX5-GPRASP2" "UBA6-AS1"       "C8orf37-AS1"   
# [250] "ALOX12P2"       "ALOX12P2"       "C8orf37-AS1"   
# [253] "C1orf204"       "RFWD2"          "KIAA1462"      
# [256] "LINC01001"      "CCDC53"         "CCDC53"        
# [259] "RRN3P3"         "TLDC1"          "C17orf51"      
# [262] "PQLC1"          "PQLC1"          "PQLC1"         
# [265] "LINC00320"      "NUP50-AS1"      "ENTPD3-AS1"    
# [268] "ARHGEF26-AS1"   "ARHGEF26-AS1"   "RPL23AP53"     
# [271] "LOC392232"      "RRN3P3"         "LOC283856"     
# [274] "MTSS1L"         "USP32P2"        "PQLC1"         
# [277] "PCBP1-AS1"      "PCBP1-AS1"      "LOC220729"     
# [280] "LOC100506746"   "LOC401320"      "RFWD2"         
# [283] "PRH1-PRR4"      "LINC00665"      "NUP50-AS1"     
# [286] "BHLHE40-AS1"    "SCOC-AS1"       "TCEB1"         
# [289] "LOC100133669"   "LOC283856"      "FLJ33534"      
# [292] "THUMPD3-AS1"    "ENTPD3-AS1"     "LINC00882"     
# [295] "LOC220729"      "NR2F1-AS1"      "APTR"          
# [298] "TCEB1"          "TP73-AS1"       "ENTPD3-AS1"    
# [301] "TCEB1"          "TMEM55A"        "TMEM72-AS1"    
# [304] "LOC441666"      "TMEM72-AS1"     "PRH1-PRR4"     
# [307] "LINC00426"      "SPG20"          "TPTE2P5"       
# [310] "LINC00282"      "LOC283856"      "USP32P2"       
# [313] "MRPL45P2"       "PQLC1"          "HPN-AS1"       
# [316] "LINC00665"      "TSSC1"          "LOC100506474"  
# [319] "PAK7"           "NUP50-AS1"      "FGD5-AS1"      
# [322] "KIAA1524"       "LOC553103"      "MIR548O2"      
# [325] "TMEM55A"        "VLDLR-AS1"      "NR2F2-AS1"     
# [328] "MTSS1L"         "AMZ2P1"         "TP53TG1"       
# [331] "MGC72080"       "RPL23AP53"      "VLDLR-AS1"     
# [334] "PLA2G16"        "LOC146880"      "TP53TG1"       
# [337] "C15orf38-AP3S2" "USP32P2"        "PCBP1-AS1"     
# [340] "LOC401320"      "FAM96A"         "VWA9"          
# [343] "FAM65C"         "ENTPD3-AS1"     "GUSBP2"        
# [346] "DLEU2"          "LOC146880"      "ALS2CR12"      
# [349] "LINC00320"      "GTF2IP1"        "GTF2IP1"       
# [352] "GTF2IP1"        "LY86-AS1"       "LINC00951"     
# [355] "TSSC1"          "NUP50-AS1"      "LYPLAL1-AS1"   
# [358] "PRKRIR"         "LINC00665"      "C2orf61"       
# [361] "LOC654342"      "GUSBP9"         "LOC100288181"  
# [364] "LOC286297"      "LOC283856"      "FLJ33534"      
# [367] "LOC100506474"   "TMEM206"        "ENTPD1-AS1"    
# [370] "CCDC53"         "RAD51L3-RFFL"   "LOC100131655"  
# [373] "HPN-AS1"        "THUMPD3-AS1"    "LINC00882"     
# [376] "SDHAP1"         "LINC01184"      "KIAA0368"      
# [379] "MRPL45P2"       "DFNB31"         "PAK7"          
# [382] "HIATL2"         "LOC653513"      "GTF2IP1"       
# [385] "LOC643339"      "TLDC1"          "AES"           
# [388] "FLJ33534"       "C8orf34-AS1"    "FAM129B"       
# [391] "SEP15"          "PRH1-PRR4"      "EFTUD1"        
# [394] "FLJ33534"       "LOC654342"      "FAM58A"        
# [397] "ZCCHC11"        "LOC441666"      "KDELC2"        
# [400] "PRH1-PRR4"      "GS1-124K5.11"   "KGFLP2"        
# [403] "LOC286297"      "CXorf36"        "ZCCHC11"       
# [406] "LARGE"          "NUP50-AS1"      "TP53TG1"       
# [409] "IKBKAP"         "TCEB2"          "MTSS1L"        
# [412] "LINC00665"      "TSSC1"          "PCBP1-AS1"     
# [415] "THUMPD3-AS1"    "LOC283856"      "KIAA1644"      
# [418] "TP53TG1"        "SEP15"          "METTL10"       
# [421] "SUV420H1"       "LOC283856"      "FLJ33534"      
# [424] "C4orf27"        "VLDLR-AS1"      "CYP4F24P"      
# [427] "TP73-AS1"       "PQLC1"          "GTF2IP1"       
# [430] "ZFHX4-AS1"      "LOC100133669"   "VLDLR-AS1"     
# [433] "CD27-AS1"       "TPTE2P5"        "LINC00461"     
# [436] "NUTM2A-AS1"     "GUSBP2"         "GTF2IP1"       
# [439] "GTF2IP1"        "LINC01278"      "LINC00882"     
# [442] "LHFP"           "NR2F2-AS1"      "LOC283922"     
# [445] "C17orf51"       "LOC440434"      "GUSBP9"        
# [448] "LOC653513"      "LOC100507377"   "LY75-CD302"    
# [451] "KGFLP2"         "KGFLP2"         "LOC100131564"  
# [454] "PRH1-PRR4"      "PCBP1-AS1"      "KIAA1644"      
# [457] "SCOC-AS1"       "MIR548O2"       "TP53TG1"       
# [460] "KIAA1462"       "LOC553103"      "CCBL1"         
# [463] "PIR-FIGF"       "C1orf204"       "TMEM254-AS1"   
# [466] "LOC392232"      "ZFHX4-AS1"      "KIAA0368"      
# [469] "CD27-AS1"       "CD27-AS1"       "HPN-AS1"       
# [472] "LOC392232"      "ZCCHC11"        "GTF2IP1"       
# [475] "LOC100131655"   "DYX1C1-CCPG1"   "FAM58A"   

# write all of the UNIQUE gene names that did not identify with a UniProt ID into a file named: "unique_gene_ids_no_uniprot_ids.csv".
write.csv(unique(gene_id_names[c(13, 23, 29, 31, 38, 39, 40, 49, 68, 83, 87, 93, 103,  110,  
                          123,  128,  131,  145,  184,  200, 201,  212,  221,  222,  
                          236,  254,  261,  266,  270,  276, 279,  286,  299,  303,  
                          315,  340,  348,  358,  359,  360, 370,  374,  378,  400,  
                          403,  412,  418,  420,  423,  436,472,  491,  492,  495,  
                          501,  506,  526,  537,  557,  575, 585,  588,  589,  590,  
                          593,  596,  602,  616,  626,  629, 630,  635, 642 , 648 , 
                          675,  677,  679,  695,  707,  708,
                          725,  754,  755,  758,  763,  764,  778,  788,  791,  796,
                          798,  804,  814,  818,  819,  821,  825,  831,  834,  847,
                          857,  858,  867,  877,  886,  897,  902,  922, 924,  943,
                          946,  956,  957,  968,  992, 1006, 1007, 1019, 1036, 1046,
                          1047, 1056, 1072, 1090, 1101, 1107, 1109, 1121, 1124, 1125,
                          1129, 1135, 1139, 1153, 1159, 1173, 1183, 1190, 1191, 1192,
                          1196, 1198, 1199, 1228, 1236, 1247, 1252, 1272, 1286, 1288,
                          1297, 1300, 1304, 1335, 1351, 1378, 1392, 1397, 1429, 1435,
                          1462, 1478, 1483, 1487, 1499, 1500, 1521, 1535, 1538, 1541,
                          1544, 1580, 1584, 1588, 1590, 1598, 1623, 1628, 1629, 1631,
                          1654, 1670, 1678, 1694, 1710, 1716, 1721, 1722, 1731, 1748,
                          1749, 1759, 1760, 1766, 1791, 1803, 1810, 1825, 1830, 1866,
                          1877, 1888, 1906, 1910, 1911, 1927, 1928, 1939, 1949, 1958,
                          1962, 1966, 1980, 1996, 2002, 2004, 2014, 2022, 2031, 2032,
                          2049, 2050, 2061, 2062, 2072, 2093, 2106, 2119, 2125, 2128,
                          2143, 2146, 2147, 2153, 2154, 2174, 2175, 2181, 2189, 2234,
                          2244, 2248, 2257, 2258, 2262, 2265, 2268, 2286, 2303, 2305, 
                          2314, 2321, 2338, 2339, 2346, 2350, 2362, 2363,2376, 2385,
                          2389, 2395, 2396, 2397, 2418, 2422, 2425, 2432,2433, 2460,
                          2463, 2500, 2504, 2505, 2509, 2512, 2520, 2521, 2536, 2540,
                          2549, 2563, 2570, 2589, 2598, 2599, 2605, 2608, 2609, 2680,
                          2697, 2719, 2721, 2727, 2730, 2740, 2751, 2757, 2770, 2787,
                          2808, 2809, 2831, 2885, 2887, 2911, 2928, 2930, 2931, 2933,
                          2961, 2966, 2972, 2979, 2983, 2984, 2995, 2997, 3020, 3030,
                          3032, 3045, 3069, 3089, 3092, 3102, 3134, 3138, 3139, 3162,
                          3163, 3164, 3167, 3171, 3174, 3180, 3196, 3198, 3203, 3216,
                          3235, 3236, 3251, 3252, 3265, 3276, 3280, 3282, 3283, 3288,
                          3289, 3290, 3321, 3323, 3337, 3340, 3349, 3353, 3380, 3384,
                          3386, 3404, 3417, 3419, 3434, 3437, 3438, 3461, 3470, 3485,
                          3500, 3509, 3511, 3525, 3528, 3534, 3539, 3560, 3568, 3572,
                          3586, 3594, 3599, 3635, 3644, 3652, 3655, 3656, 3676, 3680,
                          3683, 3702, 3717, 3720, 3739, 3743, 3745, 3749, 3756, 3757,
                          3819, 3824, 3825, 3826, 3840, 3865, 3866, 3878, 3885, 3911,
                          3912, 3917, 3920, 3921, 3929, 3958, 3964, 3970, 3975, 3981,
                          3986, 3999, 4007, 4026, 4040, 4059, 4071, 4093, 4105, 4108,
                          4109, 4110, 4120, 4124, 4146, 4164, 4187, 4190, 4191, 4193,
                          4202, 4215, 4220, 4223, 4224, 4225, 4235, 4242, 4249, 4266,
                          4281, 4282, 4294, 4300, 4318, 4327, 4330, 4335, 4342, 4347,
                          4372, 4379, 4380, 4409, 4413, 4440, 4441, 4442, 4481, 4493,
                          4497, 4520,  4522, 4550, 4563, 4566, 4571)]), "unique_gene_ids_no_uniprot_ids.csv")


# Step One Part B: Use the Alternative Gene Ids for original "no matches" to identify UniProt IDs and Entry Names

# Take the list of gene ids that did not have a UniProt ID and check them using an external database to see
# if they are protein coding or not (used GeneCards to check this; see table "no_match_category_names.xlsx" for results)
# Repeat the steps completed above with the alternative gene ids for the previous "no matches" that are protein coding
gene_ids_no_orig_matches <- read.table(file = "updated_gene_ids_for_original_no_matches.txt", header= FALSE, sep="\n", stringsAsFactors=FALSE)

# create a gene_id_names_no_orig_matches object to store only the column of the gene_id_names_no_orig_matches
gene_id_names_no_orig_matches <- gene_ids_no_orig_matches[,1]

# use the gene_id_names_no_orig_matches object in the select function to select the UniProt ID for each gene name
orig_no_matches_to_uniprot_ids <- select(uniprot_homo_sapiens,
                                         keys = gene_id_names_no_orig_matches,
                                         columns = "ID",
                                         keytype = "GENECARDS")

# write a csv file that contains the GENECARDS gene name and the UniProt ID (denoted as "ID")
write.csv(orig_no_matches_to_uniprot_ids, "gene_names_to_ids_orig_no_matches.csv")

# write another csv file that contains only the unique GENECARDS gene names and the UniProt ID (denoted as "ID")
write.csv(unique(orig_no_matches_to_uniprot_ids), "unique_gene_names_to_ids_orig_no_matches.csv")

# use the gene_id_names_no_orig_matches object in the select function to select the Entry Name for each gene name
orig_no_matches_to_uniprot_entryname <- select(uniprot_homo_sapiens,
                                               keys = gene_id_names_no_orig_matches,
                                               columns = "ENTRY-NAME",
                                               keytype = "GENECARDS")

# write a csv file that contains the GENECARDS gene name and the entry name (denoted as "ENTRY_NAME")
write.csv(orig_no_matches_to_uniprot_entryname, "gene_names_to_entrynames_orig_no_matches.csv")

# write a csv file that contains only the unique GENECARDS gene name and the entry name (denoted as "ENTRY_NAME")
write.csv(unique(orig_no_matches_to_uniprot_entryname), "unique_gene_names_to_entrynames_orig_no_matches.csv")

##################################################################################

# Step Two: Download the canonical FASTA protein sequences using the getseq() function from the bio3d package
# Use the get.seq() function to generate the FASTA files for each of the canonical sequences in the UniProt database
# get.seq() function is a part of the bio3d package
library(bio3d)
# use the file generated that contains all of the unique UniProt Ids
# this file was created through copy and paste of the following files:
# "unique_gene_names_to_ids_orig_no_matches.csv"
# "unique_gene_names_to_ids.csv"

# read in the text file that contains all of the protein coding gene ids
id_csv_file <- read.table(file = "uniprot_ids_for_canonical_download.txt", header=FALSE, sep="\n", stringsAsFactors = FALSE)
# set the ids variable to contain the V1 column from the id_csv_file object
ids <- id_csv_file$V1
# use bio3d package to download the canonical sequences for each of the unique ids and place it into an output file ("canonical_seqs.fasta")
canonical_seqs <- get.seq(ids, outfile = "canonical_seqs.fasta", db = "uniprot", verbose = FALSE)
