#!/bin/bash
set -e
set -u
set -o pipefail

esearch -db sra -query PRJNA385815  | efetch --format runinfo | cut -d ',' -f 1 | grep SRR | xargs ~/sratoolkit.2.9.2-centos_linux64/bin/fastq-dump -O ~/probiotic_gene_exp/data/rawseqs/ --split-files