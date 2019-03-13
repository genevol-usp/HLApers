library(Biostrings)
library(tidyverse)

processed_quants <- commandArgs(TRUE)[1]
out <- commandArgs(TRUE)[2]

index <- "../1-make_indices/data/hladb/hladb.fasta" %>%
    readDNAStringSet()

alleles <- read_tsv(processed_quants) %>%
    pull(allele) %>%
    unique()

writeXStringSet(index[alleles], out)
