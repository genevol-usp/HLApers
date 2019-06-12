library(hlaseqlib)
suppressPackageStartupMessages(library(tidyverse))

opts <- commandArgs(TRUE)
transcripts<- opts[1]
typings    <- opts[2]
outPrefix  <- opts[3] 

outindex <- paste0(outPrefix, "_index.fa")

index <- Biostrings::readDNAStringSet(transcripts)

typings_df <- read_tsv(typings)

alleles <- typings_df %>%
    pull(allele) %>%
    unique()

Biostrings::writeXStringSet(index[alleles], outindex)
