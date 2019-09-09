library(hlaseqlib)

opts <- commandArgs(TRUE)
transcripts<- opts[1]
typings    <- opts[2]
outPrefix  <- opts[3] 

outindex <- paste0(outPrefix, "_index.fa")

index <- Biostrings::readDNAStringSet(transcripts)

typings_df <- readr::read_tsv(typings)

alleles <- unlist(strsplit(unique(typings_df$allele), "-"))

Biostrings::writeXStringSet(index[alleles], outindex)
