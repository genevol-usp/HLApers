library(Biostrings)
library(tidyverse)

processed_quants <- commandArgs(TRUE)[1]
outdir <- commandArgs(TRUE)[2]

index <- readDNAStringSet("~/hlatx/1-make_indices/imgt_index.fa")

genos <- read_tsv(processed_quants) %>%
    distinct(subject, allele) %>%
    mutate(path = file.path(outdir, paste0("hla_", subject, ".fa"))) %>%
    split(.$subject)

map(genos, ~writeXStringSet(index[.$allele], unique(.$path)))
