library(Biostrings)
library(tidyverse)

processed_quants <- commandArgs(TRUE)[1]
outdir <- commandArgs(TRUE)[2]

index <- 
    readDNAStringSet("~/hla_rnaseq_pipeline/1-index_preparation/imgt_index.fa")

genos <- 
    read_tsv(processed_quants) %>%
    select(subject, locus, allele, est_counts) %>%
    group_by(subject, locus, allele) %>%
    summarize(est_counts = sum(est_counts)) %>%
    ungroup() %>%
    separate_rows(allele, sep = "=") %>% 
    arrange(subject, locus, desc(est_counts)) %>%
    group_by(subject, locus) %>%
    filter(est_counts >= 5L | allele == first(allele)) %>%
    ungroup() %>%
    select(subject, allele) %>%
    mutate(path = file.path(outdir, paste0("hla_", subject, ".fa"))) %>%
    split(.$subject)

map(genos, ~writeXStringSet(index[.$allele], unique(.$path)))
