library(tidyverse)

quant_file <- commandArgs(TRUE)[1]
out <- commandArgs(TRUE)[2]

index <- "../1-make_indices/data/hladb/hladb.fasta" %>%
    Biostrings::readDNAStringSet()

imgt_quants <- read_tsv(quant_file) %>%
    filter(grepl("^IMGT", Name)) %>%
    mutate(lineage = sub("^IMGT_([^:]+).*$", "\\1", Name), 
	   locus = sub("^([^\\*]+).+$", "\\1", lineage)) %>%
    select(locus, lineage, allele = Name, est_counts = NumReads, tpm = TPM)

top_alleles <- imgt_quants %>%
    group_by(locus) %>%
    top_n(5, est_counts) %>%
    ungroup() %>%
    group_by(locus, lineage) %>%
    filter(tpm/max(tpm) > 0.25) %>%
    ungroup() %>%
    pull(allele) %>%
    unique()

Biostrings::writeXStringSet(index[top_alleles], out)
