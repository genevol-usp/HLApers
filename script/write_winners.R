library(tidyverse)

top5_quants <- commandArgs(TRUE)[1]
out <- commandArgs(TRUE)[2]

quants <- top5_quants %>%
    read_tsv() %>%
    mutate(locus = sub("^IMGT_([^\\*]+).+$", "\\1", Name),
	   lineage = sub("^IMGT_([^:]+).*$", "\\1", Name)) %>%
    select(locus, lineage, allele = Name, counts = NumReads, tpm = TPM)

winner_alleles <- quants %>%
    group_by(locus) %>%
    slice(which.max(counts)) %>%
    ungroup() %>%
    distinct(allele) %>%
    pull(allele)

writeLines(winner_alleles, out)
