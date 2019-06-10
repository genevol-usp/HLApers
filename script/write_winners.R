suppressPackageStartupMessages(library(tidyverse))

opts <- commandArgs(TRUE)
top5_quants <- opts[1]
out         <- opts[2]

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
