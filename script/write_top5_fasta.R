suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

opts <- commandArgs(TRUE)
quant_file <- opts[1]
gencode    <- opts[2]
depth      <- opts[3]
out        <- opts[4]

index <- Biostrings::readDNAStringSet(gencode)

depth_df <- read_tsv(depth, col_names = FALSE) %>%
    group_by(allele = X1) %>%
    summarise(s = sum(X3), d = mean(X3 >= 10)) %>%
    ungroup()

imgt_quants <- read_tsv(quant_file) %>%
    filter(grepl("^IMGT", Name)) %>%
    mutate(lineage = sub("^IMGT_([^:]+).*$", "\\1", Name), 
	   locus = sub("^([^\\*]+).+$", "\\1", lineage)) %>%
    select(locus, lineage, allele = Name, est_counts = NumReads, tpm = TPM)

top_alleles <- imgt_quants %>% 
    group_by(locus) %>%
    top_n(5, est_counts) %>%
    ungroup() %>%
    left_join(depth_df, by = c("allele")) %>%
    group_by(locus, lineage) %>%
    filter(s/max(s) >= .95) %>%
    filter(est_counts/max(est_counts) > 0.2) %>%
    group_by(locus) %>%
    filter(d/max(d) >= .7) %>%
    ungroup() %>%
    pull(allele) %>%
    unique()

Biostrings::writeXStringSet(index[top_alleles], out)
