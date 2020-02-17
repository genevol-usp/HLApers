library(hlaseqlib)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

opts <- commandArgs(TRUE)
outPrefix   <- opts[1]

abundance_file <- file.path(outPrefix, "abundance.tsv")
outgenos <- file.path(outPrefix, "genotypes.tsv")

genos <- abundance_file %>% 
    read_tsv() %>% 
    filter(grepl("^IMGT_", target_id)) %>%
    mutate(locus = sub("^IMGT_([^\\*]+).+$", "\\1", target_id)) %>%
    select(locus, allele = target_id, counts = est_counts, tpm) %>%
    hla_genotype(th = 0.15) %>%
    filter(!is.na(allele))

write_tsv(genos, outgenos)
