devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

quants_1st <- commandArgs(TRUE)[1]
quants_2nd <- commandArgs(TRUE)[2]
outPrefix <- commandArgs(TRUE)[3] 

outgenos <- paste0(outPrefix, "_genotypes.tsv")
outindex <- paste0(outPrefix, "_index.fa")

index <- "../1-make_indices/data/hladb/hladb.fasta" %>%
    Biostrings::readDNAStringSet()

typings_1st <- quants_1st %>% 
    read_tsv() %>%
    mutate(locus = sub("^IMGT_([^\\*]+).+$", "\\1", Name)) %>%
    select(locus, allele = Name, counts = NumReads, tpm = TPM) %>%
    group_by(locus) %>%
    slice(which.max(counts)) %>%
    ungroup()

typings_2nd <- quants_2nd %>%
    read_tsv() %>%
    mutate(locus = sub("^IMGT_([^\\*]+).+$", "\\1", Name)) %>%
    select(locus, allele = Name, counts = NumReads, tpm = TPM) %>% 
    group_by(locus) %>%
    slice(which.max(counts)) %>%
    ungroup() %>%
    filter(counts > 0)

typings_df <- bind_rows(typings_1st, typings_2nd) %>%
    arrange(locus) %>%
    hla_genotype(th = 0.05)

typings_df %>% 
    select(locus, allele) %>%
    write_tsv(outgenos)

alleles <- typings_df %>%
    pull(allele) %>%
    unique()

Biostrings::writeXStringSet(index[alleles], outindex)
