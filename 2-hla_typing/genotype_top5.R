devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

top5_quants <- commandArgs(TRUE)[1]
outdir <- commandArgs(TRUE)[2]

outgenos <- file.path(outdir, "genotypes.tsv")
outwinner <- file.path(outdir, "winners.txt")

quants <- top5_quants %>%
    read_tsv() %>%
    mutate(locus = sub("^IMGT_([^\\*]+).+$", "\\1", Name),
	   lineage = sub("^IMGT_([^:]+).*$", "\\1", Name)) %>%
    select(locus, lineage, allele = Name, counts = NumReads, tpm = TPM)

lineage_typing <- quants %>%
    group_by(locus, allele = lineage) %>%
    summarise_at(vars(tpm, counts), sum) %>%
    ungroup() %>%
    hla_genotype(th = .15) %>%
    select(locus, lineage = allele) %>%
    separate_rows(lineage, sep = "-") %>%
    distinct()

quants2 <- inner_join(quants, lineage_typing)

hets <- quants2 %>%
    group_by(locus) %>%
    filter(n_distinct(lineage) > 1L) %>%
    group_by(lineage) %>%
    slice(which.max(counts)) %>%
    ungroup()

quants3 <- quants2 %>%
    group_by(locus) %>%
    filter(n_distinct(lineage) == 1L) %>%
    ungroup() %>%
    bind_rows(hets) %>%
    arrange(locus)

allele_typing <- hla_genotype(quants3, th = .15)

allele_typing %>%
    separate_rows(allele, sep = "-") %>%
    write_tsv(outgenos)

winner_alleles <- allele_typing %>%
    group_by(locus) %>%
    slice(which.max(counts)) %>%
    ungroup() %>%
    distinct(allele) %>%
    pull(allele)

writeLines(winner_alleles, outwinner)
