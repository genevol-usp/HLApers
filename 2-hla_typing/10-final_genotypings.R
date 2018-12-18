devtools::load_all("~/hlaseqlib")
library(tidyverse)

samples <- readLines("../data/sample_ids.txt")

quants_noWin <-
    file.path("./quantifications_noWinner", samples, "quant.sf") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "subject") %>%
    mutate(locus = imgt_to_gname(Name),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, gene_id, locus, allele = Name, est_counts = NumReads, tpm = TPM)

typings_2nd <- quants_noWin %>%
    group_by(subject, locus) %>%
    slice(which.max(est_counts)) %>%
    ungroup() %>%
    filter(est_counts > 0)

typings_1st <- 
    read_tsv("./quantifications_topAlleles/processed_imgt_quants.tsv") %>%
    group_by(subject, locus) %>%
    slice(which.max(est_counts)) %>%
    ungroup()

typings_df <- bind_rows(typings_1st, typings_2nd) %>%
    arrange(subject, locus) %>%
    hla_genotype_dt(th = 0.01) %>%
    as_tibble()

calls <- typings_df %>%
    select(subject, locus, allele) %>%
    mutate(subject = convert_ena_ids(subject),
	   locus = sub("HLA-", "", locus),
	   allele = sub("IMGT_", "", allele)) %>%
    filter(locus %in% pag$locus)

goldstd <- mutate(pag, allele = hla_trimnames(allele, 3))

write_tsv(calc_genotyping_accuracy(calls, goldstd), "./genotyping_concordance.tsv")

typings_df %>% 
    select(subject, locus, allele) %>%
    write_tsv("./genotype_calls.tsv")
