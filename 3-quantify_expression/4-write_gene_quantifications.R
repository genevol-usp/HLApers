devtools::load_all("~/hlaseqlib")
library(tidyverse)

quant_dir <- "./quantifications"

gencode <- gencode_chr_tx %>%
    filter(chr %in% 1:22) %>%
    select(target_id = tx_id, gene_id, gene_name)

samples <- geuvadis_info %>%
    filter(kgp_phase3 == 1, pop != "YRI") %>%
    select(name, subject = ena_id)

expression_df <- file.path(quant_dir, samples$subject, "quant.sf") %>%
    setNames(samples$subject) %>%
    map_df(read_tsv, .id = "subject") %>%
    left_join(samples, by = "subject") %>%
    select(subject = name, target_id = Name, tpm = TPM)

autosomes_set <- inner_join(expression_df, gencode, by = c("target_id")) %>%
    select(subject, gene_id, tpm)

imgt_set <- expression_df %>%
    filter(grepl("^IMGT", target_id)) %>%
    mutate(gene_name = sub("^IMGT_([^\\*]+).+$", "HLA-\\1", target_id)) %>%
    inner_join(distinct(gencode, gene_id, gene_name), by = "gene_name") %>%
    select(subject, gene_id, tpm) %>%
    complete(subject, gene_id, fill = list(tpm = 0))

gene_df <- bind_rows(autosomes_set, imgt_set) %>%
    group_by(subject, gene_id) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup()

write_tsv(gene_df, "gene_quantifications.tsv")
