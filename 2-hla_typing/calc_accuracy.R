library(tidyverse)
devtools::load_all("/home/vitor/Libraries/hlaseqlib")

sample_ids <- select(geuvadis_info, subject = name, ena_id)

samples <- readLines("../map_to_genome/samples.txt")

genos <- paste0("/home/vitor/mappings/", samples, "_genotypes.tsv") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "ena_id") %>%
    filter(locus %in% c("A", "B", "C", "DQB1", "DRB1")) %>%
    mutate(allele = gsub("IMGT_", "", allele)) %>%
    left_join(sample_ids, by = "ena_id") %>%
    select(subject, locus, allele) %>%
    arrange(subject, locus, allele)

acc <- calc_genotyping_accuracy(genos, polypheme_pag)

write_tsv(acc, "./typing_accuracy.tsv")
