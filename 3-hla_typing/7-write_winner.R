devtools::load_all("~/hlaseqlib")
library(tidyverse)

quants <- read_tsv("./quantifications_topAlleles/processed_imgt_quants.tsv")

winner_alleles <- quants %>%
    group_by(subject, locus) %>%
    slice(which.max(est_counts)) %>%
    ungroup() %>%
    select(subject, locus, allele)

winner_list <- winner_alleles %>%
    select(subject, allele) %>%
    mutate(path = file.path("./quantifications_topAlleles", subject, "winner.txt")) %>%
    split(.$subject)

map(winner_list, ~writeLines(.$allele, unique(.$path)))
