devtools::load_all("~/hlaseqlib")
library(Biostrings)
library(tidyverse)

main_loci <- c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB")

other_loci <- list.files("~/IMGTHLA/alignments/", pattern = "_nuc\\.txt") %>% 
    strsplit("_") %>% 
    map_chr(1) %>%
    .[! . %in% c(main_loci, "ClassI", "ClassII")]

loci_df <- tibble(locus = c(main_loci, other_loci)) %>%
    mutate(infer = locus %in% main_loci,
	   seqs = map2(locus, infer, hla_make_sequences, n_cores = 1)) 

seqs_df <- select(loci_df, seqs) %>% unnest()

index_df <- seqs_df %>%
    mutate(cds = hla_format_sequence(cds)) %>%
    group_by(cds) %>%
    summarize(allele = paste(allele, collapse = "/")) %>%
    ungroup() %>%
    mutate(allele3f = hla_trimnames(allele)) %>%
    group_by(allele3f) %>%
    mutate(n = n()) %>%
    ungroup() %>% 
    arrange(allele)

index <- index_df %>%
    mutate(allele = paste0("IMGT_", ifelse(n > 1L, allele, allele3f))) %>%
    select(allele, cds) %>%
    split(.$allele) %>%
    map_chr("cds") %>%
    DNAStringSet()

writeXStringSet(index, "./imgt_index.fa")	    
