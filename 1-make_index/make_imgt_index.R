devtools::load_all("~/hlaseqlib")
library(tidyverse)

locus <- commandArgs(TRUE)[1]

if (grepl("DRB\\d+", locus)) { 
    locus_nuc <- "DRB"
} else {
    locus_nuc <- locus
}

nuc_file <- paste0("~/IMGTHLA/alignments/", locus_nuc, "_nuc.txt")

hla_df <- hla_read_alignment(nuc_file) %>%
   mutate(gene = sub("^([^*]+).+$", "\\1", allele)) %>%
   filter(gene == locus) %>%
   select(-gene)

if (all(grepl("\\*", hla_df$cds))) {
    stop(paste("no complete sequence for locus", locus))
}

if (nrow(hla_df) == 1L || all(!grepl("\\*", hla_df$cds))) {
    
    final_df <- hla_df

} else {

    distmatrix <- make_dist_matrix(hla_df)

    closest_allele_df <- make_closest_allele_df(distmatrix)

    closest_allele_df$id <- closest_allele_df %>% group_indices(inc_allele)

    closest_allele_df_step2 <-
	bind_rows(select(closest_allele_df, id, allele = inc_allele),
		  select(closest_allele_df, id, allele = closest)) %>%
	distinct() %>%
	left_join(hla_df, by = "allele") %>%
	split(.$id) %>%
	map(make_dist_matrix) %>%
	map(make_closest_allele_df) %>%
	bind_rows()

    closest_within_type <- find_closest_within_type(closest_allele_df_step2)

    inferred_df <- closest_within_type %>%
	left_join(hla_df, by = c("inc_allele" = "allele")) %>%
	left_join(hla_df, by = c("closest" = "allele")) %>%
	mutate(cds = map2_chr(cds.x, cds.y, hla_attribute_seq)) %>%
	select(allele = inc_allele, cds)

    final_df <- hla_df %>%
	filter(!grepl("\\*", cds)) %>%
	bind_rows(inferred_df) %>%
	arrange(allele)
}

out_df <- final_df %>%
    mutate(cds = hla_format_sequence(cds)) %>%
    rename(transcript = cds) %>%
    mutate(allele3f = hla_trimnames(allele, 3)) %>%
    distinct(allele3f, transcript, .keep_all = TRUE) %>%
    group_by(allele3f) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    mutate(allele = ifelse(n > 1L, allele, allele3f)) %>%
    select(allele, transcript) %>%
    arrange(allele)

write_tsv(out_df, paste0("./index_", locus, ".tsv"))
