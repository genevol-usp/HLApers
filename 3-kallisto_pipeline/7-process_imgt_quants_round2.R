devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name

samples <- readLines("../data/sample_ids.txt")

imgt_quants <- read_tsv("./quantifications_2/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(target_id),
	   gene_id = gname_to_gid(locus)) %>%
    select(locus, gene_id, allele = target_id, est_counts, tpm)
  
missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

out_df <- hla_genotype_dt(imgt_quants, th = 0) %>%
    hla_apply_zigosity_threshold(th = 0.2)

write_tsv(out_df, "./quantifications_2/processed_imgt_quants.tsv") 
