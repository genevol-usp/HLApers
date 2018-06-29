devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name 

samples <- geuvadis_info %>% 
    filter(kgp_phase3 == 1L & pop != "YRI") %>%
    pull(ena_id)

imgt_quants <- read_tsv("./quantifications/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(Name),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, locus, gene_id, allele = Name,
	   est_counts = NumReads, tpm = TPM)

missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

out_df <- hla_genotype_dt(imgt_quants, th = 0) %>%
    hla_apply_zigosity_threshold(th = 0.1)

write_qsv(out_df, "./quantifications/processed_imgt_quants.tsv")
