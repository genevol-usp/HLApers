devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name 

samples <- geuvadis_info %>%
    filter(pop != "YRI", kgp_phase3 == 1L) %>%
    pull(ena_id)

imgt_quants <- read_tsv("./quantifications_MHC/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(Name),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, locus, gene_id, allele = Name, 
	   est_counts = NumReads, tpm = TPM)

top_alleles <- imgt_quants %>%
    group_by(subject, gene_id) %>%
    top_n(5, est_counts) %>%
    ungroup() %>%
    mutate(lineage = hla_trimnames(sub("IMGT_", "", allele), 1)) %>%
    group_by(subject, locus, lineage) %>%
    filter(tpm/max(tpm) > 0.25) %>%
    ungroup()

write_tsv(top_alleles, "./quantifications_MHC/imgt_quants_topAlleles.tsv")
