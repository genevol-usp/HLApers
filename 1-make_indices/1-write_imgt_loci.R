devtools::load_all("~/hlaseqlib")
library(tidyverse)

loci <-
    list.files("~/IMGTHLA/alignments/", pattern = "_nuc\\.txt") %>% 
    strsplit("_") %>% 
    map_chr(1) %>%
    .[!grepl("Class|HFE|MIC|TAP|DRB", .)]

loci <- c(loci, paste0("DRB", 1:9)) %>% sort() %>% paste0("HLA-", .)

annot <- "./gencode/gencode.v25.primary_assembly.annotation.gtf.gz" %>%
    read_tsv(comment = "##", col_names = FALSE, col_types = "c-cii-c-c",
	     progres = FALSE) %>%
    filter(X3 == "gene") %>%
    select(X9) %>%
    mutate(i = seq_len(nrow(.)), X9 = str_split(X9, "; ")) %>%
    unnest() %>%
    filter(grepl("^gene_name|^gene_type", X9)) %>%
    separate(X9, c("info", "value"), " ") %>%
    mutate(value = gsub("\"", "", value)) %>%
    spread(info, value) %>%
    filter(gene_name %in% loci, gene_type != "protein_coding")

coding_loci <- loci[! loci %in% annot$gene_name & loci != "HLA-Y"]

sub("HLA-", "", coding_loci) %>%
writeLines("./imgt_loci.txt")

gencode_pri_tx %>%
    filter(gene_name %in% coding_loci) %>%
    pull(tx_id) %>%
    writeLines("./imgt_pri_ids.txt")
