library(Biostrings)
devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_tx <-
    read_tsv("./imgt.annotations.gtf",
             comment = "##", col_names = FALSE,
             col_types = "c-cii-c-c", progress = FALSE) %>%
    filter(X3 == "exon") %>%
    select(X4, X5, X9) %>%
    mutate(i = seq_len(nrow(.)), X9 = str_split(X9, "; ")) %>%
    unnest() %>%
    filter(grepl("^gene_name|^transcript_id|^transcript_type", X9)) %>%
    separate(X9, c("info", "value"), " ") %>%
    mutate(value = gsub("\"", "", value)) %>%
    group_by(i, info, X4, X5) %>%
    summarize(value = paste(value, collapse = ";")) %>%
    ungroup() %>%
    spread(info, value) %>%
    mutate(pos = map2(X4, X5, `:`)) %>%
    unnest()

hla_coding <- hla_tx %>%
    filter(transcript_type == "protein_coding")

hla_non_coding <- hla_tx %>%
    filter(transcript_type != "protein_coding")

tx_to_keep <- hla_non_coding %>%
    group_by(transcript_id) %>%
    filter(!any(pos %in% hla_coding$pos)) %>%
    distinct(gene_name, transcript_id) %>%
    pull(transcript_id)

imgt <- readDNAStringSet("./imgt_index.fa") 

loci_in_index <- readLines("./imgt_loci.txt") %>% paste0("HLA-", .)

gencode <- 
  readDNAStringSet("./gencode/gencode.v25.PRI.uniqTranscripts.fa")

gencode_pri_imgt <- readLines("./imgt_pri_ids.txt") %>%
    .[! . %in% tx_to_keep]
 
gencode_no_imgt <- gencode[! names(gencode) %in% gencode_pri_imgt] 
  
gencode_imgt <- c(gencode_no_imgt, imgt) 

writeXStringSet(gencode_imgt, "./gencode.v25.PRI.IMGT.transcripts.fa")
writeXStringSet(gencode_no_imgt, "./gencode.v25.PRI.transcripts.noIMGT.fa")
