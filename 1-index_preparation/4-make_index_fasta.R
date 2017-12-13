library(Biostrings)
devtools::load_all("~/hlaseqlib")
library(tidyverse)

imgt <- readDNAStringSet("./imgt_index.fa") 

gencode <- readDNAStringSet("./gencode/gencode.v25.PRI.transcripts.fa")
  
names(gencode) <- sub("^([^|]+).*$", "\\1", names(gencode))  

loci_in_index <- unique(imgt_to_gname(names(imgt)))

gencode_pri_imgt <-
    gencode_pri_tx %>%
    filter(gene_name %in% loci_in_index) %>%
    arrange(gene_name, tx_id)

gencode_no_imgt <- gencode[! names(gencode) %in% gencode_pri_imgt$tx_id] 
  
gencode_imgt <- c(gencode_no_imgt, imgt) 

writeXStringSet(gencode_imgt, "./gencode.v25.PRI.IMGT.transcripts.fa")
writeXStringSet(gencode_no_imgt, "./gencode.v25.PRI.transcripts.noIMGT.fa")
