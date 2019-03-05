suppressMessages(devtools::load_all("~/Libraries/hlaseqlib"))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))

imgt_db <- "./data/IMGTHLA"

imgt_loci <-
    file.path(imgt_db, "alignments") %>%
    list.files(pattern = "_nuc\\.txt") %>% 
    strsplit("_") %>% 
    map_chr(1) %>%
    .[!. %in% c("ClassI", "ClassII", "HFE", "MICA", "MICB", "TAP1", "TAP2")]
    
imgt_genes <- imgt_loci %>%
    .[. != "DRB"] %>%
    c(paste0("DRB", 1:9)) %>% 
    paste0("HLA-", .) %>%
    sort()

annot <- gencode_chr_gene %>%
    filter(gene_name %in% imgt_genes, 
	   gene_type == "protein_coding" | grepl("transcribed", gene_type))

pers_index_loci <- sort(annot$gene_name) %>% sub("HLA-", "", .)

hladb <- tibble(locus = pers_index_loci) %>%
    mutate(data = map(locus, ~hla_compile_index(., imgt_db))) %>%
    filter(!is.na(data)) %>%
    unnest(data) %>%
    select(-locus) %>%
    mutate(allele = paste0("IMGT_", allele)) %>%
    split(.$allele) %>%
    map_chr("cds") %>%
    DNAStringSet()

writeXStringSet(hladb, "./data/hladb/hladb.fasta")

hladb_genes <- unique(imgt_to_gname(names(hladb)))

message("Processing Gencode annotations...")
hla_tx <- "./data/gencode/gencode.v25.primary_assembly.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "exon") %>%
    filter(gene_name %in% hladb_genes) %>%
    mutate(pos = map2(start, end, `:`)) %>%
    unnest()
    
hla_coding <- hla_tx %>%
    filter(tx_type == "protein_coding")

hla_non_coding <- hla_tx %>%
    filter(tx_type != "protein_coding")

hlatx_to_keep <- hla_non_coding %>%
    group_by(tx_id) %>%
    filter(!any(pos %in% hla_coding$pos)) %>%
    distinct(gene_name, tx_id) %>%
    pull(tx_id)

tx_types_rm <- c("pseudogene", "processed_pseudogene", "unprocessed_pseudogene")

filtered_annots <- gencode_pri_tx %>% filter(! tx_type %in% tx_types_rm)

gencode <- readDNAStringSet("./data/gencode/gencode.v25.PRI.transcripts.fa") %>%
    .[names(.) %in% filtered_annots$tx_id] %>%
    unique()

hladb_tx_rm <- gencode_pri_tx %>%
    filter(gene_name %in% hladb_genes, ! tx_id %in% hlatx_to_keep)
 
gencode_no_hla <- gencode[! names(gencode) %in% hladb_tx_rm$tx_id] 
  
gencode_hlasupp <- c(gencode_no_hla, hladb) 

message("Writing index files...")

writeXStringSet(gencode_hlasupp, "./data/gencode/gencode.v25.PRI.IMGT.transcripts.fa")
writeXStringSet(gencode_no_hla, "./data/gencode/gencode.v25.PRI.transcripts.noIMGT.fa")

mhc_coords <- gencode_pri_gene %>%
    filter(gene_name %in% hladb_genes) %>%
    summarise(start = min(start) -1e6, end = max(end) + 1e6)

mhc_coords %>%
    mutate(out = paste0("chr6:", start, "-", end)) %>%
    pull(out) %>%
    writeLines("mhc_coords.txt")

mhc_tx <- gencode_pri_tx %>%
    filter(start >= mhc_coords$start, end <= mhc_coords$end, 
	   !gene_name %in% hladb_genes) %>%
    pull(tx_id)

gencode_mhc <- gencode_hlasupp %>%
    .[grepl("IMGT", names(.)) | names(.) %in% mhc_tx]

writeXStringSet(gencode_mhc, "./data/gencode/gencode.v25.MHC.IMGT.transcripts.fa")
message("Done!")
