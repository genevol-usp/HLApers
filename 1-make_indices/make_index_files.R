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

annot <- gencode_all_gene %>%
    filter(gene_name %in% imgt_genes, 
	   gene_type == "protein_coding" | grepl("transcribed", gene_type)) %>%
    distinct(gene_name, .keep_all = TRUE)

pers_index_loci <- sort(annot$gene_name) %>% sub("HLA-", "", .)

hladb <- tibble(locus = pers_index_loci) %>%
    mutate(data = map(locus, ~hla_compile_index(., imgt_db))) %>%
    filter(!is.na(data)) %>%
    unnest(data) %>%
    filter(!grepl("N$", allele)) %>%
    select(-locus) %>%
    mutate(allele = paste0("IMGT_", allele)) %>%
    split(.$allele) %>%
    map_chr("cds") %>%
    DNAStringSet()

writeXStringSet(hladb, "./data/hladb/hladb.fasta")

hladb_genes <- unique(imgt_to_gname(names(hladb)))

message("Processing Gencode annotations...")
tx_types_rm <- c("pseudogene", "processed_pseudogene", "unprocessed_pseudogene")

filtered_annots <- gencode_pri_tx %>% filter(! tx_type %in% tx_types_rm)

hla_tx <- "./data/gencode/gencode.v25.primary_assembly.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "exon") %>%
    filter(gene_name %in% hladb_genes, tx_id %in% filtered_annots$tx_id) %>%
    mutate(pos = map2(start, end, `:`)) %>%
    unnest()
    
hla_coding <- hla_tx %>%
    filter(tx_type == "protein_coding")

hla_noncoding <- hla_tx %>%
    filter(tx_type != "protein_coding")

hlatx_noncoding_keep <- hla_noncoding %>%
    group_by(tx_id) %>%
    filter(!any(pos %in% hla_coding$pos)) %>%
    distinct(gene_name, tx_id) %>%
    pull(tx_id)

gencode <- readDNAStringSet("./data/gencode/gencode.v25.PRI.transcripts.fa") %>%
    .[filtered_annots$tx_id]

gencode_no_hla <- filtered_annots %>%
    filter(! gene_name %in% hladb_genes | tx_id %in% hlatx_noncoding_keep) %>%
    pull(tx_id) %>%
    gencode[.]

gencode_no_hla_uniq <-
    tibble(tx_id = names(gencode_no_hla), 
	   cds = as.character(gencode_no_hla)) %>%
    group_by(cds) %>%
    summarise(tx_id = paste(tx_id, collapse = "-")) %>%
    ungroup() %>%
    split(.$tx_id) %>%
    map_chr("cds") %>%
    DNAStringSet()

gencode_hlasupp <- c(gencode_no_hla_uniq, hladb) 

message("Writing index files...")

writeXStringSet(gencode_hlasupp, "./data/gencode/gencode.v25.PRI.IMGT.transcripts.fa")
writeXStringSet(gencode_no_hla_uniq, "./data/gencode/gencode.v25.PRI.transcripts.noIMGT.fa")

mhc_coords <- filtered_annots %>%
    filter(gene_name %in% hladb_genes) %>%
    summarise(start = min(start) -5e5, end = max(end) + 5e5)

mhc_coords %>%
    mutate(out = paste0("chr6:", start, "-", end)) %>%
    pull(out) %>%
    writeLines("mhc_coords.txt")

mhc_tx <- filtered_annots %>%
    filter(chr == 6, start >= mhc_coords$start, end <= mhc_coords$end, 
	   tx_id %in% names(gencode_no_hla)) %>%
    pull(tx_id)

gencode_mhc <-
    tibble(tx_id = names(gencode_no_hla), 
	   cds = as.character(gencode_no_hla)) %>%
    filter(tx_id %in% mhc_tx) %>%
    group_by(cds) %>%
    summarise(tx_id = paste(tx_id, collapse = "-")) %>%
    ungroup() %>%
    split(.$tx_id) %>%
    map_chr("cds") %>%
    DNAStringSet()

gencode_mhc_supp <- c(gencode_mhc, hladb)

writeXStringSet(gencode_mhc_supp, "./data/gencode/gencode.v25.MHC.IMGT.transcripts.fa")
message("Done!")
