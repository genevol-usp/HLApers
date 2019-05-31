library(hlaseqlib)
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))

# inputs
opts <- commandArgs(TRUE)
gencode_fasta <- opts[1]
gencode_annot <- opts[2]
imgt_db       <- opts[3]
out           <- opts[4]

# outputs
out_noHLA   <- file.path(out, "gencode_noHLA.fa")
out_supp    <- file.path(out, "gencode_HLAsupp.fa")
out_MHCsupp <- file.path(out, "gencode_MHC_HLAsupp.fa")

# HLA database
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

imgt_annot <- "./data/imgt_gene_types.tsv" %>%
    read_tsv(col_types = "cc") %>%
    filter(gene_type == "protein_coding" | grepl("transcribed", gene_type))

pers_index_loci <- sort(imgt_annot$gene_name) %>% sub("HLA-", "", .)

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

hladb_genes <- unique(sub("^IMGT_([^\\*]+).+$", "HLA-\\1", names(hladb)))

# Annotations
message("Processing Gencode annotations...")
g_annot <- read_tsv(gencode_annot, comment = "##", col_names = FALSE, 
		    col_types = "ccciicccc", progress = FALSE) 

tx_types_rm <- c("pseudogene", "processed_pseudogene", "unprocessed_pseudogene")

transc_annots <- g_annot %>%
    filter(X3 == "transcript") %>%
    select(chr = X1, start = X4, end = X5, X9) %>%
    mutate(i = seq_len(nrow(.)), X9 = str_split(X9, "; ")) %>%
    unnest() %>%
    filter(grepl("^gene_name|^gene_id|^transcript_id|^transcript_type", X9)) %>%
    separate(X9, c("tag", "id"), " ") %>%
    mutate(id = gsub("\"", "", id)) %>%
    spread(tag, id) %>%
    filter(! transcript_type %in% tx_types_rm) %>%
    select(-i)

mhc_coords <- transc_annots %>%
    filter(chr == "chr6", gene_name %in% hladb_genes) %>%
    summarise(start = min(start) -5e5, end = max(end) + 5e5)

mhc_coords %>%
    mutate(out = paste0("chr6:", start, "-", end)) %>%
    pull(out) %>%
    writeLines("./data/mhc_coords.txt")


# Transcript sequences
gencode <- readDNAStringSet(gencode_fasta) %>%
    `names<-`(sub("^([^\\|]+).*$", "\\1", names(.))) %>%
    .[transc_annots$transcript_id]

gencode_no_hla <- transc_annots %>%
    filter(! gene_name %in% hladb_genes) %>%
    pull(transcript_id) %>%
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

mhc_tx <- transc_annots %>%
    filter(chr == "chr6", start >= mhc_coords$start, end <= mhc_coords$end, 
	   transcript_id %in% names(gencode_no_hla)) %>%
    pull(transcript_id)

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

message("writing index files...")
writeXStringSet(gencode_hlasupp, out_supp)
writeXStringSet(gencode_no_hla_uniq, out_noHLA)
writeXStringSet(gencode_mhc_supp, out_MHCsupp)

message("Done!")
