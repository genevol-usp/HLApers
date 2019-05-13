suppressMessages(devtools::load_all("~/Libraries/hlaseqlib"))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))

###
imgt_db <- "./IMGTHLA"
hla_db <- "./hladb"
gencode_annot <- "./gencode/gencode.v25.annotation.gtf.gz"
###

gencode_fasta <- commandArgs(TRUE)[1]
gencode_annot <- commandArgs(TRUE)[2]
imgt_db <- commandArgs(TRUE)[3]
hla_db <- commandArgs(TRUE)[4]

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

writeXStringSet(hladb, file.path(hla_db, "hladb.fasta"))

###
hladb <- readDNAStringSet("./hladb/hladb.fasta")
###

hladb_genes <- unique(imgt_to_gname(names(hladb)))


# Annotations
imgt_annot <- file.path(hla_db, "imgt_gene_types.tsv") %>%
    read_tsv() %>%
    filter(gene_type == "protein_coding" | grepl("transcribed", gene_type))

g_annot <-
    read_tsv(gencode_annot, comment = "##", col_names = FALSE, 
             col_types = "c-cii-c-c", progress = FALSE) 

transc_annots <- g_annot %>%
    filter(X3 == "transcript") %>%
    select(X9) %>%
    mutate(i = seq_len(nrow(.)), X9 = str_split(X9, "; ")) %>%
    unnest() %>%
    filter(grepl("^gene_name|^gene_id|^transcript_id|^transcript_type", X9)) %>%
    separate(X9, c("tag", "id"), " ") %>%
    mutate(id = gsub("\"", "", id)) %>%
    spread(tag, id)

tx_types_rm <- c("pseudogene", "processed_pseudogene", "unprocessed_pseudogene")

transc_annots_filter <- transc_annots %>% 
    filter(! transcript_type %in% tx_types_rm)

hla_transc <- g_annot %>%
    filter(X3 == "exon") %>%
    select(chr = X1, start = X4, end = X5, X9) %>%
    mutate(i = seq_len(nrow(.)), X9 = str_split(X9, "; ")) %>%
    unnest() %>%
    filter(grepl("^gene_name|^transcript_id|^transcript_type", X9)) %>%
    separate(X9, c("tag", "id"), " ") %>%
    mutate(id = gsub("\"", "", id)) %>%
    spread(tag, id) %>%
    filter(gene_name %in% hladb_genes, 
	   transcript_id %in% transc_annots_filter$transcript_id) %>%
    mutate(pos = map2(start, end, `:`)) %>%
    unnest()
    
hla_coding <- hla_transc %>%
    filter(transcript_type == "protein_coding")

hla_noncoding <- hla_transc %>%
    filter(transcript_type != "protein_coding")

hlatx_noncoding_keep <- hla_noncoding %>%
    group_by(tx_id) %>%
    filter(!any(pos %in% hla_coding$pos)) %>%
    distinct(gene_name, tx_id) %>%
    pull(tx_id)


# Transcript sequences

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
