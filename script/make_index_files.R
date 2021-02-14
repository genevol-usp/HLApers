library(hlaseqlib)
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# inputs
opts <- commandArgs(TRUE)
transcript_fasta <- opts[1]
transcript_annot <- opts[2]
imgt_db       <- opts[3]
out           <- opts[4]

# outputs
out_noHLA   <- file.path(out, "transcripts_noHLA.fa")
out_supp    <- file.path(out, "transcripts_HLAsupp.fa")
out_MHCsupp <- file.path(out, "transcripts_MHC_HLAsupp.fa")
out_coord   <- file.path(out, "mhc_coords.txt") 

if (!file.exists(out)) dir.create(out)

# HLA database
imgt_loci <- c("A", "B", "C", "E", "F", "G", "H",
	       "DMA", "DMB", "DOA", "DOB",
	       "DPA1", "DPA2", "DPB1", "DPB2", 
	       "DQA1", "DQA2", "DQB1", 
	       "DRA", "DRB1", "DRB3", "DRB4", "DRB5")

hladb <- tibble(locus = imgt_loci) %>%
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
message("Reading transcript annotations...")
g_annot <- read_tsv(transcript_annot, comment = "#", col_names = FALSE, 
		    col_types = "ccciicccc", progress = FALSE) 

transcripts_db <- g_annot %>%
    filter(X3 == "transcript") %>%
    mutate(gene_name = sub("^.*gene_name \"([^\"]+)\";.*$", "\\1", X9),
	   gene_id = sub("^.*gene_id \"([^\"]+)\";.*$", "\\1", X9),
	   transcript_id = sub("^.*transcript_id \"([^\"]+)\";.*$", "\\1", X9)) %>%
    select(chr = X1, start = X4, end = X5, gene_name, gene_id, transcript_id)

mhc_coords <- transcripts_db %>%
    filter(chr == "chr6" | chr == 6, gene_name %in% hladb_genes) %>%
    summarise(chr = unique(chr), start = min(start) -5e5, end = max(end) + 5e5)

mhc_coords %>%
    mutate(out = paste0(chr, ":", start, "-", end)) %>%
    pull(out) %>%
    writeLines(out_coord)

# Transcript sequences
transcripts <- readDNAStringSet(transcript_fasta) %>%
    `names<-`(sub("^([^\\|]+).*$", "\\1", names(.))) 

transcripts_no_hla <- transcripts_db %>%
    filter(! gene_name %in% hladb_genes) %>%
    pull(transcript_id) %>%
    transcripts[.]

transcripts_no_hla_uniq <-
    tibble(tx_id = names(transcripts_no_hla), 
	   cds = as.character(transcripts_no_hla)) %>%
    group_by(cds) %>%
    summarise(tx_id = paste(tx_id, collapse = "-")) %>%
    ungroup() %>%
    split(.$tx_id) %>%
    map_chr("cds") %>%
    DNAStringSet()

transcripts_hlasupp <- c(transcripts_no_hla_uniq, hladb) 

mhc_transc_ids <- transcripts_db %>%
    filter(chr == "chr6" | chr == 6, start >= mhc_coords$start, end <= mhc_coords$end, 
	   transcript_id %in% names(transcripts_no_hla)) %>%
    pull(transcript_id)

transcripts_mhc <-
    tibble(tx_id = names(transcripts_no_hla), 
	   cds = as.character(transcripts_no_hla)) %>%
    filter(tx_id %in% mhc_transc_ids) %>%
    group_by(cds) %>%
    summarise(tx_id = paste(tx_id, collapse = "-")) %>%
    ungroup() %>%
    split(.$tx_id) %>%
    map_chr("cds") %>%
    DNAStringSet()

transcripts_mhc_supp <- c(transcripts_mhc, hladb)

message("writing index files...")
writeXStringSet(transcripts_hlasupp, out_supp)
writeXStringSet(transcripts_no_hla_uniq, out_noHLA)
writeXStringSet(transcripts_mhc_supp, out_MHCsupp)

message("Done!")
