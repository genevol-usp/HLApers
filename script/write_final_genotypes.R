library(hlaseqlib)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

opts <- commandArgs(TRUE)
transcripts<- opts[1]
quants_1st <- opts[2]
quants_2nd <- opts[3]
outPrefix  <- opts[4] 

outgenos <- paste0(outPrefix, "_genotypes.tsv")

index <- Biostrings::readDNAStringSet(transcripts)

typings_1st <- quants_1st %>% 
    read_tsv() %>%
    mutate(locus = sub("^IMGT_([^\\*]+).+$", "\\1", Name)) %>%
    select(locus, allele = Name, counts = NumReads, tpm = TPM) %>%
    group_by(locus) %>%
    slice(which.max(counts)) %>%
    ungroup()

if (file.exists(quants_2nd)) {

    typings_2nd <- quants_2nd %>%
	read_tsv() %>%
	mutate(locus = sub("^IMGT_([^\\*]+).+$", "\\1", Name)) %>%
	select(locus, allele = Name, counts = NumReads, tpm = TPM) %>% 
	group_by(locus) %>%
	slice(which.max(counts)) %>%
	ungroup() %>%
	filter(counts > 0)

    typings_df <- bind_rows(typings_1st, typings_2nd) %>%
	arrange(locus) %>%
	hla_genotype(th = 0.05)

} else {
    
    typings_df <- typings_1st %>%
	arrange(locus) %>%
	hla_genotype(th = 0.05)
}

typings_df %>%
    filter(!is.na(allele)) %>%
    select(locus, allele) %>%
    write_tsv(outgenos)
