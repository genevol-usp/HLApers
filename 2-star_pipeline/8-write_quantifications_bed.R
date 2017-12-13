devtools::load_all("~/hlaseqlib")
library(data.table)

setDTthreads(1)
setDT(geuvadis_info)
setDT(gencode_chr_tx)
setDT(gencode_chr_gene)

gencode <- 
    gencode_chr_tx[chr %in% 1:22, .(target_id = tx_id, gene_id, gene_name)]

expression_dt <-
    fread("zcat < ./quantifications_2/all_transcripts_quants.tsv.gz")

samples <- readLines("../data/sample_ids.txt")

target_set <- expression_dt[, .(target_id = unique(target_id))]

autosomes_set <- 
    target_set[gencode, on = .(target_id), nomatch = 0L
	     ][, .(target_id, gene_id)]

imgt_set <-
    target_set[grepl("^IMGT", target_id)
	     ][, gene_id := imgt_to_gid(target_id)
	     ][gencode, on = .(gene_id), nomatch = 0L
	     ][order(target_id), .(target_id, gene_id)]

gene_set <- rbind(autosomes_set, imgt_set)

gene_dt <- 
    expression_dt[gene_set, on = .(target_id), nomatch = 0L
		][, .(gene_tpm = sum(tpm)), by = .(subject, gene_id)]

expressedGenes <- gene_dt[, .(mean(gene_tpm > 0.1)), by = .(gene_id)][V1 >= 0.5]

gene_dt <- gene_dt[gene_id %in% expressedGenes$gene_id]

gene_dt_wide <- dcast(gene_dt, gene_id ~ subject, value.var = "gene_tpm")

gene_bed <- 
    gene_dt_wide[gencode_chr_gene, on = .(gene_id), nomatch = 0L
	       ][, `:=`(gid = gene_id, gene_name = NULL)]

setcolorder(gene_bed, 
	    c("chr", "start", "end", "gene_id", "gid", "strand", samples))

setnames(gene_bed, c("chr", "gene_id", "strand"), c("#chr", "id", "strd"))

fwrite(gene_bed, "quantifications_expressed50%.bed", sep = "\t")
