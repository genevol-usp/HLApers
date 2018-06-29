library(Biostrings)

index <- readDNAStringSet("./gencode/gencode.v25.PRI.transcripts.fa")

writeXStringSet(unique(index), "./gencode/gencode.v25.PRI.uniqTranscripts.fa")
