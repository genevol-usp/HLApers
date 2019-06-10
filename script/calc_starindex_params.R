suppressPackageStartupMessages(library(Biostrings))

opts <- commandArgs(TRUE)
index <- opts[1]
outprefix <- opts[2]

out <- file.path(outprefix, "indexparams.txt")

gen <- readDNAStringSet(index)

genlen  <- sum(width(gen))

nrefs <- length(gen)

binbits <- floor(min(18, log2(genlen/nrefs)))
saindex <- floor(min(14, log2(genlen)/2 - 1))

writeLines(as.character(c(binbits, saindex)), out)
