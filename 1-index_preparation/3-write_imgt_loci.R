library(tidyverse)

loci <-
  list.files("~/IMGTHLA/alignments/", pattern = "_nuc\\.txt") %>% 
  strsplit("_") %>% 
  map_chr(1) %>%
  .[!grepl("Class|HFE|MIC|TAP|DRB", .)]

loci <- c(loci, paste0("DRB", 1:9)) %>% sort()

writeLines(loci, "./imgt_loci.txt")
