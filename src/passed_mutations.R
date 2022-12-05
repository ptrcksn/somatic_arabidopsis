library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
mutations_map_file <- args[1]
out_file <- args[2]


passed_mutations <- read_tsv(mutations_map_file, col_names=F) %>%
    dplyr::rename(accession=X1, chr=X2, pos=X3, ref=X4, alt=X5, context=X6, depth=X7, alt_depth=X8, flag=X9) %>%
    filter(flag == "PASS" | flag == "clustered_mutation") %>%
    dplyr::select(accession:alt_depth)


write_tsv(passed_mutations, out_file)
