library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
mutations_map_file <- args[1]
out_file <- args[2]


mutations <- read_tsv(mutations_map_file)

unique_mutations <- mutations %>%
    dplyr::count(chr, pos) %>%
    filter(n == 1) %>%
    inner_join(mutations %>%
        dplyr::select(accession:context), by=c("chr", "pos")) %>%
    dplyr::select(!n) %>%
    dplyr::select(accession, chr, pos, ref, alt, context)


write_tsv(unique_mutations, out_file)
