library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
mutations_map_file <- args[1]
outliers_file <- args[2]
out_file <- args[3]


mutations <- read_tsv(mutations_map_file) %>%
    mutate(accession=as.character(accession))

outliers <- as.character(scan(outliers_file))

mutations_outlier_removed <- mutations %>%
    filter(!(accession %in% outliers)) %>%
    dplyr::select(accession, chr, pos, ref, alt, context)


write_tsv(mutations_outlier_removed, out_file)
