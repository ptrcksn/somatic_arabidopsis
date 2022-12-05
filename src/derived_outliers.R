library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
mutations_map_file <- args[1]
out_file <- args[2]


mutations <- read_tsv(mutations_map_file)

accession_mutation_count <- mutations %>%
    dplyr::count(accession)

quantiles <- quantile(accession_mutation_count$n)
iqr <- as.numeric(quantiles[4])-as.numeric(quantiles[2])
lower <- as.numeric(quantiles[2]) - 1.5*iqr
upper <- as.numeric(quantiles[4]) + 1.5*iqr

outliers <- accession_mutation_count %>%
    filter(n < lower | n > upper) %>%
    pull(accession)

write(outliers, out_file, ncolumns=1)
