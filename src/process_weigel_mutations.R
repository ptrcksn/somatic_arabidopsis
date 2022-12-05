library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
weigel_raw_mutations_file <- args[1]
out_file <- args[2]


weigel_mutations <- read_tsv(weigel_raw_mutations_file)

weigel_mutations <- weigel_mutations %>%
    inner_join(weigel_mutations %>%
        count(CHROM, POS, TYPE, src) %>%
        filter(n == 1) %>%
        dplyr::select(!n), by=c("CHROM", "POS", "TYPE", "src")) %>%
    filter(src  == "EU_som" &
        TYPE == "SNP" &
        REF %in% c("A","C","G","T") &
        ALT %in% c("A","C","G","T"))


write_tsv(weigel_mutations, out_file)
