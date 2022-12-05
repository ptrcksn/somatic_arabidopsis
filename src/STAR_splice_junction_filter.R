library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
sj_tab_file <- args[1]
out_file <- args[2]


sj_tab <- read_tsv(sj_tab_file, col_names=F) %>%
    dplyr::select(chr=X1, left=X2, right=X3) %>%
    pivot_longer(cols=left:right, names_to="end", values_to="pos")

sj_bed <- tibble(chrom=sj_tab$chr, chromStart=(sj_tab$pos-1), chromEnd=sj_tab$pos, name=".", score=".", strand=".") %>%
    arrange(chrom, chromStart)


write_tsv(sj_bed, out_file, col_names=F)

