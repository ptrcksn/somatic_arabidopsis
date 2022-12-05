library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
genome_gff_file <- args[1]
out_file <- args[2]


exon_junctions <- read_tsv(genome_gff_file, col_names=F) %>%
    dplyr::select(chr=X1, type=X3, left=X4, right=X5, strand=X7, info=X9) %>%
    filter(type == "exon") %>%
    pivot_longer(cols=left:right, names_to="end", values_to="pos")

exon_junctions_bed <- tibble(chrom=str_remove(exon_junctions$chr, "Chr"), chromStart=(exon_junctions$pos-1), chromEnd=exon_junctions$pos, name=".", score=".", strand=".") %>%
    arrange(chrom, chromStart)


write_tsv(exon_junctions_bed, out_file, col_names=F)
