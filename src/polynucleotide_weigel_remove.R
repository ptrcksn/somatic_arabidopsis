library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
weigel_mutations_file <- args[1]
polynucleotide_bed_file <- args[2]
out_file <- args[3]


weigel_mutations <- read_csv(weigel_mutations_file)

weigel_mutations_r <- GRanges(seqnames=weigel_mutations$CHROM, IRanges(start=weigel_mutations$POS, end=weigel_mutations$POS), type=weigel_mutations$TYPE, ref=weigel_mutations$REF, alt=weigel_mutations$ALT, src=weigel_mutations$src)

polynucleotide <- read_tsv(polynucleotide_bed_file, col_names=F) %>%
    dplyr::select(chrom=X1, chromStart=X2, chromEnd=X3)

polynucleotide_r <- GRanges(seqnames=polynucleotide$chrom, IRanges(start=polynucleotide$chromStart, end=(polynucleotide$chromEnd-1)))

weigel_mutations <- as_tibble(GRanges(as_tibble(weigel_mutations_r) %>%
    anti_join(as_tibble(join_overlap_intersect(polynucleotide_r, weigel_mutations_r)), by=c("seqnames", "start", "end", "ref", "alt", "src")))) %>%
    dplyr::select(CHROM=seqnames, POS=start, TYPE=type, REF=ref, ALT=alt, src=src)


write_tsv(weigel_mutations, out_file)