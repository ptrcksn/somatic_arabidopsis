library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
gff_file <- args[1]
late_ratio_file <- args[2]
early_ratio_file <- args[3]
out_file <- args[4]


gff <- read_tsv(gff_file, col_names=F) %>%
    dplyr::select(chr=X1, type=X3, left=X4, right=X5, strand=X7, info=X9) %>%
    mutate(id=info) %>%
    filter(type == "exon")

gff_r <- GRanges(seqnames=gff$chr, IRanges(start=gff$left, end=gff$right), id=gff$id, type=gff$type, strand=gff$strand, left=gff$left, right=gff$right)

late <- read_tsv(late_ratio_file, col_names=F) %>%
    dplyr::rename(chr=X1, start=X2, end=X3, signal_late=X4)

early <- read_tsv(early_ratio_file, col_names=F) %>%
    dplyr::rename(chr=X1, start=X2, end=X3, signal_early=X4)

late_early <- late %>%
    inner_join(early, by=c("chr", "start", "end")) %>%
    mutate(late_to_early=signal_late/signal_early)

early_late_r <- GRanges(seqnames=paste0("Chr", late_early$chr), IRanges(start=(late_early$start+1), end=late_early$end), signal=late_early$late_to_early)

timing <- as_tibble(join_overlap_intersect(gff_r, early_late_r)) %>%
    mutate(timing_signal=(signal*(end-start+1))) %>%
    dplyr::rename(chr=seqnames) %>%
    group_by(id, type, chr, left, right) %>%
    summarise(timing_signal=sum(timing_signal)) %>%
    ungroup()


write_tsv(timing, out_file)