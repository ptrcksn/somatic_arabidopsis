library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
l_args <- length(args)
gff_file <- args[1]
histone_mark_files <- args[2:(l_args-2)]
mark <- args[(l_args-1)]
out_file <- args[(l_args)]


gff <- read_tsv(gff_file, col_names=F) %>%
    dplyr::select(chr=X1, type=X3, left=X4, right=X5, strand=X7, info=X9) %>%
    filter(type == "exon") %>%
    mutate(id=info)

gff_r <- GRanges(seqnames=gff$chr, IRanges(start=gff$left, end=gff$right), id=gff$id, type=gff$type, strand=gff$strand, left=gff$left, right=gff$right)

signal_files <- histone_mark_files[which(str_detect(histone_mark_files, paste0(mark, "-")))]
mark_signal <- tibble(id=gff$id, type=gff$type, chr=gff$chr, left=gff$left, right=gff$right, signal=0)

for(signal_file in signal_files){

    signal <- read_tsv(signal_file, col_names=F, comment="#") %>%
    dplyr::select(chr=X1, start=X2, end=X3, raw_signal=X4)

    signal_r <- GRanges(seqnames=gsub("chr", "Chr", signal$chr), IRanges(start=(signal$start+1), end=signal$end), signal=(signal$raw_signal/max(signal$raw_signal)))

    signal_gene <- as_tibble(join_overlap_intersect(gff_r, signal_r)) %>%
        dplyr::rename(chr=seqnames) %>%
        group_by(id, type, chr, left, right) %>%
        summarise(signal=sum(signal*width)) %>%
        ungroup()

    mark_signal <- mark_signal %>%
        left_join(signal_gene, by=c("id", "type", "chr", "left", "right")) %>%
        mutate(signal.y=ifelse(is.na(signal.y), 0, signal.y)) %>%
        mutate(signal=signal.x+signal.y) %>%
        dplyr::select(id, type, chr, left, right, signal)

}

colnames(mark_signal)[6] <- paste(mark, "signal", sep="_")


write_tsv(mark_signal, file=out_file)
