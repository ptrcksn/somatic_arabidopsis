library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
gff_file <- args[1]
fasta_file <- args[2]
out_file <- args[3]


gff <- read_tsv(gff_file, col_names=F) %>%
    dplyr::select(chr=X1, type=X3, left=X4, right=X5, strand=X7, info=X9) %>%
    mutate(id=info) %>%
    filter(type == "exon")

gff_r <- GRanges(seqnames=gff$chr, IRanges(start=gff$left, end=gff$right), id=gff$id, type=gff$type, left=gff$left, right=gff$right)

fa <- readDNAStringSet(filepath=fasta_file, format="fasta")
names(fa) <- paste0("Chr", names(fa))
names(fa)[6] <- "ChrM"
names(fa)[7] <- "ChrC"

gc <- gff %>%
    mutate(GC=rowSums(oligonucleotideFrequency(getSeq(fa, gff_r), width=1)[, c("C","G")])) %>%
    dplyr::select(id, type, chr, left, right, GC)


write_tsv(gc, out_file)