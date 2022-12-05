library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
l_args <- length(args)
gff_file <- args[1]
genome_file <- args[2]
augmented_mutations_file <- args[3]
breadth_depth_file <- args[4]
GC_file <- args[5]
replication_timing_file <- args[6]
histone_mark_signal_files <- args[7:(l_args-1)]
out_file <- args[(l_args)]


gff <- read_tsv(gff_file, col_names=F) %>%
    dplyr::select(chr=X1, type=X3, left=X4, right=X5, strand=X7, info=X9) %>%
    filter(type == "exon") %>%
    mutate(id=info) %>%
    mutate(width=(right-left+1))

fa <- readDNAStringSet(filepath=genome_file, format="fasta")
names(fa) <- paste0("Chr", names(fa))
names(fa)[6] <- "ChrM"
names(fa)[7] <- "ChrC"

gff_r <- GRanges(seqnames=gff$chr, IRanges(start=gff$left, end=gff$right), id=gff$id, type=gff$type, left=gff$left, right=gff$right)

widths <- as_tibble(oligonucleotideFrequency(getSeq(fa, gff_r), width=1)) %>%
    bind_cols(gff) %>%
    pivot_longer(cols=A:T, names_to="ref", values_to="base_width") %>%
    mutate(transcribed_status=ifelse(ref == "C" | ref == "T",
        ifelse(strand == "+", "coding", "template"),
        ifelse(strand == "+", "template", "coding"))) %>%
    mutate(ref=ifelse(ref == "C" | ref == "T", ref, as.character(reverseComplement(DNAStringSet(ref)))))

variants <- read_tsv(augmented_mutations_file)

breadth_depth <- read_tsv(breadth_depth_file, col_types=list(col_character(), col_character(), col_character(), col_integer(), col_integer(), col_character(), col_integer(), col_number(), col_character(), col_character()))

GC <- read_tsv(GC_file)

replication_timing <- read_tsv(replication_timing_file)

histone_mark_signals <- histone_mark_signal_files %>%
    map(read_tsv) %>%
    purrr::reduce(left_join, by=c("id", "type", "chr", "left", "right"))

features <- gff %>%
    left_join(widths) %>%
    left_join(variants, by=c("id", "type", "chr", "left", "right", "ref", "transcribed_status")) %>%
    left_join(breadth_depth, by=c("id", "type", "chr", "left", "right", "strand", "ref", "transcribed_status")) %>%
    left_join(GC, by=c("id", "type", "chr", "left", "right")) %>%
    left_join(replication_timing, by=c("id", "type", "chr", "left", "right")) %>%
    left_join(histone_mark_signals, by=c("id", "type", "chr", "left", "right")) %>%
    mutate(breadth=ifelse(is.na(breadth), 0, breadth)) %>%
    mutate(depth=ifelse(is.na(depth), 0, depth)) %>%
    arrange(chr, left, right, strand, ref, alt, transcribed_status)


write_tsv(features, out_file)
