library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
gff_file <- args[1]
processed_mutations_file <- args[2]
processed_weigel_mutations_file <- args[3]
out_file <- args[4]


gff <- read_tsv(gff_file, col_names=F) %>%
    dplyr::select(chr=X1, type=X3, left=X4, right=X5, strand=X7, info=X9) %>%
    mutate(id=info) %>%
    filter(type == "exon")

gff_r <- GRanges(seqnames=gff$chr, IRanges(start=gff$left, end=gff$right), id=gff$id, type=gff$type, strand=gff$strand, left=gff$left, right=gff$right)

processed_weigel_mutations <- read_tsv(processed_weigel_mutations_file)

processed_mutations <- read_tsv(processed_mutations_file)

processed_mutations <- processed_mutations %>%
    bind_rows(processed_weigel_mutations)

processed_mutations_r <- GRanges(paste0("Chr", seqnames=processed_mutations$CHROM), IRanges(start=processed_mutations$POS, end=processed_mutations$POS), ref=processed_mutations$REF, alt=processed_mutations$ALT, src=processed_mutations$src)

augmented_mutations <- crossing(aggregate=paste(gff$id, gff$type, gff$chr, gff$left, gff$right, sep=","),
        ref=c("C", "T"),
        alt=c("A", "C", "G", "T"),
        transcribed_status=c("coding", "template"),
        src=c("epigenomes", "EU_som")) %>%
    filter(ref != alt) %>%
    left_join(as_tibble(join_overlap_intersect(gff_r, processed_mutations_r)) %>%
        dplyr::rename(chr=seqnames) %>%
        mutate(transcribed_status=ifelse((ref == "C" | ref == "T"),
            ifelse(strand == "+", "coding", "template"),
            ifelse(strand == "+", "template", "coding"))) %>%
        mutate(alt=ifelse((ref == "C" | ref == "T"),
            alt,
            as.character(reverseComplement(DNAStringSet(alt))))) %>%
        mutate(ref=ifelse((ref == "C" | ref == "T"),
            ref,
            as.character(reverseComplement(DNAStringSet(ref))))) %>%
        mutate(aggregate=paste(id, type, chr, left, right, sep=",")) %>%
        count(aggregate, ref, alt, transcribed_status, src), by=c("aggregate", "ref", "alt", "transcribed_status", "src")) %>%
    mutate(n=ifelse(is.na(n), 0, n)) %>%
    pivot_wider(names_from=src, values_from=c("n")) %>%
    separate(aggregate, sep=",", into=c("id", "type", "chr", "left", "right")) %>%
    mutate(across(left:right, as.numeric))


write_tsv(augmented_mutations, out_file)
