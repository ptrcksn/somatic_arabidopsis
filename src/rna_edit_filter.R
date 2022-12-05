library(tidyverse)
library(readxl)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
genome_gff_file <- args[1]
rna_edits_file <- args[2]
out_file <- args[3]


genome_gff <- read_tsv(genome_gff_file, col_names=F) %>%
    dplyr::select(chr=X1, type=X3, left=X4, right=X5, strand=X7, info=X9) %>%
    filter(type == "mRNA")

rna_edits <- read_excel(rna_edits_file)

chroms <- genome_gff$chr
names(chroms) <- str_remove(str_split_fixed(genome_gff$info, ";", 2)[,1], "ID=")

pos <- genome_gff$left
names(pos) <- str_remove(str_split_fixed(genome_gff$info, ";", 2)[,1], "ID=")

rna_edits <- rna_edits %>%
    mutate(chrom=chroms[rna_edits$`Gene model ID`]) %>%
    mutate(pos=(pos[rna_edits$`Gene model ID`]+rna_edits$`Position in gene model`-1))

rna_edits_bed <- tibble(chrom=str_remove(rna_edits$chrom, "Chr"), chromStart=(rna_edits$pos-1), chromEnd=rna_edits$pos, name=".", score=".", strand=".") %>%
    arrange(chrom, chromStart)


write_tsv(rna_edits_bed, out_file, col_names=F)
