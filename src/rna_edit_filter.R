library(tidyverse)
library(readxl)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)

args <- commandArgs(trailingOnly=TRUE)
fasta_file <- args[1]
genome_gff_file <- args[2]
rna_edits_file <- args[3]
out_file <- args[4]

fa <- readDNAStringSet(filepath=fasta_file, format="fasta")
names(fa) <- paste0("Chr", c(1:5, "M", "C"))

genome_gff <- read_tsv(genome_gff_file, col_names=F) %>%
    dplyr::select(chr=X1, type=X3, left=X4, right=X5, strand=X7, info=X9) %>%
    filter(type == "mRNA")

rna_edits <- read_excel(rna_edits_file)

chroms <- genome_gff$chr
names(chroms) <- str_remove(str_split_fixed(genome_gff$info, ";", 2)[,1], "ID=")

posi <- ifelse(genome_gff$strand == "+", genome_gff$left, genome_gff$right)
names(posi) <- str_remove(str_split_fixed(genome_gff$info, ";", 2)[,1], "ID=")

strand <- genome_gff$strand
names(strand) <- str_remove(str_split_fixed(genome_gff$info, ";", 2)[,1], "ID=")

rna_edits <- rna_edits %>%
    mutate(strand=strand[`Gene model ID`]) %>%
    mutate(chrom=chroms[`Gene model ID`]) %>%
    mutate(pos=ifelse(strand == "+",
        (posi[`Gene model ID`]+`Position in gene model`-1),
        (posi[`Gene model ID`]-`Position in gene model`+1))) %>%
    mutate(found=as.character(getSeq(fa, GRanges(seqnames=chrom, IRanges(start=pos, end=pos)))))

rna_edits_bed <- tibble(chrom=str_remove(rna_edits$chrom, "Chr"), chromStart=rna_edits$pos, chromEnd=(rna_edits$pos+1), name=".", score=".", strand=".") %>%
    arrange(chrom, chromStart)

write_tsv(rna_edits_bed, out_file, col_names=F)
