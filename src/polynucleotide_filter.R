library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
fasta_file <- args[1]
out_file <- args[2]


NEIGHBOURHOOD <- 20
CHROMOSOMES <- 1:5


fa <- readDNAStringSet(filepath=fasta_file, format="fasta")
names(fa) <- paste0("Chr", names(fa))


polynucleotide <- lapply(CHROMOSOMES, function(x){
    aaa <- gregexpr("AAAAAAA+|CCCCCCC+|GGGGGGG+|TTTTTTT+", fa[[x]])
    left <- aaa[[1]] - NEIGHBOURHOOD
    right <- aaa[[1]] + attr(aaa[[1]], "match.length") - 1 + NEIGHBOURHOOD
    tibble(chrom=x, chromStart=left, chromEnd=(right+1), name=".", score=".", strand=".")
})

polynucleotide_bed <- map_dfr(polynucleotide, bind_rows)


write_tsv(polynucleotide_bed, out_file, col_names=F)