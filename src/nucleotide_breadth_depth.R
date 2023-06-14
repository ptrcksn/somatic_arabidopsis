library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
gff_file <- args[1]
fasta_file <- args[2]
accession_file <- args[3]
bedgraph_file <- args[4]
out_file <- args[5]


gff <- read_tsv(gff_file, col_names=F) %>%
    dplyr::select(chr=X1, type=X3, left=X4, right=X5, strand=X7, info=X9) %>%
    mutate(id=info) %>%
    filter(type == "exon")

gff_r <- GRanges(seqnames=gff$chr, IRanges(start=gff$left, end=gff$right), aggregate_locus=paste(gff$id, gff$type, gff$chr, gff$left, gff$right, sep=","))

fa <- readDNAStringSet(filepath=fasta_file, format="fasta")
names(fa) <- paste0("Chr", c(1:5, "M", "C"))

accession_order <- as.character(scan(accession_file))

f <- function(x, pos){

    depth <- x %>%
        mutate(chr=paste0("Chr", X1)) %>%
        dplyr::rename(start=X2) %>%
        mutate(end=start) %>%
        dplyr::select(-X1,-X3)

    colnames(depth)[2:(1+length(accession_order))] <- accession_order

    depth_r <- GRanges(seqnames=depth$chr, IRanges(start=depth$start, end=depth$end))
    depth_gene_r <- join_overlap_intersect(gff_r, depth_r)

    ref <- as.character(getSeq(fa, depth_gene_r))

    depth_gene <- as_tibble(depth_gene_r) %>%
        dplyr::rename(chr=seqnames) %>%
        mutate(ref=factor(ref, levels=c("A","C","G","T"))) %>%
        inner_join(depth, by=c("chr", "start", "end"))

    featurewise_breadth <- depth_gene %>%
        dplyr::select(!width:aggregate_locus) %>%
        mutate(start=as.integer(start), end=as.integer(end)) %>%
        mutate(across(where(is.double), ~ifelse(.x > 0, 1, 0))) %>%
        mutate(breadth=rowSums(across(where(is.double)), na.rm=T)) %>%
        dplyr::select(chr, start, ref, breadth)

    featurewise_depth <- depth_gene %>%
        dplyr::select(!width:aggregate_locus) %>%
        mutate(start=as.integer(start), end=as.integer(end)) %>%
        mutate(depth=rowSums(across(where(is.double)), na.rm=T)) %>%
        dplyr::select(chr, start, ref, depth)

    featurewise_breadth_depth <- featurewise_breadth %>%
        left_join(featurewise_depth, by=c("chr", "start", "ref"))

}

chunked_featurewise_breadth_depth <- read_tsv_chunked(bedgraph_file, DataFrameCallback$new(f), chunk_size=100000, col_names=F)

write_tsv(chunked_featurewise_breadth_depth, out_file, col_names=F)
