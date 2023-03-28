library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)


args <- commandArgs(trailingOnly=TRUE)
gff_file <- args[1]
fasta_file <- args[2]
accession_file <- args[3]
outlier_file <- args[4]
multi_instance_file <- args[5]
bedgraph_file <- args[6]
out_file <- args[7]


gff <- read_tsv(gff_file, col_names=F) %>%
    dplyr::select(chr=X1, type=X3, left=X4, right=X5, strand=X7, info=X9) %>%
    mutate(id=info) %>%
    filter(type == "exon")

gff_r <- GRanges(seqnames=gff$chr, IRanges(start=gff$left, end=gff$right), aggregate_locus=paste(gff$id, gff$type, gff$chr, gff$left, gff$right, sep=","))

fa <- readDNAStringSet(filepath=fasta_file, format="fasta")
names(fa) <- paste0("Chr", c(1:5, "M", "C"))

accession_order <- as.character(scan(accession_file))

outlier <- as.character(scan(outlier_file))

multi_instance <- read_tsv(multi_instance_file, col_names=F) %>%
    dplyr::rename(chr=X1, start=X2) %>%
    mutate(chr=paste0("Chr", chr))

f <- function(x, pos){

    depth <- x %>%
        mutate(chr=paste0("Chr", X1)) %>%
        dplyr::rename(start=X2) %>%
        mutate(end=start) %>%
        dplyr::select(-X1,-X3) %>%
        anti_join(multi_instance, by=c("chr", "start"))

    colnames(depth)[2:(1+length(accession_order))] <- accession_order

    depth <- depth[, !(colnames(depth) %in% outlier)]

    depth_r <- GRanges(seqnames=depth$chr, IRanges(start=depth$start, end=depth$end))
    depth_gene_r <- join_overlap_intersect(gff_r, depth_r)

    ref <- as.character(getSeq(fa, depth_gene_r))

    depth_gene <- as_tibble(depth_gene_r) %>%
        dplyr::rename(chr=seqnames) %>%
        mutate(ref=factor(ref, levels=c("A","C","G","T"))) %>%
        inner_join(depth, by=c("chr", "start", "end"))

    featurewise_breadth <- depth_gene %>%
        count(aggregate_locus, ref, .drop=F) %>%
        dplyr::rename(breadth=n)

    featurewise_depth <- depth_gene %>%
        dplyr::select(!chr:strand) %>%
        mutate(depth=rowSums(across(where(is.double)), na.rm=T)) %>%
        dplyr::select(aggregate_locus, ref, depth) %>%
        group_by(aggregate_locus, ref) %>%
        summarise(depth=sum(depth, na.rm=T)) %>%
        ungroup()

    featurewise_breadth_depth <- featurewise_breadth %>%
        left_join(featurewise_depth, by=c("aggregate_locus", "ref"))

}

chunked_featurewise_breadth_depth <- read_tsv_chunked(bedgraph_file, DataFrameCallback$new(f), chunk_size=100000, col_names=F)

featurewise_breadth_depth <- chunked_featurewise_breadth_depth %>%
    mutate(ref=as.character(ref)) %>%
    group_by(aggregate_locus, ref) %>%
    summarise(breadth=sum(breadth, na.rm=T), depth=sum(depth, na.rm=T)) %>%
    ungroup() %>%
    mutate(breadth=ifelse(depth == 0, 0, breadth)) %>%
    separate(aggregate_locus, sep=",", into=c("id", "type", "chr", "left", "right")) %>%
    mutate(across(left:right, as.integer)) %>%
    filter(!is.na(ref)) %>%
    inner_join(gff %>%
        dplyr::select(id, type, chr, left, right, strand), by=c("id", "type", "chr", "left", "right")) %>% #
    mutate(transcribed_status=ifelse((ref == "C" | ref == "T"),
        ifelse(strand == "+", "coding", "template"),
        ifelse(strand == "+", "template", "coding"))) %>%
    mutate(ref=ifelse((ref == "C" | ref == "T"), ref, as.character(reverseComplement(DNAStringSet(ref)))))


write_tsv(featurewise_breadth_depth, out_file, col_names=F)