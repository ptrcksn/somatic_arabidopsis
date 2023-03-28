library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)

args <- commandArgs(trailingOnly=TRUE)
accession_file <- args[1]
multi_instance_file <- args[2]
bedgraph_file <- args[3]
out_file <- args[4]

accession_order <- as.character(scan(accession_file))

multi_instance <- read_tsv(multi_instance_file, col_names=F) %>%
    dplyr::rename(chr=X1, start=X2) %>%
    mutate(chr=as.character(chr))

f <- function(x, pos){

    depth <- x %>%
        mutate(chr=as.character(X1)) %>%
        dplyr::rename(start=X2) %>%
        mutate(end=start) %>%
        dplyr::select(-X1,-X3) %>%
        anti_join(multi_instance, by=c("chr", "start"))

    colnames(depth)[2:(1+length(accession_order))] <- accession_order

    featurewise_breadth <- depth %>%
        summarise(across(!c(chr,start,end), ~length(which(.x > 0)))) %>%
        pivot_longer(cols=everything())

    featurewise_depth <- depth %>%
        summarise(across(!c(chr, start, end), sum)) %>%
        pivot_longer(cols=everything())

    featurewise_breadth_depth <- featurewise_breadth %>%
        full_join(featurewise_depth, by=c("name")) %>%
        dplyr::rename(accession=name, breadth=value.x, depth=value.y)

}

sample_breadth_depth <- read_tsv_chunked(bedgraph_file, DataFrameCallback$new(f), chunk_size=100000, col_names=F) %>%
    group_by(accession) %>%
    summarise(across(c(breadth, depth), sum))

write_tsv(sample_breadth_depth, out_file, col_names=F)