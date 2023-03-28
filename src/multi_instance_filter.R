library(tidyverse)


args <- commandArgs(trailingOnly=TRUE)
mutations_map_file <- args[1]
out_file <- args[2]


mutations <- read_tsv(mutations_map_file, col_names=F) %>%
    dplyr::rename(accession=X1, chr=X2, pos=X3, ref=X4, alt=X5, context=X6, depth=X7, alt_depth=X8, flag=X9)

multiple_instance_mutations <- mutations %>%
    filter(flag == "PASS" | flag == "clustered_mutation") %>%
    dplyr::count(chr, pos) %>%
    filter(n > 1) %>%
    mutate(chrom=chr) %>%
    mutate(chromStart=pos) %>%
    mutate(chromEnd=(pos+1)) %>%
    mutate(name=paste(chr, pos, sep=".")) %>%
    mutate(score=".") %>%
    mutate(strand=".") %>%
    dplyr::select(chrom, chromStart, chromEnd, name, score, strand) %>%
    arrange(chrom, chromStart)


write_tsv(multiple_instance_mutations, out_file, col_names=F)