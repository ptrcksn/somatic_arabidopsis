library(tidyverse)


args <- commandArgs(trailingOnly=TRUE)
mutations_map_file <- args[1]
query_accession <- args[2]
query_flags <- args[3]
out_file <- args[4]


mutations <- read_tsv(mutations_map_file, col_names=F) %>%
    dplyr::rename(accession=X1, chr=X2, pos=X3, ref=X4, alt=X5, context=X6, depth=X7, alt_depth=X8, flag=X9)

accession_mutations <- mutations %>%
    filter(accession == query_accession)

query_flags <- str_split(query_flags, ";")[[1]]

flagged_accession_mutations <- lapply(query_flags, function(x){accession_mutations %>%
    filter(str_detect(flag, x)) %>%
    mutate(chrom=chr) %>%
    mutate(chromStart=(pos-1)) %>%
    mutate(chromEnd=pos) %>%
    mutate(name=paste(chr, pos, sep=".")) %>%
    mutate(score=".") %>%
    mutate(strand=".") %>%
    dplyr::select(chrom, chromStart, chromEnd, name, score, strand)
})

merged_flagged_accession_mutations <- map_dfr(flagged_accession_mutations, bind_rows) %>%
    arrange(chrom, chromStart)

write_tsv(merged_flagged_accession_mutations, out_file, col_names=F)
