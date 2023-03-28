library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
accession_file <- args[1]
mutations_map_file <- args[2]
sample_depth_file <- args[3]
out_file <- args[4]

#palette <- c("#F4F1DE", "#E07A5F", "#3D405B", "#468C98", "#F2CC8F", "darkgrey")

accessions_id <- as.character(scan(accession_file))
mutations_map <- read_tsv(mutations_map_file)
sample_depth <- read_tsv(sample_depth_file)

mutations_map_depth <- mutations_map %>%
    mutate(accession=factor(accession, accessions_id)) %>%
    group_by(accession) %>%
    count(accession, .drop=F) %>%
    ungroup() %>%
    full_join(sample_depth %>%
        group_by(accession) %>%
        summarise(depth=sum(depth), breadth=sum(breadth)) %>%
        mutate(accession=factor(accession)), by="accession")

mutations_map_depth <- mutations_map_depth %>%
    mutate(box=cut(depth, breaks=4)) %>%
    group_by(box) %>%
    nest() %>%
    mutate(lower=map(data, function(df) quantile(df$n)[2]-(quantile(df$n)[4]-quantile(df$n)[2]))) %>%
    mutate(upper=map(data, function(df) quantile(df$n)[4]+(quantile(df$n)[4]-quantile(df$n)[2]))) %>%
    mutate(size=map(data, function(x) nrow(x))) %>%
    ungroup() %>%
    unnest(lower) %>%
    unnest(upper) %>%
    unnest(size) %>%
    unnest(data) %>%
#    mutate(outlier=ifelse(n < lower | n > upper, 1, 0))
    mutate(outlier=ifelse((n < lower | n > upper) | size <= 1, 1, 0))

outliers <- mutations_map_depth %>%
    filter(outlier == 1) %>%
    dplyr::select(accession)

write_tsv(outliers, out_file, col_names=F)
