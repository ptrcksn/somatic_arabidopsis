set.seed(1)

library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)
library(glmnet)


args <- commandArgs(trailingOnly=TRUE)
features_file <- args[1]
study <- args[2]
ref_base <- args[3]
alt_base <- args[4]
out_file <- args[5]


fit_model <- function(m, b, X_names, fold_ids){
    linear_model <- cv.glmnet(x=X_names,
        y=m,
        offset=b,
        family=poisson(),
	foldid=fold_ids)

    beta <- coef(linear_model, s="lambda.min")

    beta
}


features <- read_tsv(features_file)

features_sub <- features %>%
    filter(ref == ref_base & alt == alt_base) %>%
    group_by(id, transcribed_status) %>%
    summarise(across(c(width, base_width, epigenomes:H4K16ac_signal), sum)) %>% # sum covariates by gene
    ungroup() %>%
    mutate(across(GC:H4K16ac_signal, ~(.x/width))) %>% # per base genic properties not dependent on expression
    rowwise() %>%
    mutate(across(c(base_width, epigenomes:H4K16ac_signal), ~ifelse(is.infinite(.x), NA, .x))) %>%
    ungroup() %>%
    drop_na() %>%
    filter(breadth != 0 & depth != 0) %>%
    mutate(norm_depth=depth/breadth) %>%
    mutate(across(timing_signal, ~log2(.x))) %>%
    rowwise() %>%
    mutate(across(timing_signal, ~ifelse(is.infinite(.x), NA, .x))) %>%
    ungroup() %>%
    drop_na() %>%
    mutate(transcribed_status=as.numeric(scale(ifelse(transcribed_status == "template", 1, 0)))) %>%
    mutate(across(GC:H4K16ac_signal, ~as.numeric(scale(.x)))) %>%
    mutate(log_norm_depth=log(norm_depth))

n_folds <- 10
fold_ids <- sample((1:nrow(features_sub) %% n_folds), size=nrow(features_sub), replace=F)
fold_ids <- ifelse(fold_ids == 0, n_folds, fold_ids)

m <- if(study == "epigenomes"){ features_sub$epigenomes }else{ features_sub$EU_som }

b <- if(study == "epigenomes"){ log(features_sub$breadth) }else{ log(features_sub$base_width) }

X_names <- if(study == "epigenomes"){
            as.matrix(tibble(
            "log_normalised_depth"=features_sub$log_norm_depth,
            "transcribed_status"=features_sub$transcribed_status,
            "GC_content"=features_sub$GC,
            "replication_timing"=features_sub$timing_signal,
            "DNA_accessibility"=features_sub$ATAC_signal,
            "DNA_methylation"=features_sub$meDIP_signal,
            "H3K14ac"=features_sub$H3K14ac_signal,
            "H3K23ac"=features_sub$H3K23ac_signal,
            "H3K27ac"=features_sub$H3K27ac_signal,
            "H3K27me1"=features_sub$H3K27me1_signal,
            "H3K27me3"=features_sub$H3K27me3_signal,
            "H3K36ac"=features_sub$H3K36ac_signal,
            "H3K36me3"=features_sub$H3K36me3_signal,
            "H3K4me1"=features_sub$H3K4me1_signal,
            "H3K4me2"=features_sub$H3K4me2_signal,
            "H3K4me3"=features_sub$H3K4me3_signal,
            "H3K56ac"=features_sub$H3K56ac_signal,
            "H3K9ac"=features_sub$H3K9ac_signal,
            "H3K9me1"=features_sub$H3K9me1_signal,
            "H3K9me2"=features_sub$H3K9me2_signal,
            "H4K16ac"=features_sub$H4K16ac_signal
        ))
    } else{
        as.matrix(tibble(
            "transcribed_status"=features_sub$transcribed_status,
            "GC_content"=features_sub$GC,
            "replication_timing"=features_sub$timing_signal,
            "DNA_accessibility"=features_sub$ATAC_signal,
            "DNA_methylation"=features_sub$meDIP_signal,
            "H3K14ac"=features_sub$H3K14ac_signal,
            "H3K23ac"=features_sub$H3K23ac_signal,
            "H3K27ac"=features_sub$H3K27ac_signal,
            "H3K27me1"=features_sub$H3K27me1_signal,
            "H3K27me3"=features_sub$H3K27me3_signal,
            "H3K36ac"=features_sub$H3K36ac_signal,
            "H3K36me3"=features_sub$H3K36me3_signal,
            "H3K4me1"=features_sub$H3K4me1_signal,
            "H3K4me2"=features_sub$H3K4me2_signal,
            "H3K4me3"=features_sub$H3K4me3_signal,
            "H3K56ac"=features_sub$H3K56ac_signal,
            "H3K9ac"=features_sub$H3K9ac_signal,
            "H3K9me1"=features_sub$H3K9me1_signal,
            "H3K9me2"=features_sub$H3K9me2_signal,
            "H4K16ac"=features_sub$H4K16ac_signal
        ))
    }

fit <- fit_model(m=m, b=b, X_names=X_names, fold_ids=fold_ids)

fit_result <- tibble(study=study, ref=ref_base, alt=alt_base, covariate=rownames(fit), estimate=as.numeric(fit))

write_tsv(fit_result, out_file, col_names=F)
