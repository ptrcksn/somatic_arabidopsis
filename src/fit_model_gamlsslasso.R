library(tidyverse)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(plyranges)
library(gamlss)
library(gamlss.lasso)


args <- commandArgs(trailingOnly=TRUE)
features_file <- args[1]
study <- args[2]
ref_base <- args[3]
alt_base <- args[4]
out_file <- args[5]


fit_model <- function(m, b, X_names, penalty_factor){
    linear_model <- gamlss(m~
            offset(b)+
            gnet(x.vars=X_names, method="IC", ICpen="AIC"),
        family=SICHEL(),
        data=features_sub,
        control=gamlss.control(n.cyc=200, c.crit=0.1))

    beta <- tail(getSmo(linear_model, "mu") ,1)[[1]]$beta[(1:length(penalty_factor))[which(penalty_factor == 1)]]

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

m <- if(study == "epigenomes"){ features_sub$epigenomes }else{ features_sub$EU_som }

b <- if(study == "epigenomes"){ log(features_sub$breadth) }else{ log(features_sub$base_width) }

X_names <- if(study == "epigenomes"){ c("log_norm_depth") }else{ c() }

X_names <- c(X_names,
        "transcribed_status",
        "GC",
        "timing_signal",
        "ATAC_signal",
        "meDIP_signal",
        "H3K14ac_signal",
        "H3K23ac_signal",
        "H3K27ac_signal",
        "H3K27me1_signal",
        "H3K27me3_signal",
        "H3K36ac_signal",
        "H3K36me3_signal",
        "H3K4me1_signal",
        "H3K4me2_signal",
        "H3K4me3_signal",
        "H3K56ac_signal",
        "H3K9ac_signal",
        "H3K9me1_signal",
        "H3K9me2_signal",
        "H4K16ac_signal")

penalty_factor <- if(study == "epigenomes"){ rep(c(0, 1), c(1, (length(X_names)-1))) }else{ rep(1, length(X_names)) }

fit <- fit_model(m=m, b=b, X_names=X_names, penalty_factor=penalty_factor)

fit_result <- tibble(study=study, ref=ref_base, alt=alt_base, covariate=names(fit), estimate=as.numeric(fit))

write_tsv(fit_result, out_file, col_names=F)
