# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 03.01.19
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

# Redirect R output
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

# Input
    q_score_filtered_F <- snakemake@input[["q_score_filtered_F"]]
    q_score_filtered_R <- snakemake@input[["q_score_filtered_R"]]

# Output
    error_profile_F <- snakemake@output[["error_profile_F"]]
    error_profile_R <- snakemake@output[["error_profile_R"]]

# Load needed libraries
    library(dada2); packageVersion("dada2")

# set.seed
    set.seed(100)

# Learn error rates
    errF <- learnErrors(q_score_filtered_F, nbases=1e8, multithread=TRUE, verbose = 1)
    errR <- learnErrors(q_score_filtered_R, nbases=1e8, multithread=TRUE, verbose = 1)

# Write these error profiles
    saveRDS(object = errF, file = error_profile_F)
    saveRDS(object = errR, file = error_profile_R)
