# Title     : DADA2 - Learn run errors profile
# Objective : Generate an error profile for each run
# Created by: valentinscherz
# Created on: 06.06.19
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

## Redirect R output
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    q_score_filtered_F <- snakemake@input[["q_score_filtered_F"]]
    q_score_filtered_R <- snakemake@input[["q_score_filtered_R"]]

## Output
    error_profile_F <- snakemake@output[["error_profile_F"]]
    error_profile_R <- snakemake@output[["error_profile_R"]]
    error_plot_F <- snakemake@output[["error_profile_F_plot"]]
    error_plot_R <- snakemake@output[["error_profile_R_plot"]]

## Load needed libraries
    library(dada2); packageVersion("dada2")

## keep only non empty fastq files
    q_score_filtered_F_non_empty = q_score_filtered_F[file.size(q_score_filtered_F) != 20]
    q_score_filtered_R_non_empty = q_score_filtered_R[file.size(q_score_filtered_R) != 20]
## Learn error rates
    errF <- learnErrors(q_score_filtered_F_non_empty, nbases=1e8, multithread = snakemake@threads, verbose = 1)
    errR <- learnErrors(q_score_filtered_R_non_empty, nbases=1e8, multithread = snakemake@threads, verbose = 1)

## Write these error profiles
    saveRDS(object = errF, file = error_profile_F)
    saveRDS(object = errR, file = error_profile_R)

## Write error_plots
    png(error_plot_F)
    plotErrors(errF, nominalQ=TRUE)
    png(error_plot_R)
    plotErrors(errR, nominalQ=TRUE)

    dev.off()
