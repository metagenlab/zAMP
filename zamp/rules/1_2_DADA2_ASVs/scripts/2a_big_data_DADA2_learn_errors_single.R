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

## Output
    error_profile_F <- snakemake@output[["error_profile_F"]]
    error_plot_F <- snakemake@output[["error_profile_F_plot"]]

## Load needed libraries
    library(dada2); packageVersion("dada2")


## Learn error rates
    errF <- learnErrors(q_score_filtered_F, nbases=1e8, multithread = snakemake@threads, verbose = 1)

## Write these error profiles
    saveRDS(object = errF, file = error_profile_F)

## Write error_plots
    png(error_plot_F)
    plotErrors(errF, nominalQ=TRUE)

    dev.off()
