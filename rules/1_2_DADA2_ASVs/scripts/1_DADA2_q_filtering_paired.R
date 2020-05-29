# Title     : DADA2 - Q-score filtering
# Objective : Filter reads based on q-score and length
# Created by: valentinscherz
# Created on: 06.06.19
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

## Redirect R output to the log file
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    fnFs <- snakemake@input[[1]]
    fnRs <- snakemake@input[[2]]

## Output
    q_score_filtered_F <- snakemake@output[["q_score_filtered_F"]]
    q_score_filtered_R <- snakemake@output[["q_score_filtered_R"]]
    q_filtering_stats_path <- snakemake@output[["q_filtering_stats"]]

## Parameters
    F_length <- snakemake@params[["F_reads_length_trim"]]
    R_length <- snakemake@params[["R_reads_length_trim"]]
    F_extected_error <- snakemake@params[["F_reads_extected_error"]]
    R_extected_error <- snakemake@params[["R_reads_extected_error"]]
    sample_name <- snakemake@params[["sample_name"]]


## Load needed libraries
    library(dada2);packageVersion("dada2")

## Filter and trim
### Reads are filtered based on the number of errors expected for the read (integration of the qscore and the length of the read). All reads with uncalled nucleotide (N) are removed too. Remaining phiX reads will be removed too. Finally, reads are cut at a length set in config.

### Filter and trim. The filtered reads are directly written while the filtering stats are being save for later compilation.
    filtering_stats <- filterAndTrim(fwd = fnFs, filt = q_score_filtered_F, rev = fnRs, filt.rev = q_score_filtered_R, truncLen=c(F_length,R_length), maxN=0, maxEE=c(F_extected_error,R_extected_error), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=snakemake@threads, verbose = TRUE)
    filtering_stats <- as.data.frame(filtering_stats)
    filtering_stats$Sample <- sample_name

### Save the stats for this sample in a R object
    saveRDS(filtering_stats, file = q_filtering_stats_path)
