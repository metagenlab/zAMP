# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 28.12.18
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

# Redirect R output
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

# Input
    q_score_filtered_F <- snakemake@input[["q_score_filtered_F"]]
    q_score_filtered_R <- snakemake@input[["q_score_filtered_R"]]
    error_profile_F <- snakemake@input[["error_profile_F"]]
    error_profile_R <- snakemake@input[["error_profile_R"]]

# Output
    infer_stats <- snakemake@output[["infer_stats"]]
    sample_seq_tab <- snakemake@output[["sample_seq_tab"]]

# Parameters
    sam <- snakemake@params[["sample_name"]]
    run <- snakemake@params[["run"]]
    x_axis_filter_column_value <- snakemake@params[["x_axis_filter_column_value"]]

# Load needed libraries
    library(dada2); packageVersion("dada2")

# set.seed
    set.seed(100)

# Create a useful function
    getN <- function(x) sum(getUniques(x))

# File renaming
    names(q_score_filtered_F) <- sam
    names(q_score_filtered_R) <- sam

# Read error rates
    errF <- readRDS(error_profile_F)
    errR <- readRDS(error_profile_R)

# Prepare named vector
    mergers <- vector("list", 1)
    names(mergers) <- sam

# Sample inference and merger of paired-end reads
    cat("Processing:", sam, "\n")
    ## Forward
        derepF <- derepFastq(q_score_filtered_F, verbose = TRUE)
        ddF <- dada(derepF, err=errF, multithread=2, verbose = 1, pool = FALSE, selfConsist = TRUE)
    ## Reverse
        derepR <- derepFastq(q_score_filtered_R, verbose = TRUE)
        ddR <- dada(derepR, err=errR, multithread=2, verbose = 1, pool = FALSE, selfConsist = TRUE)
    ## Merge
        merger <- mergePairs(dadaF = ddF, derepF= derepF, dadaR = ddR, derepR = derepR, verbose = TRUE, minOverlap = 20, maxMismatch = 0)
        mergers[[sam]] <- merger

# Save the dereplicated, corrected and merged sequences for this sample
    saveRDS(mergers, file = sample_seq_tab)

# For statistics record the number of reads
    ## Write the statistics in a dataframe
        infer <- data.frame(denoisedF = getN(ddF))
        infer$denoisedR <- getN(ddR)
        infer$merged_pairs <- getN(merger)
        infer$Sample <- sam
        infer$label <- x_axis_filter_column_value
        infer$RUN <- run

    ## Save the sequences stats for this sample
        saveRDS(infer, file = infer_stats)

