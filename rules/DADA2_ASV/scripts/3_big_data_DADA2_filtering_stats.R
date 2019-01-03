# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 03.01.19
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

## Redirect R output
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

## Input
    q_score <- snakemake@input[["q_score"]]
    run_merged_F_correct <- snakemake@input[["run_merged_F_correct"]]
    run_merged_R_correct <- snakemake@input[["run_merged_R_correct"]]
    merged <- snakemake@input[["merged"]]
    with_chim <- snakemake@input[["with_chim"]]
    no_chim <- snakemake@input[["no_chim"]]
    length_filtered <- snakemake@input[["length_filtered"]]

## Output
    filtering_stats <- snakemake@output[["filtering_stats"]]

## Load needed libraries
    library(dada2); packageVersion("dada2")

## Create a useful function
    getN <- function(x) sum(getUniques(x))

## Load the q score filtration R stats
    filtration <- do.call(rbind, lapply(q_score, readRDS))

## Load the forward and reverse reads correction stats
    dadaFs <- do.call("rbind", sapply(run_merged_F_correct, readRDS, simplify = TRUE, USE.NAMES = FALSE))
    dadaRs <- do.call("rbind", sapply(run_merged_R_correct, readRDS, simplify = TRUE, USE.NAMES = FALSE))

## Load the merge reads correction stats
    dadamerged <- do.call("rbind", sapply(merged, readRDS, simplify = TRUE, USE.NAMES = FALSE))
    dadamerged

    filtration <- cbind(filtration, dadaFs, dadaRs, dadamerged)

write.table(x = filtration, file = filtering_stats, sep = "\t", col.names = NA, row.names = TRUE)

