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
    forward_correct_seq <- snakemake@input[["forward_correct_seq"]]
    reverse_correct_seq <- snakemake@input[["reverse_correct_seq"]]
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
    dadaFs <- do.call(rbind, lapply(forward_correct_seq, readRDS))
    dadaFs <- getN(dadaFs)

    filtration <- cbind(filtration, dadaFs)

write.table(x = filtration, file = filtering_stats, sep = "\t", col.names = NA, row.names = TRUE)
