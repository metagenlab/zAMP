# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 03.01.19
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
sample_seq_table <- snakemake@input[["sample_seq_table"]]
forward_stats <- snakemake@input[["forward_stats"]]
reverse_stats <- snakemake@input[["reverse_stats"]]
merged_stats <- snakemake@input[["merged_stats"]]

## Output
run_seq_table <- snakemake@output[["run_seq_table"]]
run_forward_stats <- snakemake@output[["run_forward_stats"]]
run_reverse_stats <- snakemake@output[["run_reverse_stats"]]
run_merged_stats <- snakemake@output[["run_merged_stats"]]

## Load needed libraries
library(dada2); packageVersion("dada2")


## Merge paired-end merge reads of the run
input_M <- sapply(sample_seq_table, readRDS, simplify = TRUE, USE.NAMES = FALSE)
st.all <- makeSequenceTable(input_M)
saveRDS(object = st.all, file = run_seq_table)


# For later use in statistics, get the number of reads for each sample at the forwards and reverse stage
## Merge forward reads stats of the run
forward_correct
input_F <- do.call("rbind", lapply(forward_stats, readRDS))
input_F
saveRDS(object = input_F, file = run_forward_stats)

## Merge reverse reads stats of the run
reverse_correct
input_R <- do.call("rbind", lapply(reverse_stats, readRDS))
input_R
saveRDS(object = input_R, file = run_reverse_stats)

## Merge merged reads stats of the run
merged_stats
input_M <- do.call("rbind", lapply(merged_stats, readRDS))
input_M
saveRDS(object = input_M, file = run_merged_stats)
