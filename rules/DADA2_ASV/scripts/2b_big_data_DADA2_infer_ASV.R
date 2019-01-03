# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 28.12.18
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
q_score_filtered_Fs <- snakemake@input[["q_score_filtered_Fs"]]
q_score_filtered_Rs <- snakemake@input[["q_score_filtered_Rs"]]
error_profile_F <- snakemake@input[["error_profile_F"]]
error_profile_R <- snakemake@input[["error_profile_R"]]

## Output
forward_correct_seq <- snakemake@output[["forward_correct_seq"]]
reverse_correct_seq <- snakemake@output[["reverse_correct_seq"]]
merged_seq_tab <- snakemake@output[["merged_seq_tab"]]

## Parameters
sam <- snakemake@params[["sample_name"]]

## Load needed libraries
library(dada2); packageVersion("dada2")

set.seed(100)

# File parsing
names(q_score_filtered_Fs) <- sam
names(q_score_filtered_Rs) <- sam


# Read error rates
errF <- readRDS(error_profile_F)
errR <- readRDS(error_profile_R)

# Sample inference and merger of paired-end reads
mergers <- vector("list", 1)
names(mergers) <- sam
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(q_score_filtered_Fs)
    ddF <- dada(derepF, err=errF, multithread=TRUE, verbose = 1, pool = FALSE, selfConsist = TRUE)
    derepR <- derepFastq(q_score_filtered_Rs)
    ddR <- dada(derepR, err=errR, multithread=TRUE, verbose = 1, pool = FALSE, selfConsist = TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger


# Save the corrected forward sequences for this sample
saveRDS(ddF, file = forward_correct_seq)

# Save the corrected reverse sequences for this sample
saveRDS(ddR, file = reverse_correct_seq)

# Save the dereplicated, corrected and merged sequences for this sample
saveRDS(mergers, file = merged_seq_tab)
