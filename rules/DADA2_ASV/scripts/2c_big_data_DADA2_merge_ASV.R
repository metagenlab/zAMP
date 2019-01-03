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

## Output
#filtering_stats <- snakemake@output[["filtering_stats"]]
merged_seq_table <- snakemake@output[["merged_seq_table"]]

## Load needed libraries
library(dada2); packageVersion("dada2")

## Merge multiple runs (if necessary)
input <- sapply(sample_seq_table, readRDS, simplify = TRUE, USE.NAMES = FALSE)
st.all <- makeSequenceTable(input)
saveRDS(object = st.all, file = merged_seq_table)

