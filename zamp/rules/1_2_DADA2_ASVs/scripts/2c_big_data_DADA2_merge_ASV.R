# Title     : DADA2 - Merge samples
# Objective : Merge samples from the same run
# Created by: valentinscherz
# Created on: 06.06.19
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

## Redirect R output
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    infer_stats <- snakemake@input[["infer_stats"]]
    sample_seq_table <- snakemake@input[["sample_seq_tab"]]

## Output
    run_stats <- snakemake@output[["run_stats"]]
    run_seq_table <- snakemake@output[["run_seq_table"]]

## Load needed libraries
    library(dada2); packageVersion("dada2")

## Merge and write paired-end merge reads of the run
    input_M <- sapply(sample_seq_table, readRDS, simplify = TRUE, USE.NAMES = FALSE)
    st.all <- makeSequenceTable(input_M)
    saveRDS(object = st.all, file = run_seq_table)

## Merge and write forward reads stats of the run
    stats <- do.call("rbind", lapply(infer_stats, readRDS))
    saveRDS(object = stats, file = run_stats)

