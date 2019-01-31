# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 26.10.18

## Inspired from https://gist.github.com/erictleung/eda33bceceebb190eb9a44c93a077d32
## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[[1]]

## Output
phyloseq_filtered_object <- snakemake@output[[1]]

## Parameters
meta_column <- snakemake@params[["meta_column"]]
column_value <- snakemake@params[["column_value"]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")

## Load the phyloseq object
phyloseq_object <- readRDS(phyloseq_object)

## filter taxa
phyltered_taxa = subset_samples(phyloseq_object, get(meta_column) == as.character(column_value))


# Write the new phyloseq object
saveRDS(object = phyltered_taxa, file = phyloseq_filtered_object)
