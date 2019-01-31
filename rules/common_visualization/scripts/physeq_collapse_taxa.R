# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 26.10.18

## Inspired from https://gist.github.com/erictleung/eda33bceceebb190eb9a44c93a077d32
## Redirect R output
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[[1]]

## Output
phyloseq_filtered_object <- snakemake@output[[1]]

## Parameters
collapse_level <- snakemake@params[["collapse_level"]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")


## Load the phyloseq object
phyloseq_object <- readRDS(phyloseq_object)

rank_names(phyloseq_object)[[as.numeric(collapse_level)]]


## Collapse taxa
collapsed_physeq <- tax_glom(phyloseq_object, taxrank=rank_names(phyloseq_object)[[as.numeric(collapse_level)]], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

# Write the new phyloseq object
saveRDS(object = collapsed_physeq, file = phyloseq_filtered_object)
