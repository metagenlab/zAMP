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
phyloseq_pct <- snakemake@output[[1]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")

## Load the phyloseq phyloseq_object
phyloseq_obj <- readRDS(phyloseq_object)

### Transform counts into % of samples
trans <- transform_sample_counts(phyloseq_obj, function(x) x / sum(x))

### Save this new object
saveRDS(object = trans, file = phyloseq_pct)
