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
tax_rank <- snakemake@params[["filter_tax_rank"]]
lineage <- snakemake@params[["filter_lineage"]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")



## Load the phyloseq object
phyloseq_object <- readRDS(phyloseq_object)
phyloseq_object


## filter taxa
phyltered_taxa = subset_taxa(phyloseq_object, get(tax_rank) == as.character(lineage))

# Write the new phyloseq object
saveRDS(object = phyltered_taxa, file = phyloseq_filtered_object)
