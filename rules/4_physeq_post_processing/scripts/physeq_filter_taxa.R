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
phyloseq_object_path <- snakemake@input[["phyloseq_object"]]

## Output
phyloseq_filtered_object <- snakemake@output[["phyloseq_object"]]

## Parameters
tax_rank <- snakemake@params[["filter_tax_rank"]]
lineage <- snakemake@params[["filter_lineage"]]
filter_out_tax_rank  <- snakemake@params[["filter_out_tax_rank"]]
filter_out_lineage <- snakemake@params[["filter_out_lineage"]]


## Load needed libraries
library(phyloseq);packageVersion("phyloseq")
library(dplyr);packageVersion("dplyr")

## Load the phyloseq object
phyloseq_object <- readRDS(phyloseq_object_path)


## filter taxa
filtered_taxa <- subset_taxa(phyloseq_object, get(tax_rank) == as.character(lineage))
filtered_taxa <- subset_taxa(filtered_taxa, get(filter_out_tax_rank) != as.character(filter_out_lineage))


## Write the new phyloseq object
saveRDS(object = filtered_taxa, file = phyloseq_filtered_object)
