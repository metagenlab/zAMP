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
tax_rank <- snakemake@params[["filter_tax_rank"]]
lineage <- snakemake@params[["filter_lineage"]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")
library(dplyr);packageVersion("dplyr")

## Load the phyloseq object
phyloseq_object <- readRDS(phyloseq_object)
phyloseq_object

## filter taxa
filtered_taxa <- subset_taxa(phyloseq_object, get(tax_rank) == as.character(lineage))

## Remove already in metadata alphia diversity values
sample_data(filtered_taxa) <- select(sample_data(filtered_taxa), -c(Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher))

## Add alpha diversity indexes to metadata
alpha_div <- estimate_richness(physeq = filtered_taxa, split = TRUE)
sample_data(filtered_taxa) <- cbind(sample_data(filtered_taxa),alpha_div)

## Write the new phyloseq object
saveRDS(object = filtered_taxa, file = phyloseq_filtered_object)
