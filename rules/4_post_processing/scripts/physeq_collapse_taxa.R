# Title     : Collapse taxa
# Objective : Collapse the read for identical taxa at a given taxonomic level
# Created by: valentinscherz
# Created on: 28.05.2019


## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[[1]]

## Output
phyloseq_filtered_object <- snakemake@output[[1]]

## Parameters
collapse_level <- snakemake@params[["collapse_level"]]



## Load needed libraries
library(phyloseq);packageVersion("phyloseq")



## Import the phyloseq object
phyloseq_object <- readRDS(phyloseq_object)

## Convert the taxonomic level to its name
rank_names(phyloseq_object)[[as.numeric(collapse_level)]]

## Collapse taxa
collapsed_physeq <- tax_glom(phyloseq_object, taxrank=rank_names(phyloseq_object)[[as.numeric(collapse_level)]], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
collapsed_physeq <- prune_taxa(taxa_sums(collapsed_physeq) > 0, collapsed_physeq) ## Removes taxa not at least present in one sample


# Write the new phyloseq object
saveRDS(object = collapsed_physeq, file = phyloseq_filtered_object)
