# Title     : Transform counts to percents
# Objective : Transform counts to percents
# Created by: valentinscherz
# Created on: 28.05.2019

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



## Import the phyloseq phyloseq_object
phyloseq_obj <- readRDS(phyloseq_object)
phyloseq_obj <- prune_taxa(taxa_sums(phyloseq_obj) > 0, phyloseq_obj) ## Removes taxa not at least present in one sample

### Transform counts into % of samples
trans <- transform_sample_counts(phyloseq_obj, function(x) 100 * x / sum(x))

### Save this new object
saveRDS(object = trans, file = phyloseq_pct)
