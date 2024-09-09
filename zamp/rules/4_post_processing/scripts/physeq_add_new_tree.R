# Title     : Physeq new tree
# Objective : Add new tree to phyloseq_object
# Created by: valentinscherz
# Created on: 28.05.2019


## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]
new_tree <- snakemake@input[["new_tree"]]

## Ouput
updated_phyloseq_path <- snakemake@output[["phyloseq_object"]]



## Load libraries
library(phyloseq);packageVersion("phyloseq")



## Read the phyloseq object
phyloseq_obj <- readRDS(phyloseq_object)
phyloseq_obj <- prune_taxa(taxa_sums(phyloseq_obj) > 0, phyloseq_obj) ## Removes taxa not at least present in one sample

## Read the new tree
PHY <- read_tree(new_tree)

## Assign the new tree to the phyloseq object
phy_tree(phyloseq_obj) <- PHY

## Write the new phyloseq object
saveRDS(object = phyloseq_obj, file = updated_phyloseq_path)
