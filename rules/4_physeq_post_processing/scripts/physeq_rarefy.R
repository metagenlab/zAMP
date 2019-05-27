# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 21.11.18

# inspired from : https://rdrr.io/bioc/phyloseq/man/rarefy_even_depth.html

## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object_path <- snakemake@input[["phyloseq_object"]]

## Ouput
rarefied_phyloseq_path <- snakemake@output[["phyloseq_object"]]

## Parameters
rarefy_value <- snakemake@params[["rarefaction_value"]]
replace_empty_tax <- snakemake@params[["viz_replace_empty_tax"]]


## Load libraries
library(vegan);packageVersion("vegan")
library(dplyr);packageVersion("dplyr")
library(phyloseq);packageVersion("phyloseq")


## Set seed for reproducibility
set.seed(1)

# Read the phyloseq object
phyloseq_obj <- readRDS(phyloseq_object_path)

otu_table(phyloseq_obj) <- t(rrarefy(t(otu_table(phyloseq_obj)), sample = as.numeric(rarefy_value)))


## Add alpha diversity values

### Add alpha diversity indexes to metadata
    alpha_div <- estimate_richness(physeq = phyloseq_obj, split = TRUE)
    sample_data(phyloseq_obj) <- cbind(sample_data(phyloseq_obj),alpha_div)

### Add alpha diveristy indexes at 1% filtration threshold
    ### Keep the IDSs of the taxa above 1%
    physeqrF = filter_taxa(physeq = phyloseq_obj, function(x) mean(x) > 0.01, FALSE)
    ### Keep only those
    physeqaF <- prune_taxa(physeqrF,phyloseq_obj)
    ### Calculate new indexes
    alpha_div_1 <- estimate_richness(physeq = physeqaF, split = TRUE, measure = "Observed")
    ### Rename those
    colnames(alpha_div_1) <- paste0(colnames(alpha_div_1), ("_min_1"))
    ### Again, bind these columns
    sample_data(phyloseq_obj) <- cbind(sample_data(phyloseq_obj),alpha_div_1)

# Write the phyloseq object
saveRDS(object = phyloseq_obj, file = rarefied_phyloseq_path)
