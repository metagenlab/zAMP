# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 26.10.18

## Inspired from https://stackoverflow.com/questions/15260245/r-convert-text-field-to-function
## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[[1]]

## Output
phyloseq_filtered_object <- snakemake@output[[1]]

## Parameters
subset_formula <- snakemake@params[["subset_formula"]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")
library(plyr);packageVersion("plyr")
library(dplyr);packageVersion("dplyr")

## Load the phyloseq phyloseq_object
phyloseq_object <- readRDS(phyloseq_object)

## Convert the fonction to be used later
subset_fct <- function(x){}
body(subset_fct) <- as.quoted(subset_formula)[[1]]

## filter taxa
subset_features <- filter_taxa(phyloseq_object, subset_fct, TRUE)

## Remove already in metadata alphia diversity values
sample_data(subset_features) <- select(sample_data(subset_features), -c(Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher))

## Add alpha diversity indexes to metadata
alpha_div <- estimate_richness(physeq = subset_features, split = TRUE)
sample_data(subset_features) <- cbind(sample_data(subset_features),alpha_div)

## Write the new phyloseq object
saveRDS(object = subset_features, file = phyloseq_filtered_object)
