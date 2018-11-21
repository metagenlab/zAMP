# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 21.11.18

# inspired from : https://rdrr.io/bioc/phyloseq/man/rarefy_even_depth.html

## Redirect R output to the log file
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")


## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]
load(file =  file.path(phyloseq_object))


## Ouput
rarefied_phyloseq <- snakemake@output[["rarefied_phyloseq"]]

## Parameters
x_axis_column <- snakemake@params[["x_axis_column"]]
grouping_column <- snakemake@params[["grouping_column"]]

## Load libraries
library('phyloseq')
library('ggplot2')
library('plyr') # ldply
library('reshape2') # melt
library('vegan')


rarefy_even_depth(physeq = phyloseq_obj, sample.size = min(sample_sums(physeq)),
  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
