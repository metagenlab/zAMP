# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 16.11.18

# https://joey711.github.io/phyloseq/plot_richness-examples.html

## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


## Input
phyloseq <- snakemake@input[["phyloseq_object"]]


## Ouput

alpha_plot <- snakemake@output[["alpha_plot"]]

## Parameters
x_axis_column <- snakemake@params[["x_axis_column"]]
grouping_column <- snakemake@params[["grouping_column"]]


## Load needed libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")


## Set theme

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

p <- plot_richness(phyloseq, x = paste0('"',x_axis_column,'"'), color = paste0('"',grouping_column,'"') )


ggsave(filename = alpha_plot,  plot = p)
