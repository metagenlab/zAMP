# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 16.11.18

# https://joey711.github.io/phyloseq/plot_richness-examples.html

## Redirect R output to the log file
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")


## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]
load(file =  file.path(phyloseq_object))
Metadata_table <- snakemake@input[["Metadata_table"]]
metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE)

## Ouput
alpha_plot <- snakemake@output[["alpha_plot"]]

## Parameters
x_axis_column <- snakemake@params[["x_axis_column"]]
grouping_column <- snakemake@params[["grouping_column"]]


## Load needed libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("data.table"); packageVersion("data.table")

## Set theme
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}


## Order the x axis as in the metadata_table
sample_data(phyloseq_obj)[[x_axis_column]] = factor(sample_data(phyloseq_obj)[[x_axis_column]], levels = metadata[[x_axis_column]], ordered = TRUE)

## Plot
p <- plot_richness(phyloseq_obj, x = x_axis_column, color = grouping_column)
p <- p + theme(axis.text.x = element_text(size=5))



# Oh no, the table wasn't ordered
# library("data.table")
# newtab = data.table(p$data)
# setorder(newtab, samples)
# p$data <- newtab

## Save plot
p.width <- 7 + 0.3*nsamples(phyloseq_obj)
ggsave(filename = alpha_plot,  plot = p, width = p.width, height = 4)




