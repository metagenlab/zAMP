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
phyloseq_object <- snakemake@input[["phyloseq_object"]]
load(file =  file.path(phyloseq_object))
Metadata_table <- snakemake@input[["Metadata_table"]]
metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE)

## Ouput
alpha_plot <- snakemake@output[["alpha_plot"]]

## Parameters
x_axis_column <- snakemake@params[["x_axis_column"]]
grouping_column <- snakemake@params[["grouping_column"]]
color_column <- snakemake@params[["color_column"]]


## Load needed libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("data.table"); packageVersion("data.table")
library("RColorBrewer"); packageVersion("RColorBrewer")

## Set theme
theme_set(theme_bw())
#### BrewerColors
 getPalette = colorRampPalette(brewer.pal(n=8, "Accent"))
 ColList = unique(metadata[[grouping_column]])
 ColPalette = getPalette(length(ColList))
 names(ColPalette) = ColList
 colors_palette <- ColPalette

## Order the x axis as in the metadata_table
sample_data(phyloseq_obj)[[grouping_column]] = factor(sample_data(phyloseq_obj)[[grouping_column]], levels = unique(metadata[[grouping_column]]), ordered = TRUE)

## Plot
p <- plot_richness(phyloseq_obj, x = grouping_column, color = color_column) +
  scale_color_manual(values = colors_palette)

p <- p + theme(axis.text.x = element_text(size=5))


## Save plot
p.width <- 7 + 0.4*length(unique(metadata[[grouping_column]]))
ggsave(filename = alpha_plot,  plot = p, width = p.width, height = 4)



