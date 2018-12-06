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
output_folder <- (dirname(alpha_plot)[1])

## Parameters
x_axis_column <- snakemake@params[["x_axis_column"]]
grouping_column <- snakemake@params[["grouping_column"]]
sample_type <- snakemake@params[["sample_type"]]


## Load needed libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("data.table"); packageVersion("data.table")
library("RColorBrewer"); packageVersion("RColorBrewer")

## Order the x axis as in the metadata_table
#sample_data(phyloseq_obj)[[grouping_column]] = factor(sample_data(phyloseq_obj)[[grouping_column]], levels = unique(metadata[[grouping_column]]), ordered = TRUE)


#### BrewerColors
 getPalette = colorRampPalette(brewer.pal(n=8, "Accent"))
 ColList = unique(metadata[[sample_type]])
 ColPalette = getPalette(length(ColList))
 names(ColPalette) = ColList
 colors_palette <- ColPalette


for (g in get_variable(phyloseq_obj, grouping_column)){
    remove_idx = as.character(get_variable(phyloseq_obj, grouping_column)) == g
    g_phyloseq_obj = prune_samples(remove_idx, phyloseq_obj)

 if (x_axis_column == "Sample"){
 ## Plot
 p <- plot_richness(g_phyloseq_obj, color = sample_type) +
  scale_color_manual(values = colors_palette) +
    geom_boxplot()
 }else{
 p <- plot_richness(g_phyloseq_obj, x = x_axis_column, color = sample_type) +
  scale_color_manual(values = colors_palette) +
  geom_boxplot()
 }




p <- p + theme(axis.text.x = element_text(size=5))


## Save plot
p.width <- 7 + 0.4*length(unique(metadata[[sample_type]]))
ggsave(filename = paste0(output_folder,"/",g,"_alpha_divesity.png"),  plot = p, width = p.width, height = 4)

}

