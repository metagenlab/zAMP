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
Metadata_table <- snakemake@input[["Metadata_table"]]
metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE)

## Ouput
alpha_plot <- snakemake@output[["alpha_plot"]]
print(paste("output", alpha_plot))

output_folder <- (dirname(alpha_plot)[1])

## Parameters
x_axis_column <- snakemake@params[["x_axis_column"]]
grouping_column <- snakemake@params[["grouping_column"]]
sample_type <- snakemake@params[["sample_type"]]
grouping_column_value <- snakemake@params[["grouping_column_value"]]
print(grouping_column_value)


## Load needed libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("data.table"); packageVersion("data.table")
library("RColorBrewer"); packageVersion("RColorBrewer")

## Load the phyloseq object
phyloseq_obj <- readRDS(phyloseq_object)

## Order the x axis as in the metadata_table
sample_data(phyloseq_obj)[[sample_type]] = factor(metadata[[sample_type]], levels = unique(metadata[[sample_type]]), ordered = TRUE)
sample_data(phyloseq_obj)[[x_axis_column]] = factor(metadata[[x_axis_column]], levels = unique(metadata[[x_axis_column]]), ordered = TRUE)

### Remove sequences not assigned at the phylum level
physeq_bacteria_only <- subset_taxa(phyloseq_obj, Kingdom == "Bacteria")
physeq_no_unassigned_phylum_bact_only <- subset_taxa(physeq_bacteria_only, Phylum != "Bacteria_phy")

#### BrewerColors
 getPalette = colorRampPalette(brewer.pal(n=8, "Accent"))
 ColList = unique(metadata[[sample_type]])
 ColPalette = getPalette(length(ColList))
 names(ColPalette) = ColList
 colors_palette <- ColPalette


#for (g in get_variable(phyloseq_obj, grouping_column)){
    remove_idx = as.character(get_variable(physeq_no_unassigned_phylum_bact_only, grouping_column)) == grouping_column_value
    g_phyloseq_obj = prune_samples(remove_idx, phyloseq_obj)

 if(nsamples(g_phyloseq_obj)>0){
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
p.width <- 7 + 0.4*length(unique(metadata[[x_axis_column]]))
ggsave(filename = paste0(output_folder,"/",grouping_column_value,"_alpha_divesity.png"),  plot = p, width = p.width, height = 4)

}else{
    filename <- paste0(output_folder,"/",grouping_column_value,"_alpha_divesity.png")
    print(paste("Create empty file", filename))
    file.create(file.path(filename))
} #}
