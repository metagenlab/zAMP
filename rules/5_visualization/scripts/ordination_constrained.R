# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 19.11.18

## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]
Metadata_table <- snakemake@input[["Metadata_table"]]
metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE)

## Output
output_path<- snakemake@output[["output1"]]

## Parameters
ordination_distance = snakemake@params[["ordination_distance"]]
# sample_label <- snakemake@params[["sample_label"]]
grouping_column <- snakemake@params[["grouping_column"]]
grouping_filter_column_value <- snakemake@params[["grouping_col_value"]]
sample_type <- snakemake@params[["sample_type"]]
ordination_factor <- snakemake@params[["ordination_factor"]]
ordination_method <-   snakemake@params[["ordination_method"]]

## Load needed libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer"); packageVersion("RColorBrewer")
library("rlang"); packageVersion("rlang")

## Set seed for reproducibility
set.seed(100)

## Load the phyloseq object
phyloseq_obj <- readRDS(phyloseq_object)

### Remove sequences not assigned at the phylum level
#physeq_bacteria_only <- subset_taxa(phyloseq_obj, Kingdom == "Bacteria")
#physeq_no_unassigned_phylum_bact_only <- subset_taxa(physeq_bacteria_only, Phylum != "Bacteria_phy")

### Remove sample with abundance < 20
physeq_filtered<- prune_samples(sample_sums(phyloseq_obj)>20, phyloseq_obj)

#### BrewerColors
 getPalette = colorRampPalette(brewer.pal(n=8, "Dark2"))
 ColList = unique(metadata[[sample_type]])
 ColPalette = getPalette(length(ColList))
 names(ColPalette) = ColList
 colors_palette <- ColPalette
 ### Order the x axis as in the metadata_table
    sample_data(physeq_filtered)[[sample_type]] = factor(sample_data(physeq_filtered)[[sample_type]], levels = unique(metadata[[sample_type]]), ordered = TRUE)

### Keep only the data of the samples of interest
    remove_idx <- as.character(get_variable(physeq_filtered, grouping_column)) == grouping_filter_column_value
    g_physeq_filtered = prune_samples(remove_idx, physeq_filtered)
  
    if(nsamples(g_physeq_filtered)>3){


            # Calculate ordination
            #iMDS  <- do.call(ordinate, list(physeq = g_physeq_no_unassigned_phylum_bact_only, formula = paste0("~", ordination_factor), method = ordination_method))
            iMDS <- ordinate(physeq = g_physeq_filtered, formula = ~ get(ordination_factor), method = ordination_method)
            ## Make plot
            # Create plot, store as temp variable, p
            p <- plot_ordination(g_physeq_filtered, iMDS, color = sample_type, shape = ordination_factor) +
              scale_color_manual(values = colors_palette) +
              geom_point(size=4) + stat_ellipse(aes(group = get(sample_type), color = get(sample_type)),linetype = 2, type = "t") ## Will be needed which variable comes here. Could also be grouping_column
            # Add title to each plot
            p <- p + ggtitle(paste(ordination_method, "constrained on", ordination_factor))
            # Save the individual graph in a folder
            ggsave(plot = p, filename = output_path)

    }else{

        file.create(file.path(output_path))
    }
