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
load(file =  file.path(phyloseq_object))
Metadata_table <- snakemake@input[["Metadata_table"]]
metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE)

## Output
output_folder <- snakemake@output[["output1"]]
output_folder
output_folder <- (dirname(output_folder)[1])


## Parameters
x_axis_column <- snakemake@params[["x_axis_column"]]
grouping_column <- snakemake@params[["grouping_column"]]
sample_type <- snakemake@params[["sample_type"]]



## Load needed libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer"); packageVersion("RColorBrewer")



## Ordination

### Remove sequences not assigned at the phylum level
physeq_bacteria_only <- subset_taxa(phyloseq_obj, Kingdom == "Bacteria")
physeq_no_unassigned_phylum_bact_only <- subset_taxa(physeq_bacteria_only, Phylum != "Bacteria_phy")


#### BrewerColors
 getPalette = colorRampPalette(brewer.pal(n=8, "Accent"))
 ColList = unique(metadata[[sample_type]])
 ColPalette = getPalette(length(ColList))
 names(ColPalette) = ColList
 colors_palette <- ColPalette






  ### Order the x axis as in the metadata_table
   sample_data(physeq_no_unassigned_phylum_bact_only)[[grouping_column]] = factor(sample_data(physeq_no_unassigned_phylum_bact_only)[[grouping_column]], levels = unique(metadata[[grouping_column]]), ordered = TRUE)

  ### Create a list of all ordination methods
  dist_methods <- c("unifrac" , "wunifrac", "jsd", "bray", "jaccard", "chao")

    ### Run a loop to save in a list all plots
      ### Create a liste
      plist <- vector("list", length(dist_methods))
      ### Rename entries in the list
      names(plist) = dist_methods
      ### Loop over all methods
      for(i in dist_methods){

          # Calculate distance matrix
          iDist <- phyloseq::distance(physeq_no_unassigned_phylum_bact_only, method=i)
          # Calculate ordination
          iMDS  <- ordinate(physeq_no_unassigned_phylum_bact_only, "MDS", distance=iDist)
          ## Make plot
            # Create plot, store as temp variable, p
            p <- plot_ordination(physeq_no_unassigned_phylum_bact_only, iMDS, shape = grouping_column, color = sample_type) +
              scale_color_manual(values = colors_palette) +
              geom_point(size=4) + stat_ellipse(aes(group = get(grouping_column), color = get(grouping_column)),linetype = 2, type = "t") ## Will be needed which variable comes here. Could also be grouping_column
            # Add title to each plot
            p <- p + ggtitle(paste("MDS using distance method ", i, sep=" "))
            # Save the individual graph in a folder
            ggsave(plot = p, filename = paste0(output_folder,"/",i,".png"))

      }

