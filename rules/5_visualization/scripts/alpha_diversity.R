# Title     : Alpha-diversity
# Objective : Plot alpha-diversity based on phyloseq build-in functions
# Created by: valentinscherz
# Created on: 06.06.18

# https://joey711.github.io/phyloseq/plot_richness-examples.html

## Redirect R output to the log file
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    metadata <- snakemake@input[["metadata"]]

## Ouput
    alpha_plot <- snakemake@output[["alpha_plot"]]

## Parameters
    sample_label <- snakemake@params[["sample_label"]]
    grouping_column <- snakemake@params[["grouping_column"]]
    color_factor <- snakemake@params[["color_factor"]]
    alpha_metric <- snakemake@params[["alpha_metric"]]
    facet_plot <- snakemake@params[["facet_plot"]]
    facetting_column <- snakemake@params[["facetting_column"]]


## Load needed libraries
    library("ggplot2"); packageVersion("ggplot2")


## Load the data
    metadata_table <- read.table(file = metadata, sep = "\t", header = TRUE)

## Order the x axis as in the metadata_table
    #sample_data(phyloseq_obj)[[color_factor]] = factor(sample_data(phyloseq_obj)[[color_factor]], levels = unique(metadata[[color_factor]]), ordered = TRUE)
    #sample_data(phyloseq_obj)[[sample_label]] = factor(sample_data(phyloseq_obj)[[sample_label]], levels = unique(metadata[[sample_label]]), ordered = TRUE)


### Open pdf device
    pdf(file = alpha_plot)

### Loop for unique value in grouping_column
    #for (i in unique(metadata_table[[grouping_column]])){

        ### Keep only the data of the samples of interest
        #metadata_i <- metadata_table[metadata_table[[grouping_column]]==i,]

        p <- ggplot(metadata_table, aes(x=get(grouping_column), y=get(alpha_metric), color = get(color_factor))) +
          geom_boxplot() +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5)) +
          #facet_grid(~storage_time, drop = TRUE, scales = "free_x") +
          xlab(label = sample_label) +
          ylab(label = alpha_metric)+
          labs(color = color_factor) +
          ggtitle(paste(grouping_column, alpha_metric, "alpha diversity" ))

        if (isTRUE(facet_plot)){
            p <- p + facet_grid(~ get(facetting_column), scales = "free_x", drop = TRUE)

        }

        ## Save plot
        p.width <- 7 + 0.4*length(unique(metadata_table[[grouping_column]]))

        if (p.width >= 30){
           p.width <- 30}

         print(p)

    #    }

  dev.off()

