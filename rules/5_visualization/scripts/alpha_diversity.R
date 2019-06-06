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
    phyloseq_object <- snakemake@input[["phyloseq_object"]]
    Metadata_table <- snakemake@input[["Metadata_table"]]
    metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE)

## Ouput
    alpha_plot <- snakemake@output[["alpha_plot"]]

## Parameters
    sample_label <- snakemake@params[["sample_label"]]
    grouping_column <- snakemake@params[["grouping_column"]]
    sample_type <- snakemake@params[["sample_type"]]
    grouping_filter_column_value <- snakemake@params[["grouping_col_value"]]
    print(grouping_filter_column_value)


## Load needed libraries
    library("ggplot2"); packageVersion("ggplot2")
    library("phyloseq"); packageVersion("phyloseq")
    library("data.table"); packageVersion("data.table")
    library("RColorBrewer"); packageVersion("RColorBrewer")

## Load the phyloseq object
    phyloseq_obj <- readRDS(phyloseq_object)

## Order the x axis as in the metadata_table
    sample_data(phyloseq_obj)[[sample_type]] = factor(sample_data(phyloseq_obj)[[sample_type]], levels = unique(metadata[[sample_type]]), ordered = TRUE)
    sample_data(phyloseq_obj)[[sample_label]] = factor(sample_data(phyloseq_obj)[[sample_label]], levels = unique(metadata[[sample_label]]), ordered = TRUE)


#### BrewerColors
     getPalette = colorRampPalette(brewer.pal(n=8, "Accent"))
     ColList = unique(metadata[[sample_type]])
     ColPalette = getPalette(length(ColList))
     names(ColPalette) = ColList
     colors_palette <- ColPalette

### Open pdf device
    pdf(file = alpha_plot)

### Loop for unique value in grouping_column
    for (i in unique(get_variable(phyloseq_obj, grouping_column))){
        print(paste("Start plotting", grouping_column, i))

        ### Keep only the data of the samples of interest
            remove_idx <- as.character(get_variable(phyloseq_obj, grouping_column)) == i
            g_phyloseq_obj <- prune_samples(remove_idx, phyloseq_obj)


        if(nsamples(g_phyloseq_obj)>0){
            if (sample_label == "Sample"){
            ## Plot
            p <- plot_richness(g_phyloseq_obj, color = sample_type) +
            scale_color_manual(values = colors_palette) +
            geom_boxplot()
            }else{
            p <- plot_richness(g_phyloseq_obj, x = sample_label, color = sample_type) +
            scale_color_manual(values = colors_palette) +
            geom_boxplot()
            }


        p <- p + theme(axis.text.x = element_text(size=5))


        ## Save plot
        p.width <- 7 + 0.4*length(unique(metadata[[sample_label]]))

        if (p.width >= 30){
           p.width <- 30}

         print(p)

        #ggsave(filename = paste0(output_folder,"/",grouping_filter_column_value,"_alpha_diversity.png"),  plot = p, width = p.width, height = 4)

        }else{
            #filename <- paste0(output_folder,"/",grouping_filter_column_value,"_alpha_diversity.png")
            print(paste("Create empty file", filename))
            #file.create(file.path(filename))
            }}

  dev.off()

