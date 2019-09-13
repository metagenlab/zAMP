# Title     : Constrained PCA
# Objective : Create constrained PCA
# Created by: valentinscherz
# Created on: 06.06.19

## Redirect R output to the log file
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    phyloseq_object <- snakemake@input[["phyloseq_object"]]
    Metadata_table <- snakemake@input[["Metadata_table"]]
    metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE)

## Output
    output_path<- snakemake@output[["constrained_ordination"]]

## Parameters
    grouping_column <- snakemake@params[["grouping_column"]]
    color_factor <- snakemake@params[["color_factor"]]
    shape_factor <- snakemake@params[["shape_factor"]]
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


### Remove sample with abundance < 20
    physeq_filtered<- prune_samples(sample_sums(phyloseq_obj)>20, phyloseq_obj)

#### BrewerColors
     getPalette = colorRampPalette(brewer.pal(n=8, "Dark2"))
     ColList = unique(metadata[[color_factor]])
     ColPalette = getPalette(length(ColList))
     names(ColPalette) = ColList
     colors_palette <- ColPalette
     ### Order the x axis as in the metadata_table
     sample_data(physeq_filtered)[[color_factor]] = factor(sample_data(physeq_filtered)[[color_factor]], levels = unique(metadata[[color_factor]]), ordered = TRUE)

### Open pdf device
    pdf(file = output_path)

### Loop for unique value in grouping_column
    for (i in unique(get_variable(physeq_filtered, grouping_column))){
        print(paste("Start plotting", grouping_column, i))

    ### Keep only the data of the samples of interest
        remove_idx <- as.character(get_variable(physeq_filtered, grouping_column)) == i
        g_physeq_filtered = prune_samples(remove_idx, physeq_filtered)

        if(nsamples(g_physeq_filtered)>3){

                # Calculate ordination
                iMDS <- ordinate(physeq = g_physeq_filtered, formula = ~ get(shape_factor), method = ordination_method)
                ## Make plot
                # Create plot, store as temp variable, p
                p <- plot_ordination(g_physeq_filtered, iMDS, color = color_factor, shape = shape_factor) +
                  scale_color_manual(values = colors_palette) +
                  geom_point(size=4) + stat_ellipse(aes(group = get(color_factor), color = get(color_factor)),linetype = 2, type = "t") ## Will be needed which variable comes here. Could also be grouping_column
                # Add title to each plot
                p <- p + ggtitle(paste(i, ordination_method, "constrained on", shape_factor))
                # Save the individual graph in a folder
                #ggsave(plot = p, filename = output_path)
                print(p)
        }else{
            print("To few point to create ordination plot")
            #file.create(file.path(output_path))
    }
    }

dev.off()
