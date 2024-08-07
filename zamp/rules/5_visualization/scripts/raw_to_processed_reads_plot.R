# Title     : Raw to processed reads plot
# Objective : Create a plots of reads passing processing
# Created by: valentinscherz
# Created on: 06.06.19


## Redirect R output to the log file
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")


## Input
    raw_to_filtered_reads_stats_path <- snakemake@input[["raw_to_filtered_reads_stats"]]
    raw_to_filtered_reads_stats <- read.table(file = raw_to_filtered_reads_stats_path, sep = "\t", header = TRUE)
    metadata_path <- snakemake@input[["Metadata_table"]]
    metadata <- read.delim(file = metadata_path, sep = "\t", header = TRUE)
    multi_QC_report_path <- snakemake@input[["multi_QC_report_path"]]
    multi_QC_report <- read.table(multi_QC_report_path, header = T)

## Ouput
    reads_plot_with_filtered <- snakemake@output[["reads_plot_with_filtered"]]

## Parameters
    x_axis_column <- snakemake@params[["sample"]]

## Load needed libaries
    library("phyloseq");packageVersion("phyloseq")
    library("data.table");packageVersion("data.table")
    library("dplyr");packageVersion("dplyr")
    library("ggplot2");packageVersion("ggplot2")
    library("RColorBrewer"); packageVersion("RColorBrewer")

## Set theme
theme_set(theme_bw())
    ### Record data on the distribution of number of reads (useful later to scale plots axis)
        smin <- min(multi_QC_report$FastQC_mqc.generalstats.fastqc.total_sequences)
        print(smin)
        smean <- mean(multi_QC_report$FastQC_mqc.generalstats.fastqc.total_sequences)
        print(smean)
        smax <- max(multi_QC_report$FastQC_mqc.generalstats.fastqc.total_sequences)
        print(smax)
    ### Order the x axis as in the metadata_table
        raw_to_filtered_reads_stats[[x_axis_column]] = factor(raw_to_filtered_reads_stats[[x_axis_column]], levels = unique(metadata[[x_axis_column]]), ordered = TRUE)
        #raw_to_filtered_reads_stats[[grouping_column]] = factor(raw_to_filtered_reads_stats[[grouping_column]], levels = unique(metadata[[grouping_column]]), ordered = TRUE)

    ### Order the reads count in logical ordered
        ordered_factors <- c("Tax filtered reads", "Reads processing", "Maintained reads")
        raw_to_filtered_reads_stats$Reads <- factor(raw_to_filtered_reads_stats$Reads, ordered = TRUE, levels = ordered_factors)

    ### Generate colors colors palette
        getPalette <- colorRampPalette(brewer.pal(n=9, "Set3"))
        ColList <- ordered_factors
        ColPalette = getPalette(length(ColList))
        names(ColPalette) = ColList
        colors_palette <- ColPalette


    overall_reads_barplot <- ggplot(raw_to_filtered_reads_stats, aes(x = get(x_axis_column), y = Count, fill = Reads)) +
        geom_col() +
        scale_fill_manual(values = colors_palette) +
        labs(x= x_axis_column,  y ="Reads") +
        ggtitle(paste("Reads counts overall")) +
        scale_x_discrete(drop = TRUE) + # Keep all groups, included the ones with values. Alternative : (drop = FALSE)
        scale_y_continuous(labels = scales::comma, limits = c(0,smax)) +
        theme(axis.text.x = element_text(angle = 90, vjust=0.5, size = 7)) #+
        # guides(fill=guide_legend(title=filling_column))


    ### Save it
        w <- 7 + 0.07 * (length(unique(raw_to_filtered_reads_stats[[x_axis_column]]))) 
        ggsave(overall_reads_barplot, filename = reads_plot_with_filtered, width = w, height = 5)
