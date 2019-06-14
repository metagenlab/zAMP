# Title     : Taxonomic barplots
# Objective : Plot taxonomy barplots
# Created by: valentinscherz
# Created on: 06.06.19

## Redirect R output to the log file
    #log <- file(snakemake@log[[1]], open="wt")
    #sink(log)
    #sink(log, type="message")

## Input
    phyloseq_melted_table <- snakemake@input[["phyloseq_melted_table"]]

## Ouput
    output_path <- snakemake@output[["barplot"]]

## Parameters
    sample_label <- snakemake@params[["sample_label"]]
    grouping_column <- snakemake@params[[ "grouping_column"]]
    t_neg_PCR_sample_on_plots <-  snakemake@params[["t_neg_PCR_sample_on_plots"]]
    t_neg_PCR_group_column_value <-  snakemake@params[["t_neg_PCR_group_column_value"]]
    relative_or_absolute_filtering <-  snakemake@params[["relative_or_absolute_filtering"]]
    filtering_value <-  snakemake@params[["filtering_value"]]
    relative_or_absolute_plot <-  snakemake@params[["relative_or_absolute_plot"]]
    plotting_tax_ranks <-  snakemake@params[["barplot_tax_ranks"]]
    distinct_colors <-  snakemake@params[["distinct_colors"]]
    horizontal_plot <-  snakemake@params[["horizontal_plot"]]
    facet_plot <-  snakemake@params[["facet_plot"]]
    facetting_column <-  snakemake@params[["facetting_column"]]
    order_by_abundance <-  snakemake@params[["order_by_abundance"]]
    separated_legend <-  snakemake@params[["separated_legend"]]
    abundance_filtre_level <- snakemake@params[["abundance_filtre_level"]]

## Load needed libraries
    library(ggplot2); packageVersion("ggplot2")
    library(dplyr); packageVersion("dplyr")
    library(RColorBrewer); packageVersion("RColorBrewer")
    library(randomcoloR); packageVersion("randomcoloR")
    library(data.table); packageVersion("data.table")
    library(forcats); packageVersion("forcats")
    library(rlang); packageVersion("rlang")
    library(ggpubr); packageVersion("ggpubr")
    library(cowplot); packageVersion("cowplot")

## Load the melted phyloseq table
    melted_dataframe <- read.csv(file.path(phyloseq_melted_table), header = TRUE, sep = "\t")

## Order the x axis as in the metadata_table
    #melted_dataframe[[sample_label]] = factor(melted_dataframe[[sample_label]], levels = unique(melted_dataframe[[sample_label]]), ordered = TRUE)

################################################################################

### Create a function
    barplots_fct <- function(melted_dataframe, x_axis_column, grouping_column, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value, relative_or_absolute_filtering = c("relative", "absolute", "nofiltering"), filtering_value, relative_or_absolute_plot = c("relative", "absolute"), plotting_tax_ranks = "all", output_path, distinct_colors = TRUE, horizontal_plot = FALSE, facet_plot = FALSE, facetting_column = NULL, order_by_abundance = TRUE, separated_legend, abundance_filtre_level = c("Sample", "Group")){

        ## Select Tax ranks needed for labelling
        if(plotting_tax_ranks == "Kingdom"){
            label_ranks <- c("Kingdom", "")
        } else if(plotting_tax_ranks == "Phylum"){
            label_ranks <- c("Kingdom", "Phylum")

        } else if(plotting_tax_ranks == "Class"){
            label_ranks <- c("Phylum", "Class")

        } else if(plotting_tax_ranks == "Order"){
            label_ranks <- c("Class", "Order")

        } else if(plotting_tax_ranks == "Family"){
            label_ranks <- c("Order", "Family")

        } else if(plotting_tax_ranks == "Genus"){
            label_ranks <- c("Family", "Genus")

        } else if(plotting_tax_ranks == "Species"){
            label_ranks <- c("Phylum", "Species")

        } else if(plotting_tax_ranks == "OTU"){
            label_ranks <- c("Species", "OTU")
        }


        ## Unquote factors
        g_column <- rlang::sym(grouping_column)
        x_column <- rlang::sym(x_axis_column)
        t_column <- rlang::sym(plotting_tax_ranks)
        l_ranks  <-  rlang::syms(label_ranks)

        if (facet_plot == TRUE){
           f_column <- rlang::sym(facetting_column)
        }

        ## Define grouping level for abundance-base filtering (Sample vs group of Samples)
        if (abundance_filtre_level == "Sample"){
            abundance_sample <- rlang::sym("Sample")
        } else if (abundance_filtre_level == "Group"){
            abundance_sample <- rlang::sym(grouping_column)
        } else {
            stop('"abundance_filtre_level" must be "Sample" or "Group"')
        }

        ## Transform abundance on 100%
        physeq_subset_norm_df <- melted_dataframe  %>% # calculate % normalized Abundance
        ungroup %>%
        dplyr::group_by(Sample) %>%
        dplyr::mutate(per=as.numeric(100*Abundance/sum(Abundance))) %>%
        ungroup() %>%
        dplyr::select(-Abundance)  %>%
        dplyr::rename(Abundance = per)

        ## Choose between relative value or absolute value dataframe for plotting value
            ### Relative value plotting
            if (relative_or_absolute_plot == "relative"){
                plotting <- "Relative"
                print("Relative value plotting")
                plotted_df <- physeq_subset_norm_df

            ### Absolute value plotting
            }else if (relative_or_absolute_plot == "absolute"){
                plotting <- "Absolute"
                print("Absolute value plotting")
                plotted_df <- melted_dataframe
            }else{stop('"relative_value_plotting" must be "relative" or "absolute"')
            }


        ## Filter values at a relative or absolute quantity
            ### No filtering
            if (relative_or_absolute_filtering == "nofiltering"){
                filtering <- "Nofiltering"
                filtering_value <- "0"
                print("No filtering based on quantity")
                physeq_subset_df_filtered <- physeq_subset_norm_df

            ### Relative value filtering
            }else if (relative_or_absolute_filtering == "relative"){
                filtering <- "Relative"
                filtering_value <- filtering_value
                print("Relative value based filtering")
                physeq_subset_df_filtered <- physeq_subset_norm_df  %>%
                    dplyr::group_by(abundance_sample, !!t_column) %>% # group the dataframe by Sample and taxa
                    dplyr::mutate(sumper=as.numeric(mean(Abundance))) %>% # calculate the cumulative relative abundance of the taxa in the sample
                    dplyr::filter(sumper >= filtering_value) %>%# keep only the taxa above threshold.
                    ungroup

            ### Absolute value filtering
            }else if (relative_or_absolute_filtering == "absolute"){
                filtering <- "Absolute"
                print("Absolute value based filtering")
                physeq_subset_df_filtered <- melted_dataframe  %>%
                    dplyr::group_by(abundance_sample, !!t_column) %>% # group the dataframe by Sample and taxa
                    dplyr::mutate(sumper=as.numeric(mean(Abundance))) %>% # calculate the cumulative relative abundance of the taxa in the sample
                    dplyr::filter(sumper >= filtering_value) %>%# keep only the taxa above threshold.
                    ungroup
            }

        ## Filter the plotted DF depending on if it is above or under the threshold
        above_threshold_df <- semi_join(plotted_df, physeq_subset_df_filtered, by = c("OTU", "Sample")) ### In a dataframe, keep only the rows matching the filtered rowmanes
        under_threshold_df <- anti_join(plotted_df, physeq_subset_df_filtered, by = c("OTU", "Sample")) ### In another dataframe, keep only the lines NOT matching the filtered rowmanes

        ## Rename taxa with < percent/reads abundance, defending of the filtering applied
            ### Define the filtering tag depending of the applied filtering
                ##### Relative
                if (relative_or_absolute_filtering == "relative"){
                    filtering_tag <-paste("Relative abund. <", filtering_value , "%")

                ##### Absolute filtering
                }else if (relative_or_absolute_filtering == "absolute"){
                        filtering_tag <-paste("Absolute abund. <", filtering_value, "reads")

                ##### No filtering
                }else if (relative_or_absolute_filtering == "nofiltering"){
                        filtering_tag <-paste("nofiltering")
                }

            ### Apply the filtering tag
            if (relative_or_absolute_filtering != "nofiltering"){
                under_threshold_df$Kingdom <-filtering_tag
                under_threshold_df$Phylum <-filtering_tag
                under_threshold_df$Class <-filtering_tag
                under_threshold_df$Order <-filtering_tag
                under_threshold_df$Family <-filtering_tag
                under_threshold_df$Genus <-filtering_tag
                under_threshold_df$Species <-filtering_tag
                under_threshold_df$OTU <-filtering_tag
            }

        ## Join the two dataframes, to put back all rows together, the ones above threshold keeping their original taxonomic identifier while the others are now grouped together
            threshod_filtered_abs <- full_join(under_threshold_df,above_threshold_df)

        ## Remove the (now but needed before) Abundance = 0 rows
            threshod_filtered_abs_no_zero <- filter(threshod_filtered_abs, threshod_filtered_abs$Abundance>0)

        ## Regroup taxonomically identical rows
        if (facet_plot == TRUE){
            threshod_filtered_abs_no_zero <- threshod_filtered_abs_no_zero %>%
            group_by(!!(x_column), !!(g_column), !!(f_column), !!!(l_ranks), !!(t_column)) %>%
            summarise(Abundance = sum(Abundance)) %>%
            ungroup()
        }else{
            threshod_filtered_abs_no_zero <- threshod_filtered_abs_no_zero %>%
            group_by(!!(x_column), !!(g_column), !!!(l_ranks),  !!(t_column)) %>%
            summarise(Abundance = sum(Abundance)) %>%
            ungroup()
        }

        ## Reorder the facet factor if later used for plotting
        if (isTRUE(facet_plot)){
            threshod_filtered_abs_no_zero[[facetting_column]] <- fct_reorder(threshod_filtered_abs_no_zero[[facetting_column]], as.numeric(threshod_filtered_abs_no_zero[[x_axis_column]]))

        }else if (isFALSE(facet_plot)){print("No faceting")

        }else {stop('"facet_plot" must be TRUE or FALSE')
        }

        ## Reverse the x_axis_column column for later if using horizontal barplot
        if (isTRUE(horizontal_plot)){
            threshod_filtered_abs_no_zero[[x_axis_column]] <- fct_rev(threshod_filtered_abs_no_zero[[x_axis_column]])

        } else if (isFALSE(horizontal_plot)){print("Vertical plotting")

        } else {stop('"horizontal_plot" must be TRUE or FALSE')
        }


        ### Set colors palette
            #### Brewer colors
            if (distinct_colors == FALSE){
                getPalette <- colorRampPalette(brewer.pal(n=9, "Set1"))
                ColList <- unique(threshod_filtered_abs_no_zero[[plotting_tax_ranks]])
                ColPalette = getPalette(length(ColList))
                names(ColPalette) = ColList
                colors_palette <- ColPalette
                colors_palette[filtering_tag] <- "#d3d3d3" # Set the filtered in balck

            #### Randomcolors distinct colors
            }else if (distinct_colors == TRUE){
                set.seed(4)
                ColList <- unique(threshod_filtered_abs_no_zero[[plotting_tax_ranks]])
                print(length(unique(threshod_filtered_abs_no_zero[[plotting_tax_ranks]])))
                ColPalette <- distinctColorPalette(altCol = TRUE, k = length(unique(threshod_filtered_abs_no_zero[[plotting_tax_ranks]])))
                names(ColPalette) = ColList
                colors_palette <- ColPalette
                colors_palette[filtering_tag] <- "#d3d3d3" # Set the filtered in balck

            }else{print("distinct_colors must be TRUE or FALSE")
            }

        ## Open pdf device
        pdf(file = output_path)

        ## Loop over unique value in grouping_column
        for (i in unique(threshod_filtered_abs_no_zero[[grouping_column]])) {
            print(paste("Start plotting", grouping_column, i))
            ### filter the table for this value of the grouping columns. Depending of the used arguments, the t_neg_PCR values are kept or not on the barplots
                #### Keep t_neg_PCR rows
                if (isTRUE(t_neg_PCR_sample_on_plots) & !is.null(t_neg_PCR_group_column_value)){
                    print('Keeping t_neg_PCR values for the graphs. The "t_neg_PCR_group_column_value" must match the one in the grouping_column for this sample')
                    filtered_df_abs_i <- filter(threshod_filtered_abs_no_zero, threshod_filtered_abs_no_zero[[grouping_column]] == i | threshod_filtered_abs_no_zero[[grouping_column]] == t_neg_PCR_group_column_value)

                }else if (isTRUE(t_neg_PCR_sample_on_plots) & is.null(t_neg_PCR_group_column_value)){
                        stop('If "t_neg_PCR_sample_on_plots" is "TRUE, a "t_neg_PCR_group_column_value" indicating the value of the T neg PCR sample in the grouping column must be indicated')
                ##### Do not keep t_neg_PCR rows
                }else if (isFALSE(t_neg_PCR_sample_on_plots)){
                            filtered_df_abs_i <- filter(threshod_filtered_abs_no_zero, threshod_filtered_abs_no_zero[[grouping_column]] == i)
                }else{ stop('"t_neg_PCR_sample_on_plots" must be TRUE or FALSE')
                }

            ### Reorder by abundance
            if (isTRUE(order_by_abundance)){
                print("Ordered by abundance")
                filtered_df_abs_i <- filtered_df_abs_i %>%
                    group_by(!! t_column) %>%
                    mutate(tot = sum(Abundance)) %>%
                    ungroup() %>%
                    mutate(!! t_column := fct_reorder(!! t_column, tot)) %>%
                    arrange(desc(tot))
            }else if (isFALSE(order_by_abundance)){ print("NOT ordered by abundance")

            }else {stop('"order_by_abundance" must be TRUE or FALSE')
            }

            ### Set the filtering_tag at the top of the plot
            filtered_df_abs_i[[plotting_tax_ranks]] <- fct_relevel(filtered_df_abs_i[[plotting_tax_ranks]], filtering_tag, after = 0)

            ### Renames the values in the vector of plotted used taxrank with their related OTU name to keep them matching. Without this step, the labels do NOT match the rows
            taxalabel <- as.character(paste(toupper(substr(filtered_df_abs_i[[label_ranks[1]]],1 ,6)) , filtered_df_abs_i[[label_ranks[2]]], sep = ";"))
            names(taxalabel) <- filtered_df_abs_i[[plotting_tax_ranks]]

            ### Create the barplot
            taxrank_barplot <- filtered_df_abs_i %>%
                ggplot(aes(x = get(x_axis_column), y = Abundance, fill = get(plotting_tax_ranks))) +
                theme_bw() +
                geom_col() +
                #scale_x_discrete(labels = x_labels, drop = TRUE) + # Set false to keep empty bars
                scale_y_discrete(labels = taxalabel) +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), legend.text=element_text(size=5), plot.title = element_text(hjust = 0.5)) + # axis and title settings
                guides(fill = guide_legend(title = paste0(plotting_tax_ranks),reverse = FALSE, keywidth = 0.5, keyheight = 0.5, ncol = 1)) + # settings of the legend
                labs(x = x_axis_column,  y = paste(plotting, "abundance"), title = paste("Taxonomic composition", plotting_tax_ranks,"level", i)) + # axis and graph title
                scale_fill_manual(values = colors_palette) # colors as set previously


                #### In option, turn barplot horizontally
                if (isTRUE(horizontal_plot)){
                    taxrank_barplot <- taxrank_barplot + coord_flip() ### Reverse the order of the samples
                    #### In option, create facet view
                    if (isTRUE(facet_plot)){
                        taxrank_barplot <- taxrank_barplot + facet_grid(get(facetting_column) ~., scales = "free", space = "free")
                }

                #### In option, create facet view
                }else if(!isTRUE(horizontal_plot)){
                    #### In option, create facet view
                    if (isTRUE(facet_plot)){
                        taxrank_barplot <- taxrank_barplot + facet_grid(~ get(facetting_column), scales = "free", space = "free")
                    }
                }

                #### In option, print the legend separated or not
                if (isTRUE(separated_legend)){
                    taxrank_barplot_no_leg <- (taxrank_barplot + guides(fill = FALSE))
                    leg <- as_ggplot(get_legend(taxrank_barplot))
                    print(taxrank_barplot_no_leg)
                    print(leg)

                }else {
                    print(taxrank_barplot)
                }


        }
        dev.off()
    }
################################################################################



## Run the function
    barplots_fct(
      melted_dataframe = melted_dataframe,
        grouping_column = grouping_column,
        x_axis_column = sample_label,
        t_neg_PCR_sample_on_plots = t_neg_PCR_sample_on_plots,
        t_neg_PCR_group_column_value = t_neg_PCR_group_column_value,
        relative_or_absolute_filtering = relative_or_absolute_filtering,
        filtering_value = filtering_value,
        relative_or_absolute_plot = relative_or_absolute_plot,
        plotting_tax_ranks = plotting_tax_ranks,
        output_path = output_path,
        distinct_colors = distinct_colors,
        horizontal_plot = horizontal_plot,
        facet_plot = facet_plot,
        facetting_column = facetting_column,
        order_by_abundance = order_by_abundance,
        separated_legend = separated_legend,
        abundance_filtre_level = abundance_filtre_level)
