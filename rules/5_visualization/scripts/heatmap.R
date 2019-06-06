# Title     : Heatmaps
# Objective : Creata taxonomic heatmaps
# Created by: valentinscherz
# Created on: 06.06.19

## Redirect R output to the log file
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    phyloseq_melted_table <- snakemake@input[["phyloseq_melted_table"]]

## Ouput
    output_path <- snakemake@output[["heatmap"]]

## Parameters
    sample_label =snakemake@params[[ "sample_label"]]
    grouping_column <- snakemake@params[[ "grouping_column"]]
    t_neg_PCR_sample_on_plots = snakemake@params[["t_neg_PCR_sample_on_plots"]]
    t_neg_PCR_group_column_value = snakemake@params[["t_neg_PCR_group_column_value"]]
    relative_or_absolute_filtering = snakemake@params[["relative_or_absolute_filtering"]]
    filtering_value = snakemake@params[["filtering_value"]]
    relative_or_absolute_plot = snakemake@params[["relative_or_absolute_plot"]]
    plotting_tax_ranks = snakemake@params[["plotting_tax_ranks"]]
    horizontal_barplot = snakemake@params[["horizontal_barplot"]]
    facet_plot = snakemake@params[["facet_plot"]]
    facetting_column = snakemake@params[["facetting_column"]]
    #order_by = snakemake@params[["order_by"]]


## Load needed libraries
    library(ggplot2); packageVersion("ggplot2")
    library(dplyr); packageVersion("dplyr")
    library(phyloseq); packageVersion("phyloseq")
    library(RColorBrewer); packageVersion("RColorBrewer")
    library(data.table); packageVersion("data.table")
    library(forcats); packageVersion("forcats")
    library(rlang); packageVersion("rlang")
    library(grid); packageVersion("grid")
    library(cowplot); packageVersion("cowplot")
    library(reshape2); packageVersion("reshape2")

## Load the melted phyloseq table
    melted_dataframe<- read.csv(file.path(phyloseq_melted_table), header = TRUE, sep = "\t")

## Order the x axis as in the metadata_table
#sample_data(phyloseq_obj)[[sample_type]] = factor(sample_data(phyloseq_obj)[[sample_type]], levels = unique(metadata[[sample_type]]), ordered = TRUE)
#sample_data(phyloseq_obj)[[sample_label]] = factor(sample_data(phyloseq_obj)[[sample_label]], levels = unique(metadata[[sample_label]]), ordered = TRUE)

################################################################################
### Create a function
    heatmap_fct <- function(melted_dataframe, x_axis_column, grouping_column, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value, relative_or_absolute_filtering = c("relative", "absolute", "nofiltering"), filtering_value, relative_or_absolute_plot = c("relative", "absolute"), plotting_tax_ranks = "all", output_path, figures_leg_path, distinct_colors = TRUE, horizontal_barplot = FALSE, facet_plot = FALSE, facetting_column = NULL, order_by_abundance = TRUE, separated_legend){

        # Transform abundance on 100%
          physeq_subset_norm_df <- melted_dataframe  %>% # calculate % normalized Abundance
              ungroup %>%
              dplyr::group_by(Sample) %>%
              dplyr::mutate(per=as.numeric(100*Abundance/sum(Abundance))) %>%
              ungroup() %>%
              dplyr::select(-Abundance)  %>%
              dplyr::rename(Abundance = per)

        # Choose between relative value or absolute value dataframe for plotting value
          ## Relative value plotting
              if (relative_or_absolute_plot == "relative"){
                  plotting <- "Relative"
                  print("Relative value plotting")
                  plotted_df <- physeq_subset_norm_df
              }
          ## Absolute value plotting
              else if (relative_or_absolute_plot == "absolute"){
                  plotting <- "Absolute"
                  print("Absolute value plotting")
                  plotted_df <- melted_dataframe
              }else{stop('"relative_value_plotting" must be "relative" or "absolute"')
            }
            ## Define the taxa ranks whill will be plotted
                    tax_ranks <- plotting_tax_ranks
                    print(tax_ranks)
                    t_column <- rlang::sym(tax_ranks) # set the taxa in format compatible with dplyr

            ### Filter values at a relative or absolute quantity
                #### No filtering
                    if (relative_or_absolute_filtering == "nofiltering"){
                        filtering <- "Nofiltering"
                        filtering_value <- "0"
                        print("No filtering based on quantity")
                        physeq_subset_df_filtered <- physeq_subset_norm_df
                        }
                #### Relative value filtering
                    else if (relative_or_absolute_filtering == "relative"){
                        filtering <- "Relative"
                        filtering_value <- filtering_value
                        print("Relative value based filtering")
                        physeq_subset_df_filtered <- physeq_subset_norm_df  %>%
                            dplyr::group_by(Sample, !!t_column) %>% # group the dataframe by Sample and taxa
                            dplyr::mutate(sumper=as.numeric(sum(Abundance))) %>% # calculate the cumulative relative abundance of the taxa in the sample
                            #ungroup %>%
                            #dplyr::group_by(Sample, !!t_column) %>% # group the dataframe by Sample and taxa
                            dplyr::filter(sumper >= filtering_value) %>%# keep only the taxa above threshold. From the applied grouping, the taxa remain in all amples if over the filtering value in one sample?
                            ungroup
                        }
                #### Absolute value filtering
                    else if (relative_or_absolute_filtering == "absolute"){
                        filtering <- "Absolute"
                        print("Absolute value based filtering")
                        physeq_subset_df_filtered <- melted_dataframe  %>%
                            dplyr::group_by(Sample, !!t_column) %>% # group the dataframe by Sample and taxa
                            dplyr::mutate(sumper=as.numeric(sum(Abundance))) %>% # calculate the cumulative relative abundance of the taxa in the sample
                            #ungroup %>%
                            #dplyr::group_by(Sample, !!t_column) %>% # group the dataframe by Sample and taxa
                            dplyr::filter(sumper >= filtering_value) %>%# keep only the taxa above threshold. From the applied grouping, the taxa remain in all amples if over the filtering value in one sample?
                            ungroup
                        }


            ### Now we have a dataframe ("plotted_df") containing all values, with Abundance expressed in absolute reads number or relative quanity in %, and a dataframe containing the rows that have to be kept while the rest will be merged by tagging them with the same name
                above_threshold_df <- semi_join(plotted_df, physeq_subset_df_filtered, by = c("OTU", "Sample")) ### In a dataframe, keep only the rows matching the filtered rowmanes

                under_threshold_df <- anti_join(plotted_df, physeq_subset_df_filtered, by = c("OTU", "Sample")) ### In another dataframe, keep only the lines NOT matching the filtered rowmanes


            ### Rename taxa with < percent/reads abundance, defending of the filtering applied
                #### Define the filtering tag depending of the applied filtering
                    ##### Relative
                    if (relative_or_absolute_filtering == "relative"){
                        filtering_tag <-paste("Relative abund. <", filtering_value , "%")
                    }
                    ##### Absolute filtering
                    else if (relative_or_absolute_filtering == "absolute"){
                        filtering_tag <-paste("Absolute abund. <", filtering_value, "reads")
                    }

                    ##### No filtering
                    else if (relative_or_absolute_filtering == "nofiltering"){
                        filtering_tag <-paste("nofiltering")
                    }


                    if (relative_or_absolute_filtering != "nofiltering"){
                #### Apply the filtering tag
                    under_threshold_df$Kingdom <-filtering_tag
                    under_threshold_df$Phylum <-filtering_tag
                    under_threshold_df$Class <-filtering_tag
                    under_threshold_df$Order <-filtering_tag
                    under_threshold_df$Family <-filtering_tag
                    under_threshold_df$Genus <-filtering_tag
                    under_threshold_df$Species <-filtering_tag
                    under_threshold_df$OTU <-filtering_tag
                }

            ### Join the two dataframes, to put back all rows together, the ones above threshold keeping their original taxonomic identifier while the others are now grouped together
                threshod_filtered_abs <- full_join(under_threshold_df,above_threshold_df)


            ### Remove the (now but needed before) Abundance = 0 rows
                threshod_filtered_abs_no_zero <- filter(threshod_filtered_abs, threshod_filtered_abs$Abundance>0)


            ### Reorder the facet factor if later used for plotting
                if (isTRUE(facet_plot)){
                    #threshod_filtered_abs_no_zero[[facetting_column]] <- as.factor(threshod_filtered_abs_no_zero[[facetting_column]])
                    threshod_filtered_abs_no_zero[[facetting_column]] <- fct_reorder(threshod_filtered_abs_no_zero[[facetting_column]], as.numeric(threshod_filtered_abs_no_zero[[x_axis_column]]))

                }else if (isFALSE(facet_plot)){print("No faceting")

                }else {stop('"facet_plot" must be TRUE or FALSE')
                }

           ### Reverse the x_axis_column column for later if using horizontal barplot
                if (isTRUE(horizontal_barplot)){
                    threshod_filtered_abs_no_zero[[x_axis_column]] <- fct_rev(threshod_filtered_abs_no_zero[[x_axis_column]])

                } else if (isFALSE(horizontal_barplot)){print("Vertical plotting")

                } else {stop('"horizontal_barplot" must be TRUE or FALSE')
                }


            ### Open pdf device
            pdf(file = output_path)

            ### Loop for unique value in grouping_column
                for (i in unique(threshod_filtered_abs_no_zero[[grouping_column]])) {
                    print(paste("Start plotting", grouping_column, i))

                #### filter the table for this value of the grouping columns. Depending of the used arguments, the t_neg_PCR values are kept or not on the barplots
                    ##### Keep t_neg_PCR rows
                        if (isTRUE(t_neg_PCR_sample_on_plots) & !is.null(t_neg_PCR_group_column_value)){
                            filtered_df_abs_i <- filter(threshod_filtered_abs_no_zero, threshod_filtered_abs_no_zero[[grouping_column]] == i | threshod_filtered_abs_no_zero[[grouping_column]] == t_neg_PCR_group_column_value)

                            print('Keeping t_neg_PCR values for the graphs. The "t_neg_PCR_group_column_value" must match the one in the grouping_column for this sample')
                        }
                        else if (isTRUE(t_neg_PCR_sample_on_plots) & is.null(t_neg_PCR_group_column_value)){
                        stop('If "t_neg_PCR_sample_on_plots" is "TRUE, a "t_neg_PCR_group_column_value" indicating the value of the T neg PCR sample in the grouping column must be indicated')
                        }
                    ##### Do not keep t_neg_PCR rows
                        else if (isFALSE(t_neg_PCR_sample_on_plots)){
                            filtered_df_abs_i <- filter(threshod_filtered_abs_no_zero, threshod_filtered_abs_no_zero[[grouping_column]] == i)
                        }

                        else{ stop('"t_neg_PCR_sample_on_plots" must be TRUE or FALSE')
                        }


        ### Merge what has been merged
        if (relative_or_absolute_filtering != "nofiltering"){

            ### Merged rows that were filtered in previous step so that they are only on one line on the heatmaps
            g_column <- rlang::sym(grouping_column)
            x_column<- rlang::sym(x_axis_column)
            f_column <- rlang::sym(facetting_column)
            filtered_OTU <- filtered_df_abs_i %>%
            dplyr::filter(grepl(filtering_tag, !! t_column)) %>%
            group_by(!!(x_column), !!(g_column), !!(f_column)) %>%
            summarise(Abundance = sum(Abundance)) %>%
            filter(Abundance > 0)

            ### Rename "filtered" those rows representing filtered rows
            filtered_OTU[[tax_ranks]] <- filtering_tag

            ### Keep only the NOT filtered rows
            individual_OTU <- filtered_df_abs_i %>%
            dplyr::filter(!grepl(filtering_tag, !! t_column))

            ### Bind the filtered and individual rows
            merged_filtered_OTU <- bind_rows(filtered_OTU, individual_OTU)

        }

        if (relative_or_absolute_filtering == "nofiltering"){
            merged_filtered_OTU <- filtered_df_abs_i
        }

                #### Reorder by abundance
                    if (isTRUE(order_by_abundance)){
                        print("Ordered by abundance")
                        merged_filtered_OTU <- merged_filtered_OTU %>%
                            group_by(!! t_column) %>%
                            mutate(tot = sum(Abundance)) %>%
                            ungroup() %>%
                            mutate(!! t_column := fct_reorder(!! t_column, tot)) %>%
                            arrange(desc(tot))
                    }
                    else if (isFALSE(order_by_abundance)){ print("NOT ordered by abundance")
                    }
                    else { stop('"order_by_abundance" must be TRUE or FALSE')
                    }



    #### Set the filtering_tabl at the top of the plot
        merged_filtered_OTU[[tax_ranks]] <- fct_relevel(merged_filtered_OTU[[tax_ranks]], filtering_tag, after = 0)


      ### Generate heatmap
      heatmap <- ggplot(merged_filtered_OTU, aes(x = get(x_axis_column), y = get(tax_ranks), fill = Abundance)) +
        theme_bw() +
        geom_tile(aes(fill=Abundance), show.legend=TRUE) +
        scale_fill_gradient(na.value = "white", low="#000033", high="#CCFF66") +
        geom_text(aes(label = round(Abundance, digits = 1), color = Abundance > 15), size = 2 ) +
        scale_color_manual(guide = FALSE, values = c("white", "black")) +
        theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
        #scale_y_discrete(labels = taxalabel) +
        labs(x="Sample",  y = paste(plotting, "abundance"), title = paste("Taxonomic composition", tax_ranks,"level" ))


      ### Turn heatmap horizontally
      if (isTRUE(horizontal_barplot)){
        ### Reverse the order of the samples
        heatmap <- heatmap + coord_flip()
      }

      ### Create facet view
      if (isTRUE(facet_plot)){
        heatmap <- heatmap + facet_grid(scales = "free", space = "free", cols = vars(get(facetting_column)))
      }

      ### Save it
      ### Set the filename
      #filename_base <- file.path(output_folder, paste(sep = "_", i, relative_or_absolute_filtering, filtering_value, plotting_tax_ranks))
      ### Print the filename to follow progress
      #print(filename_base)
      ###Save the figure
      #ggsave(heatmap, filename = paste0(filename_base, "_heatmap.png"), width = 4.5, height = 7)
      print(heatmap)


      }
    dev.off()
    }


################################################################################


## Run the function
    heatmap_fct(
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
      horizontal_barplot = horizontal_barplot,
      facet_plot = facet_plot,
      facetting_column = facetting_column,
      order_by = TRUE)
