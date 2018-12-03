
### Create a function
barplots_fct <- function(melted_dataframe, x_axis_column, grouping_column, grouping_column_filtering = c(FALSE, TRUE), grouping_column_filtering_value, t_neg_PCR_sample_on_plots, t_neg_PCR_sample_grp_column_value, taxonomic_filtering = c(TRUE, FALSE), taxonomic_filtering_rank = "Kingdom" , taxonomic_filtering_value = "Bacteria" ,  quantity_filtering_type = c("relative", "absolute", "rank", "nofiltering", "absolute_and_rank"), absolute_quantity_filtering_value, relative_quantity_filtering_value, rank_quantity_filtering_value, plotting_value = c("relative", "absolute"), plotting_tax_ranks = "all", figures_save_dir, distinct_colors = TRUE, horizontal_barplot = FALSE, facet_plot = FALSE, facetting_column = NULL, order_by_abundance = TRUE, separated_legend){

# Transform values in  dplyr format
    g_column <- rlang::sym(grouping_column)
    x_column<- rlang::sym(x_axis_column)


# In option, filter for a given grouping_columns.
    ## No filter for the grouping column
        if (isFALSE(grouping_column_filtering)){
            melted_dataframe_filtered <- melted_dataframe
            print("No filtering of the grouping_column")
        }

    ## Filter for the grouping column
        else if (isTRUE(grouping_column_filtering)){
            melted_dataframe_filtered <- melted_dataframe %>% filter(get(grouping_column) == grouping_column_filtering_value)
            print(paste("Data filtered for", grouping_column, "is", grouping_column_filtering_value))

        }else{ stop("grouping_column_filtering must be TRUE or FALSE")
        }


# In option, filter at a given taxonomic rank and for a given taxonomic ID
    ## No filter
        if (isFALSE(taxonomic_filtering)) {
            physeq_subset_df <- melted_dataframe_filtered

            physeq_subset_norm_df <- physeq_subset_df  %>% # calculate % normalized Abundance
                ungroup %>%
                dplyr::group_by(Sample) %>%
                dplyr::mutate(per=as.numeric(100*Abundance/sum(Abundance))) %>%
                ungroup() %>%
                dplyr::select(-Abundance)  %>%
                dplyr::rename(Abundance = per)

            taxonomic_filtering_value <- "no_tax_rank_filtering"
            print("No filtering of any taxonomic rank")
        }

    ## A specific filter
        else if (isTRUE(taxonomic_filtering)){
            physeq_subset_df <- melted_dataframe_filtered %>%
                dplyr::filter(get(taxonomic_filtering_rank) == taxonomic_filtering_value)

            physeq_subset_norm_df <- physeq_subset_df  %>%
                dplyr::group_by(Sample) %>%
                dplyr::mutate(per=paste0((100*Abundance/sum(Abundance)))) %>%
                ungroup %>%
                dplyr::select(-Abundance)  %>%
                dplyr::rename(Abundance = per)

            physeq_subset_norm_df$Abundance <- as.numeric(physeq_subset_norm_df$Abundance)
            print(paste("Data filtered for", taxonomic_filtering_rank, "is", taxonomic_filtering_value))

        }else{stop('"taxonomic_filtering" must be TRUE or FALSE')
        }

# Choose between relative value or absolute value dataframe for plotting value
    ## Relative value plotting
        if (plotting_value == "relative"){
            plotting <- "Relative"
            print("Relative value plotting")
            plotted_df <- physeq_subset_norm_df
        }
    ## Absolute value plotting
        else if (plotting_value == "absolute"){
            plotting <- "Absolute"
            print("Absolute value plotting")
            plotted_df <- physeq_subset_df
        }else{ stop('"relative_value_plotting" must be "relative" or "absolute"')
        }

# From here on, loop for all taxonomic levels chosen for plotting, this is what changes with the alternative plotting script where this is done at the end.
    ## Define the taxa ranks whill will be plotted
        if (plotting_tax_ranks == "all"){
            tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")
        }else{
            tax_ranks <- plotting_tax_ranks
            print('Plotting at a specific taxonomic rank. Check spelling one or several of ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"). This error does NOT always mean that there was an error in spelling')
        }

    ## Run the loop
        for (t in tax_ranks){
        t_column <- rlang::sym(t) # set the taxa in format compatible with dplyr
        print(t_column)

        ### Filter values at a relative or absolute quantity
            #### No filtering
                if (quantity_filtering_type == "nofiltering"){
                    filtering <- "Nofiltering"
                    quantity_filtering_value <- "0"
                    print("No filtering based on quantity")
                    physeq_subset_df_filtered <- physeq_subset_norm_df
                }

            #### Relative value filtering
                else if (quantity_filtering_type == "relative"){
                    filtering <- "Relative"
                    quantity_filtering_value <- relative_quantity_filtering_value
                    print("Relative value based filtering")
                    physeq_subset_df_filtered <- physeq_subset_norm_df  %>%
                        dplyr::group_by(!!x_column, !!t_column) %>% # group the rows by grouping column and taxa
                        dplyr::mutate(sumper=as.numeric(sum(Abundance))) %>% # calculate the cumulative relative abundance of the taxa in the sample
                        ungroup %>%
                        dplyr::group_by(!! x_column)  %>% # group by sample
                        dplyr::mutate(totper = as.numeric(paste0(sum(Abundance)))) %>% # calculate the total abundance in the group
                        ungroup %>%
                        dplyr::mutate(ponper = as.numeric(100*sumper/totper)) %>% # divide the sum of the taxa by the sum of the group
                        dplyr::filter(ponper >= quantity_filtering_value) # keep only the taxa above threshold in the group
                }

            #### Absolute value filtering
                else if (quantity_filtering_type == "absolute"){
                    filtering <- "Absolute"
                    quantity_filtering_value <- absolute_quantity_filtering_value
                    print("Absolute value based filtering")
                    physeq_subset_df_filtered <- physeq_subset_df  %>%
                        dplyr::group_by(!! g_column, !!t_column) %>%  # group the dataframe by the grouping column
                        dplyr::mutate(sumper=as.numeric(sum(Abundance))) %>% ## calculate a absolute abundance of the taxa in the group
                        dplyr::filter(sumper >= quantity_filtering_value) %>% # keep only the taxa above threshold
                        ungroup
                }

            #### Write the table after table choice for relative of absolute value plotting , without Abundance = 0 rows to reduce its size.
                plotted_df %>%
                    filter(Abundance>0) %>%
                    write.table(file = paste0(figures_save_dir,"/quantitative_barplots/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_", "u", "_all_", x_axis_column, "_plotted_table.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)

            #### Write the table with values that passed the quantity filtering.
                physeq_subset_df_filtered %>%  ### To follow and control the process, write external .tsv tables
                    filter(Abundance>0) %>%
                    write.table(file = paste0(figures_save_dir,"/quantitative_barplots/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", quantity_filtering_value, "_all_", x_axis_column,"_",t ,"_filtration_table.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)

        ### Now we have a dataframe ("plotted_df") containing all values, with Abundance expressed in absolute reads number or relative quanity in %, and a dataframe containing the rows that have to be kept while the rest will be merged by tagging them with the same name
            above_threshold_df <- semi_join(plotted_df, physeq_subset_df_filtered, by = c("OTU", "Sample")) ### In a dataframe, keep only the rows matching the filtered rowmanes
            under_threshold_df <- anti_join(plotted_df, physeq_subset_df_filtered, by =c ("OTU", "Sample")) ### In another dataframe, keep only the lines NOT matching the filtered rowmanes

        ### In another dataframe, keep the join of the under and above tables, before masking of the filtered taxa and without Abundance = 0 rows to reduce its size.
            merged_filtered_abs <- full_join(under_threshold_df,above_threshold_df) %>%
                filter(Abundance>0)
            write.table(merged_filtered_abs, file = paste0(figures_save_dir,"/quantitative_barplots/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", quantity_filtering_value, "_all_", x_axis_column, "_merged_without_masking.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)


        ### Rename taxa with < percent/reads abundance, defending of the filtering applied
            #### Define the filtering tag depending of the applied filtering
                ##### Relative
                if (quantity_filtering_type == "relative"){
                    filtering_tag <-paste("Relative cumulated abund. <", quantity_filtering_value , "%")
                }
                ##### Absolute filtering
                else if (quantity_filtering_type == "absolute"){
                    filtering_tag <-paste("Absolute cumulated abund. <", quantity_filtering_value, "reads")
                }

            #### Apply the filtering tag
                under_threshold_df$Kingdom <-filtering_tag
                under_threshold_df$Phylum <-filtering_tag
                under_threshold_df$Class <-filtering_tag
                under_threshold_df$Order <-filtering_tag
                under_threshold_df$Family <-filtering_tag
                under_threshold_df$Genus <-filtering_tag
                under_threshold_df$Species <-filtering_tag
                under_threshold_df$OTU <-filtering_tag
  
  
        ### Join the two dataframes, to put back all rows together, the ones above threshold keeping their original taxonomic identifier while the others are now grouped together
            threshod_filtered_abs <- full_join(under_threshold_df,above_threshold_df)

        ### Remove the (now but needed before) Abundance = 0 rows
            threshod_filtered_abs_no_zero <- filter(threshod_filtered_abs, threshod_filtered_abs$Abundance>0)
  

        ### Reorder the facet factor if later used for plotting
            if (isTRUE(facet_plot)){
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
  
  
        ### Set colors palette
            #### Brewer colors
                if (distinct_colors == FALSE) {
                    getPalette <- colorRampPalette(brewer.pal(n=9, "Set1"))
                    ColList <- unique(threshod_filtered_abs_no_zero[[t]])
                    ColPalette = getPalette(length(ColList))
                    names(ColPalette) = ColList
                    colors_palette <- ColPalette
                    colors_palette[filtering_tag] <- "#d3d3d3" # Set the filtered in balck
                }

            #### Distinct colors
                else if (distinct_colors == TRUE) {
                    set.seed(2)
                    ColList <- unique(threshod_filtered_abs_no_zero[[t]])
                    ColPalette <- distinctColorPalette(altCol = TRUE, k = length(unique(threshod_filtered_abs_no_zero[[t]])))
                    names(ColPalette) = ColList
                    colors_palette <- ColPalette
                    colors_palette[filtering_tag] <- "#d3d3d3" # Set the filtered in balck

                }else{ print("distinct_colors must be TRUE or FALSE")
                }

    
        ### Loop for unique value in grouping_column
            for (i in unique(threshod_filtered_abs_no_zero[[grouping_column]])) {
                print(paste("Start plotting", grouping_column, i))

            #### filter the table for this value of the grouping columns. Depending of the used arguments, the t_neg_PCR values are kept or not on the barplots
                ##### Keep t_neg_PCR rows
                    if (isTRUE(t_neg_PCR_sample_on_plots) & !is_null(t_neg_PCR_sample_grp_column_value)){
                        filtered_df_abs_i <- filter(threshod_filtered_abs_no_zero, threshod_filtered_abs_no_zero[[grouping_column]] == i | threshod_filtered_abs_no_zero[[grouping_column]] == t_neg_PCR_sample_grp_column_value)

                        print('Keeping t_neg_PCR values for the graphs. The "t_neg_PCR_sample_grp_column_value" must match the one in the grouping_column for this sample')
                    }
                    else if (isTRUE(t_neg_PCR_sample_on_plots) & is_null(t_neg_PCR_sample_grp_column_value)){
                    stop('If "t_neg_PCR_sample_on_plots" is "TRUE, a "t_neg_PCR_sample_grp_column_value" indicating the value of the T neg PCR sample in the grouping column must be indicated')
                    }

                ##### Do not keep t_neg_PCR rows
                    else if (isFALSE(t_neg_PCR_sample_on_plots)){
                        filtered_df_abs_i <- filter(threshod_filtered_abs_no_zero, threshod_filtered_abs_no_zero[[grouping_column]] == i)
                    }

                    else{ stop('"t_neg_PCR_sample_on_plots" must be TRUE or FALSE')
                    }
      
            #### Reorder by abundance
                if (isTRUE(order_by_abundance)){
                    print("Ordered by abundance")
                    filtered_df_abs_i <- filtered_df_abs_i %>%
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
                filtered_df_abs_i[[t]] <- fct_relevel(filtered_df_abs_i[[t]], filtering_tag, after = 0)


            #### Write the content of the plot table in a external file
                write.table(filtered_df_abs_i, file = paste0(figures_save_dir,"/quantitative_barplots/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", quantity_filtering_value, "_",grouping_column, "_", i, "_", t,"_",x_axis_column, "_abundancy_table.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)


            #### Renames the values of the vector used for labeling
                x_labels <- as(filtered_df_abs_i[[x_axis_column]], "character")
                names(x_labels) <- filtered_df_abs_i[["Sample"]]


            #### Create the barplot
                taxrank_barplot <- filtered_df_abs_i %>%
                    ggplot(aes(x = Sample, y = Abundance, fill = get(t))) +
                    theme_bw() +
                    geom_col() +
                    scale_x_discrete(labels = x_labels, drop = TRUE) + # Set false to keep empty bars
                    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 1)) + # axis and title settings
                    guides(fill = guide_legend(title = paste(grouping_column, i, sep= " "),reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) + # settings of the legend
                    labs(x=x_axis_column,  y = paste(plotting, "abundance"), title = paste(t, "composition",grouping_column , i, sep= " ")) + # axis and graph title
                    scale_fill_manual(values = colors_palette) # colors as set previously


            #### In option, turn barplot horizontally
                if (isTRUE(horizontal_barplot)){
                    taxrank_barplot <- taxrank_barplot + coord_flip() ### Reverse the order of the samples
                }

            #### In option, create facet view
                if (isTRUE(facet_plot)){
                    taxrank_barplot <- taxrank_barplot + facet_grid(get(facetting_column) ~., scales = "free", space = "free")
                }
            #### In option, keep the barplot without the legend
                if (isTRUE(separated_legend)){
                    taxrank_barplot_no_leg <- (taxrank_barplot + guides(fill = FALSE))
                }else {
                    taxrank_barplot_no_leg <- taxrank_barplot
                }

            #### Save it
                ##### Set the filename
                    filename_base <- (paste0(figures_save_dir,"/quantitative_barplots/", plotting, "/",filtering,"/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", quantity_filtering_value, "_",grouping_column, "_", i, "_", x_axis_column, "_",facetting_column,"_" ,t))
                    print(filename_base) #Print the filename to follow progress
                ##### Finally, save the figure
                  ggsave(taxrank_barplot_no_leg, filename = paste0(filename_base, "_barplot.png"), width = 7, height = 7)

                #### Extract the legend and save it
                leg <- get_legend(taxrank_barplot)
                ggsave(leg, filename = paste0(filename_base, "_barplot_leg.png"), width = 5, height = 10)

    # The loop restart to create blots for all taxa
      
    }}}
