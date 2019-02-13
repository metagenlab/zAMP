# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 11.02.19

## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_melted_table <- snakemake@input[["phyloseq_melted_table"]]
melted_dataframe <- read.table(file = phyloseq_melted_table, sep = "\t", header = TRUE)

## Ouput


output_folder <- dirname(snakemake@output[["barplot"]])

## Parameters
x_axis_column =snakemake@params[[ "x_axis_column"]]
grouping_column <- snakemake@params[[ "grouping_column"]]
t_neg_PCR_sample_on_plots = snakemake@params[["t_neg_PCR_sample_on_plots"]]
t_neg_PCR_group_column_value = snakemake@params[["t_neg_PCR_group_column_value"]]
relative_or_absolute_filtering = snakemake@params[["relative_or_absolute_filtering"]]
filtering_value = snakemake@params[["filtering_value"]]
relative_or_absolute_plot = snakemake@params[["relative_or_absolute_plot"]]
plotting_tax_ranks = snakemake@params[["plotting_tax_ranks"]]
distinct_colors = snakemake@params[["distinct_colors"]]
horizontal_barplot = snakemake@params[["horizontal_barplot"]]
facet_plot = snakemake@params[["facet_plot"]]
facetting_column = snakemake@params[["facetting_column"]]
order_by_abundance = snakemake@params[["order_by_abundance"]]
separated_legend = snakemake@params[["separated_legend"]]


## Load needed libraries
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(phyloseq); packageVersion("phyloseq")
library(RColorBrewer); packageVersion("RColorBrewer")
library(randomcoloR); packageVersion("randomcoloR")
library(data.table); packageVersion("data.table")
library(forcats); packageVersion("forcats")
library(rlang); packageVersion("rlang")
library(grid); packageVersion("grid")
library(cowplot); packageVersion("cowplot")

## Load the melted phyloseq table
melted_dataframe<- read.csv(file.path(phyloseq_melted_table), header = TRUE, sep = "\t")

## Order the x axis as in the metadata_table
#sample_data(phyloseq_obj)[[sample_type]] = factor(sample_data(phyloseq_obj)[[sample_type]], levels = unique(metadata[[sample_type]]), ordered = TRUE)
#sample_data(phyloseq_obj)[[x_axis_column]] = factor(sample_data(phyloseq_obj)[[x_axis_column]], levels = unique(metadata[[x_axis_column]]), ordered = TRUE)

################################################################################

### Create a function
    barplots_fct <- function(melted_dataframe, x_axis_column, grouping_column, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value, relative_or_absolute_filtering = c("relative", "absolute", "nofiltering"), filtering_value, relative_or_absolute_plot = c("relative", "absolute"), plotting_tax_ranks = "all", output_folder, figures_leg_path, distinct_colors = TRUE, horizontal_barplot = FALSE, facet_plot = FALSE, facetting_column = NULL, order_by_abundance = TRUE, separated_legend){

        # Import melted dataframe
          physeq_subset_df <- melted_dataframe

        # Transform abundance on 100%
          physeq_subset_norm_df <- physeq_subset_df  %>% # calculate % normalized Abundance
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
                  plotted_df <- physeq_subset_df
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
                        physeq_subset_df_filtered <- physeq_subset_df  %>%
                            dplyr::group_by(Sample, !!t_column) %>% # group the dataframe by Sample and taxa
                            dplyr::mutate(sumper=as.numeric(sum(Abundance))) %>% # calculate the cumulative relative abundance of the taxa in the sample
                            #ungroup %>%
                            #dplyr::group_by(Sample, !!t_column) %>% # group the dataframe by Sample and taxa
                            dplyr::filter(sumper >= filtering_value) %>%# keep only the taxa above threshold. From the applied grouping, the taxa remain in all amples if over the filtering value in one sample?
                            ungroup
                        }

                #### Write the table after table choice for relative of absolute value plotting , without Abundance = 0 rows to reduce its size.
                    #plotted_df %>%
                        #filter(Abundance>0) %>%
                        #write.table(file = paste0(figures_save_dir,"/quantitative_barplots/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_", "u", "_all_", x_axis_column, "_plotted_table.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)


                #### Write the table with values that passed the quantity filtering.
                    #physeq_subset_df_filtered %>%  ### To follow and control the process, write external .tsv tables
                        #filter(Abundance>0) %>%
                        #write.table(file = paste0(figures_save_dir,"/quantitative_barplots/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", filtering_value, "_all_", x_axis_column,"_",t ,"_filtration_table.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)


            ### Now we have a dataframe ("plotted_df") containing all values, with Abundance expressed in absolute reads number or relative quanity in %, and a dataframe containing the rows that have to be kept while the rest will be merged by tagging them with the same name
                above_threshold_df <- semi_join(plotted_df, physeq_subset_df_filtered, by = c("OTU", "Sample")) ### In a dataframe, keep only the rows matching the filtered rowmanes

                under_threshold_df <- anti_join(plotted_df, physeq_subset_df_filtered, by = c("OTU", "Sample")) ### In another dataframe, keep only the lines NOT matching the filtered rowmanes

            ### In another dataframe, keep the join of the under and above tables, before masking of the filtered taxa and without Abundance = 0 rows to reduce its size.
                #merged_filtered_abs <- full_join(under_threshold_df,above_threshold_df) %>%
                #    filter(Abundance>0)
                #write.table(merged_filtered_abs, file = paste0(figures_save_dir,"/quantitative_barplots/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", filtering_value, "_all_", x_axis_column, "_merged_without_masking.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)

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

                print("test5")
            ### Set colors palette
                #### Brewer colors
                    if (distinct_colors == FALSE) {
                        getPalette <- colorRampPalette(brewer.pal(n=9, "Set1"))
                        ColList <- unique(threshod_filtered_abs_no_zero[[tax_ranks]])
                        ColPalette = getPalette(length(ColList))
                        names(ColPalette) = ColList
                        colors_palette <- ColPalette
                        colors_palette[filtering_tag] <- "#d3d3d3" # Set the filtered in balck
                    }
                #### Distinct colors
                    else if (distinct_colors == TRUE) {
                        set.seed(4)
                        ColList <- unique(threshod_filtered_abs_no_zero[[tax_ranks]])
                        print(length(unique(threshod_filtered_abs_no_zero[[tax_ranks]])))
                        ColPalette <- distinctColorPalette(altCol = TRUE, k = length(unique(threshod_filtered_abs_no_zero[[tax_ranks]])))
                        names(ColPalette) = ColList
                        colors_palette <- ColPalette
                        colors_palette[filtering_tag] <- "#d3d3d3" # Set the filtered in balck

                    }else{ print("distinct_colors must be TRUE or FALSE")
                    }

                print("test6")

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
                    filtered_df_abs_i[[tax_ranks]] <- fct_relevel(filtered_df_abs_i[[tax_ranks]], filtering_tag, after = 0)


                #### Write the content of the plot table in a external file
                    #write.table(filtered_df_abs_i, file = paste0(figures_save_dir,"/quantitative_barplots/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", filtering_value, "_",grouping_column, "_", i, "_", tax_ranks,"_",x_axis_column, "_abundancy_table.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)


                #### Renames the values of the vector used for labeling
                    #x_labels <- as(filtered_df_abs_i[[x_axis_column]], "character")
                    #names(x_labels) <- filtered_df_abs_i[["Sample"]]


                #### Create the barplot
                    taxrank_barplot <- filtered_df_abs_i %>%
                        ggplot(aes(x = get(x_axis_column), y = Abundance, fill = get(tax_ranks))) +
                        theme_bw() +
                        geom_col() +
                        #scale_x_discrete(labels = x_labels, drop = TRUE) + # Set false to keep empty bars
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), plot.title = element_text(hjust = 0.5)) + # axis and title settings
                        guides(fill = guide_legend(title = paste0(tax_ranks),reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) + # settings of the legend
                        labs(x="Sample",  y = paste(plotting, "abundance"), title = paste("Taxonomic composition", tax_ranks,"level" )) + # axis and graph title
                        scale_fill_manual(values = colors_palette) # colors as set previously


                #### In option, turn barplot horizontally
                    if (isTRUE(horizontal_barplot)){
                        taxrank_barplot <- taxrank_barplot + coord_flip() ### Reverse the order of the samples

                    #### In option, create facet view
                        if (isTRUE(facet_plot)){
                            taxrank_barplot <- taxrank_barplot + facet_grid(get(facetting_column) ~., scales = "free", space = "free")
                        }
                    } else {
                        #### In option, create facet view
                        if (isTRUE(facet_plot)){
                            taxrank_barplot <- taxrank_barplot + facet_grid(~ get(facetting_column), scales = "free", space = "free")
                        }
                    }


                #### In option, keep the barplot without the legend
                    if (isTRUE(separated_legend)){
                        taxrank_barplot_no_leg <- (taxrank_barplot + guides(fill = FALSE))
                    }else {
                        taxrank_barplot_no_leg <- taxrank_barplot
                    }


                    #### Save it
                        ##### Set the filename
                            filename_base <- file.path(output_folder, paste(sep = "_", i, relative_or_absolute_filtering, filtering_value, plotting_tax_ranks))
                        ##### Finally, save the figure
                          ggsave(taxrank_barplot_no_leg, filename = paste0(filename_base, "_barplot.png"), width = 7, height = 7)
                        #### Extract the legend and save it
                        leg <- get_legend(taxrank_barplot)
                        ggsave(leg, filename = paste0(filename_base, "_barplot_leg.png"), width = 5, height = 10)


        }}
################################################################################

#save.image(file = file.path(output_folder, "rdebug.RData"))


## Run the function
barplots_fct(
  melted_dataframe = melted_dataframe,
    grouping_column = grouping_column,
    x_axis_column = x_axis_column,
    t_neg_PCR_sample_on_plots = t_neg_PCR_sample_on_plots,
    t_neg_PCR_group_column_value = t_neg_PCR_group_column_value,
    relative_or_absolute_filtering = relative_or_absolute_filtering,
    filtering_value = filtering_value,
    relative_or_absolute_plot = relative_or_absolute_plot,
    plotting_tax_ranks = plotting_tax_ranks,
    output_folder = output_folder,
    distinct_colors = distinct_colors,
    horizontal_barplot = horizontal_barplot,
    facet_plot = facet_plot,
    facetting_column = facetting_column,
    order_by_abundance = order_by_abundance,
    separated_legend = separated_legend)
