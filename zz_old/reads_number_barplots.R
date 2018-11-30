# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 12.11.18




reads_barplots_fct <- function(count_table_df, figures_save_dir, grouping_column, filling_column, x_axis_column, distinct_colors = FALSE, facet_plot = FALSE, facet_value = NULL){

      ### Record data on the distribution of number of reads (useful later to scale plots axis)
      smin <- min(count_table_df$TotalReads)
      print(smin)
      smean <- mean(count_table_df$TotalReads)
      print(smean)
      smax <- max(count_table_df$TotalReads)
      print(smax)



      ### Create an overall barplot

            ### Set colors

                  #### BrewerColors
              if (distinct_colors == FALSE) {
                  getPalette = colorRampPalette(brewer.pal(n=11, "Paired"))
                  ColList = unique(count_table_df[[filling_column]])
                  ColPalette = getPalette(length(ColList))
                  names(ColPalette) = ColList
                  colors_palette <- ColPalette
              }
              if (distinct_colors == TRUE) {

                ### Distinct colors
                    set.seed(2)
                    ColList <- unique(count_table_df[[filling_column]])
                    ColPalette <- distinctColorPalette(altCol = TRUE, k = length(unique(count_table_df[[filling_column]])))
                    names(ColPalette) = ColList
                    colors_palette <- ColPalette
                }

          ### Create the barplot
          overall_reads_barplot <- ggplot(count_table_df, aes(x = get(x_axis_column), y = TotalReads, fill = get(filling_column))) +
          geom_col() +
          scale_fill_manual(values = colors_palette) +
          labs(x= grouping_column,  y ="Reads") +
          ggtitle(paste("Reads counts overall")) +
          scale_x_discrete(drop = TRUE) + # Keep all groups, included the ones with values. Alternative : (drop = FALSE)
          scale_y_continuous(labels = comma, limits = c(0,smax)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          guides(fill=guide_legend(title=filling_column))

          ### Save it
          ggsave(overall_reads_barplot, filename = paste(figures_save_dir, "/reads/barplot/overall", filling_column, "_", grouping_column,  "_reads_barplot.png",sep=""), width = 10, height = 5)

          ### Create facet view
          if (isTRUE(facet_plot)){
          facet_reads_barplot <- overall_reads_barplot +  facet_grid(.~ get(facet_value) , scales = "free", space = "free")


          ### Save it
          ggsave(facet_reads_barplot, filename = paste(figures_save_dir, "/reads/barplot/overall_facet",facet_value,"_", filling_column, "_", grouping_column,  "_reads_barplot.png",sep=""), width = 10, height = 5)
          }


      ### Create a separate plot for each of the grouping column value
      ### Run a loop for each on the values in a column
      for (i in unique(reads_counts_df[[grouping_column]])) {

        filtered_count_table_df <- filter(count_table_df, count_table_df[[grouping_column]] == i )

      ### Create the barplot
      reads_barplot <- ggplot(filtered_count_table_df, aes(x = get(x_axis_column), y = TotalReads, fill = get(filling_column))) +
      geom_col() +
      scale_fill_manual(values = colors_palette) +
      labs(x= grouping_column,  y ="Reads") +
      ggtitle(paste("Reads counts",grouping_column, i)) +
      scale_x_discrete(drop = TRUE) + # Keep all groups, included the ones with values. Alternative : (drop = FALSE)
      scale_y_continuous(labels = comma, limits = c(0,smax)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      guides(fill = FALSE)
      ### Save it
      ggsave(reads_barplot, filename = paste(figures_save_dir, "/reads/barplot/", filling_column, "_", grouping_column, i, "_reads_barplot.png",sep=""), width = 7, height = 5)
        }}

  ### Run the function in a loop to have it for each patientID
 reads_barplots_fct(count_table_df = reads_counts_df, figures_save_dir = ".", grouping_column = "patientID", filling_column = "sample_source", x_axis_column = "sample_source", distinct_colors = FALSE, facet_plot = FALSE, facet_value = "patientID")

