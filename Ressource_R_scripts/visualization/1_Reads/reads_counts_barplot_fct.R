### Create the function
reads_barplots_fct <- function(count_table_df, figures_save_dir, grouping_column, filling_column, sample_label, distinct_colors = FALSE, facet_plot = FALSE, facetting_column = NULL, zoome){
  
  
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
  overall_reads_barplot <- ggplot(count_table_df, aes(x = get(sample_label), y = TotalReads, fill = get(filling_column))) +
    theme_bw() +
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
    facet_reads_barplot <- overall_reads_barplot + facet_grid(.~get(facetting_column), scales = "free", space = "free")
    
    facet_reads_barplot <<- facet_reads_barplot
    
    ### Save it
    ggsave(facet_reads_barplot, filename = paste(figures_save_dir, "/reads/barplot/overall_facet",facetting_column,"_", filling_column, "_", grouping_column,  "_reads_barplot.png",sep=""), width = 10, height = 5) 
  }
  
  
  
  ### Create a separate plot for each of the grouping column value
  
  
  ### Run a loop for each on the values in a column
  for (i in unique(reads_counts_df[[grouping_column]])) {
    
    filtered_count_table_df <- filter(count_table_df, count_table_df[[grouping_column]] == i )

      smax_g <- max(filtered_count_table_df$TotalReads)
      print(smax_g)
    
    ### Create the barplot
    reads_barplot <- ggplot(filtered_count_table_df, aes(x = get(sample_label), y = TotalReads, fill = get(filling_column))) +
      theme_bw() +
      geom_col() +
      scale_fill_manual(values = colors_palette) + 
      labs(x= grouping_column,  y ="Reads") +
      ggtitle(paste("Reads counts",grouping_column, i)) +
      scale_x_discrete(drop = TRUE) + # Keep all groups, included the ones with values. Alternative : (drop = FALSE)
      scale_y_continuous(labels = comma, limits = c(0,smax_g)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      guides(fill=guide_legend(title=filling_column))


    ### Create facet view
    if (isTRUE(facet_plot)){
      reads_barplot <- reads_barplot + facet_grid(.~get(facetting_column), scales = "free", space = "free")

    }
    
    ### Save it
    ggsave(reads_barplot, filename = paste(figures_save_dir, "/reads/barplot/", filling_column, "_", grouping_column, i, "_reads_barplot.png",sep=""), width = 10, height = 5) 
  }}
