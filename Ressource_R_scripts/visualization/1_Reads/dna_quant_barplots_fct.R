### Create the function
dna_quant_barplots_fct <- function(count_table_df, figures_save_dir, grouping_column, filling_column, x_axis_column, distinct_colors = FALSE, facet_plot = FALSE, facetting_column = NULL, dna_quant_column){
  
  
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
  
  ### Create the barplot for DNA
  overall_dna_barplot <- ggplot(count_table_df, aes(x = get(x_axis_column), y = get(dna_quant_column), fill = get(filling_column))) + 
    geom_col() +
    scale_fill_manual(values = colors_palette) + 
    labs(x= grouping_column,  y ="Quantification [nM]") +
    ggtitle(paste("DNA quantity")) +
    scale_x_discrete(drop = TRUE) + # Keep all groups, included the ones with values. Alternative : (drop = FALSE)  
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank()) +
    guides(fill=guide_legend(title=filling_column)) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  ### Create facet view
  if (isTRUE(facet_plot)){
    overall_dna_barplot <- overall_dna_barplot + facet_grid(.~get(facetting_column), scales = "free")
    
  }
  
  
  
  
  ### Zoom on a on lower y
  zoomed_overall_dna_barplot <- overall_dna_barplot +  coord_cartesian(ylim=c(0,50)) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggtitle(paste("DNA quantity - low values")) 
  
  ### Create facet view
  if (isTRUE(facet_plot)){
    zoomed_overall_dna_barplot <- zoomed_overall_dna_barplot + facet_grid(.~get(facetting_column), scales = "free")

  }
  
  
  
  ### Create the barplot for reads
  overall_reads_barplot <- ggplot(count_table_df, aes(x = get(x_axis_column), y = TotalReads, fill = get(filling_column))) + 
    geom_col() +
    scale_fill_manual(values = colors_palette) + 
    labs(x= grouping_column,  y ="Reads") +
    ggtitle(paste("Reads counts overall")) +
    scale_x_discrete(drop = TRUE) + # Keep all groups, included the ones with values. Alternative : (drop = FALSE)  
    scale_y_continuous(labels = comma, limits = c(0,smax)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank()) +
    guides(fill=guide_legend(title=filling_column)) 
  
  
  
  
  ### Create facet view
  if (isTRUE(facet_plot)){
    overall_reads_barplot <- overall_reads_barplot + facet_grid(.~get(facetting_column), scales = "free")
    
    
    
    ### Zoom on a on lower y
    zoomed_overall_reads_barplot <- overall_reads_barplot +  coord_cartesian(ylim=c(0,200000)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      ggtitle(paste("Reads counts overal - low values")) 

    
    ### Create facet view
    if (isTRUE(facet_plot)){
      zoomed_overall_reads_barplot <- zoomed_overall_reads_barplot + facet_grid(.~get(facetting_column), scales = "free")
    }
    
 
    
  }
  
  png(paste(figures_save_dir, "/reads/barplot/DNA_Quant_", filling_column, "_", grouping_column,  "_reads_barplot.png",sep=""), width = 800, height = 900, units = "px") 
  grid.draw(rbind( ggplotGrob(zoomed_overall_dna_barplot), ggplotGrob(overall_dna_barplot), ggplotGrob(zoomed_overall_reads_barplot), ggplotGrob(overall_reads_barplot),  size = "last"))
  dev.off()
  
  
}   