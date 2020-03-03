
## qPCR reads plots

https://gist.github.com/tomhopper/faa24797bb44addeba79


```{r}

### Create the function
qPCR_dna_quant_barplots_fct <- function(count_table_df, figures_save_dir, grouping_column, filling_column, sample_label, distinct_colors = FALSE, facet_plot = FALSE, facetting_column = NULL){
  
  
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
    ColPalette <- distinctColorPalette(altCol = FALSE, k = length(unique(count_table_df[[filling_column]])))
    names(ColPalette) = ColList
    colors_palette <- ColPalette
  }
  
  
  ### Create the barplot for qPCR actin
  qPCR_16s_barplot <- ggplot(count_table_df, aes(x = get(sample_label), y = Quantity_actin, fill = get(filling_column))) +
    geom_col() +
    scale_fill_manual(values = colors_palette) + 
    labs(x= grouping_column,  y ="Quantification actin qPCR [copies/ml]") +
    ggtitle(paste("Actin quantity qPCR")) +
    scale_x_discrete(drop = TRUE) + # Keep all groups, included the ones with values. Alternative : (drop = FALSE)  
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank()) +
    guides(fill=guide_legend(title=filling_column)) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  
  ### Create the barplot for qPCR 16S
  qPCR_actin_barplot <- ggplot(count_table_df, aes(x = get(sample_label), y = Quantity_16S, fill = get(filling_column))) +
    geom_col() +
    scale_fill_manual(values = colors_palette) + 
    labs(x= grouping_column,  y ="Quantification 16S qPCR [copies/ml]") +
    ggtitle(paste("16S quantity qPCR")) +
    scale_x_discrete(drop = TRUE) + # Keep all groups, included the ones with values. Alternative : (drop = FALSE)  
    scale_y_continuous() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank()) +
    guides(fill=guide_legend(title=filling_column)) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  
  ## Create the barplot for bact to actin ratio
  bact_to_actin_ratio_barplot <- ggplot(count_table_df, aes(x = get(sample_label), y = bact_to_human_ratio, fill = get(filling_column))) +
    geom_col() +
    scale_fill_manual(values = colors_palette) + 
    labs(x= grouping_column,  y ="Quantitative DNA ratio between bacteria and actin") +
    ggtitle(paste("Quantitative ratio between bacterial and human DNA")) +
    scale_x_discrete(drop = TRUE) + # Keep all groups, included the ones with values. Alternative : (drop = FALSE)  
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank()) +
    guides(fill=guide_legend(title=filling_column)) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  
  ### Create the barplot for reads
  overall_reads_barplot <- ggplot(count_table_df, aes(x = get(sample_label), y = TotalReads, fill = get(filling_column))) +
    geom_col() +
    scale_fill_manual(values = colors_palette) + 
    labs(x= sample_label,  y ="Reads") +
    ggtitle(paste("Reads counts overall")) +
    scale_x_discrete(drop = TRUE) + # Keep all groups, included the ones with values. Alternative : (drop = FALSE)  
    scale_y_continuous(labels = comma, limits = c(0,smax)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank()) +
    guides(fill=guide_legend(title=filling_column))
  
  png(paste(figures_save_dir, "/reads/barplot/DNA_Quant_", filling_column, "_", grouping_column,  "_reads_barplot2.png",sep=""), width = 800, height = 800, units = "px") 
  grid.draw(rbind(ggplotGrob(qPCR_16s_barplot), ggplotGrob(qPCR_actin_barplot), ggplotGrob(bact_to_actin_ratio_barplot), ggplotGrob(overall_reads_barplot), size = "last"))
  dev.off()
  
  
}   

dna_quant_barplots_fct(count_table_df = reads_counts_df, figures_save_dir = "r_figures", grouping_column = "extraction_kit", filling_column = "sample_source", sample_label = "sample_label", distinct_colors = TRUE, facet_plot = FALSE)




```



