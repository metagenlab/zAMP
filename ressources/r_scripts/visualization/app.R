# Load packages
## Shiny packages
library(shiny)
library(shinythemes)
library(colourpicker)

## Specific packages
library(tidyr); packageVersion("tidyr")
library(dplyr); packageVersion("dplyr")
library(pheatmap); packageVersion("pheatmap")
library(ggplot2); packageVersion("ggplot2")
library(RColorBrewer); packageVersion("RColorBrewer")
library(randomcoloR); packageVersion("randomcoloR")
library(data.table); packageVersion("data.table")
library(forcats); packageVersion("forcats")
library(rlang); packageVersion("rlang")
library(ggpubr); packageVersion("ggpubr")
library(cowplot); packageVersion("cowplot")
library(vegan); packageVersion("vegan")
library(readxl); packageVersion("readxl")

# Create sets of functions
## Filter metadata to exclude samples
filter_metadata_fct <- function(metadata_table, filter_out, filtering_column, filter_out_value){
  

  if (isTRUE(filter_out)){
    metadata_table_filtered <- metadata_table %>% dplyr::filter(!(!! rlang::sym(filtering_column) %in% filter_out_value))
  } else {
    metadata_table_filtered <- metadata_table
  }

  metadata_table_filtered <- droplevels(metadata_table_filtered)
  
  return(metadata_table_filtered)
  
}

## Add colors to metadata for alpha diversity and NMDS
add_colors_fct <- function(metadata_table, alp_NMDS_color, alp_NMDS_palette){
  
  #### Add Brewer colors for each value in the color columns to prepare for alpha-diversity and NMDS plots
  getPalette <- colorRampPalette(brewer.pal(n=12, alp_NMDS_palette)) # Create a function to generate a palette of colors
  colors_palette <- getPalette(length(unique(metadata_table[[alp_NMDS_color]]))) # Call this function for each unique value of the column of interest
  names(colors_palette) <- unique(metadata_table[[alp_NMDS_color]]) # Name colors by the value in the column
  m <- match(metadata_table[[alp_NMDS_color]], names(colors_palette)) # match the names of the colors with the values 
  metadata_table$NMDS_alpha_colors <- colors_palette[m] # Fill the table with the colors for each value
  
  return(metadata_table)
  
}

## Prepare a melted long table from OTU, taxonomy and metadata table
filter_OTU_count_fct <- function(count_table, metadata_table, taxonomy_table, grouping_column, x_axis_column, sort_column, facetting_column, relative_or_absolute_filtering = c("relative", "absolute", "nofiltering"), filter_value, plotting_tax_ranks, distinct_colors, abundance_filtre_level = c("Sample", "Group")){
  

  ## Reformat and join count table, metadata and taxonomy
  ### Transform count table to long format
  count_table[["ASV"]] <- row.names(count_table) 
  count_table_long <- tidyr::gather(count_table, key = "Sample", value = "Count", -one_of("ASV"))
  
  ### Reformat taxonomy
  taxonomy_table_split <- taxonomy_table %>% 
    tidyr::separate(V2, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence"))
  colnames(taxonomy_table_split)[1] <- "ASV"
  
  ### Add taxonomy to long count table
  count_table_long_tax <- dplyr::left_join(count_table_long, taxonomy_table_split, by = "ASV")
  
  ## Filter taxa based on quantities
  ### Set a list with all ranks from Kingdom to the plotted rank. Used to then be able to concatenate names
  all_rank <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "ASV")
  ranks <- all_rank[1:match(plotting_tax_ranks, all_rank)]
  
  ### Absolute value filtering
  if (relative_or_absolute_filtering == "absolute"){
    ## Define a tag to fill the taxonomy of low-abundance taxa
    filtering_tag <- paste0("<", filter_value, "reads")
    
    count_table_long_tax_sum <- count_table_long_tax %>% 
      dplyr::group_by_at(c("Sample", ranks)) %>% 
      dplyr::summarise(Count = sum(Count)) %>% # Regroup counts for identical taxa within each sample
      dplyr::ungroup() %>% 
      dplyr::mutate_at(.vars = vars(ranks), .funs = list(newtax = ~ dplyr::if_else(condition = (filter_value < Count), true = unique(.), false = filtering_tag))) %>% # Each taxa represented less counts that the cut-off are replated by the tag
      dplyr::select(-one_of(ranks)) %>% # remove the orginal taxonomic columns
      dplyr::rename_at(.vars = vars(ends_with("_newtax")), 
                       .funs = ~sub(x = ., pattern = "_newtax$", replacement = "")) # Rename the columns just created by the mutate function to the original names
    
    
    ### Relative value filtering
  } else if (relative_or_absolute_filtering == "relative"){
    ## Define a tag to fill the taxonomy of low-abundance taxa
    filtering_tag <- paste0("<", filter_value, "%_of_reads")
    
    count_table_long_tax_sum <- count_table_long_tax %>% 
      dplyr::group_by_at(c("Sample", ranks)) %>% 
      dplyr::summarise(Count = sum(Count)) %>% # Regroup counts for identical taxa within each sample
      dplyr::ungroup() %>%
      dplyr::group_by_at("Sample") %>% 
      dplyr::mutate(S_sum = sum(Count)) %>% # Create a column containing the sum of counts for each samples
      dplyr::ungroup() %>%
      dplyr::group_by_at(c("Sample", ranks)) %>% 
      dplyr::mutate(Count = (100*Count)/S_sum) %>% # Compute the % for each taxa within each sample
      dplyr::mutate_at(.vars = vars(ranks), .funs = list(newtax = ~ dplyr::if_else(condition = (filter_value < Count), true = unique(.), false = filtering_tag))) %>% # Each taxa represented less counts that the cut-off are replated by the tag
      dplyr::ungroup() %>%
      dplyr::select(-one_of(ranks)) %>% # remove the orginal taxonomic columns
      dplyr::rename_at(.vars = vars(ends_with("_newtax")),
                       .funs = ~sub(x = ., pattern = "_newtax$", replacement = "")) %>% # Rename the columns just created by the mutate function to the original names
      dplyr::select(-S_sum) %>% # Remove the sum of counts per sample column
      dplyr::group_by_at(c("Sample", ranks)) %>% 
      dplyr::summarise(Count = sum(Count)) %>% # Regroup the identical taxa per sample (the filtered ones)
      dplyr::ungroup()
    
    
    
    ### No filtering
  } else if (relative_or_absolute_filtering == "nofiltering"){
    
    count_table_long_tax_sum <- count_table_long_tax
    
  }
  
  
  ## Add metadata to this long form at count table with taxonomy
  ### Keep only the used columns to keep the table small
  metadata_table_select <- select(metadata_table, c(Sample, x_axis_column, grouping_column, facetting_column, sort_column))
  
  ### Join tables
  count_table_long_tax_meta <- dplyr::right_join(count_table_long_tax_sum, metadata_table_select, by = "Sample") 
  ### Order the x_axis_column based on the column order
  count_table_long_tax_meta[[x_axis_column]] <- reorder(count_table_long_tax_meta[[x_axis_column]], count_table_long_tax_meta[[sort_column]])
  
  ## Add a color palette for the taxonomy
  #### Brewer colors
  if (distinct_colors == FALSE){
    getPalette <- colorRampPalette(brewer.pal(n=9, "Set1"))
    ColList <- unique(count_table_long_tax_meta[[plotting_tax_ranks]])
    colors_palette <- getPalette(length(ColList))
    names(colors_palette) <- ColList
    colors_palette[filtering_tag] <- "#d3d3d3" # Set the filtered in balck
    m <- match(count_table_long_tax_meta[[plotting_tax_ranks]], names(colors_palette))
    count_table_long_tax_meta$colors <- colors_palette[m]
    
    
    #### Randomcolors distinct colors
  }else if (distinct_colors == TRUE){
    set.seed(4)
    ColList <- unique(count_table_long_tax_meta[[plotting_tax_ranks]])
    colors_palette <- distinctColorPalette(altCol = FALSE, k = length(unique(count_table_long_tax_meta[[plotting_tax_ranks]])))
    names(colors_palette) <-  ColList
    colors_palette[filtering_tag] <- "#d3d3d3" # Set the filtered in balck
    m <- match(count_table_long_tax_meta[[plotting_tax_ranks]], names(colors_palette))
    count_table_long_tax_meta$colors <- colors_palette[m]
    names(count_table_long_tax_meta$colors) <- count_table_long_tax_meta[[plotting_tax_ranks]]
  }

  
  ## Return filtered table
  return(count_table_long_tax_meta)
  
}

## Create barplot
barplots_fct <- function(long_count_table, grouping_column, grouping_column_value, x_axis_column, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value, relative_or_absolute_plot = c("relative", "absolute"), plotting_tax_ranks, horizontal_plot, facet_plot, facetting_column, order_by_abundance){
  
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
    
  } else if(plotting_tax_ranks == "ASV"){
    label_ranks <- c("Species", "ASV")
  }
  
  
  ## Keep the sample matching a specified value in the groupe_column variable. In option, we integrate negative controls too.
  if(isTRUE(t_neg_PCR_sample_on_plots)){
    filtered_count_table_f <- long_count_table[long_count_table[[grouping_column]] == grouping_column_value | long_count_table[[grouping_column]] == t_neg_PCR_group_column_value,] 
  } else if (isFALSE(t_neg_PCR_sample_on_plots)){
    filtered_count_table_f <- long_count_table[long_count_table[[grouping_column]] == grouping_column_value,] 
  }
  
  ## Normalize in % the table or plot absolute abundance
  filtered_count_table_f_norm <- filtered_count_table_f  %>% # calculate % normalized Abundance
    ungroup %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Count=as.numeric(100*Count/sum(Count))) %>%
    ungroup() 
  
  ## Reorder by abundance
  ### Unquote a factors
  t_column <- rlang::sym(plotting_tax_ranks)
  ### Reorder
  if (isTRUE(order_by_abundance)){
    filtered_count_table_ord <- filtered_count_table_f_norm %>%
      group_by(!! t_column) %>%
      mutate(tot = sum(Count)) %>%
      ungroup() %>%
      mutate(!! t_column := fct_reorder(!! t_column, tot)) %>%
      arrange(desc(tot))
  }else if (isFALSE(order_by_abundance)){ print("NOT ordered by abundance")
    filtered_count_table_ord <- filtered_count_table_f_norm
  }
  
  ### Set the filtering_tag at the top of the plot
  filtered_count_table_ord[[plotting_tax_ranks]] <- fct_relevel(filtered_count_table_ord[[plotting_tax_ranks]], unique(grep(filtered_count_table_ord[[plotting_tax_ranks]], value = TRUE, pattern = "<")), after = 0)
  
  ## Generate taxonomic and color labels
  ### Renames the values in the vector of plotted used taxrank with their related OTU name to keep them matching. Without this step, the labels do NOT match the rows
  taxalabel <- as.character(paste(toupper(substr(filtered_count_table_ord[[label_ranks[1]]],1 ,6)) , filtered_count_table_ord[[label_ranks[2]]], sep = ";"))
  names(taxalabel) <- filtered_count_table_ord[[plotting_tax_ranks]]
  
  ### Colors are recovered from an column and renamed to match their taxonomic assignment
  col <- as.character(filtered_count_table_ord$colors)
  names(col) <- as.character(filtered_count_table_ord[[plotting_tax_ranks]])
  
  
  ### Create the barplot
  taxrank_barplot <- filtered_count_table_ord %>%
    ggplot(aes_string(x = x_axis_column, y = "Count", fill = plotting_tax_ranks)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), legend.text=element_text(size=10), plot.title = element_text(hjust = 0.5)) + # axis and title settings
    guides(fill = guide_legend(title = paste0(plotting_tax_ranks),reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) + # settings of the legend
    labs(x = x_axis_column,  y = paste("Counts"), title = paste(plotting_tax_ranks,"level", grouping_column)) + # axis and graph title
    scale_fill_manual(values = col, labels = taxalabel) # colors and taxa labels as set previously
  
  
  #### In option, turn barplot horizontally and create facet view
  if (isTRUE(horizontal_plot)){
    taxrank_barplot <- taxrank_barplot + coord_flip() ### Reverse the order of the samples
    if (isTRUE(facet_plot)){
      taxrank_barplot <- taxrank_barplot + facet_grid(get(facetting_column) ~., scales = "free", space = "free")
    }
    
  }else if(!isTRUE(horizontal_plot)){
    if (isTRUE(facet_plot)){
      taxrank_barplot <- taxrank_barplot + facet_grid(~ get(facetting_column), scales = "free", space = "free")
    }
  }
  
  return(taxrank_barplot)
}

## Heatmaps
heatmaps_fct <- function(long_count_table, metadata_table, grouping_column, grouping_column_value, x_axis_column, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value, relative_or_absolute_plot = c("relative", "absolute"), plotting_tax_ranks, horizontal_plot, facet_plot ,facetting_column, high_color, low_color, log_transform){
  
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
    
  } else if(plotting_tax_ranks == "ASV"){
    label_ranks <- c("Species", "ASV")
  }
  
  
  ## Keep the sample of matching a specified value in the groupe_column variable
  if(isTRUE(t_neg_PCR_sample_on_plots)){
    filtered_count_table_f <- long_count_table[long_count_table[[grouping_column]] == grouping_column_value | long_count_table[[grouping_column]] == t_neg_PCR_group_column_value,] 
  } else if (isFALSE(t_neg_PCR_sample_on_plots)){
    filtered_count_table_f <- long_count_table[long_count_table[[grouping_column]] == grouping_column_value,] 
  }
  
  ## Keep the sample of matching a specified value in the groupe_column variable
  if(isTRUE(t_neg_PCR_sample_on_plots)){
    metadata_table_f <- metadata_table[metadata_table[[grouping_column]] == grouping_column_value | metadata_table[[grouping_column]] == t_neg_PCR_group_column_value,] 
  } else if (isFALSE(t_neg_PCR_sample_on_plots)){
    metadata_table_f <- metadata_table[metadata_table[[grouping_column]] == grouping_column_value,] 
  }
  
  ## Normalize counts in %
  if (relative_or_absolute_plot == "relative"){
    ## Normalize in % the table or plot absolute abundance
    filtered_count_table_f_norm <- filtered_count_table_f  %>% # calculate % normalized Abundance
      ungroup %>%
      dplyr::group_by(Sample) %>%
      dplyr::mutate(Count=as.numeric(100*Count/sum(Count))) %>%
      ungroup()
  }else if (relative_or_absolute_plot == "absolute"){
    filtered_count_table_f_norm <- filtered_count_table_f  
  }
  
  ## Generate taxa labels
  filtered_count_table_f_norm$taxalabel <- as.character(paste(toupper(substr(filtered_count_table_f_norm[[label_ranks[1]]],1 ,6)) , filtered_count_table_f_norm[[label_ranks[2]]], sep = ";"))
  ### Renames the values in the vector of plotted used taxrank with their related OTU name to keep them matching. Without this step, the labels do NOT match the rows
  names(filtered_count_table_f_norm$taxalabel) <- filtered_count_table_f_norm[[plotting_tax_ranks]]
  
  ## Cluster taxa to order the heatmap
  ### Generate an OTU wide table from the long format
  filtered_df_abs_i_wide <- filtered_count_table_f_norm %>% 
    select(Sample, taxalabel,  Count) %>%
    spread(taxalabel, value = Count)
  
  #### Set the Sample ID as a rowname
  row.names(filtered_df_abs_i_wide) <- filtered_df_abs_i_wide[["Sample"]]
  filtered_df_abs_i_wide[["Sample"]] <- NULL
  
  ### Format the OTU table for clustering with vegdist
  #### Transpose the table
  filtered_df_abs_i_wide <- t(filtered_df_abs_i_wide)
  
  
  #### Set the most abundant taxa at the top
  filtered_df_abs_i_wide <- filtered_df_abs_i_wide[order(rowSums(filtered_df_abs_i_wide, na.rm = TRUE), decreasing=T),]
  #### Shape it into a matrix
  filtered_df_abs_i_wide <- data.matrix(filtered_df_abs_i_wide)
  
  ### Cluster samples based on Bray-Curits distance, remplacing Na by 0 for clustering only
  data.dist.g <- vegdist(t(filtered_df_abs_i_wide), na.rm = TRUE, method = "bray")
  data.dist.g[is.na(data.dist.g)] <- 0
  col.clus <- hclust(data.dist.g,  method = "ward.D2")
  
  ### Prepare annotation of the heatmap
  #### x_axis_column
  metadata_table_f <- data.frame(metadata_table_f)
  rownames(metadata_table_f) <- metadata_table_f[["Sample"]]
  
  #### in option, transform the counts
  ##### Get rid of 0 for tansformation
  filtered_df_abs_i_wide[filtered_df_abs_i_wide == 0] <- NA
  if(log_transform != "None"){
    filtered_df_abs_i_wide <- log(filtered_df_abs_i_wide, base = as.numeric(as.character(log_transform)))
  }
  
  
  if(isTRUE(facet_plot)){
    
    #### Color annotation
    anno <- data.frame(facetting_column = metadata_table_f[[facetting_column]])
    colnames(anno)[1] <- facetting_column
    rownames(anno) <- metadata_table_f[["Sample"]]
    
    ### Plot vertically
    if (isFALSE(horizontal_plot)){
      heat <- pheatmap(mat = filtered_df_abs_i_wide, labels_col = metadata_table_f[[x_axis_column]], annotation_col = anno, cluster_rows = FALSE , cluster_cols = col.clus, cellwidth = 11, cellheight = 11, color = colorRampPalette(c(low_color, high_color))(100), fontsize = 10)
      
    }
    
    ### Plot horizontally
    else if (isTRUE(horizontal_plot)){
      heat <- pheatmap(mat = t(filtered_df_abs_i_wide), labels_row = metadata_table_f[[x_axis_column]], annotation_row = anno,  cluster_rows = col.clus , cluster_cols = FALSE, cellwidth = 11, cellheight = 11, color = colorRampPalette(c(low_color, high_color))(100), fontsize = 8, angle_col = 45)
    }
    
  }else{
    
    
    ### Plot vertically
    if (isFALSE(horizontal_plot)){
      heat <- pheatmap(mat = filtered_df_abs_i_wide, labels_col = metadata_table_f[[x_axis_column]], cluster_rows = FALSE , cluster_cols = col.clus, cellwidth = 11, cellheight = 11, color = colorRampPalette(c(low_color, high_color))(100), fontsize = 10)
    }
    
    ### Plot horizontally
    else if (isTRUE(horizontal_plot)){
      heat <- pheatmap(mat = t(filtered_df_abs_i_wide), labels_row = metadata_table_f[[x_axis_column]], cluster_rows = col.clus , cluster_cols = FALSE, cellwidth = 11, cellheight = 11, color = colorRampPalette(c(low_color, high_color))(100), fontsize = 8, angle_col = 45)
    }
    
    
  }
  
  return(heat)
  
  
}

## NMDS
NMDS_fct <- function(count_table, metadata_table, grouping_column, grouping_column_value, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value, plotting_tax_ranks, distance, color_column, shape_column, color_palette){
  
  ## Keep the sample of matching a specified value in the groupe_column variable
  if(isTRUE(t_neg_PCR_sample_on_plots)){
    metadata_table_f <- metadata_table[metadata_table[[grouping_column]] == grouping_column_value | metadata_table[[grouping_column]] == t_neg_PCR_group_column_value,] 
  } else if (isFALSE(t_neg_PCR_sample_on_plots)){
    metadata_table_f <- metadata_table[metadata_table[[grouping_column]] == grouping_column_value,] 
  }
  
  ## Filter the columns of the count_table to keep only those
  count_table_col_f <- count_table[which(names(count_table) %in% metadata_table_f[["Sample"]])]
  count_table_f <- count_table_col_f[which(rowSums(count_table_col_f) != 0),]
  
  ## Compute the distance matrix
  if (distance == "jaccard"){
    beta_dist <- vegdist(t(count_table_f),method = distance, binary = TRUE) # By default, jaccard is NOT computed after presence/absence transformation, when it should
  }else{
    beta_dist <- vegdist(t(count_table_f), method = distance, binary = FALSE)
  }
  
  ## Generate the NMDS clustering
  mds <- metaMDS(beta_dist, try = 50)
  
  ## Re-add the metadata to the points
  mds_data <- as.data.frame(mds$points)
  mds_data$Sample <- rownames(mds_data)
  mds_data <- dplyr::left_join(mds_data, metadata_table_f, by = "Sample")
  
  ## Generate the clustering
  plot <- ggplot(mds_data, aes_string(x = "MDS1", y = "MDS2", color = color_column, shape = shape_column)) +
    geom_point(size=4) +
    theme_bw() +
    coord_equal() +
    scale_color_brewer(palette=color_palette) + 
    stat_ellipse(aes(group = get(color_column), color = get(color_column)),linetype = 2, type = "t")
  
  ## Create stress plot, from https://stackoverflow.com/questions/47124238/annotate-in-ggplot2-does-not-honor-newline-is-a-pasted-and-parsed-command
  ### This plot indicates the concordance between the real distance between the samples and the distances in the plots
  #### Format the table
  tib <- tibble(mds$diss, mds$dist, mds$dhat)
  colnames(tib) <- c("diss", "dist", "dhat")
  stress <- mds$stress
  coord_x <- min(tib$diss)
  coord_y <- max(tib$dist)
  nonmetric_r2 <- round(1 - stress * stress, digits = 3)
  linear_r2 <- round(summary(lm(mds$dist~mds$dhat))$adj.r.squared, 3)
  
  ### Generate the plot
  #### Generate a label
  nonmetric_label = c(paste0("Non-metric~fit~italic(R)^2 ==", nonmetric_r2),
                      paste0("Linear~fit~italic(R)^2 ==", linear_r2)) 
  #### Plot
  stress <- ggplot(tib,
                   aes(x = diss, y = dist)) +
    geom_point(color = "blue") +
    geom_step(aes(x = diss, y = dhat), color = "red") +
    annotate(
      geom = "text",
      x = coord_x,
      y = c(coord_y, 0.95*coord_y),
      hjust = 0,
      #vjust = 1,
      label = nonmetric_label, parse = TRUE) +
    labs(x = "Observed Dissimilarity",
         y = "Ordination Distance") +
    theme_classic() +
    coord_equal() +
    labs(caption=paste0("stress:", round(mds$stress, digits = 2))) +
    theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
  
  ## Return both plots in a list, that will be splitted then
  return(list(plot, stress))
  
}

## Alphy diversity
alpha_fct <- function(metadata_table, x_axis_column, grouping_column, grouping_column_value, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value,facet_plot = FALSE, facetting_column = NULL, index, color_column){
  
  ## Keep the sample of matching a specified value in the groupe_column variable
  if(isTRUE(t_neg_PCR_sample_on_plots)){
    metadata_table_f <- metadata_table[metadata_table[[grouping_column]] == grouping_column_value | metadata_table[[grouping_column]] == t_neg_PCR_group_column_value,] 
  } else if (isFALSE(t_neg_PCR_sample_on_plots)){
    metadata_table_f <- metadata_table[metadata_table[[grouping_column]] == grouping_column_value,] 
  }
  
  ## Plot alpha diversity
  ### Here, we only recover the pre-computed values. We could eventually re-compute these values on the fly 
  alpha <- ggplot(metadata_table_f, aes_string(x=x_axis_column, y=index, color = "NMDS_alpha_colors")) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab(label = x_axis_column) +
    ylab(label = index) +
    labs(color = color_column) +
    ggtitle(paste(grouping_column_value, "alpha-diversity")) +
    scale_colour_identity(guide = "legend", labels = metadata_table_f[[color_column]], breaks = metadata_table_f[["NMDS_alpha_colors"]])
  
  ### In option, generate facets
  if (isTRUE(facet_plot)){
    alpha <- alpha + facet_grid(.~get(facetting_column) , scales = "free", space = "free")
    
  }
  
  return(alpha)
  
}


# Define UI
ui <- fluidPage(theme = shinytheme("lumen"),
                
                titlePanel("Amplicon metagenomics visualization"),
                
                
                fluidRow(
                  column(3, h3("Import data"),
                         
                         ## Load data
                         ### Metadata
                         fileInput("input_metadata_table", "Load metadata", 
                                   placeholder = "metadata_table.tsv",
                                   multiple = FALSE, 
                                   accept = c("text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".tsv",
                                              ".xls",
                                              ".xlsx")),
                         
                         ### Taxonomy
                         fileInput("input_taxonomy_table", "Load taxonomy",
                                   placeholder = "dna-sequences_tax_assignments.txt",
                                   multiple = FALSE, 
                                   accept = c("text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".tsv")),
                         
                         ### Count table
                         fileInput("input_count_table", "Load count table",
                                   placeholder = "count_table.txt",
                                   multiple = FALSE, 
                                   accept = c("text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".tsv")),
                         
                         
                         ## Set filtering out of certain samples
                         ### Yes/no of exclsuion
                         checkboxInput(inputId = "filter_out", label = "Exclude samples", value = FALSE),
                         
                         ### If yes, which columns is used to select samples to exclude
                         conditionalPanel(condition = "input.filter_out == true",
                                          selectInput(inputId = "filtering_column", label = "Column for filtering column",choice = NULL)),
                         
                         ### If yes, which values in this columns
                         conditionalPanel(condition = "input.filter_out == true",
                                          selectInput(inputId = "filter_out_value", label = "Value in filtering column",choice = NULL, multiple = TRUE)),
                         
                         ## Select in which format to generate the plots
                         selectInput(inputId = "fformat", "Download plot format", choices=c("png","svg","jpeg","pdf"), selected = "png", multiple = FALSE, selectize = TRUE)
                         
                  ),
                  
                  column(3, h3("Select input columns"),
                         ## Select columm of interest
                         ### grouping column (a column for which an individual plot will be generated)
                         selectInput(inputId = "grouping_column", label = "Grouping column", choices = NULL),
                         
                         ### grouping column value (for which value in the grouping column do we want a plot)
                         selectInput(inputId ="grouping_column_value", label = "Grouping column value", choices = NULL),
                         
                         ### QC on plots (do we add the QC on each plots)
                         checkboxInput(inputId = "t_neg_PCR_sample_on_plots", label = "QC on each plot", value = FALSE),
                         
                         ### QC value in grouping column (which values have QC in the grouping column)
                         conditionalPanel(condition = "input.t_neg_PCR_sample_on_plots == true",
                                          selectInput(inputId ="t_neg_PCR_group_column_value", label = "QC value in grouping column",choices = NULL)),
                         
                         ### Facet plot 
                         checkboxInput(inputId = "facet_plot", label = "Facet plot", value = FALSE),
                         
                         ### Column value for facetting
                         conditionalPanel(condition = "input.facet_plot == true",
                                          selectInput(inputId ="facetting_column", label = "Column for facets", choices = NULL)),
                         
                         ### x axis
                         selectInput(inputId = "x_axis_column", label = "X axis column", choices = NULL),
                         
                         ### Sort column
                         ### x axis
                         selectInput(inputId = "sort_column", label = "Sorting column", choices = NULL)
                         
                  ),
                  
                  column(3, h3("Barplot and Heatmap settings"),
                         
                         h4("Both"),
                         
                         ### Relative or absolute filtering
                         selectInput(inputId = "relative_or_absolute_filtering",label = "Filtering type", c("relative", "absolute", "nofiltering")),
                         
                         ### Filter cumulated at the sample or group scale
                         #selectInput(inputId = "abundance_filtre_level",label = "Filtering scale",c("Sample", "Group")),
                         
                         ### Filtering value
                         numericInput(inputId = "filtering_value", label = "Filtering value", value = 1, min = 0),
                         
                         ### Plotting value
                         selectInput(inputId = "relative_or_absolute_plot", label = "Plotting value", c("relative", "absolute")),  
                         
                         ### Plotting tax rank
                         selectInput(inputId = "plotting_tax_ranks",label = "Plotting taxonomic rank", choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"), selected = "Genus"),
                         
                         ### Horizontal plot
                         checkboxInput(inputId = "horizontal_plot", label = "Horizontal plot", value = FALSE),
                         
                         
                         h4("Barplots only"),
                         
                         ### Order by abundance
                         checkboxInput(inputId = "order_by_abundance", label = "Order by abundance", value = FALSE),
                         
                         ### Distinct color
                         checkboxInput(inputId = "distinct_colors", label = "Random distinct color", value = TRUE),  
                         
                         ### Seprated legends
                         checkboxInput(inputId = "separated_legend", label = "Separated legend", value = TRUE),
                         
                         h4("Heatmap only"),
                         
                         colourInput(inputId = "low_color",label = "Low abundance color", "#000033"),
                         
                         colourInput(inputId = "high_color",label = "High abundance color", "#CCFF66"),
                         
                         selectInput(inputId = "log_transform",label = "Log transformation", choices = c("None", 2, 10), selected = FALSE)
                         
                         
                  ),
                  
                  column(3, h3("Alpha-diversity and NMDS settings"),
                         
                         h4("Both"),
                         ### Color column
                         selectInput(inputId ="alp_NMDS_color", label = "Colors column", choices = NULL),
                         
                         ### Color palette
                         selectInput(inputId ="alp_NMDS_palette", label = "Color palette", choices = c("Set1", "Set2", "Set3", "Pastel1", "Pastel2","Paired", "Dark2", "Accent", "Spectral"), selected = "Dark2"),
                         
                         h4("Alpha-diversity"),
                         ### Index
                         selectInput(inputId ="alpha_index", label = "Alpha-diversity index", choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson","InvSimpson", "Observed_min_1"), selected = "Chao1"),
                         
                         h4("NMDS of beta-diversity"),
                         ### Shapes column
                         selectInput(inputId ="NMDS_shape", label = "NMDS shapes", choices = NULL),
                         
                         ### Distances
                         selectInput(inputId ="NMDS_distance", label = "NMDS distance", choices = c("bray", "jaccard", "manhattan", "euclidean", "canberra","gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis"), selected = "bray")
                         
                         
                  )
                  
                ),
                
                tabsetPanel(
                  # Barplots plannel
                  tabPanel(title = "Baplots",
                           
                           downloadButton('Download_barplot', label = "Barplot"),
                           downloadButton('Download_barplot_legend', label = "Legend"),
                           splitLayout(cellWidths = c("75%", "25%"), plotOutput("b_plot"), plotOutput("legend"))
                           
                  ),
                  
                  
                  # Heatmaps
                  tabPanel(title = "Heatmaps",
                           
                           downloadButton('Download_heatmap', label = "Heatmap"),
                           plotOutput("heat_plot")
                  ),
                  
                  
                  ## Alpha diversity
                  tabPanel(title = "Alpha",
                           
                           downloadButton('Download_alpha', label = "Alpha-diversity"),
                           plotOutput("alpha_plot")
                           
                  ),
                  
                  # NMDS
                  tabPanel(title = "NMDS",
                           
                           downloadButton('Download_NMDS',label = "NMDS"),
                           splitLayout(cellWidths = c("50%", "50%"), plotOutput("NMDS", height = "400px", width = "500px"), plotOutput("NMDS_stress"))
                           
                  )
                  
                )
                
                
)

# Define server function
server <- function(input, output, session) {
  
  ## Increase default max size
  options(shiny.maxRequestSize=300*1024^2) 
  
  ## Load data
  ### Metadata table
  metadata_table <- reactive({
    req(input$input_metadata_table)
    
    if(any(grepl(x = input$input_metadata_table$datapath, pattern = ".tsv"), grepl(x = input$input_metadata_table$datapath, pattern = ".txt"))){
      m <- read.table(input$input_metadata_table$datapath, header = TRUE, sep = "\t", check.names = FALSE)
    } else if(any(grepl(x = input$input_metadata_table$datapath, pattern = ".xlsx"),grepl(x = input$input_metadata_table$datapath, pattern = ".xls"))){
      m <- read_excel(path = input$input_metadata_table$datapath, sheet = 1)
    }
    
    return(m)
  })
  
  ### Count table
  count_table <- reactive({
    req(input$input_count_table)
    c <- read.table(input$input_count_table$datapath, header = TRUE, sep = "\t", check.names = FALSE)
    return(c)
  })
  
  ### Taxonomy
  taxonomy_table <- reactive({
    req(input$input_taxonomy_table)
    t <- read.table(input$input_taxonomy_table$datapath, header = FALSE, sep = "\t", check.names = FALSE)
    return(t)
  })
  
  
  ## In option, filter out some samples from the metadata
  ### Columns used for filtering
  observeEvent({metadata_table()},{
    updateSelectInput(session = session, 'filtering_column',
                      choices=colnames(metadata_table())[!colnames(metadata_table()) %in%  c("Observed" ,"Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Observed_min_1")])
  })
  
  ### Values filtered out in this column
  observeEvent({metadata_table()
    input$filtering_column},{
      updateSelectInput(session = session, 'filter_out_value',
                        choices=unique(metadata_table()[[input$filtering_column]]))
    })
  
  ### Filter out the columns
  metadata_table_filtered <- reactive({
    req(metadata_table())
    metadata_table_filtered <- filter_metadata_fct(metadata_table = metadata_table(), 
                                                   filter_out = input$filter_out, 
                                                   filtering_column = input$filtering_column, 
                                                   filter_out_value = input$filter_out_value)
    
    return(metadata_table_filtered)
    
  })
  
  
  
  ## Add colors to the filtered metadata for later plotting of NMDS and alpha diversity
  metadata_table_filtered_col <- reactive({
    req(metadata_table_filtered())
    metadata_colored <- add_colors_fct(metadata_table = metadata_table_filtered(), 
                                       alp_NMDS_color = input$alp_NMDS_color, 
                                       alp_NMDS_palette = input$alp_NMDS_palette)
    
    return(metadata_colored)
    
  })
  
  ## Update parameters based on the content of the filtered metadata
  #### Grouping_column
  observeEvent({metadata_table_filtered()},{
    updateSelectInput(session = session, 'grouping_column',
                      choices=colnames(metadata_table_filtered())[!colnames(metadata_table_filtered()) %in%  c("Observed" ,"Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Observed_min_1")])
  })
  
  #### Grouping column value 
  observeEvent({metadata_table_filtered()
    input$grouping_column
    metadata_table_filtered()},{
      updateSelectInput(session = session, 'grouping_column_value',
                        choices=unique(metadata_table_filtered()[[input$grouping_column]]))
    })
  
  #### Grouping column value for QC
  observeEvent({metadata_table_filtered()
    input$t_neg_PCR_sample_on_plots
    input$grouping_column},{
      updateSelectInput(session = session,'t_neg_PCR_group_column_value',
                        choices=unique(metadata_table_filtered()[[input$grouping_column]]))
    })
  
  #### x_axis_column
  observeEvent({metadata_table_filtered()},{
    updateSelectInput(session = session, 'x_axis_column',
                      choices=colnames(metadata_table_filtered())[!colnames(metadata_table_filtered()) %in%  c("Observed" ,"Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Observed_min_1")])
  })
  
  #### Column for sorting of x_axis_column
  observeEvent({metadata_table_filtered()
    input$x_axis_column},{
      updateSelectInput(session = session, 'sort_column',
                        choices=colnames(metadata_table_filtered())[!colnames(metadata_table_filtered()) %in%  c("Observed" ,"Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Observed_min_1")])
    })
  
  #### Facetting column
  observeEvent({metadata_table_filtered()
    input$facet_plot},{
      updateSelectInput(session = session, 'facetting_column',
                        choices=colnames(metadata_table_filtered())[!colnames(metadata_table_filtered()) %in%  c("Observed" ,"Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Observed_min_1")])
    })
  
  
  #### Column for which we want alpha/NMDS colors
  observeEvent({metadata_table_filtered()},{
    updateSelectInput(session = session, 'alp_NMDS_color',
                      choices=colnames(metadata_table_filtered())[!colnames(metadata_table_filtered()) %in%  c("Observed" ,"Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Observed_min_1")])
  })
  
  ### NMDS shapes
  observeEvent({metadata_table_filtered()},{
    updateSelectInput(session = session, 'NMDS_shape',
                      choices=colnames(metadata_table_filtered())[!colnames(metadata_table_filtered()) %in%  c("Observed" ,"Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Observed_min_1")])
  })
  
  
  ## Generate a long table for barplot, heatmap and NMDS
  long_table <- reactive({
    req(metadata_table_filtered())
    req(count_table())
    req(taxonomy_table())
    
    long_table <- filter_OTU_count_fct(
      count_table = count_table(),
      metadata_table = metadata_table_filtered(), 
      taxonomy_table = taxonomy_table(), 
      grouping_column = input$grouping_column,
      x_axis_column = input$x_axis_column,
      sort_column = input$sort_column,
      facetting_column = input$facetting_column,
      relative_or_absolute_filtering = input$relative_or_absolute_filtering,
      filter_value = input$filtering_value, 
      plotting_tax_ranks = input$plotting_tax_ranks,
      distinct_colors = input$distinct_colors,
      abundance_filtre_level = c("Sample", "Group"))
    
    return(long_table)
    
  })
  
  
  ## Generate plots
  ### Barplots
  #### Generate the barplot
  b_plot <- reactive({
    req(long_table())
    
    barplots_fct(long_count_table = long_table(),
                 x_axis_column = input$x_axis_column,
                 grouping_column = input$grouping_column,
                 grouping_column_value = input$grouping_column_value,
                 t_neg_PCR_sample_on_plots = input$t_neg_PCR_sample_on_plots ,
                 t_neg_PCR_group_column_value = input$t_neg_PCR_group_column_value,
                 relative_or_absolute_plot = input$relative_or_absolute_plot,
                 plotting_tax_ranks = input$plotting_tax_ranks,
                 horizontal_plot = input$horizontal_plot,
                 facet_plot = input$facet_plot,
                 facetting_column =  input$facetting_column,
                 order_by_abundance = input$order_by_abundance)
  })
  
  #### Compute the widht of the barplot based on the number of values in the x_axis_column
  b_size <-  reactive({
    req(long_table())
    
    #### Filter the long table for the value of the grouping column
    table_i <- long_table()[long_table()[[input$grouping_column]] == input$grouping_column_value,]
    
    #### basic width
    w <- 50 + 30*(length(unique(table_i[[input$x_axis_column]])))
    
    ##### Add more if legend with plot
    if(isFALSE(input$separated_legend)){
      w <- w + 200
    }
    
    ### basic height
    h <- 400
    
    ### Revert widht and height if horizontal plot
    if(isTRUE(input$horizontal_plot)){
      return(list(w,h))
    }else{
      return(list(h,w))
    }
  })
  
  
  #### Remove the legend and render the plot
  output$b_plot <- renderPlot({
    req(b_plot())
    req(b_size())
    
    if (isTRUE(input$separated_legend)){
      p <- b_plot() + guides(fill = FALSE)
      
    } else if(isFALSE(input$separated_legend)){
      p <- b_plot()
    }
    
    ggsave(paste0("barplot.",input$fformat), p)
    
    return(p)
    
  }, height = function(){b_size()[[1]]}, width = function(){b_size()[[2]]}) # Recover the pre-computed length 
  
  output$legend <- renderPlot({
    req(b_plot())
    
    leg <- as_ggplot(get_legend(b_plot()))
    
    ggsave(paste0("legend.",input$fformat), leg)
    
    return(leg)
    
  })
  
  ### Heatmap
  #### Compute the widht of the barplot based on the number of values in the x_axis_column
  h_size <-  reactive({
    req(long_table())
    
    #### Filter the long table for the value of the grouping column
    table_i <- long_table()[long_table()[[input$grouping_column]] == input$grouping_column_value,]
    
    #### basic width
    w <- 300 + 11*(length(unique(table_i[[input$x_axis_column]])))
    
    ##### Add more if legend with plot
    if(isTRUE(input$facet_plot)){
      w <- w + 200
    }
    
    ### basic height
    h <- 400 + 11 *(length(unique(table_i[[input$plotting_tax_ranks]])))
    
    ### Revert widht and height if horizontal plot
    if(isTRUE(input$horizontal_plot)){
      return(list(w,h))
    }else{
      return(list(h,w))
    }
  })
  
  
  #### Generate plot
  output$heat_plot <- renderPlot({
    req(long_table())
    req(h_size())
    
    heat <- heatmaps_fct(long_count_table =  long_table(),
                         metadata_table = metadata_table_filtered(),
                         x_axis_column = input$x_axis_column,
                         grouping_column = input$grouping_column,
                         grouping_column_value = input$grouping_column_value,
                         t_neg_PCR_sample_on_plots = input$t_neg_PCR_sample_on_plots ,
                         t_neg_PCR_group_column_value = input$t_neg_PCR_group_column_value,
                         relative_or_absolute_plot = input$relative_or_absolute_plot,
                         plotting_tax_ranks = input$plotting_tax_ranks,
                         horizontal_plot = input$horizontal_plot,
                         facet_plot = input$facet_plot,
                         facetting_column = input$facetting_column,
                         high_color = input$high_color,
                         low_color = input$low_color,
                         log_transform = input$log_transform)
    
    ggsave(paste0("heatmap.",input$fformat), plot = heat)
    
    return(heat)
  }, height = function(){h_size()[[1]]}, width = function(){h_size()[[2]]}) # Recover the pre-computed length 
  
  
  ### NMDS
  #### Generate the plot
  NMDS_biplot <- reactive({
    req(metadata_table_filtered_col())
    plot <- NMDS_fct(count_table = count_table() ,
                     metadata_table = metadata_table_filtered_col(),
                     grouping_column = input$grouping_column, 
                     grouping_column_value = input$grouping_column_value, 
                     t_neg_PCR_sample_on_plots = input$t_neg_PCR_sample_on_plots, 
                     t_neg_PCR_group_column_value = input$t_neg_PCR_group_column_value, 
                     plotting_tax_ranks = input$plotting_tax_ranks, 
                     distance = input$NMDS_distance,
                     color_column = input$alp_NMDS_color,
                     shape_column = input$NMDS_shape,
                     color_palette = input$alp_NMDS_palette)
    
    print(paste0("NMDS.",input$fformat))
    
    ggsave(filename = paste0("NMDS.",input$fformat), plot = plot[[1]])
    
    return(plot)
  })
  
  #### Render the NMDS itself
  output$NMDS <- renderPlot({
    req(NMDS_biplot())
    return(NMDS_biplot()[[1]])
  })
  
  #### Render the stress plot
  output$NMDS_stress <- renderPlot({
    req(NMDS_biplot())
    return(NMDS_biplot()[[2]])
  })
  
  ### Alpha-diversity
  #### Compute the width of the table based on the number of values in the x_axis_column
  alp_width <-  reactive({
    req(long_table())
    
    table_i <- long_table()[long_table()[[input$grouping_column]] == input$grouping_column_value,]
    
    w <- 200 + 30*(length(unique(table_i[[input$x_axis_column]])))
    
    if(isTRUE(input$facet_plot)){
      w <- w + 100
      
    }
    
    return(w)
    
  })
  
  #### Generate alpha diversity plot
  output$alpha_plot <- renderPlot({
    req(metadata_table_filtered_col())
    alpha_plot <- alpha_fct(metadata_table = metadata_table_filtered_col(), 
                            grouping_column = input$grouping_column, 
                            grouping_column_value = input$grouping_column_value, 
                            t_neg_PCR_sample_on_plots = input$t_neg_PCR_sample_on_plots, 
                            t_neg_PCR_group_column_value = input$t_neg_PCR_group_column_value, 
                            x_axis_column = input$x_axis_column, 
                            facet_plot = input$facet_plot,
                            facetting_column =  input$facetting_column,
                            index = input$alpha_index, 
                            color_column = input$alp_NMDS_color)
    
    ggsave(filename = paste0("alpha.",input$fformat), plot = alpha_plot)
    
    return(alpha_plot)
    
  }, height = 400, width = function(){alp_width()})  # Recover the pre-computed length 
  
  
  
  ## Download the pre-computed plots. Here they were actually pre-saved and were here copied.
  ### Barplot
  output$Download_barplot <- downloadHandler(
    filename = function() {
      paste0("barplot.",input$fformat)
    },
    content = function(file) {
      file.copy(paste0("barplot.",input$fformat), file, overwrite=TRUE)
    }
  )
  
  ### Legend of the barplot
  output$Download_barplot_legend <- downloadHandler(
    filename = function() {
      paste0("legend.",input$fformat)
    },
    content = function(file) {
      file.copy(paste0("legend.",input$fformat), file, overwrite=TRUE)
    }
  )
  
  ### Heatmap
  output$Download_heatmap<- downloadHandler(
    filename = function() {
      paste0("heatmap.",input$fformat)
      
    },
    content = function(file) {
      file.copy(paste0("heatmap.",input$fformat), file, overwrite=TRUE)
    }
  )
  
  ### NMDS
  output$Download_NMDS<- downloadHandler(
    filename = function() {
      paste0("NMDS.",input$fformat)
    },
    content = function(file) {
      file.copy(paste0("NMDS.",input$fformat), file, overwrite=TRUE)
    }
  )
  
  ### Alpha-diversity
  output$Download_alpha<- downloadHandler(
    filename = function() {
      paste0("alpha.",input$fformat)
    },
    content = function(file) {
      file.copy(paste0("alpha.",input$fformat), file, overwrite=TRUE)
    }
  )
}


# Create Shiny object
shinyApp(ui = ui, server = server)

