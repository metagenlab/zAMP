# Load packages
library(shiny)
library(shinythemes)

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


# Load data
#melted_dataframe <- read.table("DADA2/4_physeq/qiimerdp_ezbiocloud201805.201909/nofiltering/norarefaction/no_collapse/base_with_tree_melted.tsv", header = TRUE, sep = "\t")

## Load metadata
#metadata <-  read.table("DADA2/4_physeq/qiimerdp_ezbiocloud201805.201909/nofiltering/norarefaction/no_collapse/base_export/metadata_table.tsv", header = TRUE, sep = "\t")


### Create functions
table_prep_fct <- function(melted_dataframe, metadata,x_axis_column, grouping_column, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value, relative_or_absolute_filtering = c("relative", "absolute", "nofiltering"), filtering_value, relative_or_absolute_plot = c("relative", "absolute"), plotting_tax_ranks = "all", distinct_colors = TRUE, horizontal_plot = FALSE, facet_plot = FALSE, facetting_column = NULL, order_by_abundance = TRUE, abundance_filtre_level = c("Sample", "Group")){
    
    
    #melted_dataframe[[x_axis_column]] <- factor(melted_dataframe[[x_axis_column]], levels = unique(metadata[[x_axis_column]]), ordered = TRUE)
    
    #melted_dataframe[[facetting_column]] <- factor(melted_dataframe[[facetting_column]], levels = unique(metadata[[facetting_column]]), ordered = TRUE)
    
    
    
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
            dplyr::group_by(!!abundance_sample, !!t_column) %>% # group the dataframe by Sample and taxa
            dplyr::mutate(sumper=as.numeric(mean(Abundance))) %>% # calculate the cumulative relative abundance of the taxa in the sample
            dplyr::filter(sumper >= filtering_value) %>%# keep only the taxa above threshold.
            ungroup
        
        ### Absolute value filtering
    }else if (relative_or_absolute_filtering == "absolute"){
        filtering <- "Absolute"
        print("Absolute value based filtering")
        physeq_subset_df_filtered <- melted_dataframe  %>%
            dplyr::group_by(!!abundance_sample, !!t_column) %>% # group the dataframe by Sample and taxa
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
        #threshod_filtered_abs_no_zero[[facetting_column]] <- fct_reorder(threshod_filtered_abs_no_zero[[facetting_column]], as.numeric(threshod_filtered_abs_no_zero[[x_axis_column]]))
        
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
        colors_palette <- getPalette(length(ColList))
        names(colors_palette) <- ColList
        colors_palette[filtering_tag] <- "#d3d3d3" # Set the filtered in balck
        m <- match(threshod_filtered_abs_no_zero[[plotting_tax_ranks]], names(colors_palette))
        threshod_filtered_abs_no_zero$colors <- colors_palette[m]
        
        
        #### Randomcolors distinct colors
    }else if (distinct_colors == TRUE){
        set.seed(4)
        
        ColList <- unique(threshod_filtered_abs_no_zero[[plotting_tax_ranks]])
        print(length(unique(threshod_filtered_abs_no_zero[[plotting_tax_ranks]])))
        colors_palette <- distinctColorPalette(altCol = TRUE, k = length(unique(threshod_filtered_abs_no_zero[[plotting_tax_ranks]])))
        names(colors_palette) <-  ColList
        colors_palette[filtering_tag] <- "#d3d3d3" # Set the filtered in balck
        m <- match(threshod_filtered_abs_no_zero[[plotting_tax_ranks]], names(colors_palette))
        threshod_filtered_abs_no_zero$colors <- colors_palette[m]
        names(threshod_filtered_abs_no_zero$colors) <- threshod_filtered_abs_no_zero[[plotting_tax_ranks]]
        
    }else{print("distinct_colors must be TRUE or FALSE")
    }
    
    return(threshod_filtered_abs_no_zero)
}

barplots_fct <- function(threshod_filtered_abs_no_zero, x_axis_column, grouping_column, grouping_column_value, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value, relative_or_absolute_filtering = c("relative", "absolute", "nofiltering"), filtering_value, relative_or_absolute_plot = c("relative", "absolute"), plotting_tax_ranks = "all", distinct_colors = TRUE, horizontal_plot = FALSE, facet_plot = FALSE, facetting_column = NULL, order_by_abundance = TRUE, abundance_filtre_level = c("Sample", "Group")){
    
    
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
    #g_column <- rlang::sym(grouping_column)
    #x_column <- rlang::sym(x_axis_column)
    t_column <- rlang::sym(plotting_tax_ranks)
    #l_ranks  <-  rlang::syms(label_ranks)
    
    if (facet_plot == TRUE){
        #f_column <- rlang::sym(facetting_column)
    }
    

    i <- grouping_column_value
    
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
    filtered_df_abs_i[[plotting_tax_ranks]] <- fct_relevel(filtered_df_abs_i[[plotting_tax_ranks]], unique(grep(filtered_df_abs_i[[plotting_tax_ranks]], value = TRUE, pattern = "abund")), after = 0)
    
    ### Renames the values in the vector of plotted used taxrank with their related OTU name to keep them matching. Without this step, the labels do NOT match the rows
    taxalabel <- as.character(paste(toupper(substr(filtered_df_abs_i[[label_ranks[1]]],1 ,6)) , filtered_df_abs_i[[label_ranks[2]]], sep = ";"))
    names(taxalabel) <- filtered_df_abs_i[[plotting_tax_ranks]]
    
    col <- as.character(filtered_df_abs_i$colors)
    names(col) <- as.character(filtered_df_abs_i[[plotting_tax_ranks]])
    
    ### Create the barplot
    taxrank_barplot <- filtered_df_abs_i %>%
        ggplot(aes(x = get(x_axis_column), y = Abundance, fill = get(plotting_tax_ranks))) +
        theme_bw() +
        geom_bar(stat = "identity") +
        #scale_x_discrete(labels = x_labels, drop = TRUE) + # Set false to keep empty bars
        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), legend.text=element_text(size=10), plot.title = element_text(hjust = 0.5)) + # axis and title settings
        guides(fill = guide_legend(title = paste0(plotting_tax_ranks),reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) + # settings of the legend
        labs(x = x_axis_column,  y = paste("Abundance"), title = paste(plotting_tax_ranks,"level", i)) + # axis and graph title
        scale_fill_manual(values = col, labels = taxalabel) # colors as set previously
    
    
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
    
    return(taxrank_barplot)
}

heatmaps_fct <- function(threshod_filtered_abs_no_zero, x_axis_column, grouping_column, grouping_column_value, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value, relative_or_absolute_filtering = c("relative", "absolute", "nofiltering"), filtering_value, relative_or_absolute_plot = c("relative", "absolute"), plotting_tax_ranks = "all", output_path, horizontal_plot = FALSE, facet_plot = FALSE, facetting_column = NULL, order_by_abundance = TRUE, separated_legend,abundance_filtre_level = c("Sample","Group")){
    
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
    x_column<- rlang::sym(x_axis_column)
    t_column <-  rlang::sym(plotting_tax_ranks)
    l_ranks  <-  rlang::syms(label_ranks)
    
    if (facet_plot == TRUE){
        f_column <- rlang::sym(facetting_column)
    }
    
    i <- grouping_column_value
    
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
        filtered_df_abs_i[[plotting_tax_ranks]] <- fct_relevel(filtered_df_abs_i[[plotting_tax_ranks]], unique(grep(filtered_df_abs_i[[plotting_tax_ranks]], value = TRUE, pattern = "abund")), after = 0)
        
        ### Renames the values in the vector of plotted used taxrank with their related OTU name to keep them matching. Without this step, the labels do NOT match the rows
        taxalabel <- as.character(paste(toupper(substr(filtered_df_abs_i[[label_ranks[1]]],1 ,6)) , filtered_df_abs_i[[label_ranks[2]]], sep = ";"))
        names(taxalabel) <- filtered_df_abs_i[[plotting_tax_ranks]]
        
        ### Generate heatmap

        
        filtered_df_abs_i_wide <- filtered_df_abs_i %>% 
            select(!! x_column, !! t_column, Abundance ) %>%
            spread(!! t_column, Abundance)
        
        #filtered_df_abs_i_wide[is.na(filtered_df_abs_i_wide)] <- 0 
        rownames(filtered_df_abs_i_wide) <- filtered_df_abs_i_wide[[x_axis_column]]
        filtered_df_abs_i_wide[[x_axis_column]] <- NULL
        #test <- filtered_df_abs_i_wide
        filtered_df_abs_i_wide <- t(filtered_df_abs_i_wide)
        filtered_df_abs_i_wide <- filtered_df_abs_i_wide[order(rowSums(filtered_df_abs_i_wide, na.rm = TRUE), decreasing=T),]
        filtered_df_abs_i_wide <- data.matrix(filtered_df_abs_i_wide)
        
        ## Cluster species and samples based on bray curtis distance, remplacing Na by 0
        data.dist <- vegdist(filtered_df_abs_i_wide, na.rm = TRUE, method = "jaccard")
        data.dist[is.na(data.dist)] <- 0
        row.clus <- hclust(data.dist, method = "ward.D2")
        
        data.dist.g <- vegdist(t(filtered_df_abs_i_wide), na.rm = TRUE, method = "jaccard")
        data.dist.g[is.na(data.dist.g)] <- 0
        col.clus <- hclust(data.dist.g,  method = "ward.D2")
        
        
        if (isFALSE(horizontal_plot)){
        plot <- pheatmap(filtered_df_abs_i_wide, cluster_rows = FALSE , cluster_cols = col.clus, cellwidth = 11, cellheight = 11, color = colorRampPalette(c("#000033", "#CCFF66"))(100), fontsize = 11)
        }
        
        if (isTRUE(horizontal_plot)){
            plot <- pheatmap(t(filtered_df_abs_i_wide), cluster_rows = col.clus , cluster_cols = FALSE, cellwidth = 11, cellheight = 11, color = colorRampPalette(c("#000033", "#CCFF66"))(100), fontsize = 11, angle_col = 45)
        }
        
        return(plot)
}

NMDS_fct <- function(threshod_filtered_abs_no_zero, shape_column, color_column, grouping_column, grouping_column_value, t_neg_PCR_sample_on_plots, t_neg_PCR_group_column_value, distance){
  
  ## Unquote factors
  g_column <- rlang::sym(grouping_column)
  s_column<- rlang::sym(shape_column)
  c_column<- rlang::sym(color_column)

 i <- grouping_column_value
  
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
  

 ## Generate OTU table 
    filtered_df_abs_i_wide <- filtered_df_abs_i %>% 
    select(Sample, OTU, Abundance ) %>%
    spread(OTU, Abundance)
  
  #filtered_df_abs_i_wide[is.na(filtered_df_abs_i_wide)] <- 0 
  rownames(filtered_df_abs_i_wide) <- filtered_df_abs_i_wide[["Sample"]]
  filtered_df_abs_i_wide[["Sample"]] <- NULL

  filtered_df_abs_i_wide <- t(filtered_df_abs_i_wide)
  filtered_df_abs_i_wide <- filtered_df_abs_i_wide[order(rowSums(filtered_df_abs_i_wide, na.rm = TRUE), decreasing=T),]
  filtered_df_abs_i_wide <- data.matrix(filtered_df_abs_i_wide)
  
  beta_dist <- vegdist(t(filtered_df_abs_i_wide),
                       index = distance)
  
  mds <- metaMDS(beta_dist)

  mds_data <- as.data.frame(mds$points)
  
  sample_data <- filtered_df_abs_i[!duplicated(filtered_df_abs_i$Sample),]
  
  mds_data$Sample <- rownames(mds_data)
  mds_data <- dplyr::left_join(mds_data, sample_data)
  
  plot <- ggplot(mds_data, aes_string(x = "MDS1", y = "MDS2", color = color_column, shape = shape_column)) +
    geom_point() +
    theme_bw() +
    coord_equal() 
  
  
  ## Create stress plot, from https://stackoverflow.com/questions/47124238/annotate-in-ggplot2-does-not-honor-newline-is-a-pasted-and-parsed-command
  tib <- tibble(mds$diss, mds$dist, mds$dhat)
  colnames(tib) <- c("diss", "dist", "dhat")
  stress <- mds$stress
  coord_x <- min(tib$diss)
  coord_y <- max(tib$dist)
  nonmetric_r2 <- round(1 - stress * stress, digits = 3)
  linear_r2 <- round(summary(lm(mds$dist~mds$dhat))$adj.r.squared, 3)
  
  nonmetric_label = c(paste0("Non-metric~fit~italic(R)^2 ==", nonmetric_r2),
                      paste0("Linear~fit~italic(R)^2 ==", linear_r2)) 
  
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
    labs(caption=paste0("stress:", round(mds$stress, digits = 2))) +
    theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))

  #return(list(plot, stress))
  return(list(plot, stress))
}

# Define UI
ui <- fluidPage(theme = shinytheme("lumen"),
                
                titlePanel("Amplicon metagenomics"),
     
                fluidRow(
                            column(3, h4("Import phyloseq melted dataframe"),
                                   
                        ## Load data
                            ### Phyloseq melted dataframe table
                            fileInput("input_melted", "Load phyloseq melted dataframe",
                                                          multiple = FALSE, 
                                                          accept = c("text/csv",
                                                                     "text/comma-separated-values,text/plain",
                                                                     ".tsv")),
                        
                        ### Filtering value column
                        checkboxInput(inputId = "filter_out", label = "Exclude samples", value = FALSE),
                        
                        conditionalPanel(condition = "input.filter_out == true",
                                         selectInput(inputId = "filtering_column",label = "Column for filtering column",choice = NULL)),
                        
                        conditionalPanel(condition = "input.filter_out == true",
                                         selectInput(inputId = "filter_out_value",label = "Value in filtering column",choice = NULL, multiple = TRUE)),
                        
                        
                                                    ### Metadata table
                            # fileInput("input_metadata", "Metadata",
                            #          multiple = FALSE, 
                            #          accept = c("text/csv",
                            #                     "text/comma-separated-values,text/plain",
                            #                     ".tsv")),
                            ),
                        
                        column(3, h4("Table columns settings"),
                        ## Select columm of interest
                            ### grouping column
                            selectInput(inputId = "grouping_column",label = "Grouping column",choice = NULL),
                        
                            ### grouping column value
                            selectInput(inputId ="grouping_column_value", label = "Grouping column value",choices = NULL),
                        
                            ### QC on plots
                            checkboxInput(inputId = "t_neg_PCR_sample_on_plots", label = "QC on each plot", value = FALSE),
                            
                            ### QC value in grouping column
                            conditionalPanel(condition = "input.t_neg_PCR_sample_on_plots == true",
                                selectInput(inputId ="t_neg_PCR_group_column_value", label = "QC value in grouping column",choices = NULL)),
                        
                            ### Facet plot
                            checkboxInput(inputId = "facet_plot", label = "Facet plot", value = FALSE),
                            
                            ### Column value for facetting
                            conditionalPanel(condition = "input.facet_plot == true",
                                selectInput(inputId ="facetting_column", label = "Column for facets", choices = NULL)),
                            
                            ### x axis
                            selectInput(inputId = "x_axis_column",label = "X axis column", choices = NULL)
                        ),
                        
                        column(3, h4("Plotting settings"),
                        
                        ## Select plotting parameters
                            ### Relative or absolute filtering
                            selectInput(inputId = "relative_or_absolute_filtering",label = "Filtering type", c("relative", "absolute", "nofiltering")),
                        
                            ### Filter cumulated at the sample or group scale
                            selectInput(inputId = "abundance_filtre_level",label = "Filtering scale",c("Sample", "Group")),
                        
                            ### Filtering value
                            numericInput(inputId = "filtering_value", label = "Filtering value", value = 1, min = 0),
                            
                            ### Plotting value
                            selectInput(inputId = "relative_or_absolute_plot",label = "Plotting value", c("relative", "absolute")),  
                        
                            ### Plotting tax rank
                            selectInput(inputId = "plotting_tax_ranks",label = "Plotting taxonomic rank", choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"), selected = "Genus"),
                        
                            ### Distinct color
                            checkboxInput(inputId = "distinct_colors", label = "Random distinct color", value = TRUE),        
                        
                            ### Horizontal plot
                            checkboxInput(inputId = "horizontal_plot", label = "Horizontal plot", value = FALSE),
                        
                            ### Order by abundance
                            checkboxInput(inputId = "order_by_abundance", label = "Order by abundance", value = FALSE),
                        
                            ### Seprated legends
                            checkboxInput(inputId = "separated_legend", label = "Separated legend", value = FALSE)
                        ),
                        
                        column(3, h4("NMDS settings"),
                               
                               ## Select plotting parameters
                               ### Color column
                               selectInput(inputId ="NMDS_color", label = "NMDS colors",choices = NULL),
                               
                               ### Shapes column
                               selectInput(inputId ="NMDS_shape", label = "NMDS shapes",choices = NULL),
                               
                               ### Distances
                               selectInput(inputId ="NMDS_distance", label = "NMDS distance", choices = c("bray", "jaccard", "wunifrac", "unifrac"), selected = "bray")
                        
                               
                               
                        )
                    ),
                
                
                tabsetPanel(
                    # Barplots plannel
                    tabPanel(title = "Baplots",
                             downloadButton('Download_barplot'),
                             splitLayout(cellWidths = c("75%", "25%"), plotOutput("b_plot", height = "1000px"), plotOutput("legend"))
                    ),

                    # Heatmaps
                    tabPanel(title = "Heatmaps",
                             downloadButton('Download_heatmap'),
                             plotOutput("heatmap", height = "1000px")
                    ),
                    
                    # NMDS
                    tabPanel(title = "NMDS",
                             downloadButton('Download_NMDS'),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("NMDS", height = "400px", width = "500px"), plotOutput("NMDS_stress"))

                    )
                    
                )
                
   
                    
                )

# Define server function
server <- function(input, output, session) {
    
    
    ## Increase default max size
    options(shiny.maxRequestSize=300*1024^2) 
    
    ## Load data
        ### Phyloseq melted dataframe
        melted_dataframe <- reactive({
            req(input$input_melted)
            #req(input$input_metadata)
            melted_dataframe <- read.table(input$input_melted$datapath, header = TRUE,sep = "\t")
            #    metadata <- read.table(input$input_metadata$datapath, header = TRUE,sep = "\t")
            #    melted_dataframe[[input$x_axis_column]] <- factor(x = melted_dataframe[[input$x_axis_column]], levels = metadata[[input$x_axis_column]], ordered = TRUE) 
                return(melted_dataframe)
            })
        
   
        ### Filtering column
        observeEvent({input$input_melted},{
                updateSelectInput(session = session, 'filtering_column',
                                  choices=colnames(melted_dataframe())[!colnames(melted_dataframe()) %in%  c("OTU", "Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
            })
        
        ### Filtering column value
        observeEvent({input$input_melted
            input$filtering_column},{
            updateSelectInput(session = session, 'filter_out_value',
                              choices=unique(melted_dataframe()[[input$filtering_column]]))
                              })
        
    ### Filter the melted dataframe to keep only vlaues of interest
    melted_dataframe_filtered <- reactive({
        req(melted_dataframe())
        if (input$filter_out == TRUE){
            req(input$filtering_column)
            req(input$filter_out_value)
            
            melted_dataframe_filtered <- melted_dataframe() %>% dplyr::filter(!(!! rlang::sym(input$filtering_column) %in% input$filter_out_value))
        } else {
            melted_dataframe_filtered <- melted_dataframe()
        }

        
        return(melted_dataframe_filtered)
    })   
    


    ## Update parameters
    ### Grouping column based in loaded melted dataframe
    #### React on changes in input
    observeEvent({input$input_melted
        melted_dataframe_filtered()},{
        #### Update the grouping_column input field based on the content of the uploaded table
        updateSelectInput(session = session, 'grouping_column',
                          choices=colnames(melted_dataframe_filtered())[!colnames(melted_dataframe_filtered()) %in%  c("OTU", "Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
    })
    
    ### Grouping column value based in based on grouping column
    observeEvent({input$input_melted
        input$grouping_column
        melted_dataframe_filtered()},{
            updateSelectInput(session = session, 'grouping_column_value',
                              choices=unique(melted_dataframe_filtered()[[input$grouping_column]]))
        })
    
    ### Grouping column value for QC
    observeEvent({input$input_melted
        input$t_neg_PCR_sample_on_plots
        input$grouping_column
        melted_dataframe_filtered()},{
            updateSelectInput(session = session,'t_neg_PCR_group_column_value',
                              choices=unique(melted_dataframe_filtered()[[input$grouping_column]]))
        })
    
    ### X axis columns based in loaded melted dataframe
    observeEvent({melted_dataframe_filtered()},{
        updateSelectInput(session = session, 'x_axis_column',
                          choices=colnames(melted_dataframe_filtered())[!colnames(melted_dataframe_filtered()) %in%  c("OTU", "Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
    })
    
    ### Facetting column
    observeEvent({input$input_melted
        input$facet_plot},{
            updateSelectInput(session = session, 'facetting_column',
                              choices=colnames(melted_dataframe_filtered())[!colnames(melted_dataframe_filtered()) %in%  c("OTU", "Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
        })
    
    
    ### NMDS colors
    observeEvent({input$input_melted},{
        updateSelectInput(session = session, 'NMDS_color',
                          choices=colnames(melted_dataframe_filtered())[!colnames(melted_dataframe_filtered()) %in%  c("OTU", "Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
      })
    
    ### NMDS shapes
    observeEvent({input$input_melted},{
      updateSelectInput(session = session, 'NMDS_shape',
                        choices=colnames(melted_dataframe_filtered())[!colnames(melted_dataframe_filtered()) %in%  c("OTU", "Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
    })
    
        
    ## Filter dataframe
        filtered_table <- reactive({
            req(melted_dataframe_filtered())
            ## Condition to prevent error message until no table is uploaded
            table_prep_fct(melted_dataframe =  melted_dataframe_filtered(),
                           metadata = metadata(), 
                           grouping_column = input$grouping_column,
                           x_axis_column = input$x_axis_column,
                           t_neg_PCR_sample_on_plots =input$t_neg_PCR_sample_on_plots,
                           t_neg_PCR_group_column_value = input$t_neg_PCR_group_column_value,
                           relative_or_absolute_filtering = input$relative_or_absolute_filtering,
                           filtering_value = input$filtering_value,
                           relative_or_absolute_plot = input$relative_or_absolute_plot,
                           plotting_tax_ranks = input$plotting_tax_ranks,
                           distinct_colors = input$distinct_colors,
                           horizontal_plot = input$horizontal_plot,
                           facet_plot = input$facet_plot,
                           facetting_column =  input$facetting_column,
                           order_by_abundance = input$order_by_abundance,
                           abundance_filtre_level = input$abundance_filtre_level)
        })
        
    ## Generate barplot
        ### Generate the plot
        b_plot <- reactive({
            req(filtered_table())
            barplots_fct(threshod_filtered_abs_no_zero =  filtered_table(),
                         grouping_column = input$grouping_column,
                         grouping_column_value = input$grouping_column_value,
                         x_axis_column = input$x_axis_column,
                         t_neg_PCR_sample_on_plots =input$t_neg_PCR_sample_on_plots,
                         t_neg_PCR_group_column_value = input$t_neg_PCR_group_column_value,
                         relative_or_absolute_filtering = input$relative_or_absolute_filtering,
                         filtering_value = input$filtering_value,
                         relative_or_absolute_plot = input$relative_or_absolute_plot,
                         plotting_tax_ranks = input$plotting_tax_ranks,
                         distinct_colors = input$distinct_colors,
                         horizontal_plot = input$horizontal_plot,
                         facet_plot = input$facet_plot,
                         facetting_column =  input$facetting_column,
                         order_by_abundance = input$order_by_abundance,
                         abundance_filtre_level = input$abundance_filtre_level)
            
            
        })
        
        b_width <-  reactive({
          req(filtered_table())
            table_i <- filtered_table()[filtered_table()[[input$grouping_column]] == input$grouping_column_value,]
            print(head(table_i))
            w <- 200 + 3*(unique(length(table_i[[input$x_axis_column]])))
            if(isFALSE(input$separated_legend)){
              w <- w + 100
            }
          return(w)
        })
        
        
        ### Remove the legend
        output$b_plot <- renderPlot({
            req(b_plot())
            
            if (isTRUE(input$separated_legend))
                p <- b_plot() + guides(fill = FALSE)
            
            else if(isFALSE(input$separated_legend))
                p <- b_plot()
            
            ggsave("barplot.svg", p)
            return(p)
            
        }, height = 400, width = function(){b_width()})
        
        output$legend <- renderPlot({
            req(b_plot())
            as_ggplot(get_legend(b_plot()))
            
        })
        
    
    ## Generate heatmaps
        ### Generate the plot
        output$heatmap <- renderPlot({
            req(filtered_table())
            plot <- heatmaps_fct(threshod_filtered_abs_no_zero =  filtered_table(),
                         grouping_column = input$grouping_column,
                         grouping_column_value = input$grouping_column_value,
                         x_axis_column = input$x_axis_column,
                         t_neg_PCR_sample_on_plots =input$t_neg_PCR_sample_on_plots,
                         t_neg_PCR_group_column_value = input$t_neg_PCR_group_column_value,
                         relative_or_absolute_filtering = input$relative_or_absolute_filtering,
                         filtering_value = input$filtering_value,
                         relative_or_absolute_plot = input$relative_or_absolute_plot,
                         plotting_tax_ranks = input$plotting_tax_ranks,
                         horizontal_plot = input$horizontal_plot,
                         facet_plot = input$facet_plot,
                         facetting_column =  input$facetting_column,
                         order_by_abundance = input$order_by_abundance,
                         abundance_filtre_level = input$abundance_filtre_level)
            ggsave("heatmap.svg", plot)
            return(plot)
        })
        
        
        ## Generate NMDS
        ### Generate the plot
        NMDS <- reactive({
          req(melted_dataframe_filtered())
          plot <- NMDS_fct(threshod_filtered_abs_no_zero =  melted_dataframe_filtered(),
                           color_column = input$NMDS_color,
                           shape_column = input$NMDS_shape ,
                           grouping_column = input$grouping_column,
                           grouping_column_value = input$grouping_column_value,
                           t_neg_PCR_sample_on_plots =input$t_neg_PCR_sample_on_plots,
                           t_neg_PCR_group_column_value = input$t_neg_PCR_group_column_value,
                            distance = input$NMDS_distance
                               )
          ggsave("NMDS.svg", plot[[1]])
          return(plot)
        })

        
        output$NMDS <- renderPlot({
          req(NMDS())
          return(NMDS()[[1]])
        })
        
        output$NMDS_stress <- renderPlot({
          req(NMDS())
          return(NMDS()[[2]])
        })
        
        
        output$Download_barplot <- downloadHandler(
            filename = function() {
                "barplot.svg"
            },
            content = function(file) {
                file.copy("barplot.svg", file, overwrite=TRUE)
            }
        )
        
        output$Download_heatmap<- downloadHandler(
            filename = function() {
                "heatmap.svg"
            },
            content = function(file) {
                file.copy("heatmap.svg", file, overwrite=TRUE)
            }
        )
        
        output$Download_NMDS<- downloadHandler(
          filename = function() {
            "NMDS.svg"
          },
          content = function(file) {
            file.copy("NMDS.svg", file, overwrite=TRUE)
          }
        )

}


# Create Shiny object
shinyApp(ui = ui, server = server)

