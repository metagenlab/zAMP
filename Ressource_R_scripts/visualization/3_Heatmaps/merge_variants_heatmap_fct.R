### Create a function
merge_variants_heatmap_fct <- function(melted_dataframe, x_axis_column, grouping_column, grouping_column_filtering = c(FALSE, TRUE), grouping_column_filtering_value, t_neg_PCR_sample_on_plots, t_neg_PCR_sample_grp_column_value, taxonomic_filtering = c(TRUE, FALSE), taxonomic_filtering_rank = "Kingdom" , taxonomic_filtering_value = "Bacteria" ,  quantity_filtering_type = c("relative", "absolute", "rank", "nofiltering", "absolute_and_rank"), absolute_quantity_filtering_value, relative_quantity_filtering_value, rank_quantity_filtering_value, plotting_value = c("relative", "absolute"), plotting_tax_ranks = "all", figures_save_dir, horizontal_barplot = FALSE, facet_plot = FALSE, facetting_column, order_by_abundance = TRUE){
  
  
  ### In option, filter for a given grouping_columns. In both cases keep only lines with Abundance not equal to 0 to reduce the dataframe size
  
  ### No filter for the grouping column
  if (isFALSE(grouping_column_filtering)){
    melted_dataframe_filtered <- melted_dataframe
    grouping_column_filtering_value <- "no_group_column_filtering"
    
    print("No filtering of the grouping_column")
  }
  
  ### Filter for the grouping column
  else if (isTRUE(grouping_column_filtering)){
    melted_dataframe_filtered <- melted_dataframe %>%
      filter(get(grouping_column) == grouping_column_filtering_value)
    
    print(paste("Data filtered for", grouping_column, "is", grouping_column_filtering_value))
    
  }else{ stop("grouping_column_filtering must be TRUE or FALSE")}
  
  
  ### In option, filter at a given taxonomic rank and for a given taxonomic identifier
  
  ### No filter
  if (isFALSE(taxonomic_filtering)) {
    physeq_subset_df <- melted_dataframe_filtered
    physeq_subet_norm_df <- physeq_subset_df  %>% 
      ungroup %>%
      dplyr::group_by(Sample) %>% 
      dplyr::mutate(per=as.numeric(100*Abundance/sum(Abundance))) %>% 
      dplyr::ungroup() %>%
      dplyr::select(-Abundance)  %>%
      dplyr::rename(Abundance = per)
    taxonomic_filtering_value <- "no_tax_rank_filtering"
    
    print("No filtering of any taxonomic rank")
  }
  
  ### A specific filter
  else if (isTRUE(taxonomic_filtering)){
    physeq_subset_df <- melted_dataframe_filtered %>%
      filter(get(taxonomic_filtering_rank) == taxonomic_filtering_value)
    physeq_subet_norm_df <- physeq_subset_df  %>% 
      group_by(Sample) %>%  
      mutate(per=paste0((100*Abundance/sum(Abundance)))) %>% 
      ungroup %>%
      select(-Abundance)  %>%
      dplyr::rename(Abundance = per)
    physeq_subet_norm_df$Abundance <- as.numeric(physeq_subet_norm_df$Abundance)
    
    print(paste("Data filtered for", taxonomic_filtering_rank, "is", taxonomic_filtering_value))
    
  }else{ stop('"taxonomic_filtering" must be TRUE or FALSE')
  }
  
  
  ### Choose between relative value or absolute value dataframe for plotting value
  
  ### Relative value plotting
  if (plotting_value == "relative"){
    plotting <- "Relative"
    print("Relative value plotting")
    plotted_df <- physeq_subet_norm_df
  }
  
  # Absolute value plotting
  else if (plotting_value == "absolute"){
    plotting <- "Absolute"
    print("Absolute value plotting")
    plotted_df <- physeq_subset_df
    
  }else{ stop('"relative_value_plotting" must be "relative" or "absolute"')
  }
  
  
  ### Filter values at a relative or absolute quantity threshold or at a given rank in a group ( x highest OTU in the group)
  
  ### No filtering 
  if (quantity_filtering_type == "nofiltering"){
    filtering <- "Nofiltering"
    quantity_filtering_value <- "0"
    print("No filtering based on quantity")
    physeq_subset_df_filtered <- physeq_subet_norm_df 
  }
  
  
  ### Relative value filtering
  else if (quantity_filtering_type == "relative"){
    filtering <- "Relative"
    quantity_filtering_value <- relative_quantity_filtering_value
    print("Relative value based filtering")
    physeq_subset_df_filtered <- physeq_subet_norm_df %>% filter(Abundance > quantity_filtering_value)
  }
  
  ### Absolute value filtering
  else if (quantity_filtering_type == "absolute"){
    filtering <- "Absolute"
    quantity_filtering_value <- absolute_quantity_filtering_value
    print("Absolute value based filtering")
    physeq_subset_df_filtered <- physeq_subset_df %>% filter(Abundance > quantity_filtering_value)
  }
  
  ### Abundance rank of the taxa in the group filtering
  else if (quantity_filtering_type == "rank"){
    filtering <- "Rank"
    print("Rank based filtering")
    quantity_filtering_value <- rank_quantity_filtering_value
    ### Define the variables as dplyr wants them
    g_column <- rlang::sym(grouping_column)
    x_column <- rlang::sym(x_axis_column)
    ### Generate the filtered dataframe 
    physeq_subset_df_filtered <- physeq_subset_df  %>% 
      group_by(!! g_column) %>%  ## group the dataframe by the grouping column
      mutate(per=as.numeric(paste0((100*Abundance/sum(Abundance))))) %>% ## calculate a relative abundance of the taxa in the group
      ungroup %>% 
      group_by((!! g_column), OTU)  %>% ## now group the OTU in each grouping columns
      mutate(sumper = as.numeric(paste0(sum(per)))) %>% ## Sum the relative abundance of each OTU in each group
      ungroup %>%
      group_by(!!! list(g_column, x_column)) %>%  ## Group the dataframe for each by group and filling column
      top_n(n = quantity_filtering_value*n_distinct(x_axis_column), wt = sumper) %>% ## Keep the x number of rows for each
      ungroup
  }
  
  
  ### Combined absolute and rank of the taxa in the group filtering
  else if (quantity_filtering_type == "absolute_and_rank"){
    filtering <- "Absolute_and_Rank"
    print("Absolute_and_Rank based filtering")
    ### Define the variables as dplyr wants them
    g_column <- rlang::sym(grouping_column)
    x_column<- rlang::sym(x_axis_column)
    ### Generate the filtered dataframe 
    ### Relative reads filtering
    physeq_subset_df_filtered <- physeq_subset_df  %>% 
      group_by(!! g_column) %>%  ## group the dataframe by the grouping column
      mutate(per=as.numeric(paste0((100*Abundance/sum(Abundance))))) %>% ## calculate a relative abundance of the taxa in the group
      ungroup %>% 
      group_by((!! g_column), OTU)  %>% ## now group the OTU in each grouping columns
      mutate(sumper = as.numeric(paste0(sum(per)))) %>% ## Sum the relative abundance of each OTU in each group
      ungroup %>%
      group_by(!!! list(g_column, x_column)) %>%  ## Group the dataframe for each by group and filling column
      top_n(n = rank_quantity_filtering_value*n_distinct(x_axis_column), wt = sumper) %>% ## Keep the x number of rows for each
      ungroup
    
    
    ### Absolute reads filtering
    physeq_subset_df_filtered <- physeq_subset_df_filtered %>% filter(Abundance > absolute_quantity_filtering_value)
    
    ### Write "quantity_filtering_value" for ggsave later
    quantity_filtering_value <- paste0("rank_", rank_quantity_filtering_value, "_absolute_", absolute_quantity_filtering_value)
    
  }else{ stop('quantity_filtering_type must be "nofiltering", "relative", "absolute", "rank" or "absolute_and_rank"')
  }
  
  
  ### Now we have a dataframe ("plotted_df") containing all values, with Abundance expressed in absolute reads number or relative quanity in %, and a dataframe containing the rows that have to be kept while the rest will be merged by tagging them with the same name
  
  
  ### In a dataframe, keep only the rows matching the filtered rowmanes
  above_threshold_df <- semi_join(plotted_df, physeq_subset_df_filtered, by = c("OTU", "Sample"))
  
  
  ### In another dataframe, keep only the lines NOT matching the filtered rowmanes
  under_threshold_df <- anti_join(plotted_df, physeq_subset_df_filtered, by =c ("OTU", "Sample"))
  
  
  ### To follow and control the process, write external .tsv tables
  ### Write the filtration table, without Abundance = 0 rows to reduce its size.
  physeq_subset_df_filtered %>%
    filter(Abundance>0) %>%
    write.table(file = paste0(figures_save_dir,"/heatmaps/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", quantity_filtering_value, "_all_", x_axis_column, "_filtration_table.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)
  
  ### Write the plotted table , without Abundance = 0 rows to reduce its size.
  plotted_df %>%
    filter(Abundance>0) %>%
    write.table(file = paste0(figures_save_dir,"/heatmaps/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", quantity_filtering_value, "_all_", x_axis_column, "_plotted_table.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)
  
  ### In another dataframe, keep the join of the under and above tables, before masking of the filtered taxa and without Abundance = 0 rows to reduce its size.
  merged_filtered_abs <- full_join(under_threshold_df,above_threshold_df) %>%
    filter(Abundance>0)
  write.table(merged_filtered_abs, file = paste0(figures_save_dir,"/heatmaps/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", quantity_filtering_value, "_all_", x_axis_column, "_merged_without_masking.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)
  
  
  
  ### Convert taxonomic all levels as character (needed for the next steps)
  above_threshold_df$Kingdom<-as.character(above_threshold_df$Kingdom)
  above_threshold_df$Phylum<-as.character(above_threshold_df$Phylum)
  above_threshold_df$Class<-as.character(above_threshold_df$Class)
  above_threshold_df$Order<-as.character(above_threshold_df$Order)
  above_threshold_df$Family<-as.character(above_threshold_df$Family)
  above_threshold_df$Genus<-as.character(above_threshold_df$Genus)
  above_threshold_df$Species<-as.character(above_threshold_df$Species)
  above_threshold_df$OTU<-as.character(above_threshold_df$OTU)
  
  under_threshold_df$Kingdom<-as.character(under_threshold_df$Kingdom)
  under_threshold_df$Phylum<-as.character(under_threshold_df$Phylum)
  under_threshold_df$Class<-as.character(under_threshold_df$Class)
  under_threshold_df$Order<-as.character(under_threshold_df$Order)
  under_threshold_df$Family<-as.character(under_threshold_df$Family)
  under_threshold_df$Genus<-as.character(under_threshold_df$Genus)
  under_threshold_df$Species<-as.character(under_threshold_df$Species)
  under_threshold_df$OTU<-as.character(under_threshold_df$OTU)
  
  ### Rename taxa with < percent/reads abundance, defending of the filtering applied
  ### Relative filtering
  if (quantity_filtering_type == "relative"){
    under_threshold_df$Kingdom <-paste("<_", quantity_filtering_value,"%_abund")
    under_threshold_df$Phylum <-paste("<_", quantity_filtering_value,"%_abund")
    under_threshold_df$Class <-paste("<_", quantity_filtering_value,"%_abund")
    under_threshold_df$Order <-paste("<_", quantity_filtering_value,"%_abund")
    under_threshold_df$Family <-paste("<_", quantity_filtering_value,"%_abund")
    under_threshold_df$Genus <-paste("<_", quantity_filtering_value,"%_abund")
    under_threshold_df$Species <-paste("<_", quantity_filtering_value,"%_abund")
  }
  
  ### Absolute filtering
  else if (quantity_filtering_type == "absolute"){
    under_threshold_df$Kingdom <-paste("<_", quantity_filtering_value,"reads")
    under_threshold_df$Phylum <-paste("<_", quantity_filtering_value,"reads")
    under_threshold_df$Class <-paste("<_", quantity_filtering_value,"reads")
    under_threshold_df$Order <-paste("<_", quantity_filtering_value,"reads")
    under_threshold_df$Family <-paste("<_", quantity_filtering_value,"reads")
    under_threshold_df$Genus <-paste("<_", quantity_filtering_value,"reads")
    under_threshold_df$Species <-paste("<_", quantity_filtering_value,"reads")
  }
  
  ### Abundance rank of the taxa in the group filtering
  else if (quantity_filtering_type == "rank"){
    under_threshold_df$Kingdom <-paste("<_", quantity_filtering_value,"rank")
    under_threshold_df$Phylum <-paste("<_", quantity_filtering_value,"rank")
    under_threshold_df$Class <-paste("<_", quantity_filtering_value,"rank")
    under_threshold_df$Order <-paste("<_", quantity_filtering_value,"rank")
    under_threshold_df$Family <-paste("<_", quantity_filtering_value,"rank")
    under_threshold_df$Genus <-paste("<_", quantity_filtering_value,"rank")
    under_threshold_df$Species <-paste("<_", quantity_filtering_value,"rank")
  }  
  
  ### Absolute abundance  AND rank of the taxa in the group filtering
  else if (quantity_filtering_type == "absolute_and_rank"){
    under_threshold_df$Kingdom <-paste("<_", quantity_filtering_value)
    under_threshold_df$Phylum <-paste("<_", quantity_filtering_value)
    under_threshold_df$Class <-paste("<_", quantity_filtering_value)
    under_threshold_df$Order <-paste("<_", quantity_filtering_value)
    under_threshold_df$Family <-paste("<_", quantity_filtering_value)
    under_threshold_df$Genus <-paste("<_", quantity_filtering_value)
    under_threshold_df$Species <-paste("<_", quantity_filtering_value)            }
  
  
  ### Join the two dataframes, to put back all rows together, the ones above threshold keeping their original taxonomic identifier while the others are now grouped together
  threshod_filtered_abs <- full_join(under_threshold_df,above_threshold_df)
  
  ### Remove the (now but needed before) Abundance = 0 rows
  threshod_filtered_abs_no_zero <- filter(threshod_filtered_abs, Abundance>0) 
  
  ### Reorder the facet factor if later used for plotting
  if (isTRUE(facet_plot)){
    threshod_filtered_abs_no_zero[[facetting_column]] <- fct_reorder(threshod_filtered_abs_no_zero[[facetting_column]], as.numeric(threshod_filtered_abs_no_zero[[x_axis_column]]))
  }
  else if (isFALSE(facet_plot)){
    print("No faceting")
  }
  else {
    stop('"facet_plot" must be TRUE or FALSE')
  }
  
  ### Reoder the x_label_column for later if using horizontal barplot
  if (isTRUE(horizontal_barplot)){
    threshod_filtered_abs_no_zero[[x_axis_column]] <- fct_rev(threshod_filtered_abs_no_zero[[x_axis_column]])
  }
  else if (isFALSE(horizontal_barplot)){
    print("Vertical plotting")
  }
  else {
    stop('"horizontal_barplot" must be TRUE or FALSE')
  }
  
  threshod_filtered_abs_no_zero <- threshod_filtered_abs_no_zero
  
  ### Loop for unique value in grouping_column
  for (i in unique(threshod_filtered_abs_no_zero[[grouping_column]])) {
    print(paste("Start plotting", grouping_column, i))
    
    ### filter the table for this value of the grouping columns. Depending of the used arguments, the t_neg_PCR values are kept or not on the barplots
    ### Keep t_neg_PCR rows
    if (isTRUE(t_neg_PCR_sample_on_plots) & !is_null(t_neg_PCR_sample_grp_column_value)){
      filtered_df_abs_i <- filter(threshod_filtered_abs_no_zero, get(grouping_column) == i 
                                  | get(grouping_column) == t_neg_PCR_sample_grp_column_value)
      print('Keeping t_neg_PCR values for the graphs. The "t_neg_PCR_sample_grp_column_value" must match the one in the grouping_column for this sample')
    }
    else if (isTRUE(t_neg_PCR_sample_on_plots) & is_null(t_neg_PCR_sample_grp_column_value)){
      stop('If "t_neg_PCR_sample_on_plots" is "TRUE, a "t_neg_PCR_sample_grp_column_value" indicating the value of the T neg PCR sample in the grouping column must be indicated')
    }
    
    ### Do not keep t_neg_PCR rows
    else if (isFALSE(t_neg_PCR_sample_on_plots)){
      filtered_df_abs_i <- filter(threshod_filtered_abs_no_zero, get(grouping_column) == i)
    }
    
    else{ stop('"t_neg_PCR_sample_on_plots" must be TRUE or FALSE')
    }
    

    ### Write this table in a external file
    write.table(filtered_df_abs_i, file = paste0(figures_save_dir,"/heatmaps/", plotting, "/",filtering,"/Table/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", quantity_filtering_value, "_",grouping_column, "_", i, "_", x_axis_column, "_abundancy_table.tsv"), append = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)
    
    ### No filtering 
    if (quantity_filtering_type != "nofiltering"){
  
    ### Merged rows that were filtered in previous step so that they are only on one line on the heatmaps
    g_column <- rlang::sym(grouping_column)
    x_column<- rlang::sym(x_axis_column)
    f_column <- rlang::sym(facetting_column)
    filtered_OTU <- filtered_df_abs_i %>% 
      dplyr::filter(grepl("<", Species)) %>%
      group_by(Sample, !!(x_column), !!(g_column), !!(f_column)) %>%
      summarise(Abundance = sum(Abundance)) %>% 
      filter(Abundance > 0)
    
    ### Rename "filtered" those rows representing filtered rows
    filtered_OTU$Kingdom <- "Filtered"
    filtered_OTU$Phylum <- "Filtered"
    filtered_OTU$Class <- "Filtered"
    filtered_OTU$Order <- "Filtered"
    filtered_OTU$Family <- "Filtered"
    filtered_OTU$Genus <- "Filtered"
    filtered_OTU$Species <- "Filtered"
    filtered_OTU$OTU <- "Filtered"
    
    ### Keep only the NOT filtered rows 
    individual_OTU <- filtered_df_abs_i %>% 
      dplyr::filter(!grepl("<", Species))
    
    ### Bind the filtered and individual rows
    merged_filtered_OTU <- bind_rows(filtered_OTU, individual_OTU)
    
    }

    if (quantity_filtering_type == "nofiltering"){ 
      merged_filtered_OTU <- filtered_df_abs_i
    }
    
    ### Reorder by abundance
    if (isTRUE(order_by_abundance)){
      print("Ordered by abundance")
      merged_filtered_OTU$Kingdom <- factor(merged_filtered_OTU$Kingdom, levels = unique(merged_filtered_OTU$Kingdom[order(as.numeric(merged_filtered_OTU$Abundance), decreasing = TRUE)]))
      merged_filtered_OTU$Phylum <- factor(merged_filtered_OTU$Phylum, levels = unique(merged_filtered_OTU$Phylum[order(as.numeric(merged_filtered_OTU$Abundance), decreasing = TRUE)]))
      merged_filtered_OTU$Class <- factor(merged_filtered_OTU$Class, levels = unique(merged_filtered_OTU$Class[order(as.numeric(merged_filtered_OTU$Abundance), decreasing = TRUE)]))
      merged_filtered_OTU$Order <- factor(merged_filtered_OTU$Order, levels = unique(merged_filtered_OTU$Order[order(as.numeric(merged_filtered_OTU$Abundance), decreasing = TRUE)]))
      merged_filtered_OTU$Family <- factor(merged_filtered_OTU$Family, levels = unique(merged_filtered_OTU$Family[order(as.numeric(merged_filtered_OTU$Abundance), decreasing = TRUE)]))
      merged_filtered_OTU$Genus <- factor(merged_filtered_OTU$Genus, levels = unique(merged_filtered_OTU$Genus[order(as.numeric(merged_filtered_OTU$Abundance), decreasing = TRUE)]))
      merged_filtered_OTU$Species <- factor(merged_filtered_OTU$Species, levels = unique(merged_filtered_OTU$Species[order(as.numeric(merged_filtered_OTU$Abundance), decreasing = TRUE)]))
      merged_filtered_OTU$OTU <- factor(merged_filtered_OTU$OTU, levels = unique(merged_filtered_OTU$OTU[order(as.numeric(merged_filtered_OTU$Abundance), decreasing = TRUE)]))
    }
  
    else if (isFALSE(order_by_abundance)){
      ### Set filtered at the top on the heatmap and the rest on order
      
      merged_filtered_OTU$Kingdom <- fct_rev(fct_relevel(merged_filtered_OTU$Kingdom, "Filtered", after = 0))
      merged_filtered_OTU$Phylum <- fct_rev(fct_relevel(merged_filtered_OTU$Phylum, "Filtered", after = 0))
      merged_filtered_OTU$Class <- fct_rev(fct_relevel(merged_filtered_OTU$Class, "Filtered", after = 0))
      merged_filtered_OTU$Order <- fct_rev(fct_relevel(merged_filtered_OTU$Order, "Filtered", after = 0))
      merged_filtered_OTU$Family <- fct_rev(fct_relevel(merged_filtered_OTU$Family, "Filtered", after = 0))
      merged_filtered_OTU$Genus <- fct_rev(fct_relevel(merged_filtered_OTU$Genus, "Filtered", after = 0))
      merged_filtered_OTU$Species <- fct_rev(fct_relevel(merged_filtered_OTU$Species, "Filtered", after = 0))
      merged_filtered_OTU$OTU <- fct_rev(fct_relevel(merged_filtered_OTU$OTU, "Filtered", after = 0))
      
    }
    
    
    ### Create a loop for all taxonomic ranks to be plotted
    
    if (plotting_tax_ranks == "all"){
      tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")
      
    }else{
      tax_ranks <- plotting_tax_ranks
      print('Plotting at a specific taxonomic rank. Check spelling one or several of ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"). This error does NOT always mean that there was an error in spelling')
    }
    
    for (t in tax_ranks){
      
        to_melt_df <- merged_filtered_OTU
        
        g_column <- rlang::sym(grouping_column)
        x_column<- rlang::sym(x_axis_column)
        f_column <- rlang::sym(facetting_column)
        taxa_column <- rlang::sym(t)
        merged_taxa_df <- to_melt_df %>%
          dplyr::group_by(Sample, !!(x_column), !!(g_column), !!(f_column), !!(taxa_column)) %>%
          dplyr::mutate(Abundance = sum(Abundance)) %>% 
          dplyr::ungroup() %>%
          dplyr::distinct(Sample, !!(taxa_column), .keep_all = TRUE)
        
        
      ### Renames the values in the vector of plotted used taxrank with their related OTU name to keep them matching. Without this step, the labels do NOT match the rows
      taxalabel <- as(merged_taxa_df[[t]], "character")
      names(taxalabel) <- merged_taxa_df[[t]]
      
      ### Generate heatmap
      heatmap <- ggplot(merged_taxa_df, aes(x = get(x_axis_column), y = get(t), fill = Abundance)) + 
        theme_bw() +
        geom_tile(aes(fill=Abundance), show.legend=TRUE) +
        scale_fill_gradient(na.value = "white", low="#000033", high="#CCFF66") +
        geom_text(aes(label = round(Abundance, digits = 1), color = Abundance > 15), size = 2 ) +
        scale_color_manual(guide = FALSE, values = c("white", "black")) +
        theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
        scale_y_discrete(labels = taxalabel) +
        labs(x= x_axis_column,  y = t)
      
      
      
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
      filename_base <- (paste0(figures_save_dir,"/heatmaps/", plotting, "/",filtering,"/", taxonomic_filtering_rank, "_",taxonomic_filtering_value,"_",filtering, "u", quantity_filtering_value, "_",grouping_column, "_", i, "_", x_axis_column, "_",facetting_column,"_" ,t))
      ### Print the filename to follow progress
      print(filename_base)
      ###Save the figure
      ggsave(heatmap, filename = paste0(filename_base, "_heatmap.png"), width = 10, height = 7)
      
      
      
      
    }}}

