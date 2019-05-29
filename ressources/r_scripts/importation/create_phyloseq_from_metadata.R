

# Generate a phyloseq object from what have been loaded into the R environment. 
## More on phyloseq is available here: (https://joey711.github.io/phyloseq/import-data.html)
## https://github.com/joey711/phyloseq/issues/901


### Generate a function to create a phyloseq object
create_phyloseq_fct <- function(physeq_name = "physeq", melted_df_name = "physeq_df", otu_table, phy_tree = NULL, tax_table, Metadata_table){
  
  ### Create the phyloseq object in cases where no phylogenic tree are given
  if (is.null(phy_tree)){
    physeq_obj<-phyloseq(
      otu_table(otu_table$data, taxa_are_rows = T), 
      tax_table(as.data.frame(tax_table) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
      sample_data(Metadata_table %>% as.data.frame() %>% column_to_rownames("Sample")))
  }
  else{
    ### Create the phyloseq object with a phylogenic tree
    physeq_obj<-phyloseq(
      otu_table(otu_table$data, taxa_are_rows = T), 
      phy_tree(phy_tree$data), 
      tax_table(as.data.frame(tax_table) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
      sample_data(Metadata_table %>% as.data.frame() %>% column_to_rownames("Sample")))
  }
  
  ### Assign the object into the environnement 
  assign(value = physeq_obj, x = physeq_name, envir =  .GlobalEnv)
  
  ### Melt the physeq objet into a dataframe with one row per feature.id and per sample, needed later
  physeq_df <- psmelt(physeq_obj)
  
  ### Assign the object into the environnement 
  assign(value = physeq_df, x = melted_df_name, envir =  .GlobalEnv)
  
  ### Write the phyloseq object in a file, so it doesn't have to be recalculated each time
  save(physeq_obj, file = "physeq_object")
  
  ### Write this table in a tsv file since it is ground for coming analysis and slow to compute
  write.table(x = physeq_df, file = "physeq_df.tsv", append = F, quote = F, sep = "\t", eol = "\n", row.names = F, col.names = T )
  
}