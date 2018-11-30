## Function to create a Metadata table containing the number of reads per Sample

reads_counts_df_fct <- function(Metadata_table = Metadata , out_put_name = "reads_counts_df", melted_dataframe = physeq_df ){
  
  print('"Metadata_table" column names must match the "melted_dataframe" ones')
  
  reads_counts_df <- melted_dataframe %>% group_by(Sample) %>%
    mutate(TotalReads = sum(Abundance)) %>%
    ungroup() %>%
    distinct(Sample, .keep_all = TRUE) %>%
    select(c(colnames(Metadata), TotalReads))
  
  ### Assign the object into the environnement 
  assign(x = out_put_name, value = reads_counts_df, envir =  .GlobalEnv)
  
}