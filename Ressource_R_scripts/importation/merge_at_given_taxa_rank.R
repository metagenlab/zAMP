
merge_variants_df_fct <- function(melted_dataframe, rank_sequences_merging){
  
to_melt_df <- melted_dataframe

rank_colum <- rlang::sym(rank_sequences_merging)
merged_melted_df <- to_melt_df %>% 
  dplyr::group_by(Sample, !!(rank_colum)) %>%
  dplyr::mutate(Abundance = sum(Abundance)) %>% 
  dplyr::ungroup() %>%
  dplyr::distinct(Sample, !!(rank_colum), .keep_all = TRUE)

name <- paste0(rank_sequences_merging, "_merged_df") 

assign(value = merged_melted_df, x = name, envir =  .GlobalEnv)


}