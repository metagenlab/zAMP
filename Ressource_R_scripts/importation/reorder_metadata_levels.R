
############################################################################################################

## Define levels for Metadata columns
# Define manually or automatically an order for the filling column, e.g. the sample_source column and the grouping column, e.g. the patientID column (sorting them in a way compatible with ggplot2). This way it will be shown in this given order in plots. Should be adapted depending of the desired output of plots.
#
#The level is set for two columns:
#  - the grouping column, which will be used to group samples together on individual graphs
#  - the filling column, which will be used to fill barplots of reads in different colors
#  - the sample_label, which will be used as label in plots. It should be unique as SampleID but more descriptive or duplicated (e.g.), for instance if there is one per grouping_column.
#    By default, it is ordered by the order of other columns but can be set manually
#  - the facetting_column, if it is intendend to create facets for a specific factors, creating one common graph but with multiple pannels


## Set the levels for our columns of interest

### create a function to define levels
reoder_columns_levels_fct <- function(grouping_column, grouping_col_order_vector = NULL, filling_column, filling_col_order_vector = NULL, facetting_column = NULL, facetting_col_order_vector = NULL, sample_label, x_axis_col_order_vector = NULL,  Metadata_table){
  
  ### Order the grouping_column
  ### In absence of specific order vector, order by numeric order
  if (is.null(grouping_col_order_vector)){
    Metadata_table[[grouping_column]] <- factor(Metadata_table[[grouping_column]], levels = sort(unique(Metadata_table[[grouping_column]])), ordered = TRUE)
  }else{L
    ### In option, we could set an order manually
    Metadata_table[[grouping_column]] <- factor(Metadata_table[[grouping_column]], levels = grouping_col_order_vector, ordered = TRUE)
  }
  Metadata <<- Metadata_table
  
  
  ### Order the filling_column
  ### In absence of specific order vector, order by numeric order
  if (is.null(filling_col_order_vector)){
    Metadata_table[[filling_column]] <- factor(Metadata_table[[filling_column]], levels = sort(unique(Metadata_table[[filling_column]])), ordered = TRUE)
  }else{L
    ### In option, we could set an order manually
    Metadata_table[[filling_column]] <- factor(Metadata_table[[filling_column]], levels = filling_col_order_vector, ordered = TRUE)
  }
  Metadata <<- Metadata_table
  
  
  
  ### Order the facetting_column
  ### In absence of specific order vector, order by numeric order
  if (is.null(facetting_col_order_vector)){
    Metadata_table[[facetting_column]] <- factor(Metadata_table[[facetting_column]], levels = sort(unique(Metadata_table[[facetting_column]])), ordered = TRUE)
  }else{
    ### In option, we could set an order manually
    Metadata_table[[facetting_column]] <- factor(Metadata_table[[facetting_column]], levels = facetting_col_order_vector, ordered = TRUE)
  }
  Metadata <<- Metadata_table
  
  
  
  ### By default, set the levels for the sample_label by the other
  if (is.null(x_axis_col_order_vector)){
    Metadata_table[[sample_label]] <- fct_reorder(Metadata_table[[sample_label]], as.numeric(Metadata_table[[grouping_column]]))
    Metadata_table[[sample_label]] <- fct_reorder(Metadata_table[[sample_label]], as.numeric(Metadata_table[[filling_column]]))
    Metadata_table[[sample_label]] <- fct_reorder(Metadata_table[[sample_label]], as.numeric(Metadata_table[[facetting_column]]))
  }else{
    ### In option, we could set an order manually
    Metadata_table[[sample_label]] <- factor(Metadata_table[[sample_label]], levels = x_axis_col_order_vector, ordered = TRUE)
  }
  Metadata <<- Metadata_table
  
}
