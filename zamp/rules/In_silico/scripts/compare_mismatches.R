# Title     : Mismatch histogram
# Objective : Plot an histogram of distribution of mismatches for all samples
# Created by: valentinscherz
# Created on: 06.11.2019

## Redirect R output
  log <- file(snakemake@log[[1]], open="wt")
  sink(log)
  sink(log, type="message")

## Input
  mismatch_tables_path <- snakemake@input[["mismatch_tables_path"]]

## Output
  missmatch_plot <- snakemake@output[["missmatch_plot"]]
  merged_mismatch_table <- snakemake@output[["merged_mismatch_table"]]


# Load needed libraries
  library(ggplot2); packageVersion("ggplot2")


# Merge data from multiple runs (if necessary)
  if (length(mismatch_tables_path) == 1){
    print("Unique run, no merging of seq_tabl")
    st.all <- readRDS(seq_tab)
  }else{
    print("Multiple run, merging")
    st.all <- do.call("rbind", lapply(mismatch_tables_path, read.table, header = TRUE, sep = "\t", row.names = 1))
  }

  p <-ggplot(st.all, aes(x=Mismatches)) +
    geom_histogram(color="black", fill="white", binwidth=1) +
    theme_bw() +
    scale_y_continuous(limits = c(0,100000)) +
      scale_x_continuous(limits = c(0,10))


  ggsave(missmatch_plot, plot = p)


  write.table(st.all, merged_mismatch_table, sep = "\t", row.names =  FALSE)
