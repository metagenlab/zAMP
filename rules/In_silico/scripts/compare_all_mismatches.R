# Title     : Mismatch histogram
# Objective : Plot an histogram of distribution of mismatches for all samples
# Created by: valentinscherz
# Created on: 06.11.2019

## Redirect R output
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

## Input
mismatch_tables_path <- snakemake@input[["mismatch_tables_path"]]

## Output
missmatch_plot <- snakemake@output[["missmatch_plot"]]
merged_mismatch_table <- snakemake@output[["merged_mismatch_table"]]


# Load needed libraries
library(ggplot2); packageVersion("ggplot2")



save.image(file= paste0(getwd(), "/compare_all.RData"))



# Merge data from multiple runs (if necessary)
# Merge data from multiple runs (if necessary)
if (length(mismatch_tables_path) == 1){
  print("Unique RUN, no merging of seq_tabl")
  st.all <- readRDS(seq_tab)
}else{
  print("Multiple RUN, merging")
  st.all <- do.call("rbind", lapply(mismatch_tables_path, read.table, header = TRUE, sep = "\t"))
}

p <- ggplot(st.all, aes(x=Mismatches, fill=RUN, weight = Counts)) +
  geom_histogram(aes(y=..density..), color="black", binwidth=1, position = position_dodge2(padding = 0.3, preserve = "single")) +
  theme_bw()


ggsave(missmatch_plot, plot = p)


write.table(st.all, merged_mismatch_table, sep = "\t", row.names =  FALSE)

