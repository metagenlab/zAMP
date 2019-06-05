# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 12.11.18

## Make a table with the number of reads for each sample of the phyloseq object
reads_counts_df = data.table(as(sample_data(physeq), "data.frame"), TotalReads = sample_sums(physeq), keep.rownames = TRUE)
## Rename the first column of this news dataframe -> Sample

setnames(reads_counts_df, "rn", "Sample")
