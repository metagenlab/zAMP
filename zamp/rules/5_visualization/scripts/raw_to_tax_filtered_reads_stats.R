# Title     : Raw to tax filtered reads stats
# Objective : Create stats of reads passing processing and taxonomic-based filtering
# Created by: valentinscherz
# Created on: 06.06.19

## Redirect R output to the log file
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    read_filtering <- snakemake@input[["read_filtering"]]
    taxonomic_filtering <- snakemake@input[["taxonomic_filtering"]]
    multi_QC_report_path <- snakemake@input[["multi_QC_report_path"]]

## Ouput
    raw_to_filtered_reads_stats <- snakemake@output[["raw_to_filtered_reads_stats"]]


## Load needed libaries
    library("phyloseq");packageVersion("phyloseq")
    library("data.table");packageVersion("data.table")
    library("dplyr");packageVersion("dplyr")


## Read files
    read_filtering_physeq <- readRDS(read_filtering)
    tax_filtering_physeq <- readRDS(taxonomic_filtering)
    multi_QC_report <- read.table(multi_QC_report_path, header = T)

## Create a table with the number of raw reads and filtered reads
### Raw reads
    multi_QC_report <- multi_QC_report %>% filter(grepl("R1", Sample)) %>% select(c("Sample","FastQC_mqc.generalstats.fastqc.total_sequences")) # keep only the total of raw sequences. Since it is already twice (R1, R2), keep only R1.
    multi_QC_report$Sample <- gsub(x=multi_QC_report$Sample, pattern = "_R1", replacement = "") # Remove the "_R1" in label

### Filtered reads
    filtered_reads_counts_df <- data.table(as(sample_data(read_filtering_physeq), "data.frame"), TotalReads_processing = sample_sums(read_filtering_physeq), keep.rownames = TRUE)
    setnames(filtered_reads_counts_df, "rn", "Sample") # Rename the first column of this news dataframe -> Sample

### Taxonomically filtered reads
    tax_filtered_reads_counts_df <- data.table(as(sample_data(tax_filtering_physeq), "data.frame"), TotalReads_taxonomy = sample_sums(tax_filtering_physeq), keep.rownames = TRUE)
    setnames(tax_filtered_reads_counts_df, "rn", "Sample") # Rename the first column of this news dataframe -> Sample
    tax_filtered_reads_counts_df <- select(tax_filtered_reads_counts_df, c("TotalReads_taxonomy","Sample"))

### Join the columns of interest
    merged_columns <- merge(filtered_reads_counts_df, multi_QC_report, by=c("Sample"))
    merged_columns <- merge(merged_columns, tax_filtered_reads_counts_df, by = c("Sample"))

### Calculate the differences
    merged_columns$reads_processing <- merged_columns$FastQC_mqc.generalstats.fastqc.total_sequences - merged_columns$TotalReads_processing
    merged_columns$tax_filtration <- merged_columns$TotalReads_processing - merged_columns$TotalReads_taxonomy


### Keep three versions of the table,  where the "Count" is consecutively the final number of reads, the reads lost during processing and the reads taxonomically filtered.
#### Keep the final count of reads, remove the two other column
    final_reads_count <- merged_columns %>% rename("Count" = "TotalReads_taxonomy")
    final_reads_count$Reads <- "Maintained reads"
    final_reads_count[ ,c("reads_processing", "tax_filtration")] <- list(NULL)

#### Keep the reads filtered during processing remove the two other column
    processing_reads_count <- merged_columns %>% rename("Count" = "reads_processing")
    processing_reads_count$Reads <- "Reads processing"
    processing_reads_count[ ,c("tax_filtration", "TotalReads_taxonomy")] <- list(NULL)

#### Keep the reads removed by taxonomic filtering remove the two other column
    tax_filtered_reads_count <- merged_columns %>% rename("Count" = "tax_filtration")
    tax_filtered_reads_count$Reads <- "Tax filtered reads"
    tax_filtered_reads_count[ ,c("reads_processing", "TotalReads_taxonomy")] <- list(NULL)

#### Bind the rows of the two just prepared tables
    melted_reads_count <- rbind(final_reads_count, tax_filtered_reads_count)
    melted_reads_count <- rbind(melted_reads_count, processing_reads_count)


# Write this table
    write.table(x = melted_reads_count, file = raw_to_filtered_reads_stats, sep = "\t", col.names = NA, row.names = TRUE)
