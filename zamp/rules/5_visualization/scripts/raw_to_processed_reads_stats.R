# Title     : Raw to processed reads stats
# Objective : Create stats  of reads passing processing
# Created by: valentinscherz
# Created on: 06.06.19

## Redirect R output to the log file
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    phyloseq_object <- snakemake@input[["phyloseq_object"]]
    multi_QC_report_path <- snakemake@input[["multi_QC_report_path"]]
    
## Ouput
    raw_to_filtered_reads_stats <- snakemake@output[["raw_to_filtered_reads_stats"]]

    ## Load needed libaries
    library("phyloseq");packageVersion("phyloseq")
    library("data.table");packageVersion("data.table")
    library("dplyr");packageVersion("dplyr")

## Load MultiQC report
    multi_QC_report <- read.table(multi_QC_report_path, header = T)

## Load the phyloseq object
    phyloseq_obj <- readRDS(phyloseq_object)

## Create a table with the number of raw reads and filtered reads
### Filtered reads
    reads_counts_df <- data.table(as(sample_data(phyloseq_obj), "data.frame"), TotalReads = sample_sums(phyloseq_obj), keep.rownames = TRUE)
    setnames(reads_counts_df, "rn", "sample") # Rename the first column of this news dataframe -> sample

### Raw reads
    setnames(multi_QC_report, "Sample", "sample") # rename multiQC Sample column to sample 
    multi_QC_report <- multi_QC_report %>% filter(grepl("R1|single", sample)) %>% select(c("sample","FastQC_mqc.generalstats.fastqc.total_sequences")) # keep only the total of raw sequences. Since it is already twice (R1, R2), keep only R1.
    multi_QC_report$sample <- gsub(x=multi_QC_report$sample, pattern = "_R1", replacement = "") # Remove the "_R1"
    #or
    multi_QC_report$sample <- gsub(x=multi_QC_report$sample, pattern = "_single", replacement = "") # Remove the "_single"
### Raw reads Join the two dataframes, the one with the filtered reads from phyloseq and the one with the raw reads.
    reads_counts_df_with_raw <- merge(multi_QC_report, reads_counts_df, by=c("sample"))

### Raw reads Calculate the difference of reads between the raw and the filtered reads.
    reads_counts_df_with_raw <- mutate(reads_counts_df_with_raw, filtered = FastQC_mqc.generalstats.fastqc.total_sequences - TotalReads)

### Raw reads Keep two version of the table, one with the filtered reads in "Count" and one with the difference with the raw reads in the same count column.
    reads_counts_df_raw_count_only <- reads_counts_df_with_raw %>% select(-"filtered") %>% rename("Count" = "TotalReads")
    reads_counts_df_raw_count_only$Reads <- "Maintained reads"
    reads_counts_df_raw_filtered_only <- reads_counts_df_with_raw %>% select(-"TotalReads") %>% rename("Count" = "filtered")
    reads_counts_df_raw_filtered_only$Reads <- "Reads processing"

### Raw reads Bind the rows of the two just prepared tables
    melted_reads_count <- rbind(reads_counts_df_raw_count_only, reads_counts_df_raw_filtered_only)

### Raw reads Write this table
    write.table(x = melted_reads_count, file = raw_to_filtered_reads_stats, sep = "\t", col.names = NA, row.names = TRUE)
