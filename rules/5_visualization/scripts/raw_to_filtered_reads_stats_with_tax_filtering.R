# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 29.11.18

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
### Filtered reads
filtered_reads_counts_df <- data.table(as(sample_data(read_filtering_physeq), "data.frame"), TotalReads_processing = sample_sums(read_filtering_physeq), keep.rownames = TRUE)
setnames(filtered_reads_counts_df, "rn", "Sample") # Rename the first column of this news dataframe -> Sample

### Taxonomically filtered reads
tax_filtered_reads_counts_df <- data.table(as(sample_data(tax_filtering_physeq), "data.frame"), TotalReads_taxonomy = sample_sums(tax_filtering_physeq), keep.rownames = TRUE)
setnames(tax_filtered_reads_counts_df, "rn", "Sample") # Rename the first column of this news dataframe -> Sample

### Raw reads
multi_QC_report <- multi_QC_report %>% filter(grepl("R1", Sample)) %>% select(c("Sample","FastQC_mqc.generalstats.fastqc.total_sequences")) # keep only the total of raw sequences. Since it is already twice (R1, R2), keep only R1.
multi_QC_report$Sample <- gsub(x=multi_QC_report$Sample, pattern = "_R1", replacement = "") # Remove the "_R1"


## Calculate the number of reads lost during processing
### Join the two dataframes, the one with the filtered reads from phyloseq and the one with the raw reads.
reads_filtered_during_processing <- merge(multi_QC_report, filtered_reads_counts_df, by=c("Sample"))

#### Calculate the difference of reads between the raw and the filtered reads.
reads_filtered_during_processing <- mutate(reads_filtered_during_processing, filtered = FastQC_mqc.generalstats.fastqc.total_sequences - TotalReads_processing)

## Calculate the number of reads filtered by the taxonomic filtering.
### Join with the previous table
reads_tax_filtered <- merge(reads_filtered_during_processing, tax_filtered_reads_counts_df, by=c("Sample"))

### Calculate the difference
reads_tax_filtered <- mutate(reads_tax_filtered, tax_filtered = TotalReads_processing - TotalReads_taxonomy)


# Keep three versions of the table, one per count
reads_counts_df_raw_count_only <- reads_tax_filtered %>% select(-"filtered") %>% rename("Count" = "TotalReads")
reads_counts_df_raw_count_only$Reads <- "Maintained reads"
reads_counts_df_raw_filtered_only <- reads_tax_filtered %>% select(-"TotalReads") %>% rename("Count" = "filtered")
reads_counts_df_raw_filtered_only$Reads <- "Read processing"
reads_counts_df_taxonomy <- reads_tax_filtered %>% select(-"tax_filtered") %>% rename("Count" = "tax_filtered_reads_counts_df")
reads_counts_df_taxonomy$Reads <- "Taxonomy based filtering"


# Bind the rows of the two just prepared tables
melted_reads_count <- rbind(reads_counts_df_raw_count_only, reads_counts_df_raw_filtered_only)
melted_reads_count <- rbind(melted_reads_count, reads_counts_df_taxonomy)


# Write this table
write.table(x = melted_reads_count, file = raw_to_filtered_reads_stats, sep = "\t", col.names = NA, row.names = TRUE)
