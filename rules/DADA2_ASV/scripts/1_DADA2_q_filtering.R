## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
fnFs <- snakemake@input[["R1_list"]]
fnRs <- snakemake@input[["R2_list"]]

## Output
q_score_filtered_F <- snakemake@output[["q_score_filtered_F"]]
q_score_filtered_R <- snakemake@output[["q_score_filtered_R"]]
q_filtering_stats_path <- snakemake@output[["q_filtering_stats"]]

## Parameters
F_length <- snakemake@params[["F_reads_length_trim"]]
R_length <- snakemake@params[["R_reads_length_trim"]]
F_extected_error <- snakemake@params[["F_reads_extected_error"]]
R_extected_error <- snakemake@params[["R_reads_extected_error"]]
sample_name <- snakemake@params[["sample_name"]]


## Load needed libraries
library(dada2);packageVersion("dada2")

## Filter and trim
# Reads are filtered based on the number of errors expected for the read (integration of the qscore and the length of the read). All reads with uncalled nucleotide (N) are removed too. Remaining phiX reads will be removed too. Finally, reads are cut at 280 bp for the forward and 255 for the reverse. This looks like an acceptable compromise between removing the end of the read where there would be more expected errors and allowing sufficient ovelap for merging of paired-ends reads. We cut asymetrically because the quality of the read is asymetrical.

# quality_plot <- plotQualityProfile(fnFs[1:2])

# Filter and trim
filtering_stats <- filterAndTrim(fwd = fnFs, filt = q_score_filtered_F, rev = fnRs, filt.rev = q_score_filtered_R, truncLen=c(F_length,R_length), maxN=0, maxEE=c(F_extected_error,R_extected_error), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE)
filtering_stats <- as.data.frame(filtering_stats)
filtering_stats$Sample <- sample_name

# Save the stats for this sample in a R object
saveRDS(filtering_stats, file = q_filtering_stats_path)
