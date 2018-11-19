## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
fnFs <- snakemake@input[["R1_list"]]
fnRs <- snakemake@input[["R2_list"]]

## Output
q_score_filtered_Fs <- snakemake@output[["q_score_filtered_Fs"]]
q_score_filtered_Rs <- snakemake@output[["q_score_filtered_Rs"]]
q_score_filtered_R_object <- snakemake@output[["filtering_stats"]]

## Parameters
F_length <- snakemake@params[["F_reads_length_trim"]]
R_length <- snakemake@params[["R_reads_length_trim"]]
F_extected_error <- snakemake@params[["F_extected_error"]]
R_extected_error <- snakemake@params[["R_extected_error"]]


## Load needed libraries
library(dada2);packageVersion("dada2")

## Filter and trim
# Reads are filtered based on the number of errors expected for the read (integration of the qscore and the length of the read). All reads with uncalled nucleotide (N) are removed too. Remaining phiX reads will be removed too. Finally, reads are cut at 280 bp for the forward and 255 for the reverse. This looks like an acceptable compromise between removing the end of the read where there would be more expected errors and allowing sufficient ovelap for merging of paired-ends reads. We cut asymetrically because the quality of the read is asymetrical.

# quality_plot <- plotQualityProfile(fnFs[1:2])

# Filter and trim
filtering_stats <- filterAndTrim(fwd = fnFs, filt = q_score_filtered_Fs, rev = fnRs, filt.rev = q_score_filtered_Rs, truncLen=c(F_length,R_length), maxN=0, maxEE=c(F_extected_error,R_extected_error), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

# Rename the rownames of filtering_stats by a sample name without suffix
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

rownames(filtering_stats) <- sample.names

write.table(x = filtering_stats, file = file.path(q_score_filtered_R_object), sep = "\t", col.names = NA, row.names = TRUE)

