# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 28.12.18

# FROM  https://benjjneb.github.io/dada2/bigdata_paired.html

## Redirect R output
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")


## Input
q_score_filtered_Fs <- snakemake@input[["q_score_filtered_Fs"]]
q_score_filtered_Rs <- snakemake@input[["q_score_filtered_Rs"]]
error_profile_F <- snakemake@input[["error_profile_F"]]
error_profile_R <- snakemake@input[["error_profile_R"]]

## Output
#filtering_stats <- snakemake@output[["filtering_stats"]]
seq_tab <- snakemake@output[["seq_tab"]]

## Load needed libraries
library(dada2); packageVersion("dada2")

# File parsing
sample.names <- sapply(strsplit(basename(q_score_filtered_Fs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(q_score_filtered_Rs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(q_score_filtered_Fs) <- sample.names
names(q_score_filtered_Rs) <- sample.names
set.seed(100)

# Read error rates
errF <- readRDS(error_profile_F)
errR <- readRDS(error_profile_R)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(q_score_filtered_Fs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE, verbose = 1, pool = FALSE, selfConsist = TRUE)
    derepR <- derepFastq(q_score_filtered_Rs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE, verbose = 1, pool = FALSE, selfConsist = TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(object = seqtab, file = seq_tab) # CHANGE ME to where you want sequence table saved
