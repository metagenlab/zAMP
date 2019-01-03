# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 28.12.18
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
seq_tab <- snakemake@input[["seq_tab"]]

## Output
with_chim <- snakemake@output[["with_chim"]]
no_chim <- snakemake@output[["no_chim"]]
length_filtered <- snakemake@output[["length_filtered"]]
renamed <- snakemake@output[["renamed"]]
count_table.txt <- snakemake@output[["count_table"]]
#filtering_stats <- snakemake@output[["filtering_stats"]]

## Parameters
merged_min_length <- snakemake@params[["merged_min_length"]]
merged_max_length  <- snakemake@params[["merged_max_length"]]

## Load needed libraries
library(dada2); packageVersion("dada2")

# Merge multiple runs (if necessary)
files <- seq_tab
st.all <- do.call("mergeSequenceTables", lapply(files, readRDS))
write(st.all, with_chim)

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
write(seqtab, no_chim)


# Sequences length inspection and filtration
## Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## That's the little added trick, the reason why we are using this script and not the one in Qiime2. Indeed we are here keeping only sequences between 390 and 500 bp of length after merging. This corresponds to the expected length of the V3V4 region of the 16S rRNA gene.
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(merged_min_length, merged_max_length)]
table(nchar(getSequences(seqtab2)))

## Before renaming
write(seqtab2, length_filtered)

# Export reads and count
#### We are writing in files the product of this DADA2 process. These are one .fasta file contanining the dereplicated, errors corrected, paired-end merged representative sequences and one .txt file indicating the prevalence of sequencne in each sample (this is the result of dereplication).

## giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab2)
asv_headers <- vector(dim(seqtab2)[2], mode="character")

for (i in 1:dim(seqtab2)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

## making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, renamed)

## count table:
asv_tab <- t(seqtab2)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, count_table.txt , sep="\t", quote=F)


