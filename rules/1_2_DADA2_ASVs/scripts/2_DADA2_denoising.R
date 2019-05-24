## Redirect R output
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

## Input
q_score_filtered_Fs <- snakemake@input[["q_score_filtered_Fs"]]
q_score_filtered_Fs
q_score_filtered_Rs <- snakemake@input[["q_score_filtered_Rs"]]
q_score_filtered_Rs
q_score_filtering_stats <- snakemake@input[["q_score_filtering_stats"]]

## Output
filtering_stats <- snakemake@output[["filtering_stats"]]
dna_sequences.fasta <- snakemake@output[["dna_sequences"]]
count_table.txt <- snakemake@output[["count_table"]]

## Parameters
merged_min_length <- snakemake@params[["merged_min_length"]]
merged_max_length  <- snakemake@params[["merged_max_length"]]

## Load needed libraries
library(dplyr);packageVersion("dplyr")
library(dada2);packageVersion("dada2")
library(phyloseq);packageVersion("phyloseq")


## Import samples list in R
# The samples names is inferred from the name of the .fastq
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(q_score_filtered_Fs), "_"), `[`, 1)


## Learn the Error Rates
# That the magic of DADA2... Where is learn the probabily of errors based of the reads so that it can then correct them.
errF <- learnErrors(q_score_filtered_Fs, multithread=TRUE)
errR <- learnErrors(q_score_filtered_Rs, multithread=TRUE)

# plotErrors(errF, nominalQ=TRUE)

## Dereplication
# Instead of having all, it keeps one represent of each single reads and count how many time it is in which sample
derepFs <- derepFastq(q_score_filtered_Fs, verbose=TRUE)
derepRs <- derepFastq(q_score_filtered_Rs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## Sample inference
# Here errors are corrected. We apply correction for each sample individually (pool = FALSE), to limit batch effect.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, verbose = 1, pool = FALSE, selfConsist = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, verbose = 1, pool = FALSE, selfConsist = TRUE )

dadaFs[[1]]

## Merge paired reads
# Once forward and reverse reads are corrected, they are merged together.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

## Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

## That's the little added trick, the reason why we are using this script and not the one in Qiime2. Indeed we are here keeping only sequences between 390 and 500 bp of length after merging. This corresponds to the expected length of the V3V4 region of the 16S rRNA gene.
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(merged_min_length, merged_max_length)]

table(nchar(getSequences(seqtab2)))

## Remove chimeras
# It is expected that a fair amount of reads can be chimera generated during PCR amplifications. This is what we filter here.
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


## Track reads through the pipeline
# We want to be able to follow the number of reads filtered at each step. We are here generating a table recording these numbers.

q_score_table <- read.csv(file = q_score_filtering_stats, header = TRUE, sep = "\t", row.names = "X")

getN <- function(x) sum(getUniques(x))
track <- cbind(q_score_table, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
track[["nonchim/input_%"]] <- (track$nonchim/track$input)*100
write.table(x = track, file = filtering_stats, sep = "\t", col.names = NA, row.names = TRUE)


## Export reads and count
# We are writing in files the product of this DADA2 process. These are one .fasta file contanining the dereplicated, errors corrected, paired-end merged representative sequences and one .txt file indicating the prevalence of sequencne in each sample (this is the result of dereplication).


# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

  # making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, dna_sequences.fasta)

  # count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, count_table.txt , sep="\t", quote=F)


