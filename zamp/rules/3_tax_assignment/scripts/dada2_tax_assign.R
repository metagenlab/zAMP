log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dada2); packageVersion("dada2")

seq_df <- read.table(snakemake@input[["seqs"]], sep="\t", header=FALSE)
colnames(seq_df) <- c("seq_id", "sequence")

trainset <- snakemake@input[["trainset"]]
tt <- assignTaxonomy(seq_df$sequence, trainset, verbose=TRUE, multithread=FALSE, outputBootstraps = TRUE)
tax_df <- as.data.frame(tt, stringsAsFactors=FALSE)
colnames(tax_df) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "confidence")
tax_df <- cbind(seq_id=seq_df$seq_id, tax_df)
tax_df$tax <- apply(tax_df[, c("kingdom", "phylum", "class", "order", "family", "genus", "species")], 1, paste, collapse=";")
tax_df <- tax_df[, c("seq_id", "tax", "confidence")]
write.table(tax_df, file = snakemake@output[[1]], sep = "\t", row.names = FALSE, quote = FALSE)