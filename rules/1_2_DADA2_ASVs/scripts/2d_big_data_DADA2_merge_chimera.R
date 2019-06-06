# Title     : DADA2 - merge runs, filter chimera
# Objective : Merge runs, filter chimera and based on length
# Created by: valentinscherz
# Created on: 06.06.19
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

## Redirect R output
  log <- file(snakemake@log[[1]], open="wt")
  sink(log)
  sink(log, type="message")

## Input
  seq_tab <- snakemake@input[["seq_tab"]]

## Output
  no_chim <- snakemake@output[["no_chim"]]
  length_filtered <- snakemake@output[["length_filtered"]]
  hashed_sequences <- snakemake@output[["hashed_sequences"]]
  hashed_count_table <- snakemake@output[["hashed_count_table"]]
  rep_seqs <- snakemake@output[["rep_seqs"]]
  count_table <- snakemake@output[["count_table"]]
  length_histo <- snakemake@output[["length_histo"]]

## Parameters
  merged_min_length <- snakemake@params[["merged_min_length"]]
  merged_max_length  <- snakemake@params[["merged_max_length"]]

## Load needed libraries
  library(dada2); packageVersion("dada2")

## Merge data from multiple runs (if necessary)
   if (length(seq_tab) == 1){
	print("Unique RUN, no merging of seq_tabl")
	st.all <- readRDS(seq_tab)
   }else{
	print("Multiple RUN, merging")
	st.all <- do.call("mergeSequenceTables", lapply(seq_tab, readRDS))
   }

## Remove chimeras
  print("Filter chimera")
  seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=snakemake@threads, verbose=TRUE)
  print("Chimera filtered")

## Sequences length inspection and filtration
  ### Inspect distribution of sequence lengths
  seq_length_table <- table(nchar(getSequences(seqtab)))
  png(length_histo)
  hist(seq_length_table,
      main="Unfiltered sequences length",
      xlab="Length")
  dev.off()

  ## Filter based on length
    seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(merged_min_length, merged_max_length)]
  ## Inspect distribution of sequence lengths after filtration
    table(nchar(getSequences(seqtab2)))

## Export reads and count
### We are writing in files the product of this DADA2 process. These are one .fasta file contanining the dereplicated, errors corrected, paired-end merged representative sequences and one .txt file indicating the prevalence of sequencne in each sample (this is the result of dereplication).
  ### Still hashed sequences
    asv_fasta_hashed <- c(rbind(asv_headers, asv_seqs))
    write(asv_fasta_hashed, hashed_sequences)

  ### Still hashed count table
    print("Create hashed count table")
    asv_tab_hashed <- t(seqtab2)
    #row.names(asv_tab) <- sub(">", "", asv_headers)
    write.table(asv_tab_hashed, hashed_count_table , sep="\t", quote=F)

  ### Give our seq headers more manageable names (ASV_1, ASV_2...)
    asv_seqs <- colnames(seqtab2)
    asv_headers <- vector(dim(seqtab2)[2], mode="character")
    for (i in 1:dim(seqtab2)[2]) {
      asv_headers[i] <- paste(">ASV", i, sep="_")
    }
  ## Make and write out a fasta of our final ASV seqs:
    asv_fasta <- c(rbind(asv_headers, asv_seqs))
    write(asv_fasta, rep_seqs)

  ## Count table:
    print("Create count table")
    asv_tab <- t(seqtab2)
    row.names(asv_tab) <- sub(">", "", asv_headers)
    write.table(asv_tab, count_table , sep="\t", quote=F)

  ## Write sequences objects in .rds for use in statistics
    ## Before length filtration
    saveRDS(seqtab, no_chim)
    ## After length filtration
    saveRDS(seqtab2, length_filtered)
