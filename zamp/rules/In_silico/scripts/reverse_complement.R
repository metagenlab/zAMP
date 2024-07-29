# Title     : Reverse complement sequences in wrong direction
# Objective : SimulatePCR extract in both direction. Flip the ones which needs
# Created by: valentinscherz
# Created on: 21.01.2021

## Redirect R output
  log <- file(snakemake@log[[1]], open="wt")
  sink(log)
  sink(log, type="message")

## Input
  input <- snakemake@input[["simple_fasta"]]

## Output
  output <- snakemake@output[["fasta_with_rev"]]

## Read data
library(Biostrings)
seqs <- readDNAStringSet(input)

## Split forward
forw <- seqs[grep("primer|F primer|R", names(seqs), fixed=TRUE)]
## Split reverse
reve <- seqs[grep("primer|R primer|F", names(seqs), fixed=TRUE)]
## Reverse complement reverse
reve_com <- reverseComplement(reve)
## Bind them
binded <-  c(forw, reve_com)
## Write output
writeXStringSet(binded, format = "fasta", output)
