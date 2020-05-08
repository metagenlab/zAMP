## Adapted from https://bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/ClassifySequences.pdf
# Title     : Decipher - tax assign
# Objective : Assign tax with decipher
# Created by: valentinscherz
# Created on: 07.08.19

## Redirect R output to the log file
   log <- file(snakemake@log[[1]], open="wt")
   sink(log)
   sink(log, type="message")

## Input
    trained_tax_path <- snakemake@input[["trained_tax"]]
    seq_path <- snakemake@input[["seq"]]	

## Output
    tax_path <- snakemake@output[["tax"]]
    tax_plot_path <- snakemake@output[["tax_plot"]]

## Parameters
    cpu <- snakemake@threads

## Load needed libraries
    library(DECIPHER);packageVersion("DECIPHER")
    library(dplyr);packageVersion("dplyr")
    library(tidyr);packageVersion("tidyr")

# load a training set object (trainingSet)
# see http://DECIPHER.codes/Downloads.html
train_set <- readRDS(trained_tax_path )
query_seqs <- readDNAStringSet(seq_path)


# classify the sequences
ids <- IdTaxa(query_seqs,
   train_set,
   strand="both", # or "top" if same as trainingSet
   threshold=60, # 60 (very high) or 50 (high)
   processors=cpu) # use all available processors

# look at the results
 pdf(file = tax_plot_path)
 plot(ids)
 dev.off()


## Format taxonomy and write it
taxid <- data.frame(t(sapply(ids, function(x){
  taxa <- c(unlist(x$taxon)[2:8], last(x$confidence))
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
  })))
taxid <- unite(data = taxid, 1:7, col = "taxonomy", sep = ";", remove = TRUE)
write.table(x = taxid, tax_path, row.names = TRUE, sep = "\t", col.names = FALSE, quote = FALSE)


