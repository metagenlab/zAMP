## Adapted from https://bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/ClassifySequences.pdf
# Title     : Decipher - tax prep
# Objective : Prep taxonomy for decipher
# Created by: valentinscherz
# Created on: 07.08.19

## Redirect R output to the log file
   log <- file(snakemake@log[[1]], open="wt")
   sink(log)
   sink(log, type="message")

## Input
    ref_tax <- snakemake@input[["ref_tax"]]
    ref_seqs <- snakemake@input[["ref_seqs"]]

## Output
    decipher_seqs <- snakemake@output[["decipher_seqs"]]

## Load needed libraries
    library(DECIPHER);packageVersion("DECIPHER")
    library(dplyr);packageVersion("dplyr")
    library(tidyr);packageVersion("tidyr")
    library(stringr);packageVersion("stringr")

################### Prep taxonomy tree ########################
tax_table <- read.table(file=ref_tax, sep="\t", stringsAsFactors=FALSE)
tax_table

################### Prep fasta table ########################
split_tax_table<-tax_table %>% separate(V2, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

## Format data
split_tax_table$Root <- "Root"
split_tax_table_f <- split_tax_table %>% unite(tax, c("Root", "Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = "; ", remove = TRUE)
split_tax_table_f <- split_tax_table_f %>% unite(taxonomy, c("V1", "tax"), sep = " ",remove = FALSE) %>% select(-tax)
split_tax_table_f <- split_tax_table_f[,c(2,1)]

## Define a function taken from phylotools https://github.com/cran/phylotools/blob/master/R/rename.fasta.R

rename.fasta <- function(infile = NULL, ref_table,
                         outfile = "renamed.fasta"){
  #fasta <- read.fasta(infile)
  fastaFile <- readDNAStringSet(infile)
  seq.name = names(fastaFile)
  seq.text = paste(fastaFile)
  fasta <- data.frame(seq.name, seq.text)


  ## Convert the input ref to dataframe
  ## usually the file was obtained by
  ## save the result of get.names.fasta to a csv file.
  ## edit the csv file, and provide the name to be replaced.

  colnames(ref_table)  <- c("old_name", "new_name")
  res <- merge(x = fasta, y = ref_table, by.x = "seq.name",
               by.y = "old_name", all.x = TRUE) ### Order of the sequence will change because of merge
  rename <- rep(NA, nrow(res))

  ### if the sequence was not found in the ref_table,
  ### keep the old name

  for(i in 1:nrow(res)){
    rename[i] <- ifelse(is.na(res$new_name[i]),
                        paste("old_name", "_",as.character(res$old_name[i]), sep = ""),
                        as.character(res$new_name[i]))
  }
  writeLines(paste(">", rename, "\n", as.character(res$seq.text), sep = ""), outfile )
  cat(paste(outfile, "has been saved to ", getwd(), "\n" ))
}

## Write data
 rename.fasta(infile = ref_seqs, ref_table = split_tax_table_f, outfile = decipher_seqs)
