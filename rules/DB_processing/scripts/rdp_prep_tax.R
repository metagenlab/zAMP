# Title     : RDP - tax assign
# Objective : Prep taxonomy for RDP
# Created by: valentinscherz
# Created on: 02.05.20


## Redirect R output to the log file
   log <- file(snakemake@log[[1]], open="wt")
   sink(log)
   sink(log, type="message")

## Input
    ref_seqs <- snakemake@input[["ref_seqs"]]
    ref_tax <- snakemake@input[["ref_tax"]]

## Output
    formatted_table <- snakemake@output[["formatted_table"]]

## Load needed libraries
    library(dplyr);packageVersion("dplyr")
    library(tidyr);packageVersion("tidyr")

## Read data
tax_table <- read.table(file=ref_tax, sep="\t", stringsAsFactors=FALSE)
split_tax_table<-tax_table %>% separate(V2, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

## Format data
King_to_species <- split_tax_table
colnames(King_to_species)[1] <- "Seq_ID"

## Write formatted table
write.table(King_to_species, sep = "\t", formatted_table, row.names = FALSE, quote = FALSE)

