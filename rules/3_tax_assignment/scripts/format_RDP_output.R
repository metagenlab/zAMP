# Title     : RDP - format output 
# Objective : Format RDP standard "allrank" output to match our standard, with collapsed taxonomy filtered based of confidenc cut-off
# Created by: valentinscherz
# Created on: 03.05.20

## Redirect R output to the log file
   log <- file(snakemake@log[[1]], open="wt")
   sink(log)
   sink(log, type="message")

## Input
rdp_tax <- snakemake@input["RDP_output"]

## Output
simplified_tax <-snakemake@input["formatted_output"]

## Load needed libraries
library(dplyr);packageVersion("dplyr")

## Read data
tax_table <- read.table(file=rdp_tax, sep="\t", stringsAsFactors=FALSE)

## Format table
   ### set rownames
   row.names(tax_table) <- tax_table$V1

   ### Filter unecessary columns
   tax_table[,c("V1", "V2", "V3", "V4", "V5")] <- NULL

## Set a table for confidence score filtering
   ### Keep score columns
   score_table <- tax_table[,seq(from = 3, to = 21, by = 3)]
   ### Rename columns
   colnames(score_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
   ### Create a boolean dataframe for filtering based on confidence
   score_table_bool <- score_table < 0.5


## Set a table with taxa names
   ### Keep names columns
   names_tables <- tax_table[,seq(from = 1, to = 19, by = 3)]
   ### Renames columns
   colnames(names_tables) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

## Filter taxa names based on boolean table
taxa_names_filtered <- replace(x = names_tables, list = score_table_bool, values = NA)

## Filter scores based on boolean table
score_table_filtered <- replace(x = score_table, list = score_table_bool, values = NA)

## Unite taxonomy
taxa_names_filtered_unite <- unite(data = taxa_names_filtered, col = "", sep = ";", na.rm = TRUE, remove = TRUE)

## Add confidence
taxa_names_filtered_unite$confidence <- apply(score_table_filtered, MARGIN = 1, FUN = min, na.rm = TRUE )

## Write formatted table
write.table(taxa_names_filtered_unite, sep = "\t", simplified_tax, row.names = TRUE, quote = FALSE, col.names = FALSE)

