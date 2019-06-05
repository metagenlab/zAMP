# Title     : Melt the physeq
# Objective : Melt the physeq object into a long table easily plottable
# Created by: valentinscherz
# Created on: 28.05.2019

## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[[1]]

## Output
phyloseq_melted_table <- snakemake@output[[1]]



## Load needed libraries
library(phyloseq);packageVersion("phyloseq")



## Import the physeq object
phyloseq_obj <- readRDS(phyloseq_object)

### Melt the physeq objet into a dataframe with one row per feature.id and per sample
physeq_df <- psmelt(phyloseq_obj)

### Write this table in a tsv format
write.table(x = physeq_df, file = phyloseq_melted_table, append = F, quote = F, sep = "\t", eol = "\n", row.names = F, col.names = T )
