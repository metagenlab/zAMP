# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 26.10.18

## Inspired from https://gist.github.com/erictleung/eda33bceceebb190eb9a44c93a077d32
## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]

## Output
phyloseq_melted_table <- snakemake@output[["phyloseq_melted_table"]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")

## Load the phyloseq phyloseq_object
load(phyloseq_object)

### Melt the physeq objet into a dataframe with one row per feature.id and per sample, needed later
physeq_df <- psmelt(phyloseq_obj)

### Write this table in a tsv file since it is ground for coming analysis and slow to compute
write.table(x = physeq_df, file = phyloseq_melted_table, append = F, quote = F, sep = "\t", eol = "\n", row.names = F, col.names = T )
