# Title     : Transpose tables
# Objective : Reformat transposing count and metadata tables
# Created by: valentinscherz
# Created on: 06.06.19

## Inspired from https://gist.github.com/erictleung/eda33bceceebb190eb9a44c93a077d32
## Redirect R output
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    count_table <- snakemake@input[["count_table"]]
    print(count_table)
    meta <- snakemake@input[["meta"]]

## Output
    transposed_table <- snakemake@output[["transposed_table"]]
    merged_meta <- snakemake@output[["merged_meta"]]

## Load the phyloseq phyloseq_object
    table <- read.table(count_table,  sep = "\t", header = TRUE)
    meta <- read.delim(meta,  sep = "\t", header = TRUE, na.strings = "NA")

## Tranpose table
    table_t <- t(table)

## Bind tables
    binded_tables <- cbind(meta, table_t)

### Write this table in a tsv file since it is ground for coming analysis and slow to compute
    write.table(x = table_t, file = transposed_table, append = F, quote = F, sep = "\t", eol = "\n", row.names = F, col.names = T )
    write.table(x = binded_tables, file = merged_meta, append = F, quote = F, sep = "\t", eol = "\n", row.names = F, col.names = T )
