# Title     : Vsearch count table
# Objective : Format a count table from vsearch output
# Created by: valentinscherz
# Created on: 27.11.2019

## Redirect R output
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Load needed library
    library(dplyr);packageVersion("dplyr")

## Input
    uc_path <- snakemake@input[["uc"]]

## Output
    count_table_path <- snakemake@output[["count_table"]]


uc_table <- read.table(uc_path, sep = "\t", header = FALSE)


count <- uc_table %>% filter(V1 == "C") %>% select(V9, V3)

write.table(x = count, file = count_table_path, sep="\t", col.names = FALSE, row.names = TRUE)
