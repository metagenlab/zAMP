# Title     : Create biom from count table
# Objective : Format count table into biom
# Created by: valentinscherz
# Created on: 06.06.19


## Redirect R output
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    count_table <- snakemake@input[["count_table"]]

## Output
    biom_count <- snakemake@output[["biom_count"]]

## Library
    library(biomformat);packageVersion("biomformat")
    library(phyloseq);packageVersion("phyloseq")

## Reformat
    asv_tab <- read.table(file = count_table)
    otu <- as(otu_table(asv_tab, taxa_are_rows = TRUE),"matrix")
    otu_biom <- make_biom(data=otu)

## Write
    write_biom(x = otu_biom, biom_file = biom_count)
