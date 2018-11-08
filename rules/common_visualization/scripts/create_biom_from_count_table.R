# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 08.11.18


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


#otu<-t(as(otu_table(asv_tab),"matrix")) # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE
asv_tab<-read.table(file = count_table)
otu<-as(otu_table(asv_tab, taxa_are_rows = TRUE),"matrix")
otu_biom<-make_biom(data=otu)
write_biom(x = otu_biom, biom_file = biom_count)
