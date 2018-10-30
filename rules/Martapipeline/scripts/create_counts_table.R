# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 29.10.18

## Inspired from https://gist.github.com/erictleung/eda33bceceebb190eb9a44c93a077d32
## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Load needed library
library(dplyr);packageVersion("dplyr")

## Input
full_taxonomy_table <- snakemake@input[["full_taxonomy_table"]]
count_table_samples <- snakemake@input[["count_table_samples"]]

print(full_taxonomy_table)

## Output
count_table <- snakemake@output[["count_table"]]

## Read full_full_taxonomy_table

tax_table <- read.table(file = full_taxonomy_table, header = TRUE, sep = "\t")

tax_table$OTU_ID <- as.factor(tax_table$OTU_ID)

otus_table<-as.data.frame(array(dim=c(1,3)))
colnames(otus_table) <- c("Sample", "V1", "V2")

for (xx in count_table_samples){
    sample <- gsub("_count_table.txt", "", basename(xx))
    table <- read.table(file = xx, sep="\t", as.is=T, header=F)
    table <- cbind("Sample"=sample, table)
    otus_table <- rbind(otus_table, table)
    }

otus_table <- otus_table[-1,]
colnames(otus_table) <- c("Sample", "OTU_ID", "counts")

# Add the taxonomy and save it for downstream analyses & plotting
otus_tax_table <- left_join(otus_table, tax_table)
write.table(otus_tax_table, file.path(count_table), row.names=F, col.names=T, sep="\t")
