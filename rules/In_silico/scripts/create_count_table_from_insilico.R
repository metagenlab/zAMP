# Title     : Vsearch count table
# Objective : Format a count table from vsearch output
# Created by: valentinscherz
# Created on: 06.06.1p

## Redirect R output
    #log <- file(snakemake@log[[1]], open="wt")
    #sink(log)
    #sink(log, type="message")

## Load needed library
    library(dplyr);packageVersion("dplyr")
    library(reshape2);packageVersion("reshape2")
    library(magrittr);packageVersion("magrittr")

## Input
    count_table_samples <- snakemake@input[["count_table_samples"]]

## Output
    count_table <- snakemake@output[["count_table"]]

## Reformat
    otus_table<-as.data.frame(array(dim=c(1,3)))
    colnames(otus_table) <- c("Sample", "V1", "V2")


## Reformat
otus_table<-as.data.frame(array(dim=c(1,3)))
colnames(otus_table) <- c("Sample", "V1", "V2")


for (xx in count_table_samples){
  print(paste("Processing", xx))
  sample <- gsub("_count_table.txt", "", basename(xx))
  print(paste("Sample:", sample))
  table <- as.data.frame(array(dim=c(1,2)))
  table$V1 <- "No_amp"
  table$V2 <- 0
  try(table <-read.table(file = xx, sep="\t", as.is=T, header=F))
  table <- cbind("Sample"=sample, table)
  print(table)
  otus_table <- rbind(otus_table, table)
}

otus_table <- otus_table[-1,]
colnames(otus_table) <- c("Sample", "OTU_ID", "counts")


## Transform this table to have a wide format where we have a column by sample
    transf_vsearch_table <- otus_table %>%
      group_by(Sample, OTU_ID) %>%
      summarise(counts = sum(counts)) %>%
      dcast(OTU_ID ~ Sample)

## Set OTU as rownames
    transf_vsearch_table <- set_rownames(x = transf_vsearch_table, value = transf_vsearch_table$OTU_ID)
    transf_vsearch_table[,1] <- NULL

## Write output
    write.table(x = transf_vsearch_table, file = count_table, sep="\t", quote=F, na = "0")
