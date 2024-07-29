# Title     : Vsearch count table
# Objective : Format a count table from vsearch output
# Created by: valentinscherz
# Created on: 06.06.1p

## Redirect R output
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Load needed library
    library(dplyr);packageVersion("dplyr")
    library(reshape2);packageVersion("reshape2")
    library(magrittr);packageVersion("magrittr")

## Input
    count_table_samples <- snakemake@input[["count_table_samples"]]

## Output
    count_table <- snakemake@output[["count_table"]]

## Reformat
otus_table <- data.frame(array(dim=c(0,3)))
colnames(otus_table) <- c("Sample", "OTU_ID", "counts")

### Loop over each sample file. If it is empty, then we just add the factor in the levels to have it then
for (file_path in count_table_samples){
  sample_name <- gsub("_count_table.tsv", "", basename(file_path))
  sample_otu_table <- read.table(file = file_path, sep="\t", as.is=T, check.names = F, header=T, comment.char = "",  skipNul = TRUE)
  colnames(sample_otu_table) <- c("OTU_ID", "counts")
  if (nrow(sample_otu_table)>0){
    sample_otu_table <- cbind("Sample"=sample_name, sample_otu_table)
    otus_table <- rbind(otus_table, sample_otu_table)
  }else if (nrow(sample_otu_table) == 0){
    levels(otus_table$Sample) <- c(levels(otus_table$Sample), sample_name)
  }
}



## Transform this table to have a wide format where we have a column by sample
transf_vsearch_table <- otus_table %>% 
  dplyr::group_by(Sample, OTU_ID, .drop = FALSE) %>%
  dplyr::summarise(counts = sum(counts)) %>%
  reshape2::dcast(OTU_ID ~ Sample) %>%
  dplyr::filter(!is.na(OTU_ID))

## Set OTU as rownames
transf_vsearch_table <- set_rownames(x = transf_vsearch_table, value = transf_vsearch_table$OTU_ID)
transf_vsearch_table[,1] <- NULL

## Write output
    write.table(x = transf_vsearch_table, file = count_table, sep="\t", quote=F, na = "0")
