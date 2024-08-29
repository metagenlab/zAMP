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
    # library(reshape2);packageVersion("reshape2")
    # library(magrittr);packageVersion("magrittr")
    library(tidyr)

## Input
    count_table_samples <- snakemake@input[["count_table_samples"]]

## Output
    count_table <- snakemake@output[["count_table"]]

# ## Reformat
#     otus_table<-as.data.frame(array(dim=c(1,3)))
#     colnames(otus_table) <- c("Sample", "V1", "V2")


# ## Reformat
#     otus_table<-as.data.frame(array(dim=c(1,3)))
#     colnames(otus_table) <- c("Sample", "V1", "V2")
otus_table <- data.frame()

    for (xx in count_table_samples){
        print(paste("Processing", xx))
        sample <- gsub("_count_table.tsv", "", basename(xx))
        print(paste("Sample:", sample))
        # table <- as.data.frame(array(dim=c(1,2)))
        # table$V1 <- "No_amp"
        # table$V2 <- 0
        if (file.size(xx) == 0){
            table <- cbind(Sample=sample, Seq_ID='No_amp', counts=1)
        } else {
            try(table <-read.table(file = xx, sep="\t", as.is=T, header=F))
            table <- cbind(Sample=sample, Seq_ID=table[,1] , counts=rowSums(table[,2:ncol(table), drop=FALSE]))
        }
        # table <- cbind("Sample"=sample, table)
        # print(table)
        otus_table <- rbind(otus_table, table)
    }
otus_table
    # otus_table <- otus_table[-1,]
    # colnames(otus_table) <- c("Sample", "Seq_ID", "counts")

otus_table$counts <- as.numeric(as.character(otus_table$counts))

## Transform this table to have a wide format where we have samples as columns
    # transf_vsearch_table <- otus_table %>%
    #   group_by(Sample, Seq_ID) %>%
    #   summarise(counts = sum(counts)) %>%
    #   dcast(Seq_ID ~ Sample)
transf_vsearch_table <- spread(otus_table, Sample, counts, fill=0)

## Set OTU as rownames
    rownames(transf_vsearch_table) <- transf_vsearch_table$Seq_ID
    transf_vsearch_table <- transf_vsearch_table[,-1]

## Write output
    write.table(x = transf_vsearch_table, file = count_table, sep="\t", quote=F, na = "0")
