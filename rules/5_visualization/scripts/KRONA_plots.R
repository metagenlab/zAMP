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
phyloseq_melted_table <- snakemake@input[["phyloseq_melted_table"]]

## Output
output_folder <- snakemake@output[["output"]]
output_folder <- (dirname(output_folder)[1])

## Parameters
sample_label <- snakemake@params[["sample_label"]]
grouping_column <- snakemake@params[["grouping_column"]]
grouping_filter_column_value <- snakemake@params[["grouping_col_value"]]

print(sample_label)
print(grouping_column)
print(grouping_filter_column_value)

## Load needed library
library(dplyr);packageVersion("dplyr")

## Adapted function https://rdrr.io/github/cpauvert/psadd/man/plot_krona.html.
## This adapted version goes form the melted df instead of the phyloseq object in the original version of the function

melted_dataframe<- read.csv(file.path(phyloseq_melted_table), header = TRUE, sep = "\t")


    df <- filter(melted_dataframe, melted_dataframe[[grouping_column]] == grouping_filter_column_value)
    df <- filter(df, df[["Abundance"]] != 0)

    df <- df[, c("Abundance", sample_label, "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")]
    df <- as.data.frame(unclass(df))
    df[, 2] <- gsub(" |\\(|\\)", "", df[, 2])
    df[, 2] <- as.factor(df[, 2])
    dir.create(file.path(output_folder,"/",grouping_filter_column_value))
    for (lvl in levels(df[, 2])) {
      write.table(unique(df[which(df[, 2] == lvl), -2]), file = paste0(output_folder,"/",grouping_filter_column_value, "/", lvl, "taxonomy.txt"), sep = "\t", row.names = F, col.names = F, na = "", quote = F)
    }

    krona_args <- paste0(output_folder,"/", grouping_filter_column_value, "/", levels(df[, 2]), "taxonomy.txt,", levels(df[, 2]), collapse = " ")
    output <- paste0(output_folder,"/",grouping_filter_column_value,".html")
    system(paste("ktImportText", krona_args, "-o", output, sep = " "))


