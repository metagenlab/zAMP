# Title     : Counts to missmaches
# Objective : Add counts to missmaches table
# Created by: valentinscherz
# Created on: 06.11.2019

## Redirect R output
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")


## Input
mismatch_tables_path <- snakemake@input[["mismatch_tables_path"]]
count_table_path <- snakemake@input[["count_table_path"]]

## Output
    missmatch_plot_path <- snakemake@output[["missmatch_plot_path"]]
    merged_mismatch_table_path <- snakemake@output[["merged_mismatch_table_path"]]

## Params
    RUN <- snakemake@params["RUN"]

## Libraries
    library(dplyr)
    library(ggplot2)
    library(stringr)


save.image(file= paste0(getwd(), "/myEnvironment.RData"))

mismatches_table <- read.table(mismatch_tables_path, header = TRUE, sep = "\t")

count_table <- read.table(count_table_path, header = FALSE, sep = "\t",row.names = 1)
head(count_table, 10)
names(count_table)[1] <- "ID"
names(count_table)[2] <- "Counts"


mismatches_table$Query.id <- str_remove(string = mismatches_table$Query.id, pattern = ";size=\\d*$")

merged_mismatch_table <- left_join(mismatches_table, count_table, by = c("Query.id"= "ID"))

merged_mismatch_table$RUN <- head(as.character(RUN))


write.table(x = merged_mismatch_table, file = merged_mismatch_table_path, sep = "\t", row.names = FALSE)

merged_mismatch_table$Mismatches <- factor(merged_mismatch_table$Mismatches)


p <- ggplot(merged_mismatch_table) +
  geom_col(aes(x=Mismatches, y=Counts))
  # expand_limits(x=c(0,10), y = c(0,300000))

ggsave(plot = p, filename = missmatch_plot_path)
