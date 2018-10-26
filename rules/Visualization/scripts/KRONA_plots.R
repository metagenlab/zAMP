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
phyloseq_melted_table <- snakemake@output[["phyloseq_melted_table"]]

## Output
phyloseq_melted_table <- snakemake@output[["phyloseq_melted_table"]]

## Parameters
grouping_column <- snakemake@params[["grouping_column"]]

## Adapter function https://rdrr.io/github/cpauvert/psadd/man/plot_krona.html.
## This adapted version goes form the melted df instead of the phyloseq object in the original version of the function

adapted_KRONA_fct <- function(melted_dataframe, grouping_column, r_figures, x_axis_column){

for (i in (unique(melted_dataframe[[grouping_column]]))){

  df <- filter(melted_dataframe, melted_dataframe[[grouping_column]] == i)


  df <- df[, c("Abundance", x_axis_column, "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")]
  df <- as.data.frame(unclass(df))
  df[, 2] <- gsub(" |\\(|\\)", "", df[, 2])
  df[, 2] <- as.factor(df[, 2])
  dir.create(file.path(r_figures,"/",i))
  for (lvl in levels(df[, 2])) {
    write.table(unique(df[which(df[, 2] == lvl), -2]), file = paste0(r_figures,"/",i, "/", lvl, "taxonomy.txt"), sep = "\t", row.names = F, col.names = F, na = "", quote = F)
  }

  krona_args <- paste(r_figures,"/", i, "/", levels(df[, 2]), "taxonomy.txt,", levels(df[, 2]), sep = "", collapse = " ")
  output <- paste(r_figures,"/",i,".html", sep = "")
  system(paste("ktImportText", krona_args, "-o", output, sep = " "))

}}


adapted_KRONA_fct(melted_dataframe = physeq_df, grouping_column = "extraction_kit",  r_figures = ".", x_axis_column = "sample_label")
