# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 26.10.18

## Inspired from https://gist.github.com/erictleung/eda33bceceebb190eb9a44c93a077d32
## Redirect R output
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

## Input
phyloseq_melted_table <- snakemake@input[["phyloseq_melted_table"]]

## Output
output_folder <- snakemake@output[["output"]]

output_folder <- (dirname(output_folder)[1])

# output_folder <- paste0(output_folder,"/KRONA")

print(output_folder)
dir.create(output_folder)


#grouping_column <- snakemake@params[["grouping_column"]]
x_axis_column <- snakemake@params[["x_axis_column"]]


## Load needed library

library(dplyr);packageVersion("dplyr")

## Adapter function https://rdrr.io/github/cpauvert/psadd/man/plot_krona.html.
## This adapted version goes form the melted df instead of the phyloseq object in the original version of the function


melted_dataframe<- read.csv(file.path(phyloseq_melted_table), header = TRUE, sep = "\t")



adapted_KRONA_fct <- function(melted_dataframe, r_figures, x_axis_column){

for (i in (unique(melted_dataframe[["Sample"]]))){

    df <- filter(melted_dataframe, melted_dataframe[["Sample"]] == i)


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

      }
      }


adapted_KRONA_fct(melted_dataframe = melted_dataframe, r_figures = output_folder, x_axis_column = x_axis_column)
