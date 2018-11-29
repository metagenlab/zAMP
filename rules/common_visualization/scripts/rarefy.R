# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 21.11.18

# inspired from : https://rdrr.io/bioc/phyloseq/man/rarefy_even_depth.html

## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]
load(file =  file.path(phyloseq_object))

## Ouput
rarefied_phyloseq_path <- snakemake@output[["rarefied_phyloseq"]]

## Parameters
rarefy_value <- snakemake@params[["rarefaction_value"]]

## Load libraries
library('phyloseq')


## Determine rarefy level

rarefy_value <- as.numeric(rarefy_value)


if (is.numeric(rarefy_value)){
    print(rarefy_value)
    phyloseq_obj <- rarefy_even_depth(physeq = phyloseq_obj, sample.size = rarefy_value, rngseed = TRUE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
}else{
    print("No numerical rarefaction value given, using the depth of the sample with the lowest number of reads as default")
    rarefy_value <- min(sample_sums(phyloseq_obj))
    print(rarefy_value)
    phyloseq_obj <- rarefy_even_depth(physeq = phyloseq_obj, sample.size = rarefy_value, rngseed = TRUE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
}



# Write the phyloseq object
save(x = phyloseq_obj, file = rarefied_phyloseq_path)
