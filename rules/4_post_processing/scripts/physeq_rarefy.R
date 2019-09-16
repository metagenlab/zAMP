# Title     : Rarefy count table
# Objective : Rarefy count table
# Created by: valentinscherz
# Created on: 28.05.2019


## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]

## Ouput
rarefied_phyloseq <- snakemake@output[["phyloseq_object"]]

## Parameters
rarefy_value <- snakemake@params[["rarefaction_value"]]


## Load libraries
library(vegan);packageVersion("vegan")
library(phyloseq);packageVersion("phyloseq")
library(dplyr);packageVersion("dplyr")


## Set seed for reproducibility
set.seed(1)


## Import the phyloseq object
phyloseq_obj <- readRDS(phyloseq_object)


print("start raref")
## Rarefy the count table
otu_table(phyloseq_obj) <- t(rrarefy(t(otu_table(phyloseq_obj)), sample = as.numeric(rarefy_value)))
print("stop raref)

## Compute alpha diversity indexes after this rarefaction
### Remove the previously computed values, in case
#sample_data(phyloseq_obj) <- select(sample_data(phyloseq_obj), -c(Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Observed_min_1))
drop <- c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "Observed_min_1")
sample_data(phyloseq_obj) <- sample_data(phyloseq_obj)[,!(names(sample_data(phyloseq_obj)) %in% drop)]




### Add alpha diversity indexes to metadata
alpha_div <- estimate_richness(physeq = phyloseq_obj, split = TRUE, c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson"))
sample_data(phyloseq_obj) <- cbind(sample_data(phyloseq_obj),alpha_div)

### In addition, add the Observed over 1% metric
#### Keep the IDSs of the taxa above 1%
physeqrF <- filter_taxa(physeq = phyloseq_obj, function(x) mean(x) > 0.01, FALSE)
#### Keep only those
physeqaF <- prune_taxa(physeqrF,phyloseq_obj)
#### Calculate the Observed index
alpha_div_1 <- estimate_richness(physeq = physeqaF, split = TRUE, measure = "Observed")
#### Rename this index
colnames(alpha_div_1) <- paste0(colnames(alpha_div_1), ("_min_1"))
#### Again bind this new column
sample_data(phyloseq_obj) <- cbind(sample_data(phyloseq_obj),alpha_div_1)

phyloseq_obj <- prune_taxa(taxa_sums(phyloseq_obj) > 0, phyloseq_obj) ## Removes taxa not at least present in one sample


# Write the phyloseq object
saveRDS(object = phyloseq_obj, file = rarefied_phyloseq)
