# Title     : Filter taxa
# Objective : Filter in and out taxa based on the assignement
# Created by: valentinscherz
# Created on: 28.05.2019

## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]

## Output
phyloseq_filtered_object <- snakemake@output[["phyloseq_filtered_object"]]

## Parameters
tax_rank <- snakemake@params[["filter_tax_rank"]]
lineage <- snakemake@params[["filter_lineage"]]
filter_out_tax_rank  <- snakemake@params[["filter_out_tax_rank"]]
filter_out_lineage <- snakemake@params[["filter_out_lineage"]]


## Load needed libraries
library(phyloseq);packageVersion("phyloseq")
library(dplyr);packageVersion("dplyr")


## Import the phyloseq object
phyloseq_obj <- readRDS(phyloseq_object)

## Filter taxa
filtered_taxa <- subset_taxa(phyloseq_obj, get(tax_rank) == as.character(lineage))
filtered_taxa <- subset_taxa(filtered_taxa, get(filter_out_tax_rank) != as.character(filter_out_lineage))
filtered_taxa <- prune_taxa(taxa_sums(filtered_taxa) > 0, filtered_taxa) ## Removes taxa not at least present in one sample

## Recompute alpha diversity indexes after this filtration
### Remove the previously computed values
sample_data(filtered_taxa) <- select(sample_data(filtered_taxa), -c(Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher, Observed_min_1))

### Add alpha diversity indexes to metadata
alpha_div <- estimate_richness(physeq = filtered_taxa, split = TRUE)
sample_data(filtered_taxa) <- cbind(sample_data(filtered_taxa),alpha_div)


### In addition, add the Observed over 1% metric
#### Keep the IDSs of the taxa above 1%
physeqrF <- filter_taxa(physeq = filtered_taxa, function(x) mean(x) > 0.01, FALSE)
#### Keep only those
physeqaF <- prune_taxa(physeqrF,filtered_taxa)
#### Calculate the Observed index
alpha_div_1 <- estimate_richness(physeq = physeqaF, split = TRUE, measure = "Observed")
#### Rename this index
colnames(alpha_div_1) <- paste0(colnames(alpha_div_1), ("_min_1"))
#### Again bind this new column
sample_data(filtered_taxa) <- cbind(sample_data(filtered_taxa),alpha_div_1)


## Write the new phyloseq object
saveRDS(object = filtered_taxa, file = phyloseq_filtered_object)
