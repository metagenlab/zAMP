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
phyloseq_object <- snakemake@input[[1]]

## Output
phyloseq_filtered_object <- snakemake@output[[1]]

## Parameters
tax_rank <- snakemake@params[["filter_tax_rank"]]
lineage <- snakemake@params[["filter_lineage"]]
filter_out_tax_rank  <- snakemake@params[["filter_out_tax_rank"]]
filter_out_lineage <- snakemake@params[["filter_out_lineage"]]
collapse_level <- snakemake@params[["collapse_level"]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")
library(dplyr);packageVersion("dplyr")

## Load the phyloseq object
phyloseq_object <- readRDS(phyloseq_object)
phyloseq_object

## filter taxa
filtered_taxa <- subset_taxa(phyloseq_object, get(tax_rank) == as.character(lineage))
filtered_taxa <- subset_taxa(filtered_taxa, get(filter_out_tax_rank) != as.character(filter_out_lineage))


if(collapse_level == "no_collapse"){

    ## Remove already in metadata alphia diversity values
    sample_data(filtered_taxa) <- select(sample_data(filtered_taxa), -c(Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher, Observed_min_1, Chao1_min_1, se.chao1_min_1, ACE_min_1, se.ACE_min_1, Shannon_min_1, Simpson_min_1, InvSimpson_min_1, Fisher_min_1))

    ## Add alpha diversity indexes to metadata
    alpha_div <- estimate_richness(physeq = filtered_taxa, split = TRUE)
    sample_data(filtered_taxa) <- cbind(sample_data(filtered_taxa),alpha_div)

    # Add alpha diveristy indexes at 1% filtration threshold
    ## Keep the taxa above 1%
    physeqrF = filter_taxa(filtered_taxa, function(x) mean(x) > 0.01,TRUE)
    ## IDs of taxa to be kept
    keeptaxa = taxa_names(physeqrF)
    ## All taxa
    alltaxa = taxa_names(filtered_taxa)
    ## Taxa to be kept
    myTaxa = alltaxa[alltaxa %in% keeptaxa]
    ## Keep only those
    physeqaF <- prune_taxa(myTaxa,filtered_taxa)
    ## Calculate new indexes
    alpha_div_1 <- estimate_richness(physeq = physeqaF, split = TRUE, "Observed")
    ## Rename those
    colnames(alpha_div_1) <- paste0(colnames(alpha_div_1), ("_min_1"))
    ## Again, bing these columns
    sample_data(filtered_taxa) <- cbind(sample_data(filtered_taxa),alpha_div_1)
    }

## Write the new phyloseq object
saveRDS(object = filtered_taxa, file = phyloseq_filtered_object)
