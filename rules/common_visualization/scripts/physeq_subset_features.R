# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 26.10.18

## Inspired from https://stackoverflow.com/questions/15260245/r-convert-text-field-to-function and https://github.com/joey711/phyloseq/issues/694
## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[[1]]

## Output
phyloseq_filtered_object <- snakemake@output[[1]]

## Parameters
filter_features_subset_formula <- snakemake@params[["filter_features_subset_formula"]]
filter_features_subset_relative_or_absolute <- snakemake@params[["filter_features_subset_relative_or_absolute"]]
collapse_level <- snakemake@params[["collapse_level"]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")
library(plyr);packageVersion("plyr")
library(dplyr);packageVersion("dplyr")

## Load the phyloseq phyloseq_object
phyloseq_object <- readRDS(phyloseq_object)

## Convert the fonction to be used later
subset_fct <- function(x){}
body(subset_fct) <- as.quoted(filter_features_subset_formula)[[1]]

## Generate the list of taxa to keep
    if (filter_features_subset_relative_or_absolute == "relative"){
        phyloseq_object_pct  = transform_sample_counts(phyloseq_object, function(x) 100*x / sum(x))
        subset_features_list = taxa_names(filter_taxa(physeq = phyloseq_object_pct, flist = subset_fct, prune = TRUE))
        print("relative")

    } else if(filter_features_subset_relative_or_absolute == "absolute"){
        subset_features_list <- taxa_names(filter_taxa(physeq = phyloseq_object, flist = subset_fct, prune = TRUE))

    } else{stop("filter_features_subset_relative_or_absolute must be 'absolute' or 'relative'")}


## Prune taxa to keep only the ones passing the applied filter
subset_features <- prune_taxa(x = phyloseq_object, taxa = subset_features_list)

## Remove already in metadata alphia diversity values
sample_data(subset_features) <- select(sample_data(subset_features), -c(Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher))

if(collapse_level == "no_collapse"){

    ## Remove already in metadata alphia diversity values
    sample_data(subset_features) <- select(sample_data(subset_features), -c(Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher, Observed_min_1, Chao1_min_1, se.chao1_min_1, ACE_min_1, se.ACE_min_1, Shannon_min_1, Simpson_min_1, InvSimpson_min_1, Fisher_min_1))

    ## Add alpha diversity indexes to metadata
    alpha_div <- estimate_richness(physeq = subset_features, split = TRUE)
    sample_data(subset_features) <- cbind(sample_data(subset_features),alpha_div)

    # Add alpha diveristy indexes at 1% filtration threshold
    ## Keep the taxa above 1%
    physeqrF = filter_taxa(subset_features, function(x) mean(x) > 0.01,TRUE)
    ## IDs of taxa to be kept
    keeptaxa = taxa_names(physeqrF)
    ## All taxa
    alltaxa = taxa_names(subset_features)
    ## Taxa to be kept
    myTaxa = alltaxa[alltaxa %in% keeptaxa]
    ## Keep only those
    physeqaF <- prune_taxa(myTaxa,subset_features)
    ## Calculate new indexes
    alpha_div_1 <- estimate_richness(physeq = physeqaF, split = TRUE, "Observed")
    ## Rename those
    colnames(alpha_div_1) <- paste0(colnames(alpha_div_1), ("_min_1"))
    ## Again, bing these columns
    sample_data(subset_features) <- cbind(sample_data(subset_features),alpha_div_1)
    }


## Write the new phyloseq object
saveRDS(object = subset_features, file = phyloseq_filtered_object)
