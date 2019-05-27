# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 26.10.18

## Inspired from https://stackoverflow.com/questions/15260245/r-convert-text-field-to-function and https://github.com/vmikk/metagMisc/blob/master/R/physeq_rm_na_tax.R
## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[[1]]

## Output
tree_path <- snakemake@output[["tree_path"]]
meta_path <- snakemake@output[["meta_path"]]
taxonomy_path <- snakemake@output[["taxonomy_path"]]
OTU_path  <- snakemake@output[["OTU_path"]]
filtered_rep_seq  <- snakemake@output[["filtered_rep_seq"]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")
library(dplyr);packageVersion("dplyr")
library(tidyr);packageVersion("tidyr")
library(stringr);packageVersion("stringr")
library(Biostrings);packageVersion("Biostrings")

## Load the phyloseq phyloseq_object
phyloseq_object <- readRDS(phyloseq_object)

# Write the tree of the phyloseq object
tree <- phy_tree(phyloseq_object)
ape::write.tree(tree, tree_path)

# Write the metadata table of the phyloseq object
metadata <- (sample_data(phyloseq_object))
metadata$Sample <- rownames(metadata)
metadata <- metadata %>% select(Sample, everything())
write.table(metadata, meta_path , sep="\t", quote=F, row.names = FALSE)

# Write the taxonomy table of the phyloseq object
taxa_df<- as.data.frame(tax_table(phyloseq_object), stringsAsFactors = F)
taxa_df$taxonomy <- apply(taxa_df, 1, function(x) str_c(x[!is.na(x)], collapse = ";"))
taxa_df <- taxa_df %>% select(taxonomy)
write.table(taxa_df, taxonomy_path , sep="\t", quote=F, col.names = F)

# Write the OTU table of the phyloseq object
write.table(otu_table(phyloseq_object), OTU_path , sep="\t", quote=F)

# Write the sequences
writeXStringSet(refseq(phyloseq_object), filtered_rep_seq, append=FALSE,
                compress=FALSE, format="fasta")


