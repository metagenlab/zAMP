# Title     : Export physeq object
# Objective : Export the content of physeq object
# Created by: valentinscherz
# Created on: 28.05.2019


## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]

## Output
tree_path <- snakemake@output[["tree_path"]]
meta_path <- snakemake@output[["meta_path"]]
taxonomy_path <- snakemake@output[["taxonomy_path"]]
OTU_path  <- snakemake@output[["OTU_path"]]
rep_seq_path  <- snakemake@output[["rep_seq_path"]]



## Load needed libraries
library(phyloseq);packageVersion("phyloseq")
library(dplyr);packageVersion("dplyr")
library(stringr);packageVersion("stringr")
library(Biostrings);packageVersion("Biostrings")



## Import the phyloseq phyloseq_object
phyloseq_obj <- readRDS(phyloseq_object)
phyloseq_obj <- prune_taxa(taxa_sums(phyloseq_obj) > 0, phyloseq_obj) ## Removes taxa not at least present in one sample

## Write the tree of the phyloseq object
tree <- phy_tree(phyloseq_obj)
ape::write.tree(tree, tree_path)

## Write the metadata table of the phyloseq object
metadata <- (sample_data(phyloseq_obj))
metadata$Sample <- rownames(metadata)
metadata <- metadata %>% select(Sample, everything())
write.table(metadata, meta_path , sep="\t", quote=F, row.names = FALSE)

## Write the taxonomy table of the phyloseq object
taxa_df<- as.data.frame(tax_table(phyloseq_obj), stringsAsFactors = F)
taxa_df$taxonomy <- apply(taxa_df, 1, function(x) str_c(x[!is.na(x)], collapse = ";"))
taxa_df <- taxa_df %>% select(taxonomy)
write.table(taxa_df, taxonomy_path , sep="\t", quote=F, col.names = F)

## Write the OTU table of the phyloseq object
write.table(otu_table(phyloseq_obj), OTU_path , sep="\t", quote=F)

## Write the representative sequences
writeXStringSet(refseq(phyloseq_obj), rep_seq_path, append=FALSE,
                compress=FALSE, format="fasta")


