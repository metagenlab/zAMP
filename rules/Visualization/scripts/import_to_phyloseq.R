# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 25.10.18

## Create phyloseq object and import Data into R but also write it in a file

## Function to import data in R from the output of the DADA2 pipeline
### Load the function
source("../Ressources/R_scripts/importation/dada2_to_R.R")
### Run the function
load_objects_fct(taxonomy_class_table = "../3_classified/RDP/Valentindb/dna-sequences_tax_assignments.txt", replace_empty_tax = TRUE, tree = "../4_tree/rooted-tree.qza", Metadata_table = "../Ressources/Metadata.tsv", features_counts_table =  "../2_dada2_in_R_denoised/feature-table.qza")
## Function to create phyloseq object from what has just been imported
### Load the function
source("../Ressources/R_scripts/importation/create_phyloseq_from_metadata.R")
### Run the function
create_phyloseq_fct(physeq_name = "physeq", otu_table = features_counts_table, phy_tree = tree, tax_table = taxonomy_class_table, Metadata_table = Metadata, melted_df_name = "physeq_df")

