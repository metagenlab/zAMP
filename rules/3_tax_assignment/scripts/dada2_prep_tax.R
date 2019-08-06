# Title     : DADA2 - tax assign
# Objective : Prep taxonomy for dada2
# Created by: valentinscherz
# Created on: 06.08.19


## Redirect R output to the log file
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    ref_seqs <- snakemake@input[["ref_seqs"]]
    ref_tax <- snakemake@input[["ref_tax"]]

## Output
    King_to_Species <- snakemake@output[["King_to_Species"]]
    King_to_Genus <- snakemake@output[["King_to_Genus"]]
    Genus_species <- snakemake@output[["Genus_species"]]

## Load needed libraries

    library(dada2);packageVersion("dada2")
    library(phylotools);packageVersion(phylotools)
    library(dplyr);packageVersion(dplyr)
    library(tidyr);packageVersion(tidyr)

## Read data
    tax_table <- read.table(file=ref_tax, sep="\t", stringsAsFactors=FALSE)
    split_tax_table<-tax_table %>% separate(V2, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

## Format data
    King_to_genus <- split_tax_table
    King_to_genus$Species <- NULL
    King_to_genus <- King_to_genus %>% unite(taxonomy, c("Kingdom","Phylum","Class","Order","Family","Genus"), sep = ";", remove = TRUE)

    genus_species <- split_tax_table
    genus_species <- genus_species %>% unite(taxonomy, c("V1", "Genus", "Species"), sep = " ", remove = FALSE) %>% select(-c("Kingdom","Phylum","Class","Order","Family", "Genus","Species"))
    genus_species <- genus_species[,c(2,1)]

## Write data
    rename.fasta(infile = ref_seqs, ref_table = tax_table, outfile = King_to_Species)
    rename.fasta(infile = ref_seqs, ref_table = King_to_genus, outfile = King_to_Genus)
    rename.fasta(infile = ref_seqs, ref_table = genus_species, outfile = Genus_species)
