# Title     : DADA2 - tax assign
# Objective : Assign taxonomy with DADA2 implementation of RDP
# Created by: valentinscherz
# Created on: 06.08.19


## Redirect R output to the log file
    #log <- file(snakemake@log[[1]], open="wt")
    #sink(log)
    #sink(log, type="message")

## Input
    seqs <- snakemake@input[["seqs"]]
    #King_to_Species <- snakemake@input[["King_to_Species"]]
    King_to_Genus <- snakemake@input[["King_to_Genus"]]
    Genus_species <- snakemake@input[["Genus_species"]]

## Output
    tax <- snakemake@output[["tax"]]

## Load needed libraries
    library(dada2);packageVersion("dada2")
    library(dplyr);packageVersion("dplyr")
    library(tidyr);packageVersion("tidyr")
    library(Biostrings);packageVersion("Biostrings")

## Read data
    fastaFile <- readDNAStringSet(seqs)
    seq.name = names(fastaFile)
    sequence = paste(fastaFile)
    seqs_table <- data.frame(seq.name, sequence)

## Format seqs
    seq_table <- as.character(seqs$seq.text)
    names(seq_table) <- make.unique(seqs$seq.name)

## Assign taxonomy
    taxa <- assignTaxonomy(seqs = seq_table, refFasta = King_to_Genus, taxLevels = c("Kingdom","Phylum","Class","Order","Family","Genus"), multithread=TRUE, tryRC = TRUE, minBoot = 50, verbose = TRUE)

    taxa_species <- addSpecies(taxtab = taxa, refFasta = Genus_species, verbose=TRUE, allowMultiple = TRUE, tryRC = TRUE)


## Format an write output
    taxa_table <- data.frame(cbind(Row.Names = rownames(taxa_species), taxa_species))
    taxa_table$Row.Names<-NULL
    taxa_table <- taxa_table %>% unite(taxonomy, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";", remove = TRUE)
    taxa_table <- data.frame(cbind(Row.Names = rownames(taxa_table), taxa_table))
    taxa_table$Row.Names<-NULL
    write.table(taxa_table, tax , sep="\t", quote=F, col.names = F)

