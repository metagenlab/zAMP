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
    King_to_Species <- snakemake@input[["King_to_Species"]]
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
    fasta <- data.frame(seq.name, sequence)

## Format seqs
    seq_table <- as.character(fasta$sequence)
    names(seq_table) <- as.character(fasta$seq.name)

## Assign taxonomy
    print("Assigning")
    head(seq_table, 10)
    head(King_to_Genus, 10)

    taxa <- assignTaxonomy(seqs = seq_table, refFasta = King_to_Genus, taxLevels = c("Kingdom","Phylum","Class","Order","Family","Genus"), multithread=4, tryRC = TRUE, minBoot = 50, verbose = TRUE)

    head(Genus_species, 10)

    taxa <- addSpecies(taxtab = taxa, refFasta = Genus_species, verbose=TRUE, allowMultiple = TRUE, tryRC = TRUE)


## Format an write output
    taxa_table <- data.frame(cbind(Row.Names = rownames(taxa), taxa))
    taxa_table <- taxa_table %>% unite(taxonomy, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";", remove = TRUE)
    taxa_table <- data.frame(cbind(Row.Names = rownames(taxa_table), taxa_table))
    taxa_table$Row.Names<-NULL
    print(taxa_table)

    write.table(x = taxa_table, file = tax , sep="\t", quote=F, col.names = F)

