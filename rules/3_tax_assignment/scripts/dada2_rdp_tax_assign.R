# Title     : DADA2 - tax assign
# Objective : Assign taxonomy with DADA2 implementation of RDP
# Created by: valentinscherz
# Created on: 06.08.19


## Redirect R output to the log file
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    seqs <- snakemake@input[["seqs"]]
    King_to_Species <- snakemake@input[["King_to_Species"]]
    King_to_Genus <- snakemake@input[["King_to_Genus"]]
    Genus_species <- snakemake@input[["Genus_species"]]

## Output
    tax <- snakemake@output[["tax"]]

## Parameters
    cpu <- snakemake@threads

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
    taxa <- assignTaxonomy(seqs = seq_table, refFasta = King_to_Species, taxLevels = c("Kingdom","Phylum","Class","Order","Family","Genus", "Species"), multithread=cpu, tryRC = TRUE, minBoot = 50, verbose = TRUE, outputBootstraps = TRUE)
    #taxa <- addSpecies(taxtab = taxa, refFasta = Genus_species, verbose=TRUE, allowMultiple = TRUE, tryRC = TRUE)
    # Not working for some reason

## Format an write output
    taxa_table <- data.frame(taxa)
    taxa_table <- data.frame(cbind(Row.Names = names(rownames(unlist(taxa$tax))), taxa_table))
    taxa_table <- taxa_table %>% unite(taxonomy, starts_with(match = "tax."), sep = ";", remove = TRUE)
    taxa_table <- taxa_table %>% unite(boot, starts_with(match = "boot"), sep = ";", remove = TRUE)

    #taxa_table$Row.Names<-NULL
    write.table(taxa_table, tax, row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)

