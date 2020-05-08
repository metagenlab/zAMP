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
    library(dplyr);packageVersion("dplyr")
    library(tidyr);packageVersion("tidyr")
    library(Biostrings);packageVersion("Biostrings")

## Read data
    tax_table <- read.table(file=ref_tax, sep="\t", stringsAsFactors=FALSE)
    split_tax_table<-tax_table %>% separate(V2, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

## Format data
    King_to_genus <- split_tax_table
    King_to_genus$Species <- NULL
    King_to_genus <- King_to_genus %>% unite(taxonomy, c("Kingdom","Phylum","Class","Order","Family","Genus"), sep = ";", remove = TRUE)
    King_to_genus$taxonomy <- paste0(King_to_genus$taxonomy, ";")

    King_to_species <- tax_table
    King_to_species$V2 <- paste0(King_to_species$V2, ";")

    genus_species <- split_tax_table
    genus_species <- genus_species %>% unite(taxonomy, c("V1", "Genus", "Species"), sep = " ", remove = FALSE) %>% select(-c("Kingdom","Phylum","Class","Order","Family", "Genus","Species"))
    genus_species <- genus_species[,c(2,1)]



## Define a function taken from phylotools https://github.com/cran/phylotools/blob/master/R/rename.fasta.R

rename.fasta <- function(infile = NULL, ref_table,
                         outfile = "renamed.fasta"){
    #fasta <- read.fasta(infile)
    fastaFile <- readDNAStringSet(infile)
    seq.name = names(fastaFile)
    seq.text = paste(fastaFile)
    fasta <- data.frame(seq.name, seq.text)


    ## Convert the input ref to dataframe
    ## usually the file was obtained by
    ## save the result of get.names.fasta to a csv file.
    ## edit the csv file, and provide the name to be replaced.

    colnames(ref_table)  <- c("old_name", "new_name")
    res <- merge(x = fasta, y = ref_table, by.x = "seq.name",
                 by.y = "old_name", all.x = TRUE) ### Order of the sequence will change because of merge
    rename <- rep(NA, nrow(res))

    ### if the sequence was not found in the ref_table,
    ### keep the old name

    for(i in 1:nrow(res)){
       rename[i] <- ifelse(is.na(res$new_name[i]),
                           paste("old_name", "_",as.character(res$old_name[i]), sep = ""),
                           as.character(res$new_name[i]))
    }
    writeLines(paste(">", rename, "\n", as.character(res$seq.text), sep = ""), outfile )
    cat(paste(outfile, "has been saved to ", getwd(), "\n" ))
}


## Write data
    rename.fasta(infile = ref_seqs, ref_table = King_to_species, outfile = King_to_Species)
    rename.fasta(infile = ref_seqs, ref_table = King_to_genus, outfile = King_to_Genus)
    rename.fasta(infile = ref_seqs, ref_table = genus_species, outfile = Genus_species)
