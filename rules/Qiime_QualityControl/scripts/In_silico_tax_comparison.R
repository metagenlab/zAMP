
# Title     : Import to phyloseq
# Objective : Import all elements to phyloseq
# Created by: valentinscherz
# Created on: 28.05.2019


## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")




## Load libraries
library(dplyr);packageVersion("dplyr")
library(tibble);packageVersion("tibble")
library(tidyr);packageVersion("tidyr")
library(phyloseq);packageVersion("phyloseq")
library(Biostrings);packageVersion("Biostrings")


## Input
count_table <- snakemake@input[["count_table"]]

Metadata_table <- snakemake@input[["Metadata_table"]]

taxonomy_table <- snakemake@input[["taxonomy_table"]]


## Ouput
output_table <- snakemake@output[["output_table"]]

## Parameters
replace_empty_tax <- snakemake@params[["viz_replace_empty_tax"]]


### Read count table
print("reading count table")
count_table <- read.table(file = count_table, header = TRUE, check.names=FALSE)
transposed_counts <- t(count_table)

### Read sample_data
print("reading metadata")
metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE, na.strings = "NA")


### Read and format taxonomy table
print("reading taxonomy table")
taxonomy_table<-read.table(file = taxonomy_table, header = FALSE, sep = "\t")

### Convert the table into a tabular split version
taxonomy_table<-taxonomy_table %>% as.tibble() %>% separate(V2, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

### Replace the not properly named headers into proper ones
colnames(taxonomy_table)[colnames(taxonomy_table)=="V1"] <- "Feature.ID"
colnames(taxonomy_table)[colnames(taxonomy_table)=="V3"] <- "Confidence"

    ### Convert taxonomic levels as character (needed for the next steps)
    taxonomy_table$Kingdom<-as.character(taxonomy_table$Kingdom)
    taxonomy_table$Phylum<-as.character(taxonomy_table$Phylum)
    taxonomy_table$Class<-as.character(taxonomy_table$Class)
    taxonomy_table$Order<-as.character(taxonomy_table$Order)
    taxonomy_table$Family<-as.character(taxonomy_table$Family)
    taxonomy_table$Genus<-as.character(taxonomy_table$Genus)
    taxonomy_table$Species<-as.character(taxonomy_table$Species)

    print(paste("replace empty taxonomy :", replace_empty_tax))

    if(isTRUE(replace_empty_tax)) {
      ### Replace NA by the previous order + a space_holder for each taxonomic level
      taxonomy_table$Kingdom[is.na(taxonomy_table$Kingdom)] <- (("Unkown_Kingdom")[is.na(taxonomy_table$Kingdom)])
      taxonomy_table$Phylum[is.na(taxonomy_table$Phylum)] <- ((paste(taxonomy_table$Kingdom,"_phy",sep=""))[is.na(taxonomy_table$Phylum)])
      taxonomy_table$Class[is.na(taxonomy_table$Class)] <- ((paste(taxonomy_table$Phylum,"_clas",sep=""))[is.na(taxonomy_table$Class)])
      taxonomy_table$Order[is.na(taxonomy_table$Order)] <- ((paste(taxonomy_table$Class,"_ord",sep=""))[is.na(taxonomy_table$Order)])
      taxonomy_table$Family[is.na(taxonomy_table$Family)] <- ((paste(taxonomy_table$Order,"_fam",sep=""))[is.na(taxonomy_table$Family)])
      taxonomy_table$Genus[is.na(taxonomy_table$Genus)] <- ((paste(taxonomy_table$Family,"_gen",sep=""))[is.na(taxonomy_table$Genus)])
      taxonomy_table$Species[is.na(taxonomy_table$Species)] <- ((paste(taxonomy_table$Genus,"_sp",sep=""))[is.na(taxonomy_table$Species)])

      print("table NA remplaced by spaceholders")

    }else{

      print("table NA NOT remplaced by spaceholders")
    }

    taxonomy_table$Taxonomy <- unlist(unite(data = taxonomy_table, col = "Taxonomy", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";", remove = TRUE)[2])

comparison <- metadata

for (i in comparison$AssemblyNames){
  print(i)
  transposed_counts_i <- transposed_counts[rownames(transposed_counts)==i,]
  transposed_counts_i_f <- transposed_counts_i[transposed_counts_i > 0]

  for (j in length(names(transposed_counts_i_f))){
    print(j)
    comparison[[paste0("Amplicon_", j)]][comparison$AssemblyNames==i] <- names(transposed_counts_i_f)[j]
    comparison[[paste0("Count_", j)]][comparison$AssemblyNames==i] <- transposed_counts_i_f[j]
    comparison[[paste0("Tax_", j)]][comparison$AssemblyNames==i] <- taxonomy_table$Taxonomy[taxonomy_table$Feature.ID == names(transposed_counts_i_f)[j]]

  }
}



write.table(x = comparison, file = output_table, sep="\t", quote=F, na = "0")
