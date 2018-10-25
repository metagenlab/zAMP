# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 25.10.18

# Load objects into R
# Here, we load:
# - taxonomy table
# - Metadata table
# - features (either sequences variants or OTU) counts table
# - taxonomic tree
#
# Regarding the "taxonomy table":
#  The loaded table should have NO headers and in column 1. the features ID and in column 2. seven taxonomic ranks sepatrated by ";". We have here the option to replace or not the NA by the previous taxonomic rank and a spaceholder such as "_sp" for unassigned species
#
# Regarding the "Metadata table", it should follow the template described in doc (to be written), including beeing in the .tsv format
#
# Regarding the objects comming from Qiime2 pipeline which are in .qza format :
#  They are loaded using the package available here: (https://github.com/jbisanz/qiime2R)

## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
q_score_filtered_Fs <- snakemake@input[["q_score_filtered_Fs"]]
q_score_filtered_Rs <- snakemake@input[["q_score_filtered_Rs"]]
q_score_filtering_stats <- snakemake@input[["q_score_filtering_stats"]]

## Output
filtering_stats <- snakemake@output[["filtering_stats"]]
dna_sequences.fasta <- snakemake@output[["dna_sequences"]]
count_table.txt <- snakemake@output[["count_table"]]
otu_biom.biom <- snakemake@output["otu_biom"]

## Parameters
merged_min_length <- snakemake@params[["merged_min_length"]]
merged_max_length  <- snakemake@params[["merged_max_length"]]


## Load needed libraries
library(dplyr);packageVersion("dplyr")
library(dada2);packageVersion("dada2")
library(biomformat);packageVersion("biomformat")
library(phyloseq);packageVersion("phyloseq")




### Create the function

load_objects_fct <- function(taxonomy_class_table, replace_empty_tax = TRUE, features_counts_table , Metadata_table, tree) {

  ## Taxonomy table

  ### Load into R a table with all Features ID and their taxonomic assignemnt.
  taxonomy_class_table<-read_tsv(file.path(taxonomy_class_table), col_names = FALSE)

  ### Convert the table into a tabular split version
  taxonomy_class_table<-taxonomy_class_table %>% as.tibble() %>% separate(X2, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

  ### Replace the not properly named columns headers into proper ones
  colnames(taxonomy_class_table)[colnames(taxonomy_class_table)=="X1"] <- "Feature.ID"
  colnames(taxonomy_class_table)[colnames(taxonomy_class_table)=="X3"] <- "Confidence"

  ### Convert taxonomic levels as character (needed for the next steps)
  taxonomy_class_table$Kingdom<-as.character(taxonomy_class_table$Kingdom)
  taxonomy_class_table$Phylum<-as.character(taxonomy_class_table$Phylum)
  taxonomy_class_table$Class<-as.character(taxonomy_class_table$Class)
  taxonomy_class_table$Order<-as.character(taxonomy_class_table$Order)
  taxonomy_class_table$Family<-as.character(taxonomy_class_table$Family)
  taxonomy_class_table$Genus<-as.character(taxonomy_class_table$Genus)
  taxonomy_class_table$Species<-as.character(taxonomy_class_table$Species)

  if(replace_empty_tax == TRUE) {
    ### Replace NA by the previous order + a space_holder for each taxonomic level
    taxonomy_class_table$Kingdom[is.na(taxonomy_class_table$Kingdom)] <- (("Unkown_Kingdom")[is.na(taxonomy_class_table$Kingdom)])
    taxonomy_class_table$Phylum[is.na(taxonomy_class_table$Phylum)] <- ((paste(taxonomy_class_table$Kingdom,"_phy",sep=""))[is.na(taxonomy_class_table$Phylum)])
    taxonomy_class_table$Class[is.na(taxonomy_class_table$Class)] <- ((paste(taxonomy_class_table$Phylum,"_clas",sep=""))[is.na(taxonomy_class_table$Class)])
    taxonomy_class_table$Order[is.na(taxonomy_class_table$Order)] <- ((paste(taxonomy_class_table$Class,"_ord",sep=""))[is.na(taxonomy_class_table$Order)])
    taxonomy_class_table$Family[is.na(taxonomy_class_table$Family)] <- ((paste(taxonomy_class_table$Order,"_fam",sep=""))[is.na(taxonomy_class_table$Family)])
    taxonomy_class_table$Genus[is.na(taxonomy_class_table$Genus)] <- ((paste(taxonomy_class_table$Family,"_gen",sep=""))[is.na(taxonomy_class_table$Genus)])
    taxonomy_class_table$Species[is.na(taxonomy_class_table$Species)] <- ((paste(taxonomy_class_table$Genus,"_sp",sep=""))[is.na(taxonomy_class_table$Species)])

    print("table NA remplaced by spaceholders")

  }else{

    print("table NA NOT remplaced by spaceholders")
  }

  ### Convert taxonomic levels back as factor (Not sure if really needed?)
  taxonomy_class_table$Kingdom<-as.factor(taxonomy_class_table$Kingdom)
  taxonomy_class_table$Phylum<-as.factor(taxonomy_class_table$Phylum)
  taxonomy_class_table$Class<-as.factor(taxonomy_class_table$Class)
  taxonomy_class_table$Order<-as.factor(taxonomy_class_table$Order)
  taxonomy_class_table$Family<-as.factor(taxonomy_class_table$Family)
  taxonomy_class_table$Genus<-as.factor(taxonomy_class_table$Genus)
  taxonomy_class_table$Species<-as.factor(taxonomy_class_table$Species)

  ### Save this table in global environnement
  taxonomy_class_table <<- taxonomy_class_table


  ## .tsv Metadata table

  ### Load Metadata table into R environment
  Metadata<-read_tsv(file.path(Metadata_table), col_names = TRUE)

  ### Save this table in global environnement
  Metadata <<- Metadata


  ## Qiime2 .qza features counts table

  ### Load Qiime2 .qza counts table into R environment
  features_counts_table<-read_qza(file.path(features_counts_table))

  ### Shows the unique identifier of the file
  print("features counts table")
  print(features_counts_table$uuid)

  ### Print information about the object
  names(features_counts_table)

  ### Show the first 5 samples and features
  print(features_counts_table$data[1:5,1:5])

  ### Save this table in global environnement
  features_counts_table <<- features_counts_table



  ## Qiime2 .qza tree.


  ### In case there IS a tree
  if (!is.null(tree)){

  ### Load into environment Qiime2 .qza tree
  tree<-read_qza(file.path(tree))

  ### Save this table in global environnement
  tree <<- tree

  ### Shows the unique identifier of the file
  print("Tree")
  print(tree$uuid)

  ### Shows data
  print(tree$data)

  ### Shows data
  print(tree$data$Vectors[1:5, 1:15])

  ### Save this table in global environnement
  tree <<- tree

  ### In case there is NO tree

  }else if (is.null(tree)){
    print("No tree loaded, tree is NULL")
  }


}

