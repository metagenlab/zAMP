# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 21.11.18

# inspired from : https://rdrr.io/bioc/phyloseq/man/rarefy_even_depth.html

## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]
load(file =  file.path(phyloseq_object))
Metadata_table <- snakemake@input[["Metadata_table"]]
taxonomy_table <- snakemake@input[["taxonomy_table"]]
tax_tree <- snakemake@input[["tax_tree"]]

## Ouput
rarefied_phyloseq_path <- snakemake@output[["rarefied_phyloseq"]]

## Parameters
rarefy_value <- snakemake@params[["rarefaction_value"]]
replace_empty_tax <- snakemake@params[["viz_replace_empty_tax"]]

## Load libraries
library(vegan);packageVersion("vegan")
library(dplyr);packageVersion("dplyr")
library(data.table);packageVersion("data.table")
library(tibble);packageVersion("tibble")
library(tidyr);packageVersion("tidyr")
library(readr);packageVersion("readr")
library(phyloseq);packageVersion("phyloseq")

## Rarefy
rarefy_value <- as.numeric(rarefy_value)

## Set seed for reproducibility
set.seed(1)


if (is.numeric(rarefy_value)){

    # Rarefy with vegan
    rarefied_OTU_table <- t(as.data.frame(rrarefy(t(phyloseq_obj@otu_table), sample = rarefy_value)))

    # Set the phyloseq_obj to NULL since the new one will have the same name.
    phyloseq_obj <- NULL

    ## Create a function
    load_objects_fct <- function(Metadata_table, taxonomy_table, replace_empty_tax = TRUE, tax_tree) {

        # Read count table
        # count_table <- read.table(file = features_counts_table, header = TRUE, check.names=FALSE)

        # Read sample_data
        metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE)

        # Read taxonomic tree
        PHY <- read_tree(tax_tree)

        # Read taxonomy table

            ### Load into R a table with all Features ID and their taxonomic assignemnt.
            taxonomy_table <-read.table(file = taxonomy_table, header = FALSE, sep = "\t")

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

            if(replace_empty_tax == TRUE) {
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

            taxonomy_table <<- taxonomy_table
            metadata <<- metadata
            PHY <<- PHY

            }


    ## Run the function to load all elements in R
    load_objects_fct(Metadata_table = Metadata_table, taxonomy_table = taxonomy_table, replace_empty_tax = TRUE, tax_tree = tax_tree)


    # Import all as phyloseq objects
    OTU <- otu_table(rarefied_OTU_table, taxa_are_rows = TRUE)
    TAX <- taxonomy_table %>% column_to_rownames("Feature.ID") %>% as.matrix() %>% tax_table()
    META <- metadata %>% as.data.frame() %>% column_to_rownames("Sample") %>% sample_data()


    # Sanity checks for consistent OTU names
    taxa_names(TAX)
    taxa_names(OTU)
    taxa_names(PHY)
    # Same sample names
    sample_names(OTU)
    sample_names(META)


    # Finally merge!
    phyloseq_obj <- phyloseq(OTU, TAX, META, PHY)

    # Write the phyloseq object
    save(x = phyloseq_obj, file = rarefied_phyloseq_path)


}else{
    print("No numerical rarefaction value given, using the depth of the sample with the lowest number of reads as default")
    #rarefy_value <- min(sample_sums(phyloseq_obj))
    #print(rarefy_value)
    #phyloseq_obj <- rarefy_even_depth(physeq = phyloseq_obj, sample.size = rarefy_value, rngseed = TRUE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

    # Write the phyloseq object
    #save(x = phyloseq_obj, file = rarefied_phyloseq_path)
}
