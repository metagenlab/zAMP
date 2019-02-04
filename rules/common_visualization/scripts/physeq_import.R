# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 25.10.18

## Inspired from https://gist.github.com/erictleung/eda33bceceebb190eb9a44c93a077d32
## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
features_counts_table <- snakemake@input[["count_table"]]
Metadata_table <- snakemake@input[["Metadata_table"]]
taxonomy_table <- snakemake@input[["taxonomy_table"]]
tax_tree <- snakemake@input[["tax_tree"]]

## Output
phyloseq_object <- snakemake@output[["phyloseq_object"]]

## Parameters
replace_empty_tax <- snakemake@params[["viz_replace_empty_tax"]]

## Load needed libraries
library(phyloseq);packageVersion("phyloseq")
library(dplyr);packageVersion("dplyr")
library(data.table);packageVersion("data.table")
library(tibble);packageVersion("tibble")
library(tidyr);packageVersion("tidyr")
library(readr);packageVersion("readr")

## Create a function
load_objects_fct <- function(features_counts_table, Metadata_table, taxonomy_table, replace_empty_tax = TRUE, tax_tree) {

    # Read count table
    print("reading count table")
    count_table <- read.table(file = features_counts_table, header = TRUE, check.names=FALSE)

    # Read sample_data
    print("reading metadata")
    metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE, na.strings = "NA")

    # Read taxonomic tree
    print("reading taxonomic tree")
    PHY <- read_tree(tax_tree)

    # Read taxonomy table

        ### Load into R a table with all Features ID and their taxonomic assignemnt.
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
        count_table <<- count_table
        PHY <<- PHY

        }


## Run the function to load all elements in R
load_objects_fct(features_counts_table = features_counts_table, Metadata_table = Metadata_table, taxonomy_table = taxonomy_table, replace_empty_tax = TRUE, tax_tree = tax_tree)


# Import all as phyloseq objects
OTU <- otu_table(count_table, taxa_are_rows = TRUE)
TAX <- taxonomy_table %>% column_to_rownames("Feature.ID") %>% as.matrix() %>% tax_table()
META <- metadata %>% as.data.frame() %>% column_to_rownames("SampleID") %>% sample_data()


# Sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(OTU)
taxa_names(PHY)


# Same sample names
sample_names(OTU)
sample_names(META)


# Finally merge!
phyloseq_obj <- phyloseq(OTU, TAX, META, PHY)


# Add alpha diversity indexes to metadata
alpha_div <- estimate_richness(physeq = phyloseq_obj, split = TRUE)
sample_data(phyloseq_obj) <- cbind(sample_data(phyloseq_obj),alpha_div)


# Add alpha diveristy indexes at 1% filtration threshold
## Keep the taxa above 1%
physeqrF = filter_taxa(phyloseq_obj, function(x) mean(x) > 0.01,TRUE)
## IDs of taxa to be kept
keeptaxa = taxa_names(physeqrF)
## All taxa
alltaxa = taxa_names(phyloseq_obj)
## Taxa to be kept
myTaxa = alltaxa[alltaxa %in% keeptaxa]
## Keep only those
physeqaF <- prune_taxa(myTaxa,phyloseq_obj)
## Calculate new indexes
alpha_div_1 <- estimate_richness(physeq = physeqaF, split = TRUE)
## Rename those
colnames(alpha_div_1) <- paste0(colnames(alpha_div_1), ("_min_1"))
## Again, bing these columns
sample_data(phyloseq_obj) <- cbind(sample_data(phyloseq_obj),alpha_div_1)


# Write the phyloseq object
saveRDS(object = phyloseq_obj, file = phyloseq_object)
