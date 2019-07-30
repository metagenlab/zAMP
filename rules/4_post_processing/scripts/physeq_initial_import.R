# Title     : Import to phyloseq
# Objective : Import all elements to phyloseq
# Created by: valentinscherz
# Created on: 28.05.2019


## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
count_table <- snakemake@input[["count_table"]]
Metadata_table <- snakemake@input[["Metadata_table"]]
taxonomy_table <- snakemake@input[["taxonomy_table"]]
rep_seqs <- snakemake@input[["rep_seqs"]]
tax_tree <- snakemake@input[["tax_tree"]]

## Ouput
phyloseq_object <- snakemake@output[["phyloseq_object"]]

## Parameters
rarefy_value <- snakemake@params[["rarefaction_value"]]
replace_empty_tax <- snakemake@params[["viz_replace_empty_tax"]]



## Load libraries
library(dplyr);packageVersion("dplyr")
library(tibble);packageVersion("tibble")
library(tidyr);packageVersion("tidyr")
library(phyloseq);packageVersion("phyloseq")
library(Biostrings);packageVersion("Biostrings")


## Import data

### Read count table
print("reading count table")
count_table <- read.table(file = count_table, header = TRUE, check.names=FALSE)

### Read sample_data
print("reading metadata")
metadata <- read.table(file = Metadata_table, sep = "\t", header = TRUE, na.strings = "NA")

### Read taxonomic tree
print("reading taxonomic tree")
PHY <- read_tree(tax_tree)

### Read representative sequences
print("importing representative sequences from fasta")
SEQS <- readDNAStringSet(rep_seqs)

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


## Import as physeq object where needed
OTU <- otu_table(count_table, taxa_are_rows = TRUE)
TAX <- taxonomy_table %>% column_to_rownames("Feature.ID") %>% as.matrix() %>% tax_table()
META <- metadata %>% as.data.frame() %>% column_to_rownames("Sample") %>% sample_data()


## Merge all in a phyloseq object
phyloseq_obj <- phyloseq(OTU, TAX, META, PHY, SEQS)



## Compute alpha diversity indexes after this filtration
### Remove the previously computed values, in case
#sample_data(phyloseq_obj) <- select(sample_data(phyloseq_obj), -c(Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher, Observed_min_1))

drop <- c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "Observed_min_1")
sample_data(phyloseq_obj) <- sample_data(phyloseq_obj)[,!(names(sample_data(phyloseq_obj)) %in% drop)]


### Add alpha diversity indexes to metadata
alpha_div <- estimate_richness(physeq = phyloseq_obj, split = TRUE)
sample_data(phyloseq_obj) <- cbind(sample_data(phyloseq_obj),alpha_div)

### In addition, add the Observed over 1% metric
#### Keep the IDSs of the taxa above 1%
physeqrF <- filter_taxa(physeq = phyloseq_obj, function(x) mean(x) > 0.01, FALSE)
#### Keep only those
physeqaF <- prune_taxa(physeqrF,phyloseq_obj)
#### Calculate the Observed index
alpha_div_1 <- estimate_richness(physeq = physeqaF, split = TRUE, measure = "Observed")
#### Rename this index
colnames(alpha_div_1) <- paste0(colnames(alpha_div_1), ("_min_1"))
#### Again bind this new column
sample_data(phyloseq_obj) <- cbind(sample_data(phyloseq_obj),alpha_div_1)

phyloseq_obj <- prune_taxa(taxa_sums(phyloseq_obj) > 0, phyloseq_obj) ## Removes taxa not at least present in one sample



## Write the phyloseq object
saveRDS(object = phyloseq_obj, file = phyloseq_object)
