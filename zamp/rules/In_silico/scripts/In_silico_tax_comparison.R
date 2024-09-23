
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
library(stringi);packageVersion("stringi")
library(reshape2);packageVersion("reshape2")
library(stringr);packageVersion("stringr")

## Input
count_table <- snakemake@input[["count_table"]]
asm_summary <- snakemake@input[["assembly_summary"]]
expected_taxonomy <- snakemake@input[["expected_taxonomy"]]
assigned_taxonomy <- snakemake@input[["assigned_taxonomy"]]
db_tax <- snakemake@input[["db_tax"]]

## Ouput
output_table <- snakemake@output[["output_table"]]
output_table_long <- snakemake@output[["output_table_long"]]

## Parameters
replace_empty_tax <- snakemake@params[["viz_replace_empty_tax"]]


### Read count table
print("reading count table")
count_table <- read.table(file = count_table, header = TRUE, check.names=FALSE)
transposed_counts <- t(count_table)
head(transposed_counts)

### Read assembly summary and taxonomy tables
print("reading metadata")
asm_df <- read.delim(file = asm_summary, sep = "\t", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
asm_df$path <- as.character(asm_df$path)
asm_df$assembly_name = sub("_genomic\\.fna\\.gz$","", basename(asm_df$path)) #add column with assembly prefix name
subset_asm_df <- asm_df %>% select("accession", "assembly_name")
tax_df <- read.delim(file = expected_taxonomy, sep = "\t", header = TRUE)
metadata <- merge(subset_asm_df, tax_df, by="accession")
head(metadata)

### Read pre processed DB used for zamp 
print('reading RDP taxonomy')
db <- read.table(db_tax, header=TRUE, sep = '\t', stringsAsFactors=FALSE)
in_db <- data.frame(species = metadata$species, in_DB = sapply(metadata$species, function(x) x %in% db$Species)) 
metadata <-  left_join(metadata, in_db, by="species")
head(metadata)

### Read and format taxonomy table
print("reading taxonomy table")
assigned_taxonomy <- read.table(file = assigned_taxonomy, header = FALSE, sep = "\t")
colnames(assigned_taxonomy) <- c("Feature.ID", "Taxonomy", "Confidence")
head(assigned_taxonomy)

### Convert the table into a tabular split version
taxonomy_table_split <- assigned_taxonomy %>% as_tibble() %>% 
  separate(Taxonomy, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

### Convert taxonomic levels as character (needed for the next steps)
taxonomy_table_split$Kingdom<-as.character(taxonomy_table_split$Kingdom)
taxonomy_table_split$Phylum<-as.character(taxonomy_table_split$Phylum)
taxonomy_table_split$Class<-as.character(taxonomy_table_split$Class)
taxonomy_table_split$Order<-as.character(taxonomy_table_split$Order)
taxonomy_table_split$Family<-as.character(taxonomy_table_split$Family)
taxonomy_table_split$Genus<-as.character(taxonomy_table_split$Genus)
taxonomy_table_split$Species<-as.character(taxonomy_table_split$Species)

print(paste("replace empty taxonomy :", replace_empty_tax))

if (isTRUE(replace_empty_tax)) {
  ### Replace NA by the previous order + a space_holder for each taxonomic level
  taxonomy_table_split$Kingdom[is.na(taxonomy_table_split$Kingdom)] <- (("Unkown_Kingdom")[is.na(taxonomy_table_split$Kingdom)])
  taxonomy_table_split$Phylum[is.na(taxonomy_table_split$Phylum)] <- ((paste(taxonomy_table_split$Kingdom,"_phy",sep=""))[is.na(taxonomy_table_split$Phylum)])
  taxonomy_table_split$Class[is.na(taxonomy_table_split$Class)] <- ((paste(taxonomy_table_split$Phylum,"_clas",sep=""))[is.na(taxonomy_table_split$Class)])
  taxonomy_table_split$Order[is.na(taxonomy_table_split$Order)] <- ((paste(taxonomy_table_split$Class,"_ord",sep=""))[is.na(taxonomy_table_split$Order)])
  taxonomy_table_split$Family[is.na(taxonomy_table_split$Family)] <- ((paste(taxonomy_table_split$Order,"_fam",sep=""))[is.na(taxonomy_table_split$Family)])
  taxonomy_table_split$Genus[is.na(taxonomy_table_split$Genus)] <- ((paste(taxonomy_table_split$Family,"_gen",sep=""))[is.na(taxonomy_table_split$Genus)])
  taxonomy_table_split$Species[is.na(taxonomy_table_split$Species)] <- ((paste(taxonomy_table_split$Genus,"_sp",sep=""))[is.na(taxonomy_table_split$Species)])

  print("table NA remplaced by spaceholders")

} else {

  print("table NA NOT remplaced by spaceholders")
}

## Re-unit the taxonomy after inclusion of spaceholders
taxonomy_table_split$Taxonomy <- unlist(unite(data = taxonomy_table_split, col = "Taxonomy", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";", remove = TRUE)[2])

## Create a metadata table to work with
comparison <- metadata

## Loop over the Assemblies
  for (i in comparison$assembly_name){
    ### Keep the row of the assembly
    transposed_counts_i <- transposed_counts[rownames(transposed_counts)==i,]
    ### Keep ony variants with non-null counts
    transposed_counts_i_f <- transposed_counts_i[transposed_counts_i > 0]
    ### Count the number of variants
    if('No_amp' %in% names(transposed_counts_i_f)){
      comparison[["Number_of_variants"]][comparison$assembly_name==i] <- 0
    } else {
      comparison[["Number_of_variants"]][comparison$assembly_name==i] <- length(names(transposed_counts_i_f))
    }

    for (j in 1:length(names(transposed_counts_i_f))){
      ## Create columns for each variants, count it, and give the corresponding taxonomy
      if(names(transposed_counts_i_f)[j] == 'No_amp'){
        comparison[[paste0("Amplicon_", j)]][comparison$assembly_name==i] <- 'No_amp'
        comparison[[paste0("Count_", j)]][comparison$assembly_name==i] <- 0
        comparison[[paste0("Tax_", j)]][comparison$assembly_name==i] <- "No_amplification_products"
      } else {
        comparison[[paste0("Amplicon_", j)]][comparison$assembly_name==i] <- names(transposed_counts_i_f)[j]
        comparison[[paste0("Count_", j)]][comparison$assembly_name==i] <- transposed_counts_i_f[j]
        comparison[[paste0("Tax_", j)]][comparison$assembly_name==i] <- taxonomy_table_split$Taxonomy[taxonomy_table_split$Feature.ID == names(transposed_counts_i_f)[j]]
      }
      j <- NULL
    }
    i <- NULL
  }

## Add the cumulated sum of all variants (OTUS)
# change no_amp from 1 to 0 counts
count_table_fix <- count_table 
count_table_fix['No_amp',] = 0
sum_of_count <- data.frame(Sum_of_copies = rowSums(t(count_table_fix)))
sum_of_count$assembly_name <- rownames(sum_of_count)
comparison <- left_join(comparison, sum_of_count)

#----------------------------------------------------------------#
## Create a table with the same informations but in a long format
## Create an Seq column
count_table$Feature.ID <- rownames(count_table)
count_table_tax <- left_join(count_table, assigned_taxonomy, by = 'Feature.ID')

## Melt it in long
melt_table <- melt(count_table_tax, id.vars = c("Feature.ID", "Taxonomy"), variable.name = "assembly_name") %>% filter(value>0)

compare_long <- right_join(metadata,melt_table)
write.table(x = compare_long, file = output_table_long, sep="\t",quote=F, row.names = FALSE)
#----------------------------------------------------------------#


## Check for Genus/Species assignment agreement
comparison$Genus_agreement <- NA
comparison$Species_agreement <- NA
comparison$Matching_amplicons <- NA
comparison$Discrepant_amplicons <- NA

for (assembly_ID in comparison$accession){
  print(assembly_ID)
  
  observed_taxa <- filter(compare_long, accession == assembly_ID ) %>%
    tidyr::separate(Taxonomy, sep = ";", into = c("Assigned_Kingdom","Assigned_Phylum","Assigned_Class","Assigned_Order","Assigned_Family","Assigned_Genus","Assigned_Species"))
  
  if(!(unique(observed_taxa$genus) %in% gsub("^.*__", "", observed_taxa$Assigned_Genus))){
    comparison$Genus_agreement[comparison$accession == assembly_ID] <- "Discrepant"
    comparison$Species_agreement[comparison$accession == assembly_ID] <- "Discrepant_Genus" 
    
  }else{

    comparison$Genus_agreement[comparison$accession == assembly_ID] <- "Matching"
      
    species <- unique(str_trim(str_remove(string = as.character(observed_taxa$species), as.character(unique(observed_taxa$genus)))))
    assigned_species <- gsub("[[:punct:]]","",str_remove(string = gsub("^.*__","",as.character(observed_taxa$Assigned_Species)), gsub("^.*__","",as.character(unique(observed_taxa$Assigned_Genus)))))

    if(all(grepl(pattern = species ,  x = assigned_species))) {
      comparison$Species_agreement[comparison$accession == assembly_ID] <- "Matching"
      comparison$Matching_amplicons[comparison$accession == assembly_ID] <- comparison$Number_of_variants[comparison$accession == assembly_ID]
      comparison$Discrepant_amplicons[comparison$accession == assembly_ID] <- 0

    }else  if (any(grepl(pattern = species ,  x = assigned_species))) {
      comparison$Species_agreement[comparison$accession == assembly_ID] <- "Partial_match"
      comparison$Matching_amplicons[comparison$accession == assembly_ID] <- length(grep(pattern = species ,  x = assigned_species))
      comparison$Discrepant_amplicons[comparison$accession == assembly_ID] <- comparison$Number_of_variants[comparison$accession == assembly_ID] - length(grep(pattern = species ,  x = assigned_species))
    }else{
      comparison$Species_agreement[comparison$accession == assembly_ID] <- "Discrepant"
      comparison$Discrepant_amplicons[comparison$accession == assembly_ID] <- comparison$Number_of_variants[comparison$accession == assembly_ID]
      comparison$Matching_amplicons[comparison$accession == assembly_ID] <- 0

    }
  }

}

# If the species was not in your pre-processed DB, an agreement cannot be found
comparison[comparison$in_DB == FALSE,]$Genus_agreement <- 'Species_not_in_your_DB'
comparison[comparison$in_DB == FALSE,]$Species_agreement <- 'Species_not_in_your_DB'
# If there was no amplification product:
comparison[comparison$Number_of_variants == 0,]$Genus_agreement <- 'No_PCR_amplification'
comparison[comparison$Number_of_variants == 0,]$Species_agreement <- 'No_PCR_amplification'

# relocate some columns
#(seems we are using old version of dplyr, this doesn't work)
# comparison <- comparison %>% 
#   dplyr::relocate(c("Sum_of_copies", "Genus_agreement", "Species_agreement", "Matching_amplicons", "Discrepant_amplicons"), .after = "Number_of_variants")
comparison <- comparison[,c(1:25, (ncol(comparison)-4):ncol(comparison), 26:(ncol(comparison)-5))]


print('final table head')
head(comparison)
write.table(x = comparison, file = output_table, sep="\t", quote=F, row.names = FALSE)
