
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
  Metadata_table <- snakemake@input[["Metadata_table"]]
  taxonomy_table <- snakemake@input[["taxonomy_table"]]

## Ouput
  output_table <- snakemake@output[["output_table"]]
  output_table_long <- snakemake@output[["output_table_long"]]

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
  taxonomy_table <- read.table(file = taxonomy_table, header = FALSE, sep = "\t")


### Convert the table into a tabular split version
  taxonomy_table_split<-taxonomy_table %>% as_tibble() %>% separate(V2, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

### Replace the not properly named headers into proper ones
  colnames(taxonomy_table_split)[colnames(taxonomy_table_split)=="V1"] <- "Feature.ID"
  colnames(taxonomy_table_split)[colnames(taxonomy_table_split)=="V3"] <- "Confidence"

### Convert taxonomic levels as character (needed for the next steps)
  taxonomy_table_split$Kingdom<-as.character(taxonomy_table_split$Kingdom)
  taxonomy_table_split$Phylum<-as.character(taxonomy_table_split$Phylum)
  taxonomy_table_split$Class<-as.character(taxonomy_table_split$Class)
  taxonomy_table_split$Order<-as.character(taxonomy_table_split$Order)
  taxonomy_table_split$Family<-as.character(taxonomy_table_split$Family)
  taxonomy_table_split$Genus<-as.character(taxonomy_table_split$Genus)
  taxonomy_table_split$Species<-as.character(taxonomy_table_split$Species)

  print(paste("replace empty taxonomy :", replace_empty_tax))

  if(isTRUE(replace_empty_tax)) {
    ### Replace NA by the previous order + a space_holder for each taxonomic level
    taxonomy_table_split$Kingdom[is.na(taxonomy_table_split$Kingdom)] <- (("Unkown_Kingdom")[is.na(taxonomy_table_split$Kingdom)])
    taxonomy_table_split$Phylum[is.na(taxonomy_table_split$Phylum)] <- ((paste(taxonomy_table_split$Kingdom,"_phy",sep=""))[is.na(taxonomy_table_split$Phylum)])
    taxonomy_table_split$Class[is.na(taxonomy_table_split$Class)] <- ((paste(taxonomy_table_split$Phylum,"_clas",sep=""))[is.na(taxonomy_table_split$Class)])
    taxonomy_table_split$Order[is.na(taxonomy_table_split$Order)] <- ((paste(taxonomy_table_split$Class,"_ord",sep=""))[is.na(taxonomy_table_split$Order)])
    taxonomy_table_split$Family[is.na(taxonomy_table_split$Family)] <- ((paste(taxonomy_table_split$Order,"_fam",sep=""))[is.na(taxonomy_table_split$Family)])
    taxonomy_table_split$Genus[is.na(taxonomy_table_split$Genus)] <- ((paste(taxonomy_table_split$Family,"_gen",sep=""))[is.na(taxonomy_table_split$Genus)])
    taxonomy_table_split$Species[is.na(taxonomy_table_split$Species)] <- ((paste(taxonomy_table_split$Genus,"_sp",sep=""))[is.na(taxonomy_table_split$Species)])

    print("table NA remplaced by spaceholders")

  }else{

    print("table NA NOT remplaced by spaceholders")
  }

## Re-unit the taxonomy after inclusion of spaceholders
  taxonomy_table_split$Taxonomy <- unlist(unite(data = taxonomy_table_split, col = "Taxonomy", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";", remove = TRUE)[2])

## Create a metadata table to work with
  comparison <- metadata


## Loop over the Assemblies
  for (i in comparison$AssemblyNames){
    print(i)
    ### Keep the row of the assembly
    transposed_counts_i <- transposed_counts[rownames(transposed_counts)==i,]
    ### Keep ony variants with non-null counts
    transposed_counts_i_f <- transposed_counts_i[transposed_counts_i > 0]
    ### Count the number of variants
    comparison[["Number of variants"]][comparison$AssemblyNames==i] <- length(names(transposed_counts_i_f))

    for (j in 1:length(names(transposed_counts_i_f))){
      print(j)
      ## Recover Create columns for each variants, count it, and give the corresponding taxonomy
      comparison[[paste0("Amplicon_", j)]][comparison$AssemblyNames==i] <- names(transposed_counts_i_f)[j]
      comparison[[paste0("Count_", j)]][comparison$AssemblyNames==i] <- transposed_counts_i_f[j]
      comparison[[paste0("Tax_", j)]][comparison$AssemblyNames==i] <- taxonomy_table_split$Taxonomy[taxonomy_table_split$Feature.ID == names(transposed_counts_i_f)[j]]
      j <- NULL
    }
    i <- NULL
  }

## Add the cumulated some for all variants
  sum_of_count <- data.frame(rowSums(t(count_table)))
  colnames(sum_of_count)[colnames(sum_of_count)=="rowSums.t.count_table.."] <- "Sum of copies"
  sum_of_count$AssemblyNames <- rownames(sum_of_count)
  comparison <- left_join(comparison, sum_of_count)

## Create a table with the same informations but in a long format
## Create an OTU column
  count_table$Feature.ID <- rownames(count_table)
## Match the position between the count and tax table and joined them
  m <- match(count_table$Feature.ID, taxonomy_table$V1)
  count_table_tax <- cbind(count_table, taxonomy_table$V2[m])
  colnames(count_table_tax)[colnames(count_table_tax)=="taxonomy_table$V2[m]"] <- "Class_tax"

print(1)
## Melt it in long
  melt_table <- melt(count_table_tax, id.vars = c("Feature.ID", "Class_tax"), variable.name = "AssemblyNames") %>% filter(value>0)

  compare_long <- right_join(metadata,melt_table)
  write.table(x = compare_long, file = output_table_long, sep="\t",quote=F, row.names = FALSE)

## Check for Genus/Species assignment agreement
comparison$Genus_agreement <- NA
comparison$Species_agreement <- NA

for (assembly_ID in unique(compare_long$AssemblyID)){
  
  observed_taxa <- filter(compare_long, AssemblyID == assembly_ID ) %>%
    tidyr::separate(Class_tax, sep = ";", into = c("Assigned_Kingdom","Assigned_Phylum","Assigned_Class","Assigned_Order","Assigned_Family","Assigned_Genus","Assigned_Species"))
  
  if(!(unique(observed_taxa$genus) %in% observed_taxa$Assigned_Genus)){
    comparison$Genus_agreement[comparison$AssemblyID == assembly_ID] <- "Discrepant"
    comparison$Species_agreement[comparison$AssemblyID == assembly_ID] <- "Discrepant_Genus" 
    
  }else{

    comparison$Genus_agreement[comparison$AssemblyID == assembly_ID] <- "Matching"
      
    species <- unique(str_trim(str_remove(string = as.character(observed_taxa$species), as.character(unique(observed_taxa$genus)))))
    assigned_species <- str_trim(str_remove(string = as.character(observed_taxa$Assigned_Species), as.character(unique(observed_taxa$Assigned_Genus))))

    if(all(grepl(pattern = species ,  x = assigned_species))) {
      comparison$Species_agreement[comparison$AssemblyID == assembly_ID] <- "Matching"
      
    }else  if (any(grepl(pattern = species ,  x = assigned_species))) {
      comparison$Species_agreement[comparison$AssemblyID == assembly_ID] <- "Partial match"
      
    }else{
      comparison$Species_agreement[comparison$AssemblyID == assembly_ID] <- "Discrepant"

    }
  }

}
  

write.table(x = comparison, file = output_table, sep="\t", quote=F, row.names = FALSE)

