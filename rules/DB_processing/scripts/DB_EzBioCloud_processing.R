## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
EzBio_tax1 <- snakemake@input[["tax"]]
EzBio_tax1
EzBio_uc <- snakemake@input[["uc"]]
EzBio_uc

## Output
EzBioCloud_V3V4_taxonomy <- snakemake@output[["filtrated"]]
EzBioCloud_V3V4_all_taxonomy <-  snakemake@output[["all"]]
EzBioCloud_V3V4_taxonomy_Qiime <- snakemake@output[["qiime"]]
EzBioCloud_V3V4_all_taxonomy_Qiime <- snakemake@output[["qiime_all"]]

## Parameters
# nom dans R <- snakemake@params[["nom dans params"]]

## Load needed libraries
library(dplyr);packageVersion("dplyr")
library(lattice);packageVersion("lattice")
library(stringr);packageVersion("stringr")
library(forcats);packageVersion("forcats")



## First changes to files obtained by qiime like eliminating what we do not need

# Recovery of the taxonomy table is a file in Qiime format
# to divide the string when this symbol is found ;
EzBio_tax <- read.delim(EzBio_tax1, header=FALSE, as.is = TRUE, colClasses = "character", strip.white = TRUE)
tax_split <- strsplit(EzBio_tax[,2],split = ";")
tax_table <- t(sapply(tax_split,function(x){x}))
rownames(tax_table) <- EzBio_tax[,1]

## Recovery of the uc table which is the output file of the script bash in jupyter notebook with qiime to carry out the extraction and the dereplication
# Separation of the table according to the first column (C, H and S) to obtain 2 separate tables
# We group the 2 tables, we do not take the table S because they are sequences present in duplicate in the table C
uc_table <- read.delim(EzBio_uc, header = FALSE, as.is = TRUE, colClasses = "character", strip.white = TRUE)
C_uc_table <- uc_table[uc_table[,1]=="C",]
H_uc_table <- uc_table[uc_table[,1]=="H",]
HC_uc_table <- rbind(C_uc_table,H_uc_table)

# cluster_rep_sq correspond to accessions numbers  of the table C
cluster_rep_sq <- C_uc_table[,9]
names(cluster_rep_sq) <- C_uc_table[,2]  # function names to get or set the names of an object

# list of characters (61654) containing lists with the different accessions numbers according to the identical sequences and the length
uc_split <- split(HC_uc_table[,9],HC_uc_table[,2])
cluster_size <- sapply(uc_split,length)

# a table containing only a length of 1 or a single accession number that will not be modified
# a table with a length greater than 1 which will undergo several filtering before emitting the different tables
singleton_clusters_rep_seq <- cluster_rep_sq[names(which(cluster_size==1))]
bigger_clusters <- names(which(cluster_size>1))
uc_split_bigger <- uc_split[bigger_clusters]

# These are the different numbers of accessions having a counting column to know that they are identical sequences are grouped together
bigger_clusters_table <- matrix(nrow=0,ncol=2)
for (xx in 1:length(uc_split_bigger))
{
    cl_name <- names(uc_split_bigger[xx])
    cl_tab <- cbind(cl_name,uc_split_bigger[[xx]])
    bigger_clusters_table <- rbind(bigger_clusters_table,cl_tab)
}

## Change in taxonomy, classification of the EzBioCloud database (kingdom, phylum, class, order, family and genus) by grouping bacteria on one side for example
##########

# list of reign names from a list of characters greater than 1
kingdom_split <- split(tax_table[bigger_clusters_table[,2],1],bigger_clusters_table[,1])
nr_kingdom <- sapply(kingdom_split,function(x){length(unique(x))})
table(nr_kingdom)
tax_table_cluster <- matrix(ncol=7,nrow=length(uc_split_bigger))
rownames(tax_table_cluster) <- names(kingdom_split)

#kingdom
tax_table_cluster[,1] <- sapply(kingdom_split,function(x){unique(x)})

##########
# Only accession numbers correspondig to Bacteria
Bacteria_clusters_names <- rownames(tax_table_cluster[tax_table_cluster[,1]%in%c("Bacteria"),])# "%in%" corresponds to the membership of one or more values in a vector
Bacteria_clusters_table <- bigger_clusters_table[bigger_clusters_table[,1]%in%Bacteria_clusters_names,]

# Only accession numbers correspondig to Archaea so the function said to be different from accession numbers cooresponding to Bacteria
Archaea_clusters_table <- bigger_clusters_table[!bigger_clusters_table[,1]%in%Bacteria_clusters_names,]

## genus
genus_split <- split(tax_table[Bacteria_clusters_table[,2],6],Bacteria_clusters_table[,1])

## terms of the species' names based on the number of genera and species
g2correct <- names(which(sapply(genus_split,function(x){length(unique(x))})>1))
tax_table_all <- tax_table

# c_name is the name of the local variable to the "for" function which will take a different value for each iteration
# g2correct corresponds to the vector containing the values that the take c_name
for (c_name in g2correct)
{
    # information needed to fulfill the conditions
    selected_tax_table <- tax_table[names(genus_split[[c_name]]),]
    new_genus_name <- paste(unique(selected_tax_table[,6]),collapse = "/")
    nr_genus <- length(unique(selected_tax_table[,6]))
    nr_species <- length(unique(selected_tax_table[,7]))

    # If loop when we have the number of genus superior or equal to 2 and the number of species superior to 4, we replace the name of the species by the genus and sp.
    if (nr_genus >= 2 & nr_species > 4)
    {
        new_species_name <- paste(new_genus_name,"sp.")
        new_species_name_f <- paste(new_species_name, "(", nr_species, ")", sep = "")
    }

    # If loop we have the number of species superior to 4, we replace the name of the species by the genus and sp.
    if (nr_species > 4)
    {
	      new_species_name <- paste(new_genus_name,"sp.")
	      new_species_name_f <- paste(new_species_name, "(", nr_species, ")", sep = "")
    }

    # If Loop when we have the number of species less than or equal to 4, we make a filtrer to eliminate the unknowns and we keep the same names (if there are several that are joined by "/")
    if ( nr_species <= 4)
    {
          species <- grep("[_s]$", unique(selected_tax_table[,7]),value = TRUE, invert = TRUE, ignore.case = TRUE) # filter the values of the vector directly with grep
	      new_species_name <- paste(unique(species),collapse = "/")
	      new_species_name_f <- paste(new_species_name, "(", nr_species, ")", sep = "")
	      print(new_species_name_f)
    }
    new_species_name_all <- paste(unique(selected_tax_table[,7]),collapse = "/")

    # retrieve and get new genus names and new species names
    new_entry <- c(selected_tax_table[1,1:5],new_genus_name,new_species_name_f)
    new_entry_all <- c(selected_tax_table[1,1:5],new_genus_name,new_species_name_all)
    new_selected_tax_table <- matrix(rep(new_entry,nrow(selected_tax_table)),nrow=nrow(selected_tax_table),ncol = 7,byrow = TRUE)
    new_selected_tax_table_all <- matrix(rep(new_entry_all,nrow(selected_tax_table)),nrow=nrow(selected_tax_table),ncol = 7,byrow = TRUE)

    # link the taxonomic list with the accession number
    tax_table[names(genus_split[[c_name]]),] <- new_selected_tax_table
    tax_table_all[names(genus_split[[c_name]]),] <- new_selected_tax_table_all
}

## terms of the species' names based on the presence of unknown species

# Same as « cluster name » except that "if" conditions are different
# concatenate the accession number to the species
species_split <- split(tax_table[Bacteria_clusters_table[,2],7],Bacteria_clusters_table[,1])
s2correct <- names(which(sapply(species_split,function(x){length(unique(x))})>1))

for (c_name in s2correct)
{
    selected_tax_table <- tax_table[names(species_split[[c_name]]),]
    new_genus_name <- selected_tax_table[1,6]
    nr_species <- length(unique(selected_tax_table[,7]))

    # when we do not find the symbol "_s" in the species names and that it is superior to 4, we put only the extension "sp."
    if (length(grep("[_s]$", unique(selected_tax_table[,7]),value = TRUE, invert = TRUE, ignore.case = TRUE)) > 4)
    {
        print("More than four with Species names")
        new_species_name <- paste(new_genus_name,"sp.")
        new_species_name_f <- paste(new_species_name, "(", nr_species, ")", sep = "")
        print(new_species_name_f)
    }

    # when we do not find the symbol "_s" in the species names and it is less than or equal to 4, we realize a filter then we glue the names of species one after the others separated by a "/"
    if (length(grep("[_s]$", unique(selected_tax_table[,7]),value = TRUE, invert = TRUE, ignore.case = TRUE)) <= 4)
    {
        print("Less than four with Species name")
        species <- grep("[_s]$", selected_tax_table[,7],value = TRUE, invert = TRUE, ignore.case = TRUE) # filter the values of the vector directly with grep
        new_species_name <- paste(species, collapse = "/")
        new_species_name_f <- paste(new_species_name, "(", nr_species,")", sep = "")
        print(new_species_name_f)
    }

    # allows to see if cases are not taken into account
    else
    {
        warning("Case not present in if conditions")
    }

    new_species_name_all <- paste(unique(selected_tax_table[,7]),collapse="/")

    new_entry <- c(selected_tax_table[1,1:6],new_species_name_f)
    new_entry_all <- c(selected_tax_table[1,1:6],new_species_name_all)
    new_selected_tax_table <- matrix(rep(new_entry,nrow(selected_tax_table)),nrow=nrow(selected_tax_table),ncol = 7,byrow = TRUE)
    new_selected_tax_table_all <- matrix(rep(new_entry_all,nrow(selected_tax_table)),nrow=nrow(selected_tax_table),ncol = 7,byrow = TRUE)

    tax_table[names(genus_split[[c_name]]),] <- new_selected_tax_table
    tax_table_all[names(genus_split[[c_name]]),] <- new_selected_tax_table_all

}


## grouping of tables and conditions to obtain a table with maximum 4 species
#version with four species

Bacteria_tax_table_cluster <- matrix(ncol=7,nrow=length(unique(Bacteria_clusters_table[,1])))

for(xx in 1:7)
{
    B_tax_split <- split(tax_table[Bacteria_clusters_table[,2],xx],Bacteria_clusters_table[,1])
    Bacteria_tax_table_cluster[,xx] <- sapply(B_tax_split,function(x){unique(x)})
}

# table of Bacteria
Bacteria_tax_table_cluster <- cbind(cluster_rep_sq[names(B_tax_split)],Bacteria_tax_table_cluster)


## Create a function
# Replacement of the first genus in the species column, where there are several species present but when there are two genera, the species is not modified
Double_genus <- function(Raw_table)
{
    if(grepl("/", Raw_table[7]))
    {
        print("/ in Genus")
        # keep the same name of species
        Raw_table[8] <- Raw_table[8]
    }

    else
    {
        # replacement of the name of the species by replacing the second genus with "/" to have only once the genus in the species when it has several species
        Raw_table[8] <- gsub(pattern = "/[A-Z][a-z]* ", replacement = "/", x = Raw_table[8])
    }
}
## Apply the just created function
Bacteria_tax_table_cluster[,8] <- apply(X = Bacteria_tax_table_cluster, MARGIN = 1, FUN = Double_genus)

# table of Archaea
Archaea_c_names <- cluster_rep_sq[unique(Archaea_clusters_table[,1])]
Archaea_tax_table_cluster <- cbind(Archaea_c_names,tax_table[Archaea_c_names,])


# create tĥe new table who regroup the table Bacteria, the table Archaea and the table with unique number accession
tax_table_cluster_singletons <- cbind(singleton_clusters_rep_seq,tax_table[singleton_clusters_rep_seq,])
tax_table_V3V4 <- rbind(Bacteria_tax_table_cluster,Archaea_tax_table_cluster,tax_table_cluster_singletons)
colnames(tax_table_V3V4) <- c("seq_name","kingdom","phylum","class","order","family","genus","species")
write.csv(x = tax_table_V3V4, file = EzBioCloud_V3V4_taxonomy, quote = FALSE, row.names = FALSE)

# Creation of a file type Qiime that is to say we have a space instead of a ";" after the accession number

EzBioCloud_V3V4_tax <- read.csv(file = EzBioCloud_V3V4_taxonomy, as.is = TRUE, strip.white = TRUE)

tax_collapsed <- apply(EzBioCloud_V3V4_tax[,2:8],1,paste,collapse=";")
new_tax <- cbind(EzBioCloud_V3V4_tax[,1],tax_collapsed)
write.table(x = new_tax, file = EzBioCloud_V3V4_taxonomy_Qiime ,quote = FALSE, sep="\t",row.names = FALSE ,col.names = FALSE)



## grouping tables and conditions to get a table with all the parameters
#version with all species names

# The same method as for the « version with four species » except that we use the file containing all the names and not the file sort and we do not realize the replacement of a genus in the column of the species when there are several kinds and species

Bacteria_tax_table_cluster_all <- matrix(ncol=7,nrow=length(unique(Bacteria_clusters_table[,1])))

for(xx in 1:7)
{
    B_tax_split <- split(tax_table_all[Bacteria_clusters_table[,2],xx],Bacteria_clusters_table[,1])
    Bacteria_tax_table_cluster_all[,xx] <- sapply(B_tax_split,function(x){unique(x)})
}

Bacteria_tax_table_cluster_all <- cbind(cluster_rep_sq[names(B_tax_split)],Bacteria_tax_table_cluster_all)

Archaea_c_names <- cluster_rep_sq[unique(Archaea_clusters_table[,1])]
Archaea_tax_table_cluster_all <- cbind(Archaea_c_names,tax_table_all[Archaea_c_names,])

tax_table_cluster_singletons <- cbind(singleton_clusters_rep_seq,tax_table_all[singleton_clusters_rep_seq,])
tax_table_V3V4_all <- rbind(Bacteria_tax_table_cluster_all,Archaea_tax_table_cluster_all,tax_table_cluster_singletons)

colnames(tax_table_V3V4_all) <- c("seq_name","kingdom","phylum","class","order","family","genus","species")
write.csv(x = tax_table_V3V4_all,file = EzBioCloud_V3V4_all_taxonomy,quote = FALSE, row.names = FALSE)

# Creation of a file type Qiime that is to say we have a space instead of a ";" after the accession number

EzBioCloud_V3V4_all_tax <- read.csv(file = EzBioCloud_V3V4_all_taxonomy, as.is = TRUE, strip.white = TRUE)

tax_all_collapsed <- apply(EzBioCloud_V3V4_all_tax[,2:8],1,paste,collapse=";")
new_all_tax <- cbind(EzBioCloud_V3V4_tax[,1],tax_all_collapsed)
write.table(x = new_all_tax, file = EzBioCloud_V3V4_all_taxonomy_Qiime ,quote = FALSE, sep="\t",row.names = FALSE ,col.names = FALSE)


