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
replace_empty_tax <- snakemake@params[["viz_replace_empty_tax"]]
ranks <- unlist(strsplit(snakemake@params[["ranks"]], ","))



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
metadata <- read.delim(file = Metadata_table, sep = "\t", header = TRUE, na.strings = "NA")

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
    taxonomy_table<-taxonomy_table %>% as_tibble() %>% separate(V2, sep=";", ranks)

    ### Replace the not properly named headers into proper ones
    colnames(taxonomy_table)[colnames(taxonomy_table)=="V1"] <- "Feature.ID"
    colnames(taxonomy_table)[colnames(taxonomy_table)=="V3"] <- "Confidence"

    ### Convert taxonomic levels as character (needed for the next steps)
    taxonomy_table[ranks] <- lapply(taxonomy_table[ranks], as.character)

        print(paste("replace empty taxonomy :", replace_empty_tax))

        if(isTRUE(replace_empty_tax)) {
          ### Replace NA by the previous order + a space_holder for each taxonomic level
          for (i in seq_along(ranks)) {
            current_rank <- ranks[i]
            
            if (i == 1) {
              # First rank: fill with "unknown_kingdom"
              taxonomy_table[[current_rank]][is.na(taxonomy_table[[current_rank]])] <- paste0("unknown_", current_rank)
            } else {
              parent_rank <- ranks[i - 1]
              
              # Build default name using parent + suffix
              suffix <- paste0("_", substr(current_rank, 1, 4))  # e.g., "_clas", "_orde"
              
              # Fill NAs in the current rank
              taxonomy_table[[current_rank]][is.na(taxonomy_table[[current_rank]])] <-
                paste0(
                  taxonomy_table[[parent_rank]][is.na(taxonomy_table[[current_rank]])],
                  suffix
                )
            }
          }

          print("table NA remplaced by spaceholders")

          }else{

          print("table NA NOT remplaced by spaceholders")
          }


## Import as physeq object where needed
OTU <- otu_table(count_table, taxa_are_rows = TRUE)
TAX <- taxonomy_table %>% column_to_rownames("Feature.ID") %>% as.matrix() %>% tax_table()
META <- metadata %>% as.data.frame() %>% column_to_rownames("sample") %>% sample_data()


## Merge all in a phyloseq object
phyloseq_obj <- phyloseq(OTU, TAX, META, PHY, SEQS)



## Compute alpha diversity indexes after this filtration
### Remove the previously computed values, in case
#sample_data(phyloseq_obj) <- select(sample_data(phyloseq_obj), -c(Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher, Observed_min_1))

drop <- c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "Observed_min_1")
sample_data(phyloseq_obj) <- sample_data(phyloseq_obj)[,!(names(sample_data(phyloseq_obj)) %in% drop)]


### Add alpha diversity indexes to metadata
alpha_div <- estimate_richness(physeq = phyloseq_obj, split = TRUE, c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson"))
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
