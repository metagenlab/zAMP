## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
DB_tax_path <- snakemake@input[["tax"]]
DB_uc_path <- snakemake@input[["uc"]]

## Output
formatted_tax_path <- snakemake@output[["formatted_tax"]]
all_tax_path <- snakemake@output[["all"]]
problematic_tax_path <- snakemake@output[["problematic"]]


## Parameters
n_species <- snakemake@params[["numbers_species"]]
n_genus <- snakemake@params[["numbers_genus"]]

## Load needed libraries
library(tidyr);packageVersion("tidyr")
library(dplyr);packageVersion("dplyr")
library(lattice);packageVersion("lattice")
library(stringr);packageVersion("stringr")
library(forcats);packageVersion("forcats")


## Import data
  ### Taxonomy
  DB_tax <- read.delim(DB_tax_path, header=FALSE, colClasses = "character", as.is = TRUE, strip.white = TRUE)
  ## Dereplication table
  uc_table <- read.delim(DB_uc_path, header = FALSE, as.is = TRUE, colClasses = "character", strip.white = TRUE)



## Pre-process taxonomy
  ### Format in case Silva database
  DB_tax$V2 <- str_remove_all(DB_tax$V2, pattern = "'")
  DB_tax$V2 <- str_remove_all(DB_tax$V2, pattern = "D_\\d{1,2}__")

  ### Split taxonomy based on presence of ";"
  tax_split <- strsplit(DB_tax[,2],split = ";")

  ### Fill in a dataframe with one column per tax level
  tax_table <- t(sapply(tax_split,function(x){x}))
  colnames(tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  ## Rename the row with the seq ids
  tax_table <- data.frame(tax_table)
  tax_table$id <- DB_tax[,1]



## Pre-process dereplication table
  ### Keep clusters, i.e. deredeplicated representative sequences
  C_uc_table <- uc_table[uc_table[,1]=="C",]
  ### Keep replicated sequeces, i.e. sequences corresponding to a given cluster
  H_uc_table <- uc_table[uc_table[,1]=="H",]
  ### Bind all sequences to have all ids
  HC_uc_table <- rbind(C_uc_table,H_uc_table)
  ### Format this table
  HC_uc_table <- data.frame(HC_uc_table)
  HC_uc_table <- HC_uc_table %>%
    rename(
      "clust_id" = V2,
      "seq_id" = V9,
      "clust_rep" = V10
      )
  ### Keep only usefull columns
  HC_uc_table_selec <- select(HC_uc_table, c("clust_id", "seq_id", "clust_rep"))
  ### Set set cluster rep. sequence as its own id (instead of *)
  HC_uc_table_selec$clust_rep[HC_uc_table_selec$clust_rep== "*"] <- HC_uc_table_selec$seq_id[HC_uc_table_selec$clust_rep== "*"]


## Merge taxonomy to clusters and format
  merged_tax_uc <- left_join(HC_uc_table_selec, tax_table, by = c("seq_id"="id"))
  ### Put singletons aside (they don't need to be formatted)
  single_tax <- merged_tax_uc %>% group_by(clust_id) %>% filter(n_distinct(seq_id)==1) %>% ungroup()
  multiple_tax <- merged_tax_uc %>% group_by(clust_id) %>% filter(n_distinct(seq_id)>1) %>% ungroup()

  ### Renomve redundant Genus name in Species ID, keeping only the last word of the vector
  multiple_tax$Species <- sapply(strsplit(as.character(multiple_tax$Species), " "), tail, 1)


## Reformat the tax id helped by dplyr. Here we first group by clust_id and clust_rep to handle each cluster individually. Then using summarise, we aggregate the taxonomy.
    #### From Kingdom to Family levels, we keep only single taxonomy IDs. If not single then, then we put a spaceholder indicating the discrepancy. This discrepancy is passed to all
    #### levels below the discrepant taxonomic level.
    ####
    #### For Genus and species levels, we merge the tax ids if there are an acceptable number of different names (limit set by user). If there more taxonomic names than the limit,
    #### then we use a spaceholder.


### Formatted
formatted_tax <- multiple_tax %>%
  group_by(clust_id, clust_rep) %>%
  summarise(new_Kingdom =
              case_when(# If there is more than one Kindgom, indicate it with a spaceholder
                        n_distinct(Kingdom) > 1 ~ paste0("Disc.King_)", paste(unique(Kingdom), collapse = "/"), "(", n_distinct(Kingdom), ")"),
                        # Else Keep only one kingdom for all sequence of the cluster
                        TRUE ~ paste(unique(Kingdom), collapse = "/")
                        )
            ,
            new_Phylum =
              case_when(# If there is more than one Kindgom, indicate it with a spaceholder here too. Same for the Phylum
                        n_distinct(Kingdom) > 1 ~ paste0("Disc.King_)", paste(unique(Kingdom), collapse = "/"), "(", n_distinct(Kingdom), ")"),
                        n_distinct(Phylum) > 1 ~ paste0("Disc.Phy_)", paste(unique(Phylum), collapse = "/"), "(", n_distinct(Phylum), ")"),
                        TRUE ~ paste(unique(Phylum), collapse = "/"),
                        )
            ,
            new_Class =
              case_when(n_distinct(Kingdom) > 1 ~ paste0("Disc.King_)", paste(unique(Kingdom), collapse = "/"), "(", n_distinct(Kingdom), ")"),
                        n_distinct(Phylum) > 1 ~ paste0("Disc.Phy_)", paste(unique(Phylum), collapse = "/"), "(", n_distinct(Phylum), ")"),
                        n_distinct(Class) > 1 ~ paste0("Disc.Class_)", paste(unique(Class), collapse = "/"), "(", n_distinct(Class), ")"),
                        TRUE ~ paste(unique(Class), collapse = "/"),
                        )
            ,
            new_Order =
              case_when(n_distinct(Kingdom) > 1 ~ paste0("Disc.King_)", paste(unique(Kingdom), collapse = "/"), "(", n_distinct(Kingdom), ")"),
                        n_distinct(Phylum) > 1 ~ paste0("Disc.Phy_)", paste(unique(Phylum), collapse = "/"), "(", n_distinct(Phylum), ")"),
                        n_distinct(Class) > 1 ~ paste0("Disc.Class_)", paste(unique(Class), collapse = "/"), "(", n_distinct(Class), ")"),
                        n_distinct(Order) > 1 ~ paste0("Disc.Ord_", paste(unique(Order), collapse = "/"), "(", n_distinct(Order), ")"),

                        TRUE ~ paste(unique(Order), collapse = "/"),
                        )
            ,
            new_Family =
              case_when(n_distinct(Kingdom) > 1 ~ paste0("Disc.King_)", paste(unique(Kingdom), collapse = "/"), "(", n_distinct(Kingdom), ")"),
                        n_distinct(Phylum) > 1 ~ paste0("Disc.Phy_)", paste(unique(Phylum), collapse = "/"), "(", n_distinct(Phylum), ")"),
                        n_distinct(Class) > 1 ~ paste0("Disc.Class_)", paste(unique(Class), collapse = "/"), "(", n_distinct(Class), ")"),
                        n_distinct(Order) > 1 ~ paste0("Disc.Ord_", paste(unique(Order), collapse = "/"), "(", n_distinct(Order), ")"),

                        TRUE ~ paste(unique(Family), collapse = "/"),
                        )
            ,
            new_Genus =
              case_when(# Same as for above levels
                        n_distinct(Kingdom) > 1 ~ paste0("Disc.King_)", paste(unique(Kingdom), collapse = "/"), "(", n_distinct(Kingdom), ")"),
                        n_distinct(Phylum) > 1 ~ paste0("Disc.Phy_)", paste(unique(Phylum), collapse = "/"), "(", n_distinct(Phylum), ")"),
                        n_distinct(Class) > 1 ~ paste0("Disc.Class_)", paste(unique(Class), collapse = "/"), "(", n_distinct(Class), ")"),
                        n_distinct(Order) > 1 ~ paste0("Disc.Ord_", paste(unique(Order), collapse = "/"), "(", n_distinct(Order), ")"),
                        n_distinct(Family) > 1 ~ paste0("Disc.Fam_", paste(unique(Family), collapse = "/"), "(", n_distinct(Family), ")"),


                        TRUE ~ case_when(# Here we tolerate more than one Genus.
                                         n_distinct(Genus) == 1 ~ paste0(unique(Genus), collapse = ""),
                                         # If there are an acceptable number of, then we collapse them.
                                         n_distinct(Genus) > 1 & n_distinct(Genus) <= n_genus ~ paste0(paste(unique(Genus), collapse = "/"), "(", n_distinct(Genus), ")"),
                                         # Otherwise, we use a spaceholder
                                         n_distinct(Genus) > n_genus ~ paste0("gen.(", n_distinct(Genus),")"),
                                         )
              ),

            new_Species =
              case_when(n_distinct(Kingdom) > 1 ~ paste0("Disc.King_)", paste(unique(Kingdom), collapse = "/"), "(", n_distinct(Kingdom), ")"),
                        n_distinct(Phylum) > 1 ~ paste0("Disc.Phy_)", paste(unique(Phylum), collapse = "/"), "(", n_distinct(Phylum), ")"),
                        n_distinct(Class) > 1 ~ paste0("Disc.Class_)", paste(unique(Class), collapse = "/"), "(", n_distinct(Class), ")"),
                        n_distinct(Order) > 1 ~ paste0("Disc.Ord_", paste(unique(Order), collapse = "/"), "(", n_distinct(Order), ")"),
                        n_distinct(Family) > 1 ~ paste0("Disc.Fam_", paste(unique(Family), collapse = "/"), "(", n_distinct(Family), ")"),

                        TRUE ~ case_when(# Same as for Genus, only that here we also indicate the Genus in Species.
                                        n_distinct(Genus) > 1 ~ case_when(n_distinct(Species) <= n_species ~ paste0(paste(unique(Genus), unique(Species), collapse = "/"), "(", n_distinct(Species), ")"),
                                                                           n_distinct(Species) > n_species ~ paste0("sp.(", n_distinct(Species),")"),
                                                                            ),
                                         TRUE ~ case_when(n_distinct(Species) <= n_species ~ paste0(paste(unique(Genus), collapse = ""), " ", paste(unique(Species), collapse = "/"), "(", n_distinct(Species), ")"),
                                                                                            n_distinct(Species) > n_species ~ paste0("sp.(", n_distinct(Species),")")
                                                                                            )
                                                                            )
                ),

            problematic_tax = case_when(n_distinct(Kingdom) > 1 ~ paste0("Disc.King_)", paste(unique(Kingdom), collapse = "/"), "(", n_distinct(Kingdom), ")"),
                                        n_distinct(Phylum) > 1 ~ paste0("Disc.Phy_)", paste(unique(Phylum), collapse = "/"), "(", n_distinct(Phylum), ")"),
                                        n_distinct(Class) > 1 ~ paste0("Disc.Class_)", paste(unique(Class), collapse = "/"), "(", n_distinct(Class), ")"),
                                        n_distinct(Order) > 1 ~ paste0("Disc.Ord_", paste(unique(Order), collapse = "/"), "(", n_distinct(Order), ")"),
                                        n_distinct(Family) > 1 ~ paste0("Disc.Fam_", paste(unique(Family), collapse = "/"), "(", n_distinct(Family), ")")
                                        )
                ) %>%
  ungroup()


### Version without collapse
raw_tax <- multiple_tax %>%
  group_by(clust_id, clust_rep) %>%
  summarise(new_Kingdom = paste0(unique(Kingdom), collapse = "/"),
            new_Phylum =  paste0(unique(Phylum), collapse = "/"),
            new_Class = paste0(unique(Class), collapse = "/"),
            new_Order = paste0(unique(Order), collapse = "/"),
            new_Family = paste0(unique(Family), collapse = "/"),
            new_Genus = paste0(unique(Genus), collapse = "/"),
            new_Species = paste0(paste(unique(Genus), collapse = ""), " ", paste(unique(Species), collapse = "/"), "(", n_distinct(Species), ")"),
            problematic_tax = case_when(n_distinct(Kingdom) > 1 ~ paste0("Disc.King_)", paste(unique(Kingdom), collapse = "/"), "(", n_distinct(Kingdom), ")"),
                                        n_distinct(Phylum) > 1 ~ paste0("Disc.Phy_)", paste(unique(Phylum), collapse = "/"), "(", n_distinct(Phylum), ")"),
                                        n_distinct(Class) > 1 ~ paste0("Disc.Class_)", paste(unique(Class), collapse = "/"), "(", n_distinct(Class), ")"),
                                        n_distinct(Order) > 1 ~ paste0("Disc.Ord_", paste(unique(Order), collapse = "/"), "(", n_distinct(Order), ")"),
                                        n_distinct(Family) > 1 ~ paste0("Disc.Fam_", paste(unique(Family), collapse = "/"), "(", n_distinct(Family), ")")
            )
            ) %>%
  ungroup()






## Generate output
### Problematic taxonomy
prob_tax <- raw_tax %>% filter(!is.na(problematic_tax))
write.csv(x = prob_tax, file = problematic_tax_path, quote = FALSE, row.names = FALSE, sep = ";")

## Singeltons
taxa_col <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
single_tax_collapsed <- single_tax %>%
  unite(remove = TRUE,  col = "taxonomy", taxa_col, sep = ";") %>%
  select(-c(seq_id, clust_id))

### Formatted taxonomy
taxa_col <- c("new_Kingdom","new_Phylum","new_Class","new_Order","new_Family","new_Genus","new_Species")
formatted_tax_collapsed <- formatted_tax %>%
  unite(remove = TRUE,  col = "taxonomy", taxa_col, sep = ";") %>%
  select(-c(problematic_tax, clust_id))
formatted_tax_collapsed_w_s <- rbind(single_tax_collapsed, formatted_tax_collapsed)
write.table(x = formatted_tax_collapsed_w_s, file = formatted_tax_path ,quote = FALSE, sep="\t",row.names = FALSE ,col.names = FALSE)

### Unformatted taxonomy
raw_tax_collapsed <- raw_tax %>%
  unite(remove = TRUE,  col = "taxonomy", taxa_col, sep = ";") %>%
  select(-c(problematic_tax, clust_id))
raw_tax_collapsed_w_s <- rbind(single_tax_collapsed, raw_tax_collapsed)
write.table(x = raw_tax_collapsed_w_s, file = all_tax_path ,quote = FALSE, sep="\t",row.names = FALSE ,col.names = FALSE)

