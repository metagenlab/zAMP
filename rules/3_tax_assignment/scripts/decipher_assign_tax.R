## Adapted from https://bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/ClassifySequences.pdf
# Title     : Decipher - tax assign
# Objective : Assign tax with decipher
# Created by: valentinscherz
# Created on: 07.08.19

## Redirect R output to the log file
   log <- file(snakemake@log[[1]], open="wt")
   sink(log)
   sink(log, type="message")

## Input
    decipher_seqs <- snakemake@input[["decipher_seqs"]]

## Output
    trained_tax <- snakemake@output[["trained_tax"]]
    training_plot <- snakemake@output[["training_plot"]]

## Load needed libraries
    library(DECIPHER);packageVersion("DECIPHER")
    library(dplyr);packageVersion("dplyr")
    library(tidyr);packageVersion("tidyr")
    library(stringr);packageVersion("stringr")

## Adapted from

library(DECIPHER)
library(dplyr)
library(tidyr)
library(stringr)

ref_tax <- "EzBioCloud_V3V4_taxonomy.txt"

## Read taxonomy data
tax_table <- read.table(file=ref_tax, sep="\t", stringsAsFactors=FALSE)
split_tax_table<-tax_table %>% separate(V2, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

## Extract and format taxa
taxa_names <- colnames(split_tax_table[,-1])
taxa_letters <- paste0(str_to_lower(substr(taxa, start = 1, stop = 1)), "__")

for (i in colnames(split_tax_table[,-1])){
  print(i)
  l <- str_to_lower(substr(i, start = 1, stop = 1))
  split_tax_table[[i]] <- paste(l, split_tax_table[[i]], sep = "__")
}

united_tax_table <- split_tax_table %>% unite(taxonomy, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";", remove = TRUE) %>% select(-V1)
ranks <- unlist(united_tax_table)

ranks <-strsplit(ranks, ";", fix=T)
taxa <- setNames(taxa_names,taxa_letters)

count <- 1L
groups <- "Root"
index <- -1L
level <- 0L
rank <- "rootrank"
pBar <- txtProgressBar(style=3)
for (i in seq_along(ranks)) {
  for (j in seq_along(ranks[[i]])) {
    rank_level <- taxa[substring(ranks[[i]][j], 1, 3)]
    group <- substring(ranks[[i]][j], 4)
    w <- which(groups==group & rank==rank_level)
    if (length(w) > 0) {
      parent <- match(substring(ranks[[i]][j - 1], 4),
                      groups)
      if (j==1 || any((parent - 1L)==index[w]))
        next # already included
    }
    count <- count + 1L
    groups <- c(groups, group)
    if (j==1) {
      15
      index <- c(index, 0)
    } else {
      parent <- match(substring(ranks[[i]][j - 1], 4),
                      groups)
      index <- c(index,
                 parent - 1L)
    }
    level <- c(level, j)
    rank <- c(rank, taxa[j])
  }
  setTxtProgressBar(pBar, i/length(ranks))
}

groups <- gsub("^[ ]+", "", groups)
groups <- gsub("[ ]+$", "", groups)
taxid <- paste(0:(length(index) - 1L), groups, index, level, rank, sep="*")
head(taxid, n=10)

writeLines(taxid, con="decipher_taxa_file")


########Format fasta

################### Prep fasta table ########################
split_tax_table<-tax_table %>% separate(V2, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

## Format data
split_tax_table$Root <- "Root"
split_tax_table_f <- split_tax_table %>% unite(tax, c("Root", "Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = "; ", remove = TRUE)
split_tax_table_f <- split_tax_table_f %>% unite(taxonomy, c("V1", "tax"), sep = " ",remove = FALSE) %>% select(-tax)
split_tax_table_f <- split_tax_table_f[,c(2,1)]

## Define a function taken from phylotools https://github.com/cran/phylotools/blob/master/R/rename.fasta.R

rename.fasta <- function(infile = NULL, ref_table,
                         outfile = "renamed.fasta"){
  #fasta <- read.fasta(infile)
  fastaFile <- readDNAStringSet(infile)
  seq.name = names(fastaFile)
  seq.text = paste(fastaFile)
  fasta <- data.frame(seq.name, seq.text)


  ## Convert the input ref to dataframe
  ## usually the file was obtained by
  ## save the result of get.names.fasta to a csv file.
  ## edit the csv file, and provide the name to be replaced.

  colnames(ref_table)  <- c("old_name", "new_name")
  res <- merge(x = fasta, y = ref_table, by.x = "seq.name",
               by.y = "old_name", all.x = TRUE) ### Order of the sequence will change because of merge
  rename <- rep(NA, nrow(res))

  ### if the sequence was not found in the ref_table,
  ### keep the old name

  for(i in 1:nrow(res)){
    rename[i] <- ifelse(is.na(res$new_name[i]),
                        paste("old_name", "_",as.character(res$old_name[i]), sep = ""),
                        as.character(res$new_name[i]))
  }
  writeLines(paste(">", rename, "\n", as.character(res$seq.text), sep = ""), outfile )
  cat(paste(outfile, "has been saved to ", getwd(), "\n" ))
}


## Write data
rename.fasta(infile = "EzBioCloud_V3V4.fasta", ref_table = split_tax_table_f, outfile = "EzBioCloud_V3V4_decipher.fasta")




### Train dataset

# specify the path to your file of training sequences:
seqs_path <- "EzBioCloud_V3V4_decipher.fasta"
# read the sequences into memory
seqs <- readDNAStringSet(seqs_path)
# NOTE: use readRNAStringSet for RNA sequences

# (optionally) specify a path to the taxid file:
rank_path <- "decipher_taxa_file"
taxid <- read.table(rank_path,
                      header=FALSE,
                      col.names=c('Index', 'Name', 'Parent', 'Level', 'Rank'),
                      sep="*", # asterisks delimited
                      quote="", # preserve quotes
                     stringsAsFactors=FALSE)

taxid <-NULL

seqs <- RemoveGaps(seqs)

# Pruning the training set
### Not really usefull here since we have already prune to have only once each sequence..

#obtain the taxonomic assignments
groups <- names(seqs) # sequence names
# assume the taxonomy begins with 'Root;'
groups <- gsub("(.*)(Root;)", "\\2", groups) # extract the group label
groupCounts <- table(groups)
u_groups <- names(groupCounts) # unique groups
length(u_groups) # number of groups

maxGroupSize <- 10 # max sequences per label (>= 1)
remove <- logical(length(seqs))
for (i in which(groupCounts > maxGroupSize)) {
 index <- which(groups==u_groups[i])
 keep <- sample(length(index),
                maxGroupSize)
 remove[index[-keep]] <- TRUE
}
sum(remove) # number of sequences eliminated


## Train dataset

 maxIterations <- 3 # must be >= 1
 allowGroupRemoval <- FALSE
 probSeqsPrev <- integer() # suspected problem sequences from prior iteration
 for (i in seq_len(maxIterations)) {
  cat("Training iteration: ", i, "\n", sep="")
  # train the classifier
  trainingSet <- LearnTaxa(seqs[!remove],
                           names(seqs)[!remove],
                           taxid)
  # look for problem sequences
  probSeqs <- trainingSet$problemSequences$Index
  if (length(probSeqs)==0) {
    cat("No problem sequences remaining.\n")
    break
  } else if (length(probSeqs)==length(probSeqsPrev) &&
             all(probSeqsPrev==probSeqs)) {
    cat("Iterations converged.\n")
    break
  }
  if (i==maxIterations)
    break
  probSeqsPrev <- probSeqs
  # remove any problem sequences
  index <- which(!remove)[probSeqs]
  remove[index] <- TRUE # remove all problem sequences
  if (!allowGroupRemoval) {
    # replace any removed groups
    missing <- !(u_groups %in% groups[!remove])
    missing <- u_groups[missing]
    if (length(missing) > 0) {
      index <- index[groups[index] %in% missing]
      remove[index] <- FALSE # don't remove
    }
  }
}
 sum(remove) # total number of sequences eliminated
 length(probSeqs) # number of remaining problem sequences

 input_seq <- readDNAStringSet("dna-sequences.fasta")

 test_taxa <- IdTaxa(test = input_seq, trainingSet = trainingSet, type = "collapsed", strand = "both", threshold = 60, bootstraps = 100, processors = 2, verbose = TRUE, minDescend = 0.98)
 ?IdTaxa

 taxa_table <- data.frame(cbind(Row.Names = names(rownames(test_taxa)), test_taxa))

 gsub("\\[.\\]", "", taxa_table$test_taxa)

 str_remove_all(taxa_table$test_taxa, "\\[.*\\]")
