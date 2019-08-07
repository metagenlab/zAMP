## Adapted from https://bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/ClassifySequences.pdf
# Title     : Decipher - tax prep
# Objective : Prep taxonomy for decipher
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


### Train dataset
    seqs_path <- "EzBioCloud_V3V4_decipher.fasta"
# read the sequences into memory
    seqs <- readDNAStringSet(decipher_seqs)

# specify a path to the taxid file:
    taxid <- NULL

# Pruning the training set
### Not really usefull here since we have already prune to have only once each sequence..But here in case

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

     maxIterations <- 5 # must be >= 1
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


## Summary of the training
#We can view summary properties of the training set (trainingSet) by printing it:
     trainingSet
     pdf(file = training_plot)
     plot(trainingSet)
     dev.off()

    ## Save trainingset object
     saveRDS(object = trainingSet, file = trained_tax)
