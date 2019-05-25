# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 16.11.18

# From https://f1000research.com/articles/5-1492/v2 and https://rstudio-pubs-static.s3.amazonaws.com/345955_fba1ccbdcd8f424aa5505c15bfd75bf7.html

## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


## Input
seqs_path <- snakemake@input[["seqs_path"]]
Metadata_table <- snakemake@input[["Metadata_table"]]

## Ouput
tree_path <- snakemake@output[["tree_path"]]


## Parameters

## Load needed libraries
library(DECIPHER)
library(dada2)
library(phangorn)
library(Biostrings)


## Read reprsentative sequences
seqs <- readDNAStringSet(seqs_path)
## Rename the sequences
names(seqs) <- seqs
## Align sequences with DECIPHER
alignment <- AlignSeqs(seqs)

## Construct phylogenetic tree with phangorn
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

bs <- bootstrap.pml(fitGTR, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

## Write tree
write.tree(bs, file=tree_path)
