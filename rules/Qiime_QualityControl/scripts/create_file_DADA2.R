## Redirect R output
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

## Input
    DADA2_seq <- snakemake@input[["DADA2_seq"]]

## Output
    rep_seqs <- snakemake@output[["rep_seqs"]]

# Parameters
    sample <- snakemake@params[["sample"]]
    print(sample)

# Load needed libraries
    library(dada2); packageVersion("dada2")
    library(Rcpp); packageVersion("Rcpp")


# Reading the file
    DADA2 <- readRDS(DADA2_seq)

# Export reads
    seqs <- DADA2[[sample]][["sequence"]]

    ## giving our seq headers more manageable names (ASV_1, ASV_2...)
    asv_headers <- vector(length(seqs), mode="character")
    for (i in 1:length(seqs)) {
      asv_headers[i] <- paste(">ASV", i, sep="_")
    }

    ## making and writing out a fasta of our final ASV seqs
    asv_fasta <- c(rbind(asv_headers, seqs))
    write(asv_fasta, rep_seqs)


