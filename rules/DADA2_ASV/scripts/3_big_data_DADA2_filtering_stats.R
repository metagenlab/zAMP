# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 03.01.19
# Modified from :https://benjjneb.github.io/dada2/bigdata_paired.html

# Redirect R output
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

# Input
    q_filtering_stats <- snakemake@input[["q_filtering_stats"]]
    run_stats <- snakemake@input[["run_stats"]]
    no_chim <- snakemake@input[["no_chim"]]
    l_filtered <- snakemake@input[["length_filtered"]]

# Output
    filtering_stats <- snakemake@output[["filtering_stats"]]

# Load needed libraries
    library(dada2); packageVersion("dada2")
    library(dplyr); packageVersion("dplyr")

# Load the q score filtration R stats
    filtration <- do.call("rbind", lapply(q_filtering_stats, readRDS))

# Load the forward, reverse and paired reads correction stats for each run
    dada_infer <- do.call("rbind", lapply(run_stats, readRDS))

# Load reads with filtered chimera
    no_chimera <- do.call("rowSums", lapply(no_chim, readRDS))
    no_chimera <- as.data.frame(no_chimera)
    no_chimera$Sample <- row.names(no_chimera)

# Load length_filtered sequences
    length_filtered <- do.call("rowSums", lapply(l_filtered, readRDS))
    length_filtered <- as.data.frame(length_filtered)
    length_filtered$Sample <- row.names(length_filtered)

# Merge all dataframe together
    all_stats <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Sample", all.x = TRUE),
       list(filtration, dada_infer, no_chimera, length_filtered))

# set RUN at the very end of the table
    all_stats <- all_stats%>%select(-RUN,everything())

# Write the  stat table
write.table(x = all_stats, file = filtering_stats, sep = "\t", row.names = FALSE)



