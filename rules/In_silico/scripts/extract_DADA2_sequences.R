# Title     : Extract fasta from DADA2 merged sequences count table
# Objective : Format a count table from vsearch output
# Created by: valentinscherz
# Created on: 06.06.1p

## Redirect R output
    #log <- file(snakemake@log[[1]], open="wt")
    #sink(log)
    #sink(log, type="message")


## Input
    seq_tab <- snakemake@input[["dada2_rds"]]

## Output
    no_chim <- snakemake@output[["no_chim"]]
    length_filtered <- snakemake@output[["length_filtered"]]
    rep_seqs <- snakemake@output[["rep_seqs"]]
    count_table <- snakemake@output[["count_table"]]
    length_histo <- snakemake@output[["length_histo"]]

# Parameters
  merged_min_length <- snakemake@params[["merged_min_length"]]
  merged_max_length  <- snakemake@params[["merged_max_length"]]

# Load needed libraries
  library(dada2); packageVersion("dada2")
  library(ggplot2); packageVersion("ggplot2")


## Import
    st.all <- makeSequenceTable(readRDS(seq_tab))



#save.image(file= paste0(getwd(), "/myEnvironment.RData"))

# Remove chimeras
    print("Filter chimera")
    seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=snakemake@threads, verbose=TRUE)
    print("Chimera filtered")

# Sequences length inspection and filtration
#### That's the little added trick, the reason why we are using this script and not the one in Qiime2. Indeed we typically are here keeping only sequences between 390 and 500 bp of length after merging. This tcorresponds to the expected length of the V3V4 region of the 16S rRNA gene.
## Inspect distribution of sequence lengths
    unfiltered_length_table <- tapply(colSums(seqtab), factor(nchar(getSequences(seqtab))), sum)
    unfiltered_length_table
    unfiltered_length_table_df <- data.frame(length=as.numeric(names(unfiltered_length_table)), count=unfiltered_length_table)
    png(length_histo)
    p <- ggplot(data=unfiltered_length_table_df, aes(x=length, y=count)) + geom_col(binwidth=20)
    print(p)
    dev.off()

# Filter based on length
    seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(merged_min_length, merged_max_length)]
# Inspect distribution of sequence lengths after filtration
    table(nchar(getSequences(seqtab2)))

# Export reads and count
#### We are writing in files the product of this DADA2 process. These are one .fasta file contanining the dereplicated, errors corrected, paired-end merged representative sequences and one .txt file indicating the prevalence of sequencne in each sample (this is the result of dereplication).
    asv_seqs <- names(seqtab2)

## giving our seq headers more manageable names (ASV_1, ASV_2...)
    asv_headers <- vector(length(names(seqtab2)), mode="character")
    for (i in 1:length(names(seqtab2))) {
      asv_headers[i] <- paste(">ASV", i, sep="_")
    }

## making and writing out a fasta of our final ASV seqs:
    asv_fasta <- c(rbind(asv_headers, asv_seqs))
    write(asv_fasta, rep_seqs)

## count table:
    print("Create count table")
    asv_tab <- seqtab2

## Renamed

    names(asv_tab) <- sub(">", "", asv_headers)
    asv_tab <- data.frame(names(asv_tab), asv_tab)

    write.table(asv_tab, count_table , sep="\t", col.names = FALSE)

## Write sequences objects in .rds for use in statistics
    ## Before length filtration
    saveRDS(seqtab, no_chim)
    ## After length filtration
    saveRDS(seqtab2, length_filtered)


