# Title     : Vsearch count table
# Objective : Format a count table from vsearch output
# Created by: valentinscherz
# Created on: 06.06.1p

## Redirect R output
    #log <- file(snakemake@log[[1]], open="wt")
    #sink(log)
    #sink(log, type="message")












# Merge data from multiple runs (if necessary)
   if (length(seq_tab) == 1){
	print("Unique RUN, no merging of seq_tabl")
	st.all <- readRDS(seq_tab)
   }else{
	print("Multiple RUN, merging")
	st.all <- do.call("mergeSequenceTables", lapply(seq_tab, readRDS))
   }

    do.call("rbind", l)
