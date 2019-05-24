# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 29.10.18

## Inspired from https://gist.github.com/erictleung/eda33bceceebb190eb9a44c93a077d32
## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
taxonomy_table <- snakemake@input[["tax_table"]]

## Output
full_taxonomy_table <- snakemake@output[["full_taxonomy_table"]]

#Create OTU table with taxonomic classification from EzbioCloud database
rdp_table<-read.table(taxonomy_table,as.is=T,sep="\t")
#creating tax table
get_tax<-function(x) {
  tax_split<-strsplit(x,";")[[1]]
  N_tax_fields<-length(tax_split)
  tax_vector<-c(tax_split,rep(NA,7-N_tax_fields))
  return(tax_vector)
}
tax_table<-matrix(ncol=7,nrow=nrow(rdp_table))
for (xx in 1:nrow(rdp_table)){
  tax_table[xx,]<-get_tax(rdp_table[xx,2])
}

colnames(tax_table)<-c("kingdom","phylum","class","order","family","genus","species")
tax_table <- as.data.frame(cbind("OTU_ID"=rdp_table$V1, tax_table))


# and now save the resulting table for further use
write.table(x = tax_table, file = full_taxonomy_table, append = F, quote = F, sep = "\t", eol = "\n", row.names = F, col.names = T)

