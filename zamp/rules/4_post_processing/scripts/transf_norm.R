# Title     : Normalize and transform counts
# Objective : Normalize and transform counts
# Created by: valentinscherz
# Created on: 16.09.2019

## Redirect R output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]

## Output
phyloseq_norm <- snakemake@output[["phyloseq_norm"]]

## Parameters
#lower all normalization to ease matching
normalization <- tolower(snakemake@params[["normalization"]])
min_abundance <- as.numeric(snakemake@params[["min_abundance_value"]])
min_prevalence <- as.numeric(snakemake@params[["min_prevalence_value"]])
print(paste("normalization setting is", normalization))

## Load needed libraries
library("phyloseq");packageVersion("phyloseq")
library("vegan");packageVersion("vegan")
library("edgeR");packageVersion("edgeR")
library("metagenomeSeq");packageVersion("metagenomeSeq")

## Import phyloseq
physeq <- readRDS(phyloseq_object)

## Recover OTU table
OTU <- data.frame(otu_table(physeq), check.names = FALSE)

## List methods computed by vegan
vegan_methods <- c("total", "max", "freq", "normalize", "range", "pa", "chi.square", "hellinger" ,"log")

## CLR must we computed sample-wise, on rows. In Phyloseq we have taxa_are_rows = TRUE. Thus, must be computed on the t() if OTU table
if (normalization == "clr"){
  print("clr normalization")
  OTU1 <- OTU + 1
  reads_trfs <- t(chemometrics::clr(t(OTU1)))

## Vegan decostand normalizations (with double t() to take in account the transposed structure of the counts used here, modified from https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/decostand)
}else if (normalization %in% vegan_methods){
  print(paste(normalization, "normalization by vegan"))
  reads_trfs <- t(vegan::decostand(t(OTU), method = normalization))

## Percent normalization with basic R
}else if (normalization == "pct"){
  print("Percent normalization with base R")
  reads_trfs <- sapply(OTU, function(x) 100*x/sum(x))
  rownames(reads_trfs) <- rownames(OTU)

## CSS normalization with metagenomeseq (modified form https://bioconductor.org/packages/release/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf )
}else if (normalization == "css"){
  print("css normalization by metagenomeseq")
  OTU_mtgs <- metagenomeSeq::newMRexperiment(counts = OTU, phenoData = NULL, featureData = NULL, libSize = NULL, normFactors = NULL)
  p <- metagenomeSeq::cumNormStat(OTU_mtgs)
  OTU_norm <-  metagenomeSeq::cumNorm(OTU_mtgs, p = p)
  reads_trfs = as.data.frame(metagenomeSeq::MRcounts(OTU_norm, norm = TRUE, log = FALSE))

## TMM normalization with edgeR (modified from https://joey711.github.io/phyloseq-extensions/edgeR.html)
} else if (normalization == "tmm"){
  print("TMM normalization by edgeR")
  OTU1 <- OTU + 1
  group <- rownames(get_variable(physeq))
  tax <- tax_table(physeq, errorIfNULL=FALSE)
  y <- edgeR::DGEList(counts=OTU1, group=group, genes=tax, remove.zeros = TRUE)
  z = edgeR::calcNormFactors(y, method="TMM")
  reads_trfs <- data.frame(z$counts, check.names = FALSE)
  #reads_trfs[reads_trfs==0.5] <- 0

} else if (normalization == "none"){
  print('No normalization (normalization = "none"')
  reads_trfs <- OTU

}else{
  stop(paste('"normalization" was', normalization, '.Must be one of : "clr", "pct", "css", "tmm", "total", "max", "freq", "normalize", "range", "pa", "chi.square", "hellinger" or "log", ')) 
}
  
## Filter based on counts and prevalence
### Filter based on the cumulated abundance
if(is.numeric(min_abundance)){
  print(paste("filter based on cumulated min_abundance:", min_abundance))
  filtered_data <- reads_trfs[,colSums(reads_trfs) > min_abundance, drop = FALSE]
}else if (min_abundance == "none"){
  print("NO filtering based on cumulated abundance")
  filtered_data <- reads_trfs
}else{
  stop(paste(min_abundance, 'must be a numeric value or "none')) 
}

### Filter based on prevalence
if(is.numeric(min_prevalence)){
  print(paste("filter based on prevalence:",min_prevalence))
  prevalence_cutoff <- (min_prevalence/100) * nsamples(physeq)
  print(paste("Features must be found in more than", prevalence_cutoff, "samples to be retained"))
  filtered_data <- filtered_data[,colSums(filtered_data != 0) > prevalence_cutoff, drop = FALSE]
}else if (min_prevalence == "none"){
  print("NO filtering based on cumulated abundance")
}else{
  stop(paste(min_prevalence, 'must be a numeric value or "none'))
}

## Repopulate physeq object (identical script than for taxa filter)

physeq_filtered <- physeq
sample_names(physeq_filtered)
head(filtered_data)
otu_table(physeq_filtered) <- otu_table(round(filtered_data, digits = 0), taxa_are_rows = TRUE)
filtered_taxa <- prune_taxa(taxa_sums(physeq_filtered) > 0, physeq_filtered) ## Removes taxa not at least present in one sample


## Recompute alpha diversity indexes after this filtration
### Remove the previously computed values
#sample_data(filtered_taxa) <- select(sample_data(filtered_taxa), -c(Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher, Observed_min_1))
drop <- c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "Observed_min_1")
sample_data(filtered_taxa) <- sample_data(filtered_taxa)[,!(names(sample_data(filtered_taxa)) %in% drop)]


### Add alpha diversity indexes to metadata
alpha_div <- estimate_richness(physeq = filtered_taxa, split = TRUE, c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson"))
sample_data(filtered_taxa) <- cbind(sample_data(filtered_taxa),alpha_div)


### In addition, add the Observed over 1% metric
#### Keep the IDSs of the taxa above 1%
physeqrF <- filter_taxa(physeq = filtered_taxa, function(x) mean(x) > 0.01, FALSE)
#### Keep only those
physeqaF <- prune_taxa(physeqrF,filtered_taxa)
#### Calculate the Observed index
alpha_div_1 <- estimate_richness(physeq = physeqaF, split = TRUE, measure = "Observed")
#### Rename this index
colnames(alpha_div_1) <- paste0(colnames(alpha_div_1), ("_min_1"))
#### Again bind this new column
sample_data(filtered_taxa) <- cbind(sample_data(filtered_taxa),alpha_div_1)


## Write the new phyloseq object
saveRDS(object = filtered_taxa, file = phyloseq_norm)
