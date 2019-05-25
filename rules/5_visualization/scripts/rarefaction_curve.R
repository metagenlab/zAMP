# Title     : TODO
# Objective : TODO
# Created by: valentinscherz
# Created on: 21.11.18


# inspired from : https://github.com/joey711/phyloseq/issues/143

## Redirect R output to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


## Input
phyloseq_object <- snakemake@input[["phyloseq_object"]]

## Ouput
rarefaction_curve <- snakemake@output[["rarefaction_curve"]]

## Parameters
sample_label <- snakemake@params[["sample_label"]]
sample_type <- snakemake@params[["sample_type"]]

## Load libraries
library('phyloseq')
library('ggplot2')
library('plyr') # ldply
library('reshape2') # melt
library('vegan')

## Load the phyloseq object
phyloseq_obj <- readRDS(phyloseq_object)


## Create a modified version of estimate_richness function of phyloseq to accept sample names containing numbers only

modified_estimated_richness <- function (physeq, split = TRUE, measures = NULL)
{
  if (!any(otu_table(physeq) == 1)) {
    warning("The data you have provided does not have\n",
            "any singletons. This is highly suspicious. Results of richness\n",
            "estimates (for example) are probably unreliable, or wrong, if you have already\n",
            "trimmed low-abundance taxa from the data.\n", "\n",
            "We recommended that you find the un-trimmed data and retry.")
  }
  if (!split) {
    OTU <- taxa_sums(physeq)
  }
  else if (split) {
    OTU <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(physeq)) {
      OTU <- t(OTU)
    }
  }
  renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
                "InvSimpson", "Fisher")
  names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon",
                        "simpson", "invsimpson", "fisher")
  if (is.null(measures)) {
    measures = as.character(renamevec)
  }
  if (any(measures %in% names(renamevec))) {
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in%
                                                            measures]
  }
  if (!any(measures %in% renamevec)) {
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  outlist = vector("list")
  estimRmeas = c("Chao1", "Observed", "ACE")
  if (any(estimRmeas %in% measures)) {
    outlist <- c(outlist, list(t(data.frame(estimateR(OTU), check.names = FALSE))))
  }
  if ("Shannon" %in% measures) {
    outlist <- c(outlist, list(shannon = diversity(OTU, index = "shannon")))
  }
  if ("Simpson" %in% measures) {
    outlist <- c(outlist, list(simpson = diversity(OTU, index = "simpson")))
  }
  if ("InvSimpson" %in% measures) {
    outlist <- c(outlist, list(invsimpson = diversity(OTU, index = "invsimpson")))
  }
  if ("Fisher" %in% measures) {
    fisher = tryCatch(fisher.alpha(OTU, se = TRUE), warning = function(w) {
      warning("phyloseq::estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
      suppressWarnings(fisher.alpha(OTU, se = TRUE)[, c("alpha",
                                                        "se")])
    })
    if (!is.null(dim(fisher))) {
      colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
      outlist <- c(outlist, list(fisher))
    }
    else {
      outlist <- c(outlist, Fisher = list(fisher))
    }
  }


  out = do.call("cbind", outlist)
  namechange = intersect(colnames(out), names(renamevec))
  colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
  colkeep = sapply(paste0("(se\\.){0,}", measures), grep, colnames(out),
                   ignore.case = TRUE)
  out = out[, sort(unique(unlist(colkeep))), drop = FALSE]
  out <- as.data.frame(out)
  return(out)
}



# Replicate curve function from : https://github.com/joey711/phyloseq/issues/143

set.seed(42)

calculate_rarefaction_curves <- function(psdata, measures, depths) {

  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)

    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)

    alpha_diversity <- modified_estimated_richness(rarified_psdata, measures = measures)

    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')

    molten_alpha_diversity
  }

  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))

  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]

  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(psdata = phyloseq_obj, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"), depth = rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
summary(rarefaction_curve_data)


####

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))


####

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(phyloseq_obj)), by.x = 'Sample', by.y = 'row.names')

rarefaction_curve_data_summary_verbose

p <- ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = sample_type,
    group = get("Sample"))) +
      labs(col = "Sample type") +
      geom_line(size = 0.2) +
      geom_pointrange(size = 0.2) +
      facet_wrap(facets = ~ Measure, scales = 'free_y')


ggsave(filename = rarefaction_curve,  plot = p, width = 10)
