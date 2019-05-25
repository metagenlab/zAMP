
### Create the function, from "plot_richness" of phyloseq but with modifications to allows facetting by values but also for a factor in data
modif_plot_richness_fct <- function (physeq, x = "samples", color = NULL, shape = NULL, 
                                     title = NULL, scales = "free_y", nrow = NULL, shsi = NULL, measures = NULL, 
                                     sortby = NULL, facetting_column = NULL, r_figures , boxplot = TRUE, jitter = TRUE) 
{
  erDF = estimate_richness(physeq, split = TRUE, measures = measures)
  measures = colnames(erDF)
  ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
  measures = measures[!measures %in% ses]
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    DF <- data.frame(erDF, sample_data(physeq))
  }
  else {
    DF <- data.frame(erDF)
  }
  if (!"samples" %in% colnames(DF)) {
    DF$samples <- sample_names(physeq)
  }
  if (!is.null(x)) {
    if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
      x <- "samples"
    }
  }
  else {
    x <- "samples"
  }
  mdf = reshape2::melt(DF, measure.vars = measures)
  mdf$se <- NA_integer_
  if (length(ses) > 0) {
    selabs = ses
    names(selabs) <- substr(selabs, 4, 100)
    substr(names(selabs), 1, 1) <- toupper(substr(names(selabs), 
                                                  1, 1))
    mdf$wse <- sapply(as.character(mdf$variable), function(i, 
                                                           selabs) {
      selabs[i]
    }, selabs)
    for (i in 1:nrow(mdf)) {
      if (!is.na(mdf[i, "wse"])) {
        mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      }
    }
    mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
  }
  if (!is.null(measures)) {
    if (any(measures %in% as.character(mdf$variable))) {
      mdf <- mdf[as.character(mdf$variable) %in% measures, 
                 ]
    }
    else {
      warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
    }
  }
  if (!is.null(shsi)) {
    warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
  }
  if (!is.null(sortby)) {
    if (!all(sortby %in% levels(mdf$variable))) {
      warning("`sortby` argument not among `measures`. Ignored.")
    }
    if (!is.discrete(mdf[, x])) {
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if (all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[, 
                                                                x])) {
      wh.sortby = which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x], levels = names(sort(tapply(X = mdf[wh.sortby, 
                                                                      "value"], INDEX = mdf[wh.sortby, x], mean, na.rm = TRUE, 
                                                              simplify = TRUE))))
    }
  }
  richness_map <<- aes_string(x = x, y = "value", colour = color, 
                            shape = shape)
  mdf <<- mdf
  
  p <- ggplot(mdf, richness_map) + geom_point(na.rm = TRUE) 
  if (any(!is.na(mdf[, "se"]))) {
    p = p + geom_errorbar(aes(ymax = value + se, ymin = value - 
                                se), width = 0.1)
  }
  
  if (isTRUE(boxplot)) {
    p = p + geom_boxplot()
  }
  if (isTRUE(jitter)) {
    p = p + geom_boxplot() + geom_jitter(width = 0.1)
  }
  
  
  p = p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, 
                                           hjust = 0))
  p = p + ylab("Alpha Diversity Measure")
  
  if (!is.null(facetting_column)){
    p = p + facet_grid(cols = vars(get(facetting_column)), rows = (vars(variable)), scales = scales)
  }
  else if (is.null(facetting_column)){
    p = p + facet_grid(rows = (vars(variable)), scales = scales)
  }
  
  if (!is.null(title)) {
    p = p + ggtitle(title)

    
  }
  ggsave(p, filename = paste0(r_figures ,"/Diversity/overall_richness.png"), width = 5, height = 10)
  return(p)
  
}





