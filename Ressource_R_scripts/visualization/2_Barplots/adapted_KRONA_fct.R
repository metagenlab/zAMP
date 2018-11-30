
## Adapter function https://rdrr.io/github/cpauvert/psadd/man/plot_krona.html.
## This adapted version goes form the melted df instead of the phyloseq object in the original version of the function


adapted_KRONA_fct <- function(melted_dataframe, grouping_column, r_figures, x_axis_column){
  
  df0 <- melted_dataframe

for (i in (unique(df0[[grouping_column]]))){
  
  df <- filter(df0, df0[[grouping_column]] == i)
  
  
  df <- df[, c("Abundance", x_axis_column, "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")]
  df <- as.data.frame(unclass(df))
  df[, 2] <- gsub(" |\\(|\\)", "", df[, 2])
  df[, 2] <- as.factor(df[, 2])
  dir.create(file.path(r_figures,"/",i))
  for (lvl in levels(df[, 2])) {
    write.table(unique(df[which(df[, 2] == lvl), -2]), file = paste0(r_figures,"/",i, "/", lvl, "taxonomy.txt"), sep = "\t", row.names = F, col.names = F, na = "", quote = F)
  }
  
  krona_args <- paste(r_figures,"/", i, "/", levels(df[, 2]), "taxonomy.txt,", levels(df[, 2]), sep = "", collapse = " ")
  output <- paste(r_figures,"/",i,".html", sep = "")
  system(paste("ktImportText", krona_args, "-o", output, sep = " "))
  
}}


