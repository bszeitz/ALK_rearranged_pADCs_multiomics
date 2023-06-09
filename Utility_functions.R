##############################
#
# Utility functions
# Author: Bea Szeitz
#
##############################



##############################
# Read in Excel sheets formatted as 1.3 gct files
# 
# Arguments for this function:
# - excel_file: the Excel file name with extension
# - sheet: the sheet number or sheet name within the Excel file
#
# This function returns a list contaning 3 separate tables:
# (1) sample annotations, (2) expression table, (3) row annotations
#
# Packages required: 
# - openxlsx
#
read_gct_13 <- function(excel_file, sheetnr){
  raw <- read.xlsx(excel_file, sheet = sheetnr, colNames = FALSE)
  
  size_info <- na.omit(as.numeric(raw[2,]))
  col_annot_size <- size_info[4]
  row_annot_size <- size_info[3]
  
  column_names <- as.character(raw[3,])[-c(1:((row_annot_size)+1))]
  
  expr_mat <- raw[-c(1: (3+col_annot_size) ),-c(1:((row_annot_size)+1))]
  colnames(expr_mat) <- column_names
  
  if (any(duplicated(raw[-c(1: (3+col_annot_size) ),1]))){
    row.names(expr_mat) <- paste(raw[-c(1: (3+col_annot_size) ),1], seq(1,nrow(raw[-c(1: (3+col_annot_size) ),]),1),
                             sep="_")
  } else {
    row.names(expr_mat) <- raw[-c(1: (3+col_annot_size) ),1]
  }
  
  if (sapply(expr_mat, class)[1] !="numeric"){
    expr_mat <- as.matrix(expr_mat)
    expr_mat <- as.data.frame(apply(expr_mat,c(1,2), as.numeric))
  }
  
  sample_annot <- as.data.frame(t(raw[c(3: (3+col_annot_size) ),-c(1:((row_annot_size)+1))]))
  colnames(sample_annot) <- c("Sample",as.character(raw[,1])[c(4:(3+col_annot_size))])
  row.names(sample_annot) <- sample_annot$Sample
  
  row_annot <- as.data.frame(raw[-c(1: (3+col_annot_size) ),c(1:((row_annot_size)+1))])
  colnames(row_annot) <- c("Row_name",as.character(raw[3,])[c(2:((row_annot_size)+1))])
  if (any(duplicated(row_annot$Row_name))){
    row_annot$Row_name <- paste(row_annot$Row_name, seq(1,nrow(row_annot),1),
                                sep="_")
  }
  row.names(row_annot) <- row_annot$Row_name
  
  return(list(Sample_annot = sample_annot,
              Row_annot = row_annot,
              Expr = expr_mat))
}


##############################
# Filter for valid values based on percentage
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
# - percentage: the minimum percentage of valid values required for 
#   each protein
#
# This function returns a protein expression table that is filtered
# for valid values.
#
# Packages required: none
# 
filter_missingvalues <- function (df, percentage) {
  df <- as.data.frame(df)
  perc <- ceiling(ncol(df)*(percentage/100))
  vv <- rowSums(!is.na(df))
  to_be_kept <- unlist(lapply(vv, function(x){ 
    if (x < perc) { y <- F}
    else y <- T}))
  return(as.data.frame(df[to_be_kept,]))
}


##############################
# Calculate singScores with the singscore algorithm
# 
# Arguments for this function:
# - ranked_df: a DGEList object that was passed through the rankGenes function  
# - geneset_list: parsed XML files for gene sets with getBroadSets
# - min_overlap: what is the min. number of genes that should be common between
#   the gene set and our expression table
#
# This function returns a list contaning 3 separate tables:
# (1) singScores, (2) dispersions
#
# Packages required: 
# - singscore
#
calculate_singScores <- function(ranked_df, geneset_list, min_overlap){
  sample_order <- colnames(ranked_df)
  
  GS <- 1
  scoredfs <- lapply(seq(1, length(geneset_list)), function(GS){
    overlapping_genes <- geneset_list[[GS]]@geneIds
    overlapping_genes <- overlapping_genes[overlapping_genes %in% row.names(ranked_df)]
    # only keep genes that are present in our data
    
    if (length(overlapping_genes) < min_overlap){
      ss <- data.frame(Col1 = rep(NA, length(sample_order)), 
                       Col2 = rep(NA, length(sample_order)), 
                       row.names = sample_order)
    } else {
      ss <- simpleScore(ranked_df, upSet = overlapping_genes, knownDirection = TRUE, centerScore =TRUE)
      ss <- ss[sample_order,]
    }
    colnames(ss) <- c(names(geneset_list)[GS], paste(names(geneset_list)[GS],"TotalDisp", sep="_"))
    return(ss)
  })
  
  scoredfs_merged <- do.call(cbind,scoredfs)
  
  Disp_tab <- scoredfs_merged[,grepl("_TotalDisp", colnames(scoredfs_merged))]
  Scores <- scoredfs_merged[,!grepl("_TotalDisp", colnames(scoredfs_merged))]
  return(list(Scores = Scores,
              Dispersions = Disp_tab))
}


##############################
# Create a long (melted) table
# 
# Arguments for this function:
# - table_to_melt: the table to be transformed  
# - row_annot: row annotation that we want to merge it with
# - col_annot: column annotation that we want to merge it with
#
# This function returns a a long (melted) table
#
# Packages required:
# - reshape2
# - varhandle
#
create_long_table <- function(table_to_melt, row_annot, col_annot){
  newtab <- melt(as.matrix(table_to_melt))
  colnames(newtab) <- c("Row_name", "Sample", "Expr")
  newtab$Sample <- unfactor(newtab$Sample)
  newtab$Row_name <- unfactor(newtab$Row_name)
  newtab <- merge(newtab, col_annot, by="Sample")
  newtab <- merge(newtab, row_annot, by="Row_name")
  return(newtab)
}



##############################
# Prepare table in long format for ggplot
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
#
# This function returns the protein expression table in a long
# format with columns Sample, Protein and Intensity.
#
# Packages required:
# - reshape2
# - varhandle
#
prepare_forggplot <- function (df) {
  df <- t(as.matrix(df))
  df <- cbind(row.names(df),df)
  df <- reshape2::melt(df, id.vars="V1")
  df <- varhandle::unfactor(df)
  df$value <- suppressWarnings(as.numeric(as.character(df$value)))
  df <- df[df[,2] !="",]
  colnames(df) <- c("Sample", "Protein","Intensity")
  return(df)
}


##############################
# Visualize protein intensity distributions across samples
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
# - name: the title of the plot
#
# This function returns a series of boxplots to visualize the
# protein intensity distributions in each sample.
#
# Packages required:
# - ggplot2
# 
create_sample_boxplots <- function(df, name){
  sampleorder <- colnames(df)
  df <- prepare_forggplot(df)
  df.expression.gg_noNA <- na.omit(df)
  df.expression.gg_noNA$Sample <- factor(df.expression.gg_noNA$Sample, levels = sampleorder)
  median.value <- median(df.expression.gg_noNA$Intensity)
  g_cell_raw <- ggplot(df.expression.gg_noNA, aes(x=Sample, y=Intensity))
  g_cell_raw + geom_violin() + geom_hline(yintercept = median.value)+ 
    theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(text= element_text(size=12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.0),
          legend.position = "none")+
    geom_boxplot(width=0.1) +
    stat_summary(fun.y=median, geom="point", size=2, color="red") +
    labs(title=name, 
         y="Expression", x="Samples")
  
}



##############################
# Perform median normalization (centering around the global median)
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows

# This function returns the a new protein expression table where 
# samples are median-normalized (centered around the global median)
#
# Packages required: none
#
normalize_median <- function (df) {
  df2 <- df
  med <- median(as.matrix(df), na.rm=T)
  for (i in 1:ncol(df)){
    df2[,i] <- (df[,i] - median(df[,i], na.rm = T) + med)
  }
  return(df2)
}



##############################
# Perform ORA
# 
# Arguments for this function:
# - pathways: the list of gene sets to check
# - genes: the vector of genes to be tested
# - universe: a universe from which genes were selected
# - minSize: min. size of a gene set to test
# - maxSize: max. size of a gene set to test
#
# This function returns a summary list with ORA results.
#
# Packages required:
# - fgsea
#
perform_fora <- function(pathways,genes, universe, minSize = 1, maxSize = Inf){
  
  y <- fgsea::fora(pathways = pathways, genes = genes, 
            universe=universe,
            minSize = minSize, 
            maxSize = maxSize)
  
  for (i in 1:nrow(y)){
    y$genes[i] <- paste(y$overlapGenes[[i]], collapse = "/")
  }
  y$overlapGenes <- NULL
  
  return(as.data.frame(y))
}


##############################
# Calculate mean expression within patients
# 
# Arguments for this function:
# - df: the gene/protein expression table
# - patient_vector: the vector of patient IDs
#
# This function returns an expression table
# where ROI-level expressions are averaged over patients.
#
# Packages required: none
#
mean_over_patient <- function(df, patient_vector){
  newdf <- matrix(nrow=nrow(df), ncol=length(patient_vector))
  colnames(newdf) <- patient_vector
  row.names(newdf) <- row.names(df)
  
  for (i in 1:ncol(newdf)){
    for (j in 1:nrow(newdf)){
      newdf[j,i] <- mean(as.numeric(df[j,grep(colnames(newdf)[i], colnames(df))]))
    }
  }
  return(newdf)
}
