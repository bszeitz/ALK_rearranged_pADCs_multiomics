##############################
#
# Load packages
# Author: Bea Szeitz
#
##############################

packages = c("ggplot2", "readxl", "openxlsx", "ggpubr", "ggbiplot", 
             "ComplexHeatmap", "circlize", 
             "clusterProfiler", 
             "ggrepel", "edgeR",
             "PerformanceAnalytics", "varhandle", "reshape2", "ggvenn", 
             "gghighlight","RColorBrewer",
             "stringr", "singscore","imputeLCMD",
             "preprocessCore", "ggbeeswarm","GSEABase","qusage", 
             "cowplot","ggfortify","fgsea", "ConsensusClusterPlus", "gridExtra",
             "survival","survminer","glmmSeq", "ggpmisc", "SpatialDecon")

package.load <- lapply(packages, function(x) {
    library(x, character.only = T)
})

