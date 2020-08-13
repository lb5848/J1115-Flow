rm(list = ls())
# Load packages
library(ggrastr)
library(rstudioapi)
library(devtools)
library("flowCore")
library("flowWorkspace")
library(cytofCore)
library(FlowSOM)
# library(cytofkit)
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(flowViz)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(uwot)
library(CATALYST)
library(Rphenograph)
library(diffcyt)
library(SummarizedExperiment)
library(stringr)
library(ggcyto)
library(SingleCellExperiment)
library(scran)
library(scater)
library(flowVS)
library(readxl)
library(flowStats)
library(FlowSOMworkshop)
library(flowAI)
library(CytoNorm)
library(PeacoQC)
library(CytoML)
options(java.parameters="-Xmx60G")
library(tidyverse)
library(data.table)
library(scuttle)
library(iMUBAC)



# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory
# Define workingDirectory
wdName <- "Working_DirectoryFCS"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")

setwd(workingDirectory)

# Load workspace and SCEobject

sce <- readRDS("SCE_part2_DR.rds")
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "median")
plotAbundances(sce, k = "meta8", by = "cluster_id", group_by = "condition")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta8", by = "cluster_id", scale = "last", bars = TRUE, perc = TRUE)



FDR_cutoff <- 0.05
ei <- sce@metadata$experiment_info

# DA using GLMM

# (da_formula1 <- createFormula(ei, 
#                               cols_fixed = "condition",
#                               cols_random = c("patient_id")))
# contrast2 <- createContrast(c(0,1))
# da_res2 <- diffcyt(sce,
#                    formula = da_formula1, contrast = contrast2,
#                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
#                    clustering_to_use = "meta10", verbose = TRUE, subsampling = TRUE,
#                    transform = FALSE, normalize = FALSE) 
# da2 <- rowData(da_res2$res)
# 
# plotDiffHeatmap(sce, da2, top_n = 12, all = TRUE, fdr = FDR_cutoff)

# DA using edgeR
design <- createDesignMatrix(ei,
                             cols_design = c("condition"))
contrast <- createContrast(c(0, 1))

nrow(contrast) == ncol(design)
out_DA <- diffcyt(sce,
                  experiment_info = ei, design = design, contrast = contrast,
                  analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "meta10", verbose = TRUE, subsampling = TRUE,
                  transform = FALSE, normalize = FALSE)

da <- rowData(out_DA$res)
plotDiffHeatmap(sce, da, top_n = 10, all = TRUE, fdr = FDR_cutoff)


# Save current workspace
save(list = ls(), file = "workspaceSCEDA.rds")
