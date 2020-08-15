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
# library(scuttle)
# library(iMUBAC)
library(ggpubr)


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

outputDirectory <- getwd()
outputDirectory <- paste(outputDirectory, "output", sep = "/")
dir.create(outputDirectory)
setwd(outputDirectory)
load("workspaceFinal.rds")
sce
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "median")

getwd()
plotFolder <- paste(outputDirectory, "plots", sep = "/")
dir.create(plotFolder)
setwd(plotFolder)

# plotAbundances w/ stats

# library(ggpubr)
stat.test <- as_tibble(da)
p.adj.signif <- c("**", "*", rep("ns", 6))
# y.position <- c(90, 15, 90)
group1 <- (rep("relapse",8))
group2 <- (rep("responder", 8))
y.position <- c(9, 9, 40, 75, 40, 40, 20, 1.1)
stat.test <- cbind(stat.test, group1, group2, p.adj.signif, y.position)
bxp <- plotAbundances(sce, k = "meta8", by = "cluster_id", group_by = "condition")
bxp <- bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 2.5)




display.brewer.all(colorblindFriendly = TRUE)
png(filename = "colorblindFriendly.png", bg = "white")
display.brewer.all(colorblindFriendly = TRUE)
dev.off()
png(filename = "brewerAll.png", bg = "white")
display.brewer.all(colorblindFriendly = FALSE)
dev.off()
# MDS plot
tiff(filename = "MDSplot.tiff", compression = "lzw", bg = "white")
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition")
dev.off()
svg(filename = "MDSplot.svg", bg = "white")
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition")
dev.off()
# boxplot abundances per cluster
tiff(filename = "boxplot_clusters.tiff", compression = "lzw", bg = "white")
plotAbundances(sce, k = "meta8", by = "cluster_id", group_by = "condition")
dev.off()
svg(filename = "boxplot_clusters.svg", bg = "white")
plotAbundances(sce, k = "meta8", by = "cluster_id", group_by = "condition")
dev.off()

# Expression Heatmap - clusters
tiff(filename = "ExpressionHeatmap.tiff", compression = "lzw", bg = "white")
plotExprHeatmap(sce, features = type_markers(sce), k = "meta8", by = "cluster_id", 
                scale = "last", bars = TRUE, perc = TRUE, hm_pal = rev(brewer.pal(11, "RdYlBu")))
dev.off()
svg(filename = "ExpressionHeatmap.svg", bg = "white")
plotExprHeatmap(sce, features = type_markers(sce), k = "meta8", by = "cluster_id")
dev.off()

# Expression - single cluster per marker
tiff(filename = "Cluster_exprs.tiff", compression = "lzw", bg = "white")
plotClusterExprs(sce, k = "meta8", features = "type")
dev.off()
svg(filename = "Cluster_exprs.svg", bg = "white")
plotClusterExprs(sce, k = "meta8", features = "type")
dev.off()

# UMAP color_by = Clusters facet_by = condition
tiff(filename = "UMAP_clusters_condition.tiff", compression = "lzw", bg = "white")
CATALYST::plotDR(sce, dr = "UMAP", color_by = "Clusters", facet_by = "condition") + scale_color_brewer(palette = "Dark2") + 
  geom_density2d(binwidth = 0.006, colour = "black")
dev.off()
svg(filename = "UMAP_clusters_condition.svg", bg = "white")
CATALYST::plotDR(sce, dr = "UMAP", color_by = "Clusters", facet_by = "condition") + scale_color_brewer(palette = "Dark2") +
  geom_density2d(binwidth = 0.006, colour = "black")
dev.off()

# UMAP color_by = Clusters facet_by = sample_id
tiff(filename = "UMAP_clusters_sample.tiff", compression = "lzw", bg = "white")
plot <- CATALYST::plotDR(sce, dr = "UMAP", color_by = "Clusters", facet_by = "sample_id") + scale_color_brewer(palette = "Dark2") +
  geom_density2d(binwidth = 0.006, colour = "black")
plot$facet$params$ncol <- 3
plot
dev.off()
svg(filename = "UMAP_clusters_sample.svg", bg = "white")
plot <- CATALYST::plotDR(sce, dr = "UMAP", color_by = "Clusters", facet_by = "sample_id") + scale_color_brewer(palette = "Dark2") +
  geom_density2d(binwidth = 0.006, colour = "black")
plot$facet$params$ncol <- 3
plot
dev.off()

# UMAP color_by = CD27 and DNAM1 facet_by = condition
tiff(filename = "UMAP_CD27_DNAM1.tiff", compression = "lzw", bg = "white")
CATALYST::plotDR(sce, dr = "UMAP", color_by = c("CD27", "DNAM1", "CD69"), facet_by = "condition")
dev.off()
svg(filename = "UMAP_CD27_DNAM1.svg", bg = "white")
plotDR(sce, dr = "UMAP", color_by = c("CD27", "DNAM1"), facet_by = "condition")
dev.off()

# UMAP color_by = Clusters
tiff(filename = "UMAP_clusters_allPts.tiff", compression = "lzw", bg = "white")
plotDR(sce, dr = "UMAP", color_by = "Clusters") + scale_color_brewer(palette = "Dark2")
dev.off()
svg(filename = "UMAP_clusters_allPts.svg", bg = "white")
plotDR(sce, dr = "UMAP", color_by = "Clusters") + scale_color_brewer(palette = "Dark2")
dev.off()

# plotCounts
tiff(filename = "plotCounts.tiff", bg = "white", compression = "lzw")
plotCounts(sce, group_by = "sample_id", color_by = "condition")
dev.off()
svg(filename = "plotCounts.svg", bg = "white")
plotCounts(sce, group_by = "sample_id", color_by = "condition")
dev.off()

cellCounts <- n_cells(sce)
cellCounts <- as.data.frame(cellCounts)
colnames(cellCounts) <- c("sample_id", "cell #")

# plotExprs
tiff(filename = "plotExprs.tiff", bg = "white", compression = "lzw")
p <- plotExprs(sce, features = NULL, color_by = "condition")
p$facet$params$ncol <- 4
p
dev.off()
svg(filename = "plotExprs.svg", bg = "white")
p <- plotExprs(sce, features = NULL, color_by = "condition")
p$facet$params$ncol <- 4
p
dev.off()
svg(filename = "plotNRS.svg", bg = "white")
plotNRS(sce, features = type_markers(sce), color_by = "condition")
dev.off()

svg(filename = "plotExprHeatmap.svg", bg = "white")
plotExprHeatmap(sce, features = type_markers(sce), k= "meta8", bin_anno = TRUE, row_anno = TRUE)
dev.off()

svg(filename = "delta_area.svg", bg = "white")
delta_area(sce)
dev.off()

svg(filename = "tSNE_clusters_condition.svg", bg = "white")
plotDR(sce, dr = "TSNE", color_by = "Clusters", facet_by = "condition") + scale_color_brewer(palette = "Dark2")
dev.off()
svg(filename = "DiffMap_clusters_condition.svg", bg = "white")
plotDR(sce, dr = "DiffusionMap", color_by = "Clusters", facet_by = "condition") + scale_color_brewer(palette = "Dark2")
dev.off()

svg(filename = "DiffHeatmap_meta8.svg", bg = "white")
plotDiffHeatmap(sce, da2, top_n = 12, all = TRUE, fdr = FDR_cutoff)
dev.off()