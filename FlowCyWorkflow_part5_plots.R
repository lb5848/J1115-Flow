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
library(scales)
library(gplots)

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

stat.test <- as_tibble(da)
stat.test$cluster_id <- paste0("C", stat.test$cluster_id)
p.adj.signif <- c("**", "*", rep("ns", 6))
# y.position <- c(90, 15, 90)
group1 <- (rep("relapse",8))
group2 <- (rep("responder", 8))
y.position <- c(9, 9, 40, 75, 40, 40, 20, 1.1)
stat.test <- cbind(stat.test, group1, group2, p.adj.signif, y.position)
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp <- bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 2.5)
bxp

df <- as.data.frame(colData(sce))
summary(df)

cond_cluster_df <- df %>% select(condition, cluster_annotation)
freq_table <- table(cond_cluster_df$condition, cond_cluster_df$cluster_annotation)
tot_rel <- sum(freq_table["relapse", ])
tot_res <- sum(freq_table["responder", ])

freq_table["relapse",] <- freq_table["relapse", ]/tot_rel*100
freq_table["responder",] <- freq_table["responder", ]/tot_res*100
freq_table <- t(freq_table)

plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", by = "cluster_id",  fun = "mean",
                scale = "first", bars = TRUE, perc = TRUE)
freqdf <- cbind(rownames(freq_table), freq_table[, 1], freq_table[, 2])
colnames(freqdf) <- c("Cluster", "relapse", "responder")
write.csv(freqdf, file = "freq_table.csv", row.names = FALSE)
balloon <- read.table(file.choose(), header = TRUE, sep = ",", stringsAsFactors = FALSE)


balloon$Cluster <- factor(balloon$Cluster, levels = rev(unique(balloon$Cluster)))
balloon_melted <- melt(balloon, sort = FALSE)


p <- ggplot(balloon_melted, aes(x = variable, y = Cluster))

pp <- p + geom_point(aes(size = value), colour = "grey34") + theme(panel.background = element_blank()) +
  scale_size_area(max_size = 12)

pp

matrix_r <- read.table(file.choose(), sep = ",", stringsAsFactors = FALSE, header = TRUE)
matrix_r <- matrix_r[-c(9:11), ]
colnames(matrix_r)
matrix_r <- matrix_r[, c(2:15 )]
colnames(matrix_r)
matrix_r[, grep("Freq", colnames(matrix_r))]
freqs <- matrix_r %>% select(grep("Freq", colnames(matrix_r)))
MFI <- matrix_r %>% select(grep("Mean", colnames(matrix_r)))
iMFI <- freqs*MFI
colnames(iMFI) <- c("CD27", "CD57", "CD69", "DNAM1", "PD1", "TIGIT", "TIM3")
scaled_iMFI <- iMFI
for(i in c(1:ncol(iMFI))){
  scaled_iMFI[, i] <- rescale(iMFI[, i], to = c(0, 100))
}

heat.colors<-colorRampPalette(c("navy","blue4","blue","skyblue","khaki1","lightgoldenrod1","goldenrod1","orange"))(100)
##### Make a heatmap  ##change parameters if needed (type help(heatmap.2) to get a full description of the options)
matrix_r <- scaled_iMFI
scale.data <- as.matrix((t(matrix_r)-apply(t(matrix_r),1,mean))/apply(t(matrix_r),1,sd))
colnames(scale.data) <- rownames(freq_table)
##### Plot heatmap 
heatmap.2(as.matrix(t(scale.data)),
          dendrogram="both", scale="none",  na.color="grey",
          col = heat.colors, trace = "none", labRow = rownames(t(scale.data)), key = TRUE, keysize = 1, cexCol=1,
          density.info = "none", symkey = FALSE, margins = c(7, 10), main = "iMFI", 
          xlab = "Markers", ylab = "CLUSTERS")

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
svg(filename = "UMAP_clusters_condition.svg", bg = "white")
CATALYST::plotDR(sce, dr = "UMAP", color_by = "cluster_annotation", facet_by = "condition") + scale_color_brewer(palette = "Dark2") +
  geom_density2d(binwidth = 0.006, colour = "black")
dev.off()

# UMAP color_by = Clusters facet_by = sample_id
svg(filename = "UMAP_clusters_sample.svg", bg = "white")
plot <- CATALYST::plotDR(sce, dr = "UMAP", color_by = "cluster_annotation", facet_by = "sample_id") + scale_color_brewer(palette = "Dark2") +
  geom_density2d(binwidth = 0.006, colour = "black")
plot$facet$params$ncol <- 3
plot
dev.off()

# UMAP color_by = CD27 and DNAM1 facet_by = condition
svg(filename = "UMAP_CD27_DNAM1.svg", bg = "white")
CATALYST::plotDR(sce, dr = "UMAP", color_by = c("CD27", "DNAM1"), facet_by = "condition")
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