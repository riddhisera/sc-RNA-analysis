#Package Management and Setup
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("slingshot", force = TRUE)
BiocManager::install("SingleCellExperiment", force = TRUE)
BiocManager::install("Seurat", force = TRUE)

library(slingshot)
library(SingleCellExperiment)
library(Seurat)
library(scater)
library(scran)
library(igraph)
library(leiden)

#Data Preparation
slingdata <- read.csv("Downloads/data.csv")
slingdata <- t(as.matrix(slingdata))

#Creating Seurat Object
slingseu <- CreateSeuratObject(counts = slingdata)
slingseu

#Data Normalization and Identification of Variable Features
all.genes <- rownames(slingseu)
slingseu <- NormalizeData(slingseu, normalization.method = "LogNormalize")
slingseu <- FindVariableFeatures(slingseu, selection.method = "vst")
top10 <- head(VariableFeatures(slingseu), 10)
slingseu <- ScaleData(slingseu, features = all.genes)

#Principal Component Analysis (PCA)
slingseu <- RunPCA(slingseu, features = VariableFeatures(object = slingseu))

#Clustering Analysis
print(slingseu[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(slingseu, dims = 1:2, reduction = "pca")

DimPlot(slingseu, reduction = "pca") + NoLegend()
DimHeatmap(slingseu, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(slingseu)

slingseu <- FindNeighbors(slingseu, dims = 1:10)
slingseu <- FindClusters(slingseu, resolution = 0.5)

head(Idents(slingseu), 5)

#UMAP
slingseu <- RunUMAP(slingseu, dims = 1:10)
DimPlot(slingseu, reduction = "umap", label = TRUE)

#Saving file for future
saveRDS(slingseu, file = "slingseu.rds")

#Marker Gene Analysis
cluster2.markers <- FindMarkers(slingseu, ident.1 = 2)
head(cluster2.markers, n = 5)

VlnPlot(slingseu, features = c("HEXB", "MAFA"))

FeaturePlot(slingseu, features = c("GRIN3A", "CHGA", "APOC1", "HOPX", "FEV", "RPS17"))
FeaturePlot(slingseu, features = c("GCK", "LDHA", "HK1", "INS", "GAD2", "FXYD2"))

slingseu.markers <- FindAllMarkers(slingseu, only.pos = TRUE)

#Heatmap Visualization
slingseu.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3
DoHeatmap(slingseu, features = top3$gene) + NoLegend()

## yet to do - Assigning cell type identity to clusters

#Trajectory Analysis with Slingshot
sce <- as.SingleCellExperiment(slingseu)
reducedDims(sce)$UMAP <- Embeddings(slingseu, "umap")
colData(sce)$cluster <- Idents(slingseu)

sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP')

lnes <- getLineages(reducedDim(sce,"UMAP"), sce$ident)
lnes@lineages

sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = "UMAP",
                 allow.breaks = FALSE, start.clus="5")

lnes <- getLineages(reducedDim(sce,"UMAP"), sce$ident, start.clus = "5")
lnes@metadata

######################################
library(Polychrome)
library(ggbeeswarm)
library(ggthemes)

# this define the cluster color. You can change it with different color scheme.
my_color <- createPalette(length(levels(sce$ident)), c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(sce$ident))

plot(reducedDims(sce)$UMAP, col = my_color[as.character(sce$ident)], pch=16, asp = 1)
legend("bottomleft",legend = names(my_color[levels(sce$ident)]),  
       fill = my_color[levels(sce$ident)])
lines(SlingshotDataSet(lnes), lwd=2, type = 'lineages', col = c("black"))

#####################################
col_data <- colData(sce)
num_rows <- nrow(col_data)
num_cols <- length(col_data)

# Initialize the data frame with the correct dimensions
slingshot_df <- data.frame(matrix(nrow = num_rows, ncol = num_cols))

# Populate the data frame
for (i in seq_along(col_data)) {
  slingshot_df[[i]] <- col_data[[i]]
}

# Set the column names
colnames(slingshot_df) <- names(col_data)

# re-order y-axis for better figure: This should be tailored with your own cluster names
slingshot_df$ident = factor(slingshot_df$ident, levels=c(5,7,9,8,2,10,14,13,0))

ggplot(slingshot_df, aes(x = slingPseudotime_1, y = ident, colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime") + 
  scale_colour_manual(values = my_color)

# Assuming slingshot_df is already filtered to remove NAs
ggplot(slingshot_df, aes(x = slingPseudotime_1_2, y = as.factor(ident), colour = as.factor(ident))) +
  geom_quasirandom(groupOnX = TRUE, alpha = 0.7, size = 1.5) + # adjust alpha for transparency and size for point size
  theme_classic() +
  theme(
    legend.position = "right", # adjust legend position
    axis.text.x = element_text(angle = 45, hjust = 1), # angle x axis text for better readability
    axis.title = element_text(size = 12), # adjust title size
    legend.title = element_text(size = 10), # adjust legend title size
    legend.text = element_text(size = 8) # adjust legend text size
  ) +
  scale_color_brewer(palette = "Set3") + # change to a color brewer palette
  xlab("First Slingshot pseudotime") +
  ylab("Cell Type") +
  ggtitle("Cells ordered by Slingshot Pseudotime") +
  labs(colour = "Identity") # rename the legend title

######################################




