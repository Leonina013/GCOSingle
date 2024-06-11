setwd("C:/SCPS/Kidney/GSE111360/BII/BII/integrated")

#____________________-CLEANER_______________(NEW TAB)
# Read the CSV file into a data frame
your_data <- read.csv("T6_CD45_NEG.csv")

# Identify and remove duplicate rows based on all columns
your_data_no_duplicates <- your_data[!duplicated(your_data[, 1]), ]


# If you want to identify and remove duplicates based on specific columns, you can use the subset function
# For example, if you want to consider only the "column_name" for identifying duplicates
# your_data_no_duplicates <- your_data[!duplicated(your_data$column_name), ]

# Write the cleaned data frame to a new CSV file if needed
write.csv(your_data_no_duplicates, "cleaned_T6_CD45_NEG.csv", row.names = FALSE)

#_____________________________________DONE________________________________



install.packages("SeuratObject")
install.packages("devtools")
install.packages("hdf5r")
install.packages("Matrix")
install.packages("htmltools")
# Install devtools package if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install BiocManager package if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install SingleR from Bioconductor
BiocManager::install("SingleR")
BiocManager::install("celldex")
remotes::install_github('satijalab/seurat-wrappers')

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

install.packages("devtools")
devtools::install_version("dbplyr", version = "2.3.4")

install.packages("singleCellHaystack")

library(BiocParallel)
#remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE) #to install seurat
register(MulticoreParam(workers = 8, progressbar = TRUE))
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(hdf5r)
library(SingleR)
library(celldex)
library(pheatmap)
library(gridExtra)
library(harmony)
library(monocle3)
library(SeuratWrappers)
library(singleCellHaystack)

# 1. read raw counts/expression matrix ------
#merged_seurat_filtered.sub <- Read10X_h5("PAAD_CRA001160_expression.h5")
merged_seurat_filtered.sub <- read.csv("OMT3_CD45_NEG.csv", row.names = 1, check.names = FALSE)
str(merged_seurat_filtered.sub)

#___________________________________________________________________________
#DO YOU WANT TO INTEGRATE? THEN KEEP FOLLOWING,ELSE SCROLL DOWN TO THE merged_seurat_filtered 
#CCA Integration (NEW TAB)
options(Seurat.object.assay.version = "v5")
# get data location
dirs <- list.files(path = "C:/SCPS/Kidney/GSE111360/BII/BII/integrated")

for(x in dirs){
  name <- gsub('','', x)
  
  cts <- read.csv(name, row.names = 1, check.names = FALSE)
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}


# merge datasets

merged_seurat <- merge(N1_NORMAL_OVARY.csv, y = c(OMT1_CD45_NEG.csv, OMT1_CD45_POS.csv, OMT3_CD45_NEG.csv, OMT3_CD45_POS.csv, T1_TUMOR_CAN.csv, T6_CD45_NEG.csv, T6_TUMOR_CAN.csv),
                       add.cell.ids = dirs,
                       project = 'OVA')

#If there are breaks in the index then run this
# Assuming ls() returns a character vector with the specified names

merged_seurat

# QC & filtering -----------------------

View(merged_seurat@meta.data)
# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Location', 'Type', 'CD45State/CAN', 'Barcode'), 
                                    sep = '_')

#______________________________________________________________________
#HarmonyIntegration (NEW TAB)

merged_seurat$mito.percent <- PercentageFeatureSet(merged_seurat, pattern = '^MT-')
View(merged_seurat@meta.data)
# explore QC

# filter
merged_seurat
merged_seurat.filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 & 
                          mito.percent < 5)

# standard workflow steps
merged_seurat.filtered <- NormalizeData(merged_seurat.filtered)
merged_seurat.filtered <- FindVariableFeatures(merged_seurat.filtered)
merged_seurat.filtered <- ScaleData(merged_seurat.filtered)
merged_seurat.filtered <- RunPCA(merged_seurat.filtered)
ElbowPlot(merged_seurat.filtered)
merged_seurat.filtered <- RunUMAP(merged_seurat.filtered, dims = 1:20, reduction = 'pca')

before <- DimPlot(merged_seurat.filtered, reduction = 'umap', group.by = 'CD45State/CAN')


# run Harmony -----------
merged_seurat.harmony <- merged_seurat.filtered %>%
  RunHarmony(group.by.vars = 'CD45State/CAN', plot_convergence = FALSE)

merged_seurat.harmony@reductions

merged_seurat.harmony.embed <- Embeddings(merged_seurat.harmony, "harmony")
merged_seurat.harmony.embed[1:10,1:10]



# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
merged_seurat.harmony <- merged_seurat.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize 
after <- DimPlot(merged_seurat.harmony, reduction = 'umap', group.by = 'CD45State/CAN')
before|after
#______________________________________________________________________

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

# explore QC
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)


# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 500 &
                                   nFeature_RNA > 200 &
                                   mitoPercent < 5)

view(merged_seurat_filtered@meta.data)

merged_seurat




# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)

top10 <- head(VariableFeatures(merged_seurat_filtered), 10)
top10

# ___plot variable features with and without labels ---------
plot3 <- VariableFeaturePlot(merged_seurat_filtered)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Location')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))
P6 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'CD45State/CAN')
grid.arrange(p1, p2, P6, ncol = 2, nrow = 2)


# perform integration to correct for batch effects ------
#for V5
merged_seurat_filtered <- IntegrateLayers(
  object = merged_seurat_filtered, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

#FOR V3
#obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Location')
#for(i in 1:length(obj.list)){
#  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
##  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
#}


# select integration features
#features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
#anchors <- FindIntegrationAnchors(object.list = obj.list,
                                 # anchor.features = features)

# Check the number of dimensions available in the integrated.cca reduction

# integrate data
#seurat.integrated <- IntegrateData(anchorset = anchors)
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, reduction = "integrated.cca", dims = 1:20)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 2, cluster.name = "cca_clusters")
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")
p5 <- DimPlot(
  merged_seurat_filtered,
  reduction = "umap.cca",
  group.by = c("Location", "Type", "CD45State/CAN"),
  combine = FALSE, label.size = 2
)



# Scale data, run PCA and UMAP and visualize integrated data
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


p3 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Location')
p4 <- DimPlot(merged_seurat_filtered, reduction = 'umap.cca', group.by = 'Type',
              cols = c('red','green','blue'))


grid.arrange(p1, p2, p3, p4, ncol = 3, nrow = 3)
p1+p2+p3+p4+p5
p3+p4
p4+p5
p8<-DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'seurat_clusters')
p7<-DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'seurat_clusters_manual')
p3|p7|p8
#______________________________________________________
#Manual Annotation (IGNORE THIS)

# Create a data frame with the cluster mapping
cluster_mapping <- data.frame(
  Cluster = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29),
  Desc = c(
    "Epithelial_Cells", "T_Cells", "T_Cells", "T_Cells", "T_Cells", "NK Cells", 
    "T_Cells", "T_Cells", "T_Cells", "T_Cells", "NK Cells", 
    "NK Cells", "T_Cells", "Epithelial_Cells", "T_Cells", 
    "Monocyte", "T_Cells", "T_Cells", "B_Cells", "Monocyte", 
    "Monocyte", "Epithelial_Cells", "T_Cells", 
    "Epithelial_Cells", "Monocyte", "Fibroblast", 
    "B_Cells", "Epithelial_Cells", "Epithelial_Cells", 
    "Monocyte"
  )
)

# Merge the Seurat object with the cluster mapping
merged_seurat_filtered$seurat_clusters_manual <- factor(
  merged_seurat_filtered$seurat_clusters,
  levels = cluster_mapping$Cluster,
  labels = cluster_mapping$Desc
)

# Plot the updated clusters
DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'seurat_clusters_manual')
#________________________________________________

VlnPlot(merged_seurat_filtered, features = c("FAP"), group.by = 'seurat_clusters_manual' )
view(merged_seurat_filtered@meta.data)
#______________________________________________________

merged_seurat_filtered <- JackStraw(merged_seurat_filtered, num.replicate = 100)
merged_seurat_filtered <- ScoreJackStraw(merged_seurat_filtered, dims = 1:20)

JackStrawPlot(merged_seurat_filtered, dims = 1:20)


#_______________________________________________________________________________________
#before going for clusters 

merged_seurat_filtered <- JoinLayers(merged_seurat_filtered)
cluster1.markers <- FindMarkers(merged_seurat_filtered, ident.1 = 2, min.pct = 0.25)

#Find all markers  #RUN ONE OF THE TWO PLEASE
merged_seurat_filtered.markers <- FindAllMarkers(merged_seurat_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markers <- merged_seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
print(markers, n = Inf)
# Assuming 'markers.csv' is the desired output file name
write.csv(markers, file = 'newmarkers.csv', row.names = FALSE)

#Find all markers
merged_seurat_filtered.markers <- FindAllMarkers(merged_seurat_filtered,
               logfc.threshold = 0.25,
               min.pct = 0.25,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')

features <- merged_seurat_filtered@commands$RunPCA.RNA$features
VlnPlot(merged_seurat_filtered, features = c("TMIGD3"), group.by = 'seurat_clusters_manual' )

# ___VlnPlot() - you can plot raw counts as well ---------
VlnPlot(merged_seurat_filtered, features = c("MUC16", "CD14"), slot = "counts", log = TRUE)

plot <- FeaturePlot(merged_seurat_filtered, features = c("CD45RA"))
HoverLocator(plot = plot, information = FetchData(merged_seurat_filtered, vars = c("ident", "PC_1", "nFeature_RNA")))


FeaturePlot(merged_seurat_filtered, features = c("MUC16","CD27"), blend = TRUE)
top10 <- merged_seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(merged_seurat_filtered, features = top10$gene)

markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
DotPlot(merged_seurat_filtered, features = markers.to.plot, cols = c("blue", "red", "green", "yellow", "purple"), dot.scale = 8, split.by = "Location") +
  RotatedAxis()
DimPlot(merged_seurat_filtered, reduction = "umap", split.by = "seurat_clusters_manual")

FeaturePlot(merged_seurat_filtered, features = c("VEGF", "CD5L", "PPARG"), split.by = "seurat_clusters_manual", max.cutoff = 3, cols = c("grey",
                                                                                                    "red"), reduction = "umap")
merged_seurat_filtered

FeaturePlot(merged_seurat_filtered, features = c("GZMK", "CD19", "MUC1", "FAP", "CD14", "NCR1"), split.by = "seurat_clusters_manual", max.cutoff = 3, cols = c("grey",
                                                                                                                      "red"), reduction = "umap")
#___________________________________________________________________________
##THIS IS NON INTEGRATED (NEW TAB)
# ___create Seurat Object with count data --------
# include only genes that are are expressed in 3 or more cells and cells with complexity of 200 genes or more
merged_seurat_filtered <- CreateSeuratObject(counts = merged_seurat_filtered.sub, project = "BII", min.cells = 3, min.features = 200)
str(merged_seurat_filtered)

# count matrix
merged_seurat_filtered@assays$RNA@counts[1:10,1:10]

# 2. QC --------
merged_seurat_filtered[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_filtered, pattern = "^MT-")
str(merged_seurat_filtered)

# Show QC metrics for the first 5 cells
head(merged_seurat_filtered@meta.data, 5)

# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts

# ___Visualize QC metrics as a violin plot -------
VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# ___feature-feature or gene-gene relationship --------
plot1 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

merged_seurat_filtered <- subset(merged_seurat_filtered, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
str(merged_seurat_filtered)

unique(merged_seurat_filtered@meta.data$singleR.labels)
##LINE BELOW NOT REQUIRED UNLESS YOU ARE DOING TRAJECTORY
Idents(merged_seurat_filtered) <- merged_seurat_filtered$singleR.labels

# 3. Normalization ----------
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
str(merged_seurat_filtered)

# ___identification of highly variable features ---------
merged_seurat_filtered <- FindVariableFeatures(merged_seurat_filtered, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_seurat_filtered), 10)
top10

# ___plot variable features with and without labels ---------
plot3 <- VariableFeaturePlot(merged_seurat_filtered)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

# 4. scaling the data (performed prior to linear dim reduction) ---------
all.genes <- rownames(merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(merged_seurat_filtered, features = all.genes)

str(merged_seurat_filtered)


# 5. Linear Dimensionality Reduction ----------
merged_seurat_filtered <- RunPCA(merged_seurat_filtered, features = VariableFeatures(object = merged_seurat_filtered))

# ___Examine and visualize PCA results a few different ways -------
print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)

# ___plot-1 --------
VizDimLoadings(merged_seurat_filtered, dims = 1:2, reduction = "pca")

# ___plot-2 --------
DimPlot(merged_seurat_filtered, reduction = "pca")

# ___plot-3 heatmap -------
# allows for easy exploration of the primary sources of heterogeneity in a dataset
# and can be useful when trying to decide which PCs to include for further downstream analyses
DimHeatmap(merged_seurat_filtered, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(merged_seurat_filtered, dims = 1:5, cells = 500, balanced = TRUE)


# ___to dertermine "dimensionality" of the dataset -------
# essentially determine how many PCs to consider - we would ideally want to consider PCs that show maximum variations

# JackStraw Procedure!
# identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time

merged_seurat_filtered <- JackStraw(merged_seurat_filtered, num.replicate = 100)
merged_seurat_filtered <- ScoreJackStraw(merged_seurat_filtered, dims = 1:20)

JackStrawPlot(merged_seurat_filtered, dims = 1:20)
# The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). 
# ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line).


# An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function).
ElbowPlot(merged_seurat_filtered)
# from the plot, it looks like majority of true signal is captured in the first 15 PCs.
# PCs to consider = 15

# 6. Cluster cells --------
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)

# The FindClusters() function contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. 
# We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. 
# Optimal resolution often increases for larger datasets. 
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(merged_seurat_filtered), 5)



# 7. Run non-linear dimensional reduction (UMAP/tSNE) ---------
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(merged_seurat_filtered, reduction = "umap")
b1 <- DimPlot(merged_seurat_filtered, reduction = "umap.cca", label = TRUE)
b2 <- DimPlot(merged_seurat_filtered, reduction = "umap.cca", group.by = "Location", label = TRUE)
b3 <- DimPlot(merged_seurat_filtered, reduction = "umap.cca", group.by = "singleR.labels", label = TRUE)
b1|b2|b3
a1<-DimPlot(merged_seurat_filtered, reduction = "umap")

#_____________________________________________________________________
#singleCellHaystack WIP (IGNORE)

str(merged_seurat_filtered)
# Run haystack on the first 20 principal components
set.seed(123)
res.pc20 <- haystack(x = merged_seurat_filtered@reductions$pca@cell.embeddings[, 1:20], expression = merged_seurat_filtered@features)

# Get the top 1000 DEGs
sorted.table <- show_result_haystack(res.haystack = res.pc20, n = 1000)
gene.subset <- row.names(sorted.table)

# Cluster the genes by their expression pattern in the input space (first 20 PCs)
res.hc <- hclust_haystack(merged_seurat_filtered@reductions$pca@cell.embeddings[, 1:20], merged_seurat_filtered@data[gene.subset, ], grid.coordinates = res.pc20$info$grid.coordinates)

# Cut the hierarchical clustering tree into 5 clusters
res.hc.clusters <- cutree(res.hc, k = 5)

# Calculate the average expression of genes in each cluster
for (cluster in unique(res.hc.clusters)) {
  merged_seurat_filtered[[paste0("cluster_", cluster)]] <- colMeans(merged_seurat_filtered@data[names(which(res.hc.clusters == cluster)), ])
}

# Plot the averaged detection pattern for each cluster
plots <- lapply(paste0("cluster_", 1:5), function(cluster) {
  ggplot(merged_seurat_filtered, aes(reduction.umap_1, reduction.umap_2, color = .data[[cluster]])) +
    geom_point() +
    scale_color_distiller(palette = "Spectral")
})

# Display the plots
wrap_plots(plots)
#_____________________________________________________________________________________________
#AUTO ANNOTATION (NEW TAB)
# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()

View(as.data.frame(colData(ref)))
View(merged_seurat_filtered@meta.data)

# expression values are log counts (log normalized counts)


# run SingleR (default mode) ---------
# default for SingleR is to perform annotation of each individual cell in the test dataset

pbmc_counts <- GetAssayData(merged_seurat_filtered, layer = 'counts')

pred <- SingleR(test = pbmc_counts,
                ref = ref,
                labels = ref$label.main)

pred

merged_seurat_filtered$singleR.labels <- pred$labels[match(rownames(merged_seurat_filtered@meta.data), rownames(pred))]
a2 <- DimPlot(merged_seurat_filtered, reduction = 'umap.cca', group.by = 'singleR.labels')
a3 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'singleR.labels')

P6+a2
p1+a3
a2|a3
a3
# Annotation diagnostics ----------


# ...Based on the scores within cells -----------
pred
pred$scores

plotScoreHeatmap(pred)


# ...Based on deltas across cells ----------

plotDeltaDistribution(pred)




# ...Comparing to unsupervised clustering ------------

tab <- table(Assigned=pred$labels, Clusters=merged_seurat_filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(20))
#_____________________________________________________________________OR______
#____________________________________________________________________________--
#___________________________________________________________________________

# run SingleR with multiple reference datasets (default mode) ---------

# for pbmc data, we will use two datasets
hpca <- celldex::HumanPrimaryCellAtlasData()
dice <- celldex::DatabaseImmuneCellExpressionData()

# ...1. Strategy 1: Using reference-specific labels ----------
hpca$label.main
dice$label.main

# adding ref info to labels
hpca$label.main <- paste0('HPCA.', hpca$label.main)
dice$label.main <- paste0('DICE.', dice$label.main)

# create a combined ref based on shared genes
shared <- intersect(rownames(hpca), rownames(dice))
combined <- cbind(hpca[shared,], dice[shared,])
combined
combined$label.main

# run singleR using combined ref
# savings counts into a separate object
pbmc_counts <- GetAssayData(merged_seurat_filtered, layer = 'counts')

com.res1 <- SingleR(test = pbmc_counts, ref = combined, labels = combined$label.main)
table(com.res1$labels)

merged_seurat_filtered$com.res1.labels <- com.res1[match(rownames(merged_seurat_filtered@meta.data), rownames(com.res1)), 'labels']
View(merged_seurat_filtered@meta.data)

DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'com.res1.labels', label = TRUE)

# ...2. Strategy 2: Comparing scores across references ----------

hpca$label.main
dice$label.main
hpca$label.main <- gsub('HPCA\\.','', hpca$label.main)
dice$label.main <- gsub('DICE\\.','', dice$label.main)

com.res2 <- SingleR(test = pbmc_counts, 
                    ref = list(HPCA = hpca, DICE = dice),
                    labels = list(hpca$label.main, dice$label.main))

# Check the final label from the combined assignment.
table(com.res2$labels)

# which reference scored best for which label?
grouping <- paste0(com.res2$labels,'.', com.res2$reference)
best_ref <- as.data.frame(split(com.res2, grouping))
view(best_ref)

# get de. genes from each individual references
metadata(com.res2$orig.results$HPCA)$de.genes
metadata(com.res2$orig.results$DICE)$de.genes

merged_seurat_filtered$com.res2.labels <- com.res2[match(rownames(merged_seurat_filtered@meta.data), rownames(com.res2)), 'labels']
View(merged_seurat_filtered@meta.data)

# Combined diagnostics
plotScoreHeatmap(com.res2)
DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'com.res2.labels', label = TRUE)

# ...3. Strategy 3: Using Harmonized Labels ----------

hpca.ont <- celldex::HumanPrimaryCellAtlasData(cell.ont = 'nonna')
dice.ont <- celldex::DatabaseImmuneCellExpressionData(cell.ont = 'nonna')

# Using the same sets of genes:
shared <- intersect(rownames(hpca.ont), rownames(dice.ont))
hpca.ont <- hpca.ont[shared,]
dice.ont <- dice.ont[shared,]

# Showing the top 10 most frequent terms:
tail(sort(table(hpca.ont$label.ont)),10)
tail(sort(table(dice.ont$label.ont)), 10)

# using label.ont instead on label.main while running SingleR

com.res3 <- SingleR(test = pbmc_counts,
                    ref = list(HPCA = hpca.ont, DICE = dice.ont),
                    labels = list(hpca.ont$label.ont, dice.ont$label.ont))


table(com.res3$labels)

merged_seurat_filtered$com.res3.labels <- com.res3[match(rownames(merged_seurat_filtered@meta.data), rownames(com.res3)), 'labels']
View(merged_seurat_filtered@meta.data)
DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'com.res3.labels', label = TRUE)


# How to map cell ontology terms? ----------------

colData(hpca.ont)
colData(dice.ont)

hpca.fle <- system.file("mapping","hpca.tsv", package = "celldex")
hpca.mapping <- read.delim(hpca.fle, header = F)

DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'com.res3.labels', label = TRUE)


#___________________________________________________________________________

# 8. Finding differentially expressed features (cluster biomarkers) ---------
# Seurat can help you find markers that define clusters via differential expression. 

# ___find all markers of cluster 1 --------
cluster1.markers <- FindMarkers(merged_seurat_filtered, ident.1 = 2, min.pct = 0.25)
head (cluster1.markers, n = 5)

cluster1.markers_sorted <- cluster1.markers %>%
  arrange(desc(avg_log2FC))
head (cluster1.markers_sorted, n = 5)

# ___find all markers distinguishing cluster 5 from clusters 0 and 3 --------
cluster5.markers <- FindMarkers(merged_seurat_filtered, ident.1 = 1, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

cluster5.markers_sorted <- cluster5.markers %>%
  arrange(desc(avg_log2FC))
head (cluster5.markers_sorted, n = 20)

# ___find markers for every cluster compared to all remaining cells, report only the positive ones ---------
merged_seurat_filtered.markers <- FindAllMarkers(merged_seurat_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
merged_seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

#Find all markers
FindAllMarkers(merged_seurat_filtered,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')

#Find Conserved Marker

markers_cluster3 <- FindConservedMarkers(merged_seurat_filtered,
                                         ident.1 = 1,
                                         grouping.var = 'Location')

head(markers_cluster3)
write.csv(markers_cluster3, file = 'conservedmarkers.csv', row.names = TRUE)

#_________________________________________________________

#________________________________________________________
# Assuming 'singleR.labels' is the column you're using for grouping
# Assuming 'CD5L' is the feature you want to plot

# Count the number of cells in each group
group_counts <- table(merged_seurat_filtered$singleR.labels)

# Identify groups with at least 10 cells
valid_groups <- names(group_counts[group_counts >= 40])

# Filter Seurat object to include only valid groups
filtered_seurat <- subset(merged_seurat_filtered, subset = singleR.labels %in% valid_groups)

# Create the violin plot
VlnPlot(filtered_seurat, features = c("CD5L"), group.by = 'singleR.labels')
#__________________________________________________________
# 9. Visualization ---------
# VlnPlot() (shows expression probability distributions across clusters)
# FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. 
# RidgePlot(), CellScatter(), and DotPlot() as additional methods to view your dataset.

str(merged_seurat_filtered)
features <- merged_seurat_filtered@commands$RunPCA.RNA$features
VlnPlot(merged_seurat_filtered, features = c("CD5L"), group.by = 'seurat_clusters_manual')
VlnPlot(merged_seurat_filtered, features = c("COL1A1", "IGFBP7","FN1","LAT2","C1QA","CSTA","S100A9","LYZ"), group.by = 'Location')

# ___VlnPlot() - you can plot raw counts as well ---------
VlnPlot(merged_seurat_filtered, features = c("COL1A1", "IGFBP7"), slot = "counts", log = TRUE)


# ___FeaturePlot()- visualize feature expression in low-dimensional space ---------
FeaturePlot(merged_seurat_filtered, features = features[1:5])
FeaturePlot(merged_seurat_filtered, features = c("CCR7","CD27"), blend = TRUE)

# Visualize co-expression of two features simultaneously
FeaturePlot(merged_seurat_filtered, features = c("GZMK", "CCR7"), blend = TRUE)

# ___interactive plots --------
# Include additional data to display alongside cell names by passing in a data frame of
# information Works well when using FetchData
# works only with one feature
plot <- FeaturePlot(merged_seurat_filtered, features = c("CD5L"))
HoverLocator(plot = plot, information = FetchData(merged_seurat_filtered, vars = c("ident", "PC_1", "nFeature_RNA", "seurat_clusters_manual")))

# ___doHeatmap() --------
top10 <- merged_seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(merged_seurat_filtered, features = top10$gene) + NoLegend()

# ___RidgePlot() - Visualize single cell expression distribution in each cluster -------
RidgePlot(merged_seurat_filtered, features = features[1:5], ncol=2)
RidgePlot(merged_seurat_filtered, features = c("COL1A1", "IGFBP7","GZMK"), ncol=2)
# ___Dot plots - the size of the dot corresponds to the percentage of cells expressing the feature --------
# in each cluster. The color represents the average expression level
DotPlot(merged_seurat_filtered, features = features[1:5]) + RotatedAxis()

# ___Single cell heatmap of feature expression -------
DoHeatmap(subset(merged_seurat_filtered, downsample = 100), features = features[1:5], size = 3)

# assigning cell type identity to clusters
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
#                      "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## let's visualize top features
#FeaturePlot(ifnb_harmony, features = c('FCGR3A'), min.cutoff = 'q10')


# rename cluster 3 ident
#Idents(ifnb_harmony)
#ifnb_harmony <- RenameIdents(ifnb_harmony, `3` = 'CD16 Mono')

#DimPlot(ifnb_harmony, reduction = 'umap', label = T)

# cells already have annotations provided in the metadata
#View(ifnb_harmony@meta.data)
#______________________________________________________________________
#Differential Expression (NEW TAB)

# Normalize the data
Idents(merged_seurat_filtered) <- merged_seurat_filtered$Location

# Normalize data
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)

# Find DE features between
T1.de.markers <- FindMarkers(merged_seurat_filtered, ident.1 = "T1", ident.2 = "N1")

T1.de.markers <- T1.de.markers[order(-T1.de.markers$avg_log2FC), ]
T1.de.markers.asc <- T1.de.markers[order(+T1.de.markers$avg_log2FC), ]
# View results
head(T1.de.markers)
head(T1.de.markers.asc)
v1 <- VlnPlot(merged_seurat_filtered, features = c('SPP1','CCL7','SAA1','CCL8','MT1H'), idents = c("T1", "N1"), group.by = "Location") 
#change gene names as required
v2 <- VlnPlot(merged_seurat_filtered, features = c('IGKC','FGFBP2','IGLC2','HBB','IGLC3'), idents = c("T1", "N1"), group.by = "Location") 

v1|v2
#______________________________________________________________________
##Pseudobulk ((NEW TAB))

pseudo_merged_seurat_filtered <- AggregateExpression(merged_seurat_filtered, assays = "RNA", return.seurat = T, group.by = c("Location", "seurat_clusters_manual"))

tail(Cells(pseudo_merged_seurat_filtered))

Idents(pseudo_merged_seurat_filtered) <- pseudo_merged_seurat_filtered$Location

bulk.T1.de <- FindMarkers(object = pseudo_merged_seurat_filtered, 
                            ident.1 = "T1", 
                            ident.2 = "N1",
                            test.use = "DESeq2")

bulk.T1.de <- bulk.T1.de[order(-bulk.T1.de$avg_log2FC), ]
bulk.T1.de.asc <- bulk.T1.de[order(+bulk.T1.de$avg_log2FC), ]

head(bulk.T1.de, n = 15)
head(bulk.T1.de.asc, n = 15)

vb1 <- VlnPlot(merged_seurat_filtered, features = c('SPP1','SAA1','CCL7','MMP7','CCL8'), idents = c("T1", "N1"), group.by = "Location") 
#change gene names as required
vb2 <- VlnPlot(merged_seurat_filtered, features = c('IGKC','IGLC2','FGFBP2','HBB','IGLC3'), idents = c("T1", "N1"), group.by = "Location") 
vb1|vb2
#______________________________________________________________________
#Compare between SCDE and PSCDE ((NEW TAB))

names(bulk.T1.de) <- paste0(names(bulk.T1.de), ".bulk")
bulk.T1.de$gene <- rownames(bulk.T1.de)

names(T1.de.markers) <- paste0(names(T1.de.markers), ".sc")
T1.de.markers$gene <- rownames(T1.de.markers)

merge_dat <- merge(T1.de.markers, bulk.T1.de, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                 merge_dat$p_val.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val.bulk > 0.05 & 
                                  merge_dat$p_val.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                    merge_dat$p_val.sc > 0.05)]
print(paste0('# Common: ',length(common)))
print(paste0('# DEG Only present single-cell: ',length(only_sc)))
print(paste0('# DEG Only in bulk: ',length(only_bulk)))

# create a new column to annotate sample-condition-celltype in the single-cell dataset
merged_seurat_filtered$donor_id.stim <- paste0(merged_seurat_filtered$Location, "-", merged_seurat_filtered$seurat_clusters_manual)

# generate violin plot 
Idents(merged_seurat_filtered) <- merged_seurat_filtered$Location
print(merge_dat[merge_dat$gene%in%common[1:2],c('gene','p_val.sc','p_val.bulk')])

VlnPlot(merged_seurat_filtered, features = common[1:5], idents = c("T1", "N1"), group.by = "Location") 

VlnPlot(merged_seurat_filtered, features = common[1:2], idents = c("T1", "N1"), group.by = "donor_id.stim", ncol = 1) 

VlnPlot(merged_seurat_filtered, features = only_sc[1:5], idents = c("T1", "N1"), group.by = "Location") 

print(merge_dat[merge_dat$gene%in%c('KLF2','CCND3'),c('gene','p_val.sc','p_val.bulk')])

VlnPlot(merged_seurat_filtered, features <- c('SIX1','AC022509.1'), idents = c("T1", "N1"), group.by = "Location") 
VlnPlot(merged_seurat_filtered, features <- c('KLF2','CCND3'), idents = c("T1", "N1"), group.by = "donor_id.stim", ncol = 1) 
#_________________________________________________________________
##Monocle3 Workflow (NEW TAB)

unique(merged_seurat_filtered@meta.data$singleR.labels)

Idents(merged_seurat_filtered) <- merged_seurat_filtered$singleR.labels
e.seu <- subset(merged_seurat_filtered, idents = "Epithelial_cells")
e.seu
unique(e.seu@meta.data$singleR.labels)

# pre-processing using seurat
e.seu <- NormalizeData(e.seu)
e.seu <- FindVariableFeatures(e.seu)
e.seu <- ScaleData(e.seu)
e.seu <- RunPCA(e.seu)
e.seu <- FindNeighbors(e.seu, dims = 1:30)
e.seu <- FindClusters(e.seu, resolution = 0.9)
e.seu <- RunUMAP(e.seu, dims = 1:30, n.neighbors = 50)

a1 <- DimPlot(e.seu, reduction = 'umap', group.by = 'singleR.labels', label = T)
a2 <- DimPlot(e.seu, reduction = 'umap', group.by = 'seurat_clusters', label = T)
a1|a2
# ...1 Convert to cell_data_set object ------------------------

cds <- as.cell_data_set(e.seu)
cds

# to get cell metadata
colData(cds)
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)



# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- e.seu@active.ident
cds@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- e.seu@reductions$umap@cell.embeddings



# plot

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "singleR.labels",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names



# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'singleR.labels',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)


# ...4. Order the cells in pseudotime -------------------

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 0]))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
##CHANGE THE seurat_clusters if you have actual clusters
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) +
  geom_boxplot()




# ...5. Finding genes that change as a function of pseudotime --------------------
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

FeaturePlot(e.seu, features = c('MALAT1', 'CCNB2', 'RPS4X', 'RPS14', 'RPLP1', 'RPS2'))


# visualizing pseudotime in seurat

e.seu$pseudotime <- pseudotime(cds)
Idents(e.seu) <- e.seu$seurat_clusters
FeaturePlot(e.seu, features = "pseudotime", label = T)
#_________________________________________________________________
#SAVE (NEW TAB)

# 10. saving processed data ------
saveRDS(merged_seurat_filtered, file = "C:/SCPs/Kidney/GSE111360/BII/BII/integrated/processed.rds")
