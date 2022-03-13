library(reticulate)
library(dplyr)
library(tidyverse)
library(readxl)
library(stringr)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(scales)
library(plotly)
library(RcppAnnoy)
library(circlize)
library(ComplexHeatmap)
library(zoo)
library(Cairo)
library(Seurat)

ISC_data <- readRDS("Mouse_Data_Regev_SmartSeq.log2.TPM.Capitalized.RDS")
ISC_data <- NormalizeData(ISC_data) %>% ScaleData(features = rownames(ISC_data))
ISC_data <- FindVariableFeatures(ISC_data, selection.method = "vst", nfeatures = 2000)

#Setup for Seurat object
ISC_counts <- Read10X("/labs/longaker/USR/KE_Bauer-Rowe/kemgfolder/205 Project/Atlas1")
ISC_seurat <- CreateSeuratObject(counts = ISC_counts, min.features = 100)
ISC_seurat[["percent.mt"]] <- PercentageFeatureSet(ISC_seurat, pattern = "^mt-");  
#VlnPlot(ISC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3); plot4=plot3;
#plot1 <- FeatureScatter(ISC_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(ISC_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ISC_seurat <- subset(ISC_seurat, subset = nFeature_RNA > 0 & nFeature_RNA < 3000 & nCount_RNA < 30000 & percent.mt < 15)
ISC_seurat <- NormalizeData(ISC_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
ISC_seurat <- FindVariableFeatures(ISC_seurat, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(ISC_seurat)
ISC_seurat <- ScaleData(ISC_seurat, features = all.genes)
ISC_seurat <- RunPCA(ISC_seurat, features = VariableFeatures(object = ISC_seurat))
#DimPlot(ISC_seurat, reduction = "pca")
#DimHeatmap(ISC_seurat, dims = 1, cells = 500, balanced = TRUE)
ISC_seurat <- FindNeighbors(ISC_seurat, dims = 1:20)
ISC_seurat <- FindClusters(ISC_seurat, resolution = 0.50)
ISC_seurat <- RunUMAP(ISC_seurat, dims = 1:20)
DimPlot(ISC_seurat, reduction = "umap")
VlnPlot(ISC_seurat, features = c("Tmem238"))
FeaturePlot(ISC_seurat, features = c("Tmem238"))
ISC.markers <- FindAllMarkers(ISC_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ISC.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(ISC_seurat, features = top10$gene, size = 3.0) + NoLegend() + theme(text = element_text(size = 5)) + ggtitle("Mesentery Global Clusters")

#Rename the clusters 
new.cluster.ids <- c("ISC", "TA", "EP Early", "EP Late", "Goblet", "EEC",
                     "Paneth", "Enterocyte", "Tuft")
names(new.cluster.ids) <- levels(ISC_seurat)
ISC_seurat <- RenameIdents(ISC_seurat, new.cluster.ids)
ISC_seurat$new.cluster.ids <- ISC_seurat@active.ident


#ISC_data <- FindVariableFeatures(ISC_data, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(ISC_data)
#ISC_data <- ScaleData(ISC_data, features = all.genes)
#ISC_data <- RunUMAP(ISC_data, dims = 1:30)
#DimPlot(ISC_data, reduction = "umap",group.by = "phenotype")

order = c("ISC", "TA", "EP Early", "EP Late", "Enterocyte", "EEC", "Goblet", "Paneth", "Tuft")

#Import table 
m.table <- read.table("mm10_ordered_v3.txt", sep = " ", header = T) # read table
m.table <- m.table[!duplicated(m.table$gene),] # remove any duplicate genes
rownames(m.table) <- m.table$gene # set the rownames to genes
#head(m.table, n=20) # check to make sure table looks good

ChromE(ISC_seurat, m.table, annotation = "new.cluster.ids", order = order, assay = "RNA", k = 10, plot.range = c(0,6))
ISC_seurat_small <- subset(ISC_seurat, select = c("ISC", "TA", "EP Early", "EP Late", "Enterocye" ))
#Boxplot

ISC_split <- SplitObject(ISC_seurat, split.by = "new.cluster.ids")
boxplot <- reg_count_boxplot(ISC_seurat, order = order, split.by = "new.cluster.ids", reg.cutoff = 0, frac = F)
boxplot


