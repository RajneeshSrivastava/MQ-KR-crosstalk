## Project: *Identification of a physiologic vasculogenic fibroblast state to achieve tissue repair*

### Upload libraries
```
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(sctransform)
library(monocle3)
library(SingleCellExperiment)
```

### Create the main seurat object; including HADF, MI day 1, MI day 3, MI day 5, MI day 7 
```
# Store path to directories for each sample
MI_Day_1 = "/MI-Treated_HDAF_D1/filtered_feature_bc_matrix"
MI_Day_3 = "/MI-Treated_HDAF_D3/filtered_feature_bc_matrix"
MI_Day_5 = "/MI-Treated_HDAF_D5/filtered_feature_bc_matrix"
MI_Day_7 = "/MI-Treated_HDAF_D7/filtered_feature_bc_matrix"
HDAF = "/Untreated_HDAF/filtered_feature_bc_matrix"

# Read the count matrix, barcodes and gene names files for each sample
MI_Day_1 = Read10X(MI_Day_1)
MI_Day_3 = Read10X(MI_Day_3)
MI_Day_5 = Read10X(MI_Day_5)
MI_Day_7 = Read10X(MI_Day_7)
HDAF = Read10X(HDAF)

# Create individual Seurat object for each sample
MI_Day_1 = CreateSeuratObject(counts = MI_Day_1, min.cells = 3, min.features = 200, project = "MI Day 1" )
MI_Day_3 = CreateSeuratObject(counts = MI_Day_3, min.cells = 3, min.features = 200, project = "MI Day 3" )
MI_Day_5 = CreateSeuratObject(counts = MI_Day_5, min.cells = 3, min.features = 200, project = "MI Day 5" )
MI_Day_7 = CreateSeuratObject(counts = MI_Day_7, min.cells = 3, min.features = 200, project = "MI Day 7" )
HDAF = CreateSeuratObject(counts = HDAF, min.cells = 3, min.features = 200, project = "HDAF" )

# Calculate mitochondrial gene expression percentage per cell
MI_Day_1[["percent.mt"]] <- PercentageFeatureSet(MI_Day_1, pattern = "^MT-")
MI_Day_3[["percent.mt"]] <- PercentageFeatureSet(MI_Day_3, pattern = "^MT-")
MI_Day_5[["percent.mt"]] <- PercentageFeatureSet(MI_Day_5, pattern = "^MT-")
MI_Day_7[["percent.mt"]] <- PercentageFeatureSet(MI_Day_7, pattern = "^MT-")
HDAF[["percent.mt"]] <- PercentageFeatureSet(HDAF, pattern = "^MT-")

# Filter out low quality cells
MI_Day_1 <- subset(MI_Day_1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA < 60000 & nCount_RNA > 500)
MI_Day_3 <- subset(MI_Day_3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA < 60000 & nCount_RNA > 500)
MI_Day_5 <- subset(MI_Day_5, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA < 60000 & nCount_RNA > 500)
MI_Day_7 <- subset(MI_Day_7, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA < 60000 & nCount_RNA > 500)
HDAF <- subset(HDAF, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA < 60000 & nCount_RNA > 500)

# Perform log normalization
MI_Day_1 <- NormalizeData(MI_Day_1, normalization.method = "LogNormalize", scale.factor = 10000)
MI_Day_3 <- NormalizeData(MI_Day_3, normalization.method = "LogNormalize", scale.factor = 10000)
MI_Day_5 <- NormalizeData(MI_Day_5, normalization.method = "LogNormalize", scale.factor = 10000)
MI_Day_7 <- NormalizeData(MI_Day_7, normalization.method = "LogNormalize", scale.factor = 10000)
HDAF <- NormalizeData(HDAF, normalization.method = "LogNormalize", scale.factor = 10000)

# Find top 2,000 variable genes for each sample
MI_Day_1 <- FindVariableFeatures(MI_Day_1, selection.method = "vst", nfeatures = 2000)
MI_Day_3 <- FindVariableFeatures(MI_Day_3, selection.method = "vst", nfeatures = 2000)
MI_Day_5 <- FindVariableFeatures(MI_Day_5, selection.method = "vst", nfeatures = 2000)
MI_Day_7 <- FindVariableFeatures(MI_Day_7, selection.method = "vst", nfeatures = 2000)
HDAF <- FindVariableFeatures(HDAF, selection.method = "vst", nfeatures = 2000)

# Integrate the individual Seurat objects into 1 main object
anchors <- FindIntegrationAnchors(object.list = list(HDAF, MI_Day_1, MI_Day_3, MI_Day_5, MI_Day_7), dims = 1:20)
main <- IntegrateData(anchorset = anchors, dims = 1:20)

# Perform scaling and PCA
main <- ScaleData(main, verbose = TRUE)
main <- RunPCA(object = main, npcs = 50, verbose = FALSE)
ElbowPlot(main, ndims = 40)

# Run tSNE
main = RunTSNE(main, dims = 1:15)

# Louvain Clustering
main <- FindNeighbors(object = main, dims = 1:20)
main <- FindClusters(main, resolution = 0.1)

# Isolate cluster 1 cells
cl1 = subset(main, subset = seurat_clusters == 1)
DefaultAssay(cl1)='RNA'
cl1 <- FindVariableFeatures(object = cl1, selection.method = "vst", nfeatures = 2000)
cl1 <- ScaleData(cl1, verbose = TRUE)
cl1 <- RunPCA(object = cl1, npcs = 50, verbose = FALSE)
ElbowPlot(cl1, ndims = 40)
cl1 = RunTSNE(cl1, dims = 1:20)
cl1 <- FindNeighbors(object = cl1, dims = 1:20)
cl1 <- FindClusters(cl1, resolution = 0.25)
```
### Trajectory inference for cluster 1 cells
```
# Get high variable genes among cluster 1 cells
genes = cl1@assays$RNA@var.features

# Select 6,000 cells as downsampling
object.downsample = subset(cl1, cells = sample(Cells(cl1), 6000))

# Trajectory inference using Monocle2
library(monocle)
pd <- new('AnnotatedDataFrame', data = object.downsample@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds = estimateSizeFactors(monocle_cds)
monocle_cds = estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(monocle_cds),
                                    num_cells_expressed >= 10))
pData(monocle_cds)$Total_mRNAs <- Matrix::colSums(exprs(monocle_cds))
monocle_cds <- monocle_cds[,pData(monocle_cds)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(monocle_cds)$Total_mRNAs)) +
                     2*sd(log10(pData(monocle_cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(monocle_cds)$Total_mRNAs)) -
                     2*sd(log10(pData(monocle_cds)$Total_mRNAs)))

monocle_cds <- monocle_cds[,pData(monocle_cds)$Total_mRNAs > lower_bound &
                             pData(monocle_cds)$Total_mRNAs < upper_bound]
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)

L <- log(exprs(monocle_cds[expressed_genes,]))

library(reshape)
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

monocle_cds <- setOrderingFilter(monocle_cds, genes)

plot_ordering_genes(monocle_cds)

monocle_cds <- reduceDimension(monocle_cds, max_components = i,
                               method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "orig.ident")+facet_wrap(~State, nrow = 1)
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("PIK3R1"))+facet_wrap(~State,nrow = 2)

# Cluster 0 subclustering
cl0 = subset(main, subset=seurat_clusters==0)
cl0 = FindVariableFeatures(cl0, selection.method = "vst", nfeatures = 2000)
cl0 = ScaleData(cl0)
cl0 = RunPCA(cl0, npcs = 50)
cl0 = RunTSNE(cl0, dims =1:15)
cl0 = FindNeighbors(cl0, dims = 1:15)
cl0 = FindClusters(cl0, resolution = 0.3)

# Trajectory inference using all clusters
library(monocle3)
data <- main@assays$RNA@data # Use the log normalized counts
cell_metadata <- new('AnnotatedDataFrame', data = main@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
gene_metadata <- new('AnnotatedDataFrame', data = fData)

gene_metadata <- as(gene_metadata, "data.frame")
cell_metadata <- as(cell_metadata, "data.frame")

#Construct monocle object
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata
                         )

#Convert seurat clusters column into categorical column
cds@colData$seurat_clusters = as.factor(cds@colData$seurat_clusters)

cds <- preprocess_cds(cds, num_dim = 20)
cds = reduce_dimension(cds, umap.min_dist = 1)
cds = cluster_cells(cds)
cds = learn_graph(cds, use_partition = F)
cds = order_cells(cds)
```
### Main figures
```
# Main Fig 1B
DimPlot(main, split.by = 'orig.ident',group.by = 'seurat_clusters')

# Main Fig 1C
Idents(main) = 'seurat_clusters'
VlnPlot(main, 'THY1')
VlnPlot(main, 'S100A4')
VlnPlot(main, 'VIM')
VlnPlot(main, 'FAP')

# Main Fig 1D
genes = c('COL1A1','COL3A1','COL1A2','TAGLN','FKBP1A','USP53','PTX3','IGFBP2','HIST1H4C','CENPF','HMGB2')
DotPlot(main, features = genes)+theme(axis.text.x = element_text(angle = 45, hjust=1))

# Main Fig 1E
genes = c('COL1A1','COL5A2','COL3A1','VEGFB','VEGFC','NRP1')
cl1 = subset(main, subset = seurat_clusters==1)
Idents(cl1) = 'orig.ident' # To plot per sample (time point)
for (i in genes)
{
  VlnPlot(main, i)
  VlnPlot(cl1, i)
}

# Main Fig 3A
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")

# Main Fig 3B
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")+facet_wrap(~orig.ident, nrow = 1)
plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = FALSE)+facet_wrap(vars(orig.ident))+geom_density_2d(bins = 8, n = 30, h = 5, alpha = 0.7)

# Main Fig 3C
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("S100A4"))+facet_wrap(~orig.ident,nrow = 2)
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("VEGFB"))+facet_wrap(~orig.ident,nrow = 2)

# Main Fig 3E
plot_cell_trajectory(monocle_cds, color_by = "State")

# Main Fig 3F
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")+facet_wrap(~seurat_clusters, nrow = 1)

# Main Fig 4A
DimPlot(main, group.by = 'seurat_clusters', cols = c('grey','#817fcc','grey','grey'))
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")

# Main Fig 4B
DimPlot(cl1, group.by = 'seurat_clusters')

# Main Fig 4D
DimPlot(cl1, group.by = 'seurat_clusters', split.by = 'orig.ident')
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")+facet_wrap(~orig.ident, nrow = 1)
```
### Extended Data Figures
```
# Extended data Fig 1F
genes = c('CD151','COL14A1','PLOD3','LOXL1','CTSL','SERPINH1','COLGALT1','CTSB','PPP1R14B','DST','MMP1','COL1A1','COL3A1','COL1A2','COL5A2','COL4A2','LOX','P4HA1','COL4A1','P4HA2','COL6A2','COL5A3','COL5A2','COL6A1','COL6A3','PLEC')
Idents(main) = 'orig.ident'
main = ScaleData(main, features = genes)
DoHeatmap(main,genes)

# Extended data Fig 2A and B
genes = c('VEGFA','NRP2','FLT1','KDR','NOS3','AKT1','MAP2K1','MAP2K2','RAF1','KRAS','PIK3CA','PTEN')
for (i in genes)
{
  VlnPlot(main, i)
  VlnPlot(cl1, i)
}

# Extended data Fig 2C
# Create an object using HADF sample, day 1 CI (control inhibitor), day 1 MI (miR-200b inhibitor), day 7 CI and day 7 MI
hadf = Read10X('Path to HADF sample folder')
ci_1 = Read10X('Path to CI day 1 sample folder')
mi_1 = Read10X('Path to MI day 1 sample folder')
ci_7 = Read10X('Path to CI day 7 sample folder')
mi_7 = Read10X('Path to MI day 7 sample folder')

hadf = CreateSeuratObject(counts = hadf, min.cells = 3, min.features = 0, project = "HADF")
ci_1 = CreateSeuratObject(counts = ci_1, min.cells = 3, min.features = 0, project = "CI day 1")
mi_1 = CreateSeuratObject(counts = mi_1, min.cells = 3, min.features = 0, project = "MI day 1")
ci_7 = CreateSeuratObject(counts = ci_7, min.cells = 3, min.features = 0, project = "CI day 7")
mi_7 = CreateSeuratObject(counts = mi_7, min.cells = 3, min.features = 0, project = "MI day 7")

objects.list = c(hadf,ci_1,mi_1,ci_7,mi_7)

# Calculate mitochondrial gene expression percentage
for (i in 1:length(objects.list))
{
  objects.list[[i]][["percent.mt"]] = PercentageFeatureSet(objects.list[[i]], pattern = "^MT-")
}

# Log normalization
for (i in 1:length(objects.list))
{
  objects.list[[i]] <- NormalizeData(objects.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}

# Remove low quality cells
for (i in 1:length(objects.list))
{
  objects.list[[i]] <- subset(objects.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA < 60000 & nCount_RNA > 500)
}

# Integrate samples using CCA
anchors <- FindIntegrationAnchors(object.list = objects.list, dims = 1:20)
merged <- IntegrateData(anchorset = anchors, dims = 1:20)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 50)
merged <- FindNeighbors(object = merged, dims = 1:20)
merged <- FindClusters(merged, resolution = 0.1)

DefaultAssay(merged) = 'RNA'
VlnPlot(merged, 'VIM')
VlnPlot(merged, 'S100A4')
VlnPlot(merged, 'VEGFB')
VlnPlot(merged, 'GLUL')

# Extended data Fig 2D
merged_cluster1 = subset(merged, subset = seurat_clusters == 1)
VlnPlot(merged_cluster1, 'VIM')
VlnPlot(merged_cluster1, 'S100A4')
VlnPlot(merged_cluster1, 'VEGFB')
VlnPlot(merged_cluster1, 'GLUL')

# Extended data Fig 3A
genes = c('LRP1','TNFRSF1A','ITGB1','IL6ST','IGF2R','TSPO','B2M','ANTXR1','KTN1','FCGRT','ITGA2','SFRP1','TNFRSF11B','DNER','TNFRSF12A','PLXNA2','SIGIRR','HMMR','ITGA5')
Idents(main) = 'seurat_clusters'
DotPlot(main, features = genes)+theme(axis.text.x = element_text(angle = 45, hjust=1))

# Extended data Fig 3D
# Cell Cycle Scoring
main = CellCycleScoring(main, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
clusters0_1 = subset(main, subset = seurat_clusters==0 | seurat_clusters==1)
plot=ggplot(clusters0_1@meta.data,aes(clusters0_1@meta.data$S.Score, clusters0_1@meta.data$G2M.Score, color=clusters0_1@meta.data$Phase)) + 
  geom_point(size = 0.1) +
  facet_wrap(~seurat_clusters)+
  scale_color_brewer(palette="Paired")+
  geom_rug(size=0.1)+
  stat_ellipse()+
  labs(title="Cell cycle scoring",
       x="S score ", y = "G2M score")+
  theme_classic()

# Extended data Fig 3E
# Heatmap for cell cycle genes
cycling_genes = append(cc.genes$s.genes,cc.genes$g2m.genes)
cl0 = subset(main, subset = seurat_clusters == 0)
cl1 = subset(main, subset = seurat_clusters == 1)
DoHeatmap(cl0, features = test, group.by = 'orig.ident')
DoHeatmap(cl1, features = test, group.by = 'orig.ident')

# Extended data Fig 4A
DimPlot(main,group.by = 'seurat_clusters')
DimPlot(object.downsample,group.by = 'seurat_clusters')
DimPlot(main,group.by = 'seurat_clusters', cols = c('grey','#817fcc','grey','grey'))
DimPlot(object.downsample,group.by = 'seurat_clusters', cols = c('grey','#817fcc','grey','grey'))
DimPlot(cl1,group.by = 'seurat_clusters', cols = c('grey','#817fcc','grey','grey'))
DimPlot(cl1.downsample,group.by = 'seurat_clusters', cols = c('grey','#817fcc','grey','grey'))

# Extended data Fig 4B
# Upload genes provided in supplementary files
Idents(cl1) = 'seurat_clusters'
cl1 = ScaleData(cl1, features = genes)
DoHeatmap(cl1, genes)

# Extended data Fig 4C
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("COL1A2"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("COL6A1"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("COL6A3"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("VEGFB"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("CD147"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("HMGB3"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("SERPINF1"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("SPARC"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("ATP2B4"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("CDC42"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("FLNA"))
plot_cell_trajectory(monocle_cds, use_color_gradient = TRUE, markers = c("ITGB1"))

# Extended data Fig 4D
DimPlot(main, cols = c('red','grey','grey','grey'))
DimPlot(cl0, group.by = 'seurat_clusters')

# Extended data Fig 4E
DimPlot(cl0, group.by = 'seurat_clusters', split.by = 'orig.ident')

# Extended data Fig 4F
jpeg(paste0("Trajectory Inference (by pseudotime).jpg"), res=300,width = 1700,height = 1200);print(plot_cells(cds, color_cells_by = "pseudotime", label_leaves = F,label_cell_groups = F,label_roots = F, label_groups_by_cluster = F, label_branch_points = F, trajectory_graph_segment_size = 0.5));dev.off()

# Extended data Fig 4G
jpeg(paste0("Trajectory Inference (by main cluster).jpg"), res=300,width = 1800,height = 1200);print(plot_cells(cds, color_cells_by = "seurat_clusters",label_cell_groups = F,label_roots = F, label_leaves = F, label_groups_by_cluster = F, label_branch_points = F, trajectory_graph_segment_size = 0.5));dev.off()

# Extended data Fig 4H
jpeg(paste0("Trajectory Inference (by subclusters).jpg"), res=300,width = 1800,height = 1200);print(plot_cells(cds, color_cells_by = "subclusters_both", label_leaves = F,label_cell_groups = F,label_roots = F, label_groups_by_cluster = F, label_branch_points = F, trajectory_graph_segment_size = 0.5));dev.off()

# Extended data Fig 4I
genes = c('CD151','COL14A1','PLOD3','LOXL1','CTSL','SERPINH1','COLGALT1','CTSB','PPP1R14B','DST','MMP1','COL1A1','COL3A1','COL1A2','COL5A1','COL4A2','LOX','P4HA1','COL4A1','P4HA2','COL6A2','COL5A3','COL5A2','COL6A1','COL6A3','PLEC')
custom.subsets = subset(main, subset = subclusters_both == '0B'|subclusters_both == '0C'|subclusters_both == '1A'|subclusters_both == '1B')
custom.subsets = ScaleData(custom.subsets, features = genes)
DoHeatmap(custom.subsets, genes)

# Extended data Fig 4J
# Read the gene signature of fibrosis "fibro" and mechanotransduction "mechano"
#Add module score function
main = AddModuleScore(main, features = list(fibro), name = "fibro")
main = AddModuleScore(main, features = list(mechano), name = "mechano")

main@meta.data %>%
  ggplot( aes(x=fibro1)) +
  geom_histogram(binwidth=0.04, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Fibrosis score") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
plot_cells(cds, color_cells_by = "fibro1", label_leaves = F,label_cell_groups = F,label_roots = F, label_groups_by_cluster = F, label_branch_points = F, trajectory_graph_segment_size = 0.3)+scale_colour_gradient2(limits=c(-0.1, 0.3),low = "green",high = "red", midpoint = 0, space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "colour")

# Extended data Fig 4K
main@meta.data %>%
  ggplot( aes(x=mechano1)) +
  geom_histogram(binwidth=0.04, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Mechanotransduction score") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
plot_cells(cds, color_cells_by = "mechano1", label_leaves = F,label_cell_groups = F,label_roots = F, label_groups_by_cluster = F, label_branch_points = F, trajectory_graph_segment_size = 0.3)+scale_colour_gradient2(limits=c(-0.1, 0.3),low = "green",high = "red", midpoint = 0, space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "colour")

```
#### Extended data Fig. 5
```
# UPLOAD DATA from GSE165816 (from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165816)]
DFU.dir = "path/GSE165816_dataset/"

#Create seurat objects
DFU.list=list("G6", "G9", "G33", "G34", "G39", "G7", "G8", "G49", "G42", "G45", "G23", "G4", "G2", "G15")

for (file in DFU.list)
  {
  DFU_data <- read.csv(file = paste0(DFU.dir, file,".csv.gz"), header = TRUE, row.names = 1)
  DFU_obj <- CreateSeuratObject(counts = 
                                DFU_data,
                              min.cells=3,
                              min.features = 200,
                              project = file)
                              assign(file, DFU_obj)
  }

sample.list=list(G6,G9,G33,G34,G39,G7,G8,G49,G42,G45,G23,G4,G2,G15)

# Quality Filtering and SCTransformation
for (i in 1:length(sample.list)) {
  sample.list[[i]][["percent.mt"]] <- PercentageFeatureSet(sample.list[[i]], pattern = "^MT-")
  sample.list[[i]] <- subset(sample.list[[i]], 
                             subset = nFeature_RNA > 200 & 
                               nFeature_RNA < 5000 & percent.mt < 15 & 
                               nCount_RNA < 25000 & nCount_RNA > 2000)}


#SCT transformation
sample.list[[i]] <- SCTransform(sample.list[[i]],
                                vars.to.regress = "percent.mt", 
                                verbose = FALSE)
  
# Data integration
sample.features <- SelectIntegrationFeatures(
  object.list = sample.list,    
  nfeatures = 3000)

sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- RunPCA(x, features = sample.features)
})
sample.list <- PrepSCTIntegration(object.list = sample.list, 
                                  anchor.features = sample.features,
                                  verbose = FALSE)           
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, 
                                         normalization.method = "SCT", 
                                         anchor.features = sample.features, 
                                         reference = c(1, 2), 
                                         reduction = "rpca", 
                                         verbose = TRUE)
sample.integrated <- IntegrateData(anchorset = sample.anchors, 
                                   normalization.method = "SCT",
                                   verbose = TRUE)

# Dimentionality reduction and clustering
sample.integrated <- RunPCA(object = sample.integrated, verbose = FALSE) 

#RUN tSNE/UMAP
sample.integrated = RunUMAP(sample.integrated, dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated)

#tweak-in sample code as NHDHU or HDHU
meta=read.table("metadata.txt",sep="\t", header=T)
head(meta,5)

GSM=sample.integrated@meta.data
GSM$group=1
for(i in 1:nrow(GSM)){
  for (j in 1:nrow(meta)){
    if (GSM[i,1] == meta[j,1]){
      GSM[i,7] = meta[j,2] }
  }
}

sample.integrated@meta.data=GSM

## Find clusters with defined resolution
sample25 <- FindClusters(sample.integrated, resolution = 0.25)
sample25$group <- factor(sample25$group, levels = c("NHDFU", "HDFU"))

# Extended data Fig 5A
DimPlot(sample25)

#Find all markers
sample25_markers=FindAllMarkers(sample25,only.pos = T,
                                logfc.threshold = 0.30, min.pct = 0.10)
write.table(sample25_markers,"All_markers_sample25.txt")

# Explore the table and identify cell types using published articles and PanglaoDb (https://panglaodb.se/))
Fb_markers=c("COL1A1","COL1A2","DCN")
#FeaturePlot(sample25, features=Fb_markers, cols=c("grey","red","black"))

# Extended data Fig 5B
# Isolate Fb cluster and do re-clustering
#DefaultAssay(sample25)="integrated"
Fb_0=subset(sample25,subset=seurat_clusters==0)
Fb_re <- RunPCA(object = Fb_0, verbose = FALSE)
Fb_re <- RunUMAP(Fb_re, dims = 1:30)
Fb_re <- FindNeighbors(Fb_re)
Fb_re10 <- FindClusters(Fb_re, resolution = 0.10)
DimPlot(Fb_re10) #Extended Data Fig. 5b

# Extended data Fig 5C
## Find average expression across samples
Idents(Fb_0)= "orig.ident"
AvgExpFb_0 = AverageExpression(Fb_0, return.seurat = FALSE, verbose = TRUE)
#write.table(AvgExpFb_0$SCT,"Fb_0_avgExp_SCT.txt") # Explore the table for VF signatures (refer Fig. 6a)

Idents(Fb_0)= "seurat_clusters"
VFS=list(c("EPAS1", "FOSL1", "GLUL", "HIF1A", "MMP14", "NRP2", "VEGFC", "TGFBI")) #significant VF signatures
DoHeatmap(Fb_0,features=VFS) #Extended Data Fig. 5c

# Extended data Fig 5D
#Compute module score for VF like cells based on above mentioned significant markers
scored.object <- AddModuleScore(object = Fb_0, features = VFS, name = "VFS")
FeaturePlot(object = scored.object, features = "VFS1", cols=c("grey","red","black"))

#Sort VF and non-VF population across the subclusters:
p <- VlnPlot(object = scored.object, features =c("VFS1"))
table(p$data$VF_score1>0.0,p$data$split) #Complute proportions

#Idents(object)="VFS1"
VF=subset(object,subset=VFS1>=0)
non-VF=subset(object,subset=VFS1<0)

# Extended data Fig 5E
Idents(VF)="seurat_clusters"
DimPlot(VF)     #Extended Data Fig. 5e

# Extended data Fig 5F
Idents(non-VF)="seurat_clusters"
DimPlot(non-VF) #Extended Data Fig. 5f

# Extended data Fig 5G-J
#####Cell-Chat-Analysis
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
future::plan("multiprocess", workers = 4)

#Preparation of the cellchat objects
#replacing original fibroblast cluster (#0) with VF and non-VF cells

Idents(sample25)="seurat_clusters"
DFU_0=subset(sample25,subset=seurat_clusters=="0")
DFU_N0=subset(sample25,subset=seurat_clusters!="0")
#saveRDS(DFU_N0,"DFU_N0.rds")

#merging seurat objects
#Merging two seurat objects (with existing normalization)
DFU_VF=merge(DFUN0,VF, merge.data = TRUE)
DFU_non-VF=merge(DFUN0,non-VF, merge.data = TRUE)

#DefaultAssay(DFU_VF)="SCT"
#DefaultAssay(DFU_non-VF)="SCT"

#Create CellChat objects
DFU_VF_cellchat=createCellChat(object=DFU_VF,meta=DFU_VF@meta.data,group.by="integrated_snn_res.0.25")
DFU_non-VF_cellchat=createCellChat(object=DFU_non-VF,meta=DFU_non-VF@meta.data,group.by="integrated_snn_res.0.25")

# Assign Idents
DFU_VF_cellchat <-setIdent(DFU_VFS_HDFU_cellchat, ident.use="integrated_snn_res.0.25")
DFU_non-VF_cellchat <-setIdent(DFU_VFS_NHDFU_cellchat, ident.use="integrated_snn_res.0.25")

#Database usage
CellChatDB<- CellChatDB.human
#showDatabaseCategory(CellChatDB)
CellChatDB.use<- CellChatDB
DFU_VF_cellchat@DB <-CellChatDB.use
DFU_non-VF_cellchat@DB <-CellChatDB.use

#subset the dataset for enhanced performance
DFU_VF_cellchat <-subsetData(DFU_VF_cellchat)
DFU_non-VF_cellchat <-subsetData(DFU_non-VF_cellchat)

# Individual Analysis
#VF included cellchat
DFU_VF_cellchat <- identifyOverExpressedGenes(DFU_VF_cellchat)
DFU_VF_cellchat <- identifyOverExpressedInteractions(DFU_VF_cellchat)
DFU_VF_cellchat <- projectData(DFU_VF_cellchat, PPI.human) #(optional)
DFU_VF_cellchat <- computeCommunProb(DFU_VF_cellchat, raw.use = TRUE)
DFU_VF_cellchat <- filterCommunication(DFU_VF_cellchat,min.cells = 10)
DFU_VF_cellchat <- computeCommunProbPathway(DFU_VF_cellchat)
DFU_VF_cellchat <- netAnalysis_computeCentrality(DFU_VF_cellchat, slot.name = "netP")
DFU_VF_cellchat <- aggregateNet(DFU_VF_cellchat)

#non-VF included cellchat
DFU_non-VF_cellchat <- identifyOverExpressedGenes(DFU_non-VF_cellchat)
DFU_non-VF_cellchat <- identifyOverExpressedInteractions(DFU_non-VF_cellchat)
DFU_non-VF_cellchat <- projectData(DFU_non-VF_cellchat, PPI.human) #(optional)
DFU_non-VF_cellchat <- computeCommunProb(DFU_non-VF_cellchat, raw.use = TRUE)
DFU_non-VF_cellchat <- filterCommunication(DFU_non-VF_cellchat,min.cells = 10)
DFU_non-VF_cellchat <- computeCommunProbPathway(DFU_non-VF_cellchat)
DFU_non-VF_cellchat <- netAnalysis_computeCentrality(DFU_non-VF_cellchat, slot.name = "netP")
DFU_non-VF_cellchat <- aggregateNet(DFU_non-VF_cellchat)

#Merging cell chats
cellchat.list=c(VF=DFU_VF_cellchat,non-VF=DFU_non-VF_cellchat)
merged_cellchat <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))

#cell_chat_objects
####saveRDS(DFU_VF_cellchat,"DFU_VF_cellchat.rds")
####saveRDS(DFU_non-VF_cellchat,"DFU_non-VF_cellchat.rds")
####saveRDS(merged_cellchat,"VF_non-VF_merged_cellchat.rds")

# Extended data Fig 5G
#Comparison of cell chat objects
#Interactions
weight.max <- getMaxWeight(cellchat.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchat.list)) {
  netVisual_circle(cellchat.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(DFU_nonVFS_chat.list)[[i]]))
}

# Extended data Fig 5H
g1=compareInteractions(merged_cellchat, show.legend = F, group = c(1,2))
g2=compareInteractions(merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
g1+g2 

#Relative Information_flow
rankNet(merged_cellchat, comparison=c(1,2),stacked = T,do.stat = TRUE)

# Extended data Fig 5I
#DE_network
netVisual_diffInteraction(merged_cellchat, weight.scale = T,comparison = c(1,2), top=0.10)

# Extended data Fig 5J
#Pathway specific analysis 
#top 4 pathways significantly enriched in DFU_VF cell chat as identified in information flow module above.
netVisual_aggregate(DFU_Fb0_HDFU_cellchat, signaling = "ncWNT",
                    layout="circle")
```
### Thank you
