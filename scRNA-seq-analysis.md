## Project: *The Mitochondrial Bioenergetics of Functional Wound Closure is Dependent on Macrophage-Keratinocyte Exosomal Crosstalk*


##### Upload libraries
```
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(hdf5r)
library(future)
library(sctransform)
library(CellChat)
library(patchwork)

plan("multisession", workers = 4)
options(future.globals.maxSize = 15000 * 1024^2)
options(stringsAsFactors = FALSE)
```
### Single-cell data analysis using seurat
##### Upload single-cell RNA-seq data obtained from GSE165816 (from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165816)]
###### Create seurat objects
```
DFU.dir = "path/GSE165816_dataset/"
DFU.list=list("G6", "G9", "G33", "G34", "G39", "G7", "G8", "G49", "G42", "G45", "G23", "G4", "G2", "G15")

for (file in DFU.list)
              {DFU_data <- read.csv(file = paste0(DFU.dir, file,".csv.gz"), header = TRUE, row.names = 1)
               DFU_obj <- CreateSeuratObject(counts = DFU_data,
                                              min.cells=3,
                                              min.features = 200,
                                              project = file)
                                              assign(file, DFU_obj)
               }

sample.list=list(G6,G9,G33,G34,G39,G7,G8,G49,G42,G45,G23,G4,G2,G15)
```
###### Quality Filtering and normalization using SCTransform
```
for (i in 1:length(sample.list)) {
  sample.list[[i]][["percent.mt"]] <- PercentageFeatureSet(sample.list[[i]], pattern = "^MT-")
  sample.list[[i]] <- subset(sample.list[[i]], 
                             subset = nFeature_RNA > 200 & 
                               nFeature_RNA < 5000 & percent.mt < 15 & 
                               nCount_RNA < 25000 & nCount_RNA > 2000)
#SCT transformation
  sample.list[[i]] <- SCTransform(sample.list[[i]],
                                vars.to.regress = "percent.mt", 
                                verbose = FALSE)
                                 }
```
###### Data integration
```
sample.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)
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
```
###### Dimensionality reduction and clustering
```
sample.integrated <- RunPCA(object = sample.integrated, verbose = FALSE) 
sample.integrated = RunUMAP(sample.integrated, dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated)

```
###### tweak-in to assign groups; HDHU and NHDHU
```
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
```
###### Find clusters with defined resolution
```
sample25 <- FindClusters(sample.integrated, resolution = 0.25)
sample25$group <- factor(sample25$group, levels = c("NHDFU", "HDFU"))

```
##### Figure S2A
```
DimPlot(sample25)
#Find all markers
sample25_markers=FindAllMarkers(sample25,only.pos = T,
                                logfc.threshold = 0.30, min.pct = 0.10, recorrect_umi=FALSE)
#write.csv(sample25_markers,"All_markers_sample25.csv")
##### Explore the table and identify cell types using signatures from published articles and PanglaoDb (https://panglaodb.se/))
```
##### Figure S2B
###### Isolate keratinocytes and myeloid cluster and do re-clustering analysis
```
#DefaultAssay(sample25)="integrated"
MK=subset(sample25,subset=seurat_clusters==2|seurat_clusters==5)
MK_re <- RunPCA(object = MK, verbose = FALSE)
MK_re <- RunUMAP(MK_re, dims = 1:30)
MK_re <- FindNeighbors(MK_re)
MK_re20 <- FindClusters(MK_re, resolution = 0.20)
DimPlot(MK_re20)
```
##### Figure S2C
```
VlnPlot(MK_re20, features=c("KRT1", "KRT10", "KRT14", "LYZ"), flip=T,stack=T)
```
##### Figure S2D
###### Identify top 5 signatures for subsets
```
DefaultAssay(MK_re20)="SCT"
MK_re20_markers=FindAllMarkers(MK_re20,only.pos = T,logfc.threshold = 0.30, min.pct = 0.10,recorrect_umi=FALSE)
write.csv(MK_re20_markers,"MK_re20_markers_SCT.csv")

top5 <- MK_re20_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5genes=unique(top5$gene)
DotPlot(MK_re20, features = top10genes) + 
    theme(axis.text.x=element_text(size=8.5)) + 
    theme(axis.text.y=element_text(size=10.0))+
    scale_colour_gradient2(low = "grey", mid = "red", high = "#000000")
```
### Cell-chat analysis
###### Prepare cellchat objects and perform groupwise interactome analysis
```
# Database usage
CellChatDB<- CellChatDB.human
# showDatabaseCategory(CellChatDB)
CellChatDB.use<- CellChatDB

group.list = list("HDFU", "NHDFU")
cellchat.list = c()

for(i in 1:length(group.list)){
  Idents(group.list[[i]])="ID"
  temp=subset(MK_re20,subset=group==group.list[[i]])
  group_cellchat=createCellChat(object=temp, meta=temp@meta.data, group.by="ID")
  group_cellchat <- setIdent(group_cellchat,ident.use="ID")
  group_cellchat@DB <-CellChatDB.use
  group_cellchat <-subsetData(group_cellchat)
  group_cellchat <- identifyOverExpressedGenes(group_cellchat)
  group_cellchat <- identifyOverExpressedInteractions(group_cellchat)
  group_cellchat <- projectData(group_cellchat, PPI.human) #(optional)
  group_cellchat <- computeCommunProb( object = group_cellchat, raw.use = TRUE)
  group_cellchat <- filterCommunication(group_cellchat, min.cells = 10)
  group_cellchat <- computeCommunProbPathway(group_cellchat)
  group_cellchat <- netAnalysis_computeCentrality(group_cellchat, slot.name = "netP")
  group_cellchat <- aggregateNet(group_cellchat)
  cellchat.list <- append(cellchat.list, group_cellchat)
                             } 

# Name the elements of cellchat.list
names(cellchat.list) <- c("HDFU", "NHDFU")  
# saveRDS(cellchat.list,"MK_re20_cellchat.list.rds")
merged_cellchat <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))
saveRDS(merged_cellchat,"MK_re20_merged_cellchat.rds")
```
##### Figure S2E
###### Comparison of cell chat objects
###### Interactions
```
#g1=compareInteractions(merged_cellchat, show.legend = F, group = c(1,2))
g2=compareInteractions(merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
g2 
```
##### Figure S2F
###### subclusters interactome
```
weight.max <- getMaxWeight(cellchat.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchat.list)) {
            netVisual_circle(cellchat.list[[i]]@net$count,
            weight.scale = T, label.edge= F, edge.weight.max = weight.max[2],
            edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cellchat.list)[[i]]))
                                   }
```
###### Relative Information_flow
```
rankNet(merged_cellchat, comparison=c(1,2),stacked = T,do.stat = TRUE)
```
##### Figure S2G-H
```
# Pathway specific analysis 
pathways.show=c("EGF", "IGF")
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(cellcaht.list)){
   for (x in pathways.show){
     netVisual_aggregate(cellcaht.list[[i]], 
                         signaling = x, 
                         layout = "circle")
                           }      }
```
##### Figure S2I-K
```
plotGeneExpression(merged_cellchat, signaling = "EGF", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(merged_cellchat, signaling = "IGF", split.by = "datasets", colors.ggplot = T)
VlnPlot(MK_re20, "TOMM70")
```
### Thank you
