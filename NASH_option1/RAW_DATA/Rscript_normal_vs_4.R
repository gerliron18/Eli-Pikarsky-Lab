library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(rlang)
library(stringr)
library(patchwork)
library(Matrix)

rawdata <- read.delim("C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NASH_option1/RAW_DATA/GSE162694_raw_counts_normal_vs_4.csv", sep=",")

#rownames(rawdata) <- rawdata[,1]
rownames(rawdata) = make.names(rawdata[,1], unique=TRUE)
rawdata <- rawdata[,-1]

#converted_names <- map.ids(rawdata, org=hsa, from="ENSEMBL", to="ENTREZID" )

SingleCellSeq <- CreateSeuratObject(counts = rawdata,project = "SingleCellSeq", assay = "RNA", min.cells = 3, min.features = 200)

SingleCellSeq[["percent.mt"]] <- PercentageFeatureSet(object = SingleCellSeq, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(SingleCellSeq, features = c("nCount_RNA"), ncol = 3)

FinalData <- NormalizeData(SingleCellSeq, normalization.method = "LogNormalize", scale.factor = 10000)

FinalData <- FindVariableFeatures(FinalData, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(FinalData), 10)
print(top10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(FinalData)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2

# Scaling the data
all.genes <- rownames(FinalData)
FinalData <- ScaleData(FinalData, features = all.genes)

# Perform linear dimensional reduction
FinalData <- RunPCA(FinalData,assay = NULL,features = NULL,npcs = 30,
                    rev.pca = FALSE,weight.by.var = TRUE,verbose = TRUE,
                    ndims.print = 1:5,nfeatures.print = 30, approx=FALSE,
                    reduction.name = "pca",reduction.key = "PC_",seed.use = 42)

# Examine and visualize PCA results a few different ways
print(FinalData[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(FinalData, dims = 1:5, reduction = "pca")
DimPlot(FinalData, reduction = "pca")
DimHeatmap(FinalData, dims = 1:5, cells = 42, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
FinalData <- JackStraw(FinalData, num.replicate = 100)
FinalData <- ScoreJackStraw(FinalData, dims = 1:5)
JackStrawPlot(FinalData, dims = 1:5)
ElbowPlot(FinalData)

# Cluster the cells
FinalData <- FindNeighbors(FinalData,reduction = "pca",dims = 1:10,k.param = 20,
                           prune.SNN = 1/15,nn.method = "annoy",n.trees = 50,
                           annoy.metric = "euclidean",nn.eps = 0,verbose = TRUE)
FinalData <- FindClusters(FinalData,modularity.fxn = 1,resolution = 0.8,
                          method = "matrix",algorithm = 1,n.start = 10,
                          n.iter = 10,random.seed = 0,group.singletons = TRUE,
                          verbose = TRUE)

head(Idents(FinalData), 42)
RidgePlot(FinalData, features = top10, ncol = 2)
RidgePlot(FinalData, features = "ENSG00000120738", ncol = 2)

# Run non-linear dimensional reduction
FinalData <- RunUMAP(FinalData,reduction = "pca",
                     assay = DefaultAssay(object = FinalData),slot = "data",
                     umap.method = "uwot",n.neighbors = 30L,n.components = 2L,
                     metric = "cosine",learning.rate = 1,min.dist = 0.3,
                     spread = 1,set.op.mix.ratio = 1,local.connectivity = 1L,
                     repulsion.strength = 1,negative.sample.rate = 5L,
                     seed.use = 42L,dens.lambda = 2,dens.frac = 0.3,
                     dens.var.shift = 0.1,verbose = TRUE,
                     reduction.name = "umap",reduction.key = "UMAP_", dims=1:5)

DimPlot(FinalData, reduction = "umap", label = TRUE, pt.size = 2.5)
umap <- DimPlot(FinalData, reduction = "umap", label = TRUE, pt.size = 2.5)
LabelPoints(plot = umap, points = TopCells(object = FinalData), repel = TRUE)

#LabelPoints(plot = umap, points = cells_vector, repel = TRUE)
#DimPlot(FinalData, reduction = "pca", cells = cells_vector, 
#        cols = type_colors_vector, pt.size = 2.5)

# Finding differentially expressed features (cluster biomarkers)
FinalData.markers <- FindAllMarkers(FinalData, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
FinalData.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

cells_vector = c('548nash100' = 'normal','548nash103' = 'normal','548nash108' = 'normal',
                 '548nash109' = 'normal','548nash1119' = 'stage 4','548nash112' = 'normal',
                 '548nash1120' = 'stage 4','548nash1121' = 'stage 4','548nash1122' = 'stage 4',
                 '548nash1123' = 'stage 4','548nash1124' = 'stage 4','548nash1125' = 'stage 4',
                 '548nash1131' = 'stage 4','548nash114' = 'normal','548nash116' = 'stage 4',
                 '548nash144' = 'stage 4','548nash26' = 'normal','548nash3' = 'stage 4',
                 '548nash32' = 'normal','548nash36' = 'normal','548nash38' = 'normal',
                 '548nash45' = 'normal','548nash49' = 'normal','548nash54' = 'normal',
                 '548nash57' = 'normal','548nash59' = 'normal','548nash60' = 'normal',
                 '548nash62' = 'normal','548nash64' = 'normal','548nash65' = 'normal',
                 '548nash69' = 'normal','548nash7' = 'normal','548nash71' = 'normal',
                 '548nash72' = 'normal','548nash73' = 'normal','548nash75' = 'normal',
                 '548nash76' = 'normal','548nash79' = 'normal','548nash86' = 'stage 4',
                 '548nash88' = 'normal','548nash90' = 'normal','548nash99' = 'normal')


#Idents(object = FinalData) <- cells_vector
#head(x = Idents(object = FinalData))


FinalData[["Identity"]]<- cells_vector

# ALL
VlnPlot(FinalData, features = top10, slot = "counts", log = TRUE, group.by = 'Identity')
FeaturePlot(FinalData, features = top10)

# EGR1
VlnPlot(FinalData, features = c("ENSG00000120738"), slot = "counts", log = TRUE, group.by = 'Identity')
FeaturePlot(FinalData, features = c("ENSG00000120738"), pt.size = 2.5)

# RNU1-1
VlnPlot(FinalData, features = c("ENSG00000206652"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000206652"), pt.size = 2.5)

# RNU1
VlnPlot(FinalData, features = c("ENSG00000207418"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000207418"), pt.size = 2.5)

# PCSK5
VlnPlot(FinalData, features = c("ENSG00000099139"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000099139"), pt.size = 2.5)

# MT-ND5
VlnPlot(FinalData, features = c("ENSG00000198786"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000198786"), pt.size = 2.5)

# C7
VlnPlot(FinalData, features = c("ENSG00000112936"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000112936"), pt.size = 2.5)

# THBS1
VlnPlot(FinalData, features = c("ENSG00000137801"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000137801"), pt.size = 2.5)

# FASN
VlnPlot(FinalData, features = c("ENSG00000169710"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000169710"), pt.size = 2.5)

# TP53
VlnPlot(FinalData, features = c("ENSG00000272060"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000272060"), pt.size = 2.5)

