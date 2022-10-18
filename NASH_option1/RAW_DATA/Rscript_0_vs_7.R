library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(rlang)
library(stringr)
library(patchwork)
library(Matrix)

rawdata <- read.delim("C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NASH_option1/RAW_DATA/GSE162694_raw_counts_0_vs_7.csv", sep=",")

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
FinalData <- RunPCA(FinalData,assay = NULL,features = NULL,npcs = 50,
                    rev.pca = FALSE,weight.by.var = TRUE,verbose = TRUE,
                    ndims.print = 1:5,nfeatures.print = 30, approx=FALSE,
                    reduction.name = "pca",reduction.key = "PC_",seed.use = 42)

# Examine and visualize PCA results a few different ways
print(FinalData[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(FinalData, dims = 1:5, reduction = "pca")
DimPlot(FinalData, reduction = "pca")
DimHeatmap(FinalData, dims = 1:5, cells = 143, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
FinalData <- JackStraw(FinalData, num.replicate = 100)
FinalData <- ScoreJackStraw(FinalData, dims = 1:5)
JackStrawPlot(FinalData, dims = 1:5)
ElbowPlot(FinalData)

# Cluster the cells
FinalData <- FindNeighbors(FinalData,reduction = "pca",dims = 1:10,k.param = 10,
                           prune.SNN = 1/15,nn.method = "annoy",n.trees = 50,
                           annoy.metric = "euclidean",nn.eps = 0,verbose = TRUE)
FinalData <- FindClusters(FinalData,modularity.fxn = 1,resolution = 0.5,
                          method = "matrix",algorithm = 1,n.start = 10,
                          n.iter = 10,random.seed = 0,group.singletons = TRUE,
                          verbose = TRUE)

head(Idents(FinalData), 143)
RidgePlot(FinalData, features = top10, ncol = 2)

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

cells_vector = c('548nash100' = '0','548nash101' = '7','548nash103' = '0',
                 '548nash108' = '0','548nash109' = '0','548nash110' = '6',
                 '548nash112' = '0','548nash1126' = '7','548nash1128' = '7',
                 '548nash1129' = '7','548nash1132' = '6','548nash1133' = '6',
                 '548nash1135' = '7','548nash1136' = '7','548nash1139' = '7',
                 '548nash114' = '0','548nash116' = '6','548nash117' = '6',
                 '548nash12' = '6','548nash24' = '1','548nash26' = '0',
                 '548nash28' = '6','548nash3' = '0','548nash30' = '7',
                 '548nash32' = '0','548nash33' = '1','548nash34' = '7',
                 '548nash36' = '0','548nash38' = '0','548nash39' = '1',
                 '548nash40' = '0','548nash41' = '6','548nash45' = '0',
                 '548nash49' = '0','548nash51' = '1','548nash52' = '1',
                 '548nash53' = '1','548nash54' = '0','548nash57' = '0',
                 '548nash59' = '0','548nash6' = '6','548nash60' = '0',
                 '548nash62' = '0','548nash64' = '0','548nash65' = '0',
                 '548nash66' = '1','548nash69' = '0','548nash7' ='0',
                 '548nash71' ='0','548nash72' = '0','548nash73' = '0',
                 '548nash75' = '0','548nash76' = '0','548nash77' = '1',
                 '548nash78' = '1','548nash79' = '0','548nash85' = '1',
                 '548nash88' = '0','548nash89' = '6','548nash9' = '1',
                 '548nash90' = '0','548nash91' = '6','548nash93' = '1',
                 '548nash95' = '6','548nash99' = '0')


#Idents(object = FinalData) <- cells_vector
#head(x = Idents(object = FinalData))


FinalData[["Identity"]]<- cells_vector

VlnPlot(object = FinalData, features = 'ENSG00000120738', group.by = 'Identity')

# ALL
VlnPlot(FinalData, features = top10, slot = "counts", log = TRUE, group.by = 'Identity')
FeaturePlot(FinalData, features = top10)

# EGR1
VlnPlot(FinalData, features = c("ENSG00000120738"), slot = "counts", log = TRUE)
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

