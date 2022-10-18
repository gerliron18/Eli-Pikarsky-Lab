library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(rlang)
library(stringr)
library(patchwork)
library(Matrix)

rawdata <- read.delim("C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/RAW_DATA/arranged_control_nash_genes.csv", sep=",")

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
FinalData <- RunPCA(FinalData,assay = NULL,features = NULL,npcs = 29,
                    rev.pca = FALSE,weight.by.var = TRUE,verbose = TRUE,
                    ndims.print = 1:7,nfeatures.print = 30, approx=FALSE,
                    reduction.name = "pca",reduction.key = "PC_",seed.use = 42)

# Examine and visualize PCA results a few different ways
print(FinalData[["pca"]], dims = 1:7, nfeatures = 5)
VizDimLoadings(FinalData, dims = 1:7, reduction = "pca")
DimPlot(FinalData, reduction = "pca")
DimHeatmap(FinalData, dims = 1:7, cells = 29, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
FinalData <- JackStraw(FinalData, num.replicate = 100)
FinalData <- ScoreJackStraw(FinalData, dims = 1:7)
JackStrawPlot(FinalData, dims = 1:7)
ElbowPlot(FinalData)

# Cluster the cells
FinalData <- FindNeighbors(FinalData,reduction = "pca",dims = 1:7,k.param = 10,
                           prune.SNN = 1/15,nn.method = "annoy",n.trees = 50,
                           annoy.metric = "euclidean",nn.eps = 0,verbose = TRUE)
FinalData <- FindClusters(FinalData,modularity.fxn = 1,resolution = 0.8,
                          method = "matrix",algorithm = 1,n.start = 10,
                          n.iter = 10,random.seed = 0,group.singletons = TRUE,
                          verbose = TRUE)

head(Idents(FinalData), 29)
RidgePlot(FinalData, features = top10, ncol = 2)

# Run non-linear dimensional reduction
FinalData <- RunUMAP(FinalData,reduction = "pca",
                     assay = DefaultAssay(object = FinalData),slot = "data",
                     umap.method = "uwot",n.neighbors = 10L,n.components = 2L,
                     metric = "cosine",learning.rate = 1,min.dist = 0.3,
                     spread = 1,set.op.mix.ratio = 1,local.connectivity = 1L,
                     repulsion.strength = 1,negative.sample.rate = 5L,
                     seed.use = 42L,dens.lambda = 2,dens.frac = 0.3,
                     dens.var.shift = 0.1,verbose = TRUE,
                     reduction.name = "umap",reduction.key = "UMAP_", dims=1:4)

DimPlot(FinalData, reduction = "umap", label = TRUE, pt.size = 2.5)
umap <- DimPlot(FinalData, reduction = "umap", label = TRUE, pt.size = 2.5)
LabelPoints(plot = umap, points = cells_vector, repel = TRUE)

cells_vector <- c("GSM1178970", "GSM1178971", "GSM1178972", "GSM1178973", 
                  "GSM1178974", "GSM1178975", "GSM1178977", "GSM1178978", 
                  "GSM1178979", "GSM1178980", "GSM1178987", "GSM1178995", 
                  "GSM1178996", "GSM1178998", "GSM1179001", "GSM1179002", 
                  "GSM1179003", "GSM1179004", "GSM1179006", "GSM1179008", 
                  "GSM1179010", "GSM1179015", "GSM1179017", "GSM1179018", 
                  "GSM1179024", "GSM1179026", "GSM1179033", "GSM1179035", 
                  "GSM1179036")

colors_vector = c('GSM1178970' = 'cyan4', 'GSM1178971' = 'cyan4', 'GSM1178972' = 'cyan4', 
                  'GSM1178973' = 'cyan4', 'GSM1178974' = 'cyan4', 'GSM1178975' = 'salmon1',
                  'GSM1178976' = 'gold', 'GSM1178977' = 'cyan4', 'GSM1178978' = 'cyan4',
                  'GSM1178979' = 'cyan4', 'GSM1178980' = 'salmon1', 'GSM1178981' = 'gold',
                  'GSM1178982' = 'gold', 'GSM1178983' = 'gold', 'GSM1178984' = 'gold', 
                  'GSM1178985' = 'gold', 'GSM1178986' = 'hotpink3', 'GSM1178987' = 'salmon1',
                  'GSM1178988' = 'hotpink3', 'GSM1178989' = 'hotpink3', 'GSM1178990' = 'gold',
                  'GSM1178991' = 'gold', 'GSM1178992' = 'gold', 'GSM1178993' = 'hotpink3',
                  'GSM1178994' = 'gold', 'GSM1178995' = 'salmon1', 'GSM1178996' = 'salmon1',
                  'GSM1178998' = 'cyan4', 'GSM1178999' = 'hotpink3', 'GSM1179001' = 'salmon1',
                  'GSM1179002' = 'salmon1', 'GSM1179003' = 'salmon1', 'GSM1179004' = 'salmon1', 
                  'GSM1179005' = 'gold', 'GSM1179006' = 'salmon1', 'GSM1179007' = 'gold',
                  'GSM1179008' = 'salmon1', 'GSM1179010' = 'cyan4', 'GSM1179012' = 'gold',
                  'GSM1179015' = 'salmon1', 'GSM1179016' = 'gold', 'GSM1179017' = 'salmon1',
                  'GSM1179018' = 'cyan4', 'GSM1179021' = 'hotpink3', 'GSM1179024' = 'cyan4',
                  'GSM1179025' = 'hotpink3', 'GSM1179026' = 'salmon1', 'GSM1179027' = 'hotpink3',
                  'GSM1179030' = 'gold', 'GSM1179033' = 'salmon1', 'GSM1179035' = 'salmon1',
                  'GSM1179036' = 'salmon1', 'GSM1179037' = 'hotpink3', 'GSM1179038' = 'gold')

type_colors_vector <- c("Control" = 'cyan4', "Nash" = 'salmon1', "Healthy obese" = 'gold', 
                        "Steatosis" = 'hotpink3')

LabelPoints(plot = umap, points = cells_vector, repel = TRUE)
DimPlot(FinalData, reduction = "pca", cells = cells_vector, 
        cols = type_colors_vector, pt.size = 2.5)

# Finding differentially expressed features (cluster biomarkers)
FinalData.markers <- FindAllMarkers(FinalData, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
FinalData.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# AKR1B10
VlnPlot(FinalData, features = c("ENSG00000198074"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000198074"), pt.size = 2.5)

# EGR1
VlnPlot(FinalData, features = c("ENSG00000120738"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000120738"), pt.size = 2.5)
FeaturePlot(FinalData, features = c("ENSG00000120738"), pt.size = 2.5, min.cutoff = 0.5, max.cutoff = 0.8)

# TP53
VlnPlot(FinalData, features = c("ENSG00000141510"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000141510"), pt.size = 2.5)

# RARa
VlnPlot(FinalData, features = c("ENSG00000131759"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000131759"), pt.size = 2.5)

# RARb
VlnPlot(FinalData, features = c("ENSG00000077092"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000077092"), pt.size = 2.5)

# RARg
VlnPlot(FinalData, features = c("ENSG00000172819"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000172819"), pt.size = 2.5)

# RXRa
VlnPlot(FinalData, features = c("ENSG00000186350"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000186350"), pt.size = 2.5)

# RXRb
VlnPlot(FinalData, features = c("ENSG00000204231"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000204231"), pt.size = 2.5)

# RARRES2
VlnPlot(FinalData, features = c("ENSG00000106538"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000106538"), pt.size = 2.5)

# ALL
VlnPlot(FinalData, features = c("ENSG00000198074", "ENSG00000120738", 
                                "ENSG00000141510", "ENSG00000131759", 
                                "ENSG00000077092", "ENSG00000172819",
                                "ENSG00000186350", "ENSG00000204231",
                                "ENSG00000204231"), 
        slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000198074", "ENSG00000120738", 
                                    "ENSG00000141510", "ENSG00000131759", 
                                    "ENSG00000077092", "ENSG00000172819",
                                    "ENSG00000186350", "ENSG00000204231",
                                    "ENSG00000204231"))

# According to PCA
# MOGAT1
VlnPlot(FinalData, features = c("ENSG00000124003"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000124003"), pt.size = 2.5)

# RNA5SP218
VlnPlot(FinalData, features = c("ENSG00000200058"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000200058"), pt.size = 2.5)

# ANXA2
VlnPlot(FinalData, features = c("ENSG00000182718"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000182718"), pt.size = 2.5)

# RND3
VlnPlot(FinalData, features = c("ENSG00000115963"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("ENSG00000115963"), pt.size = 2.5)
