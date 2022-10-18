library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(rlang)
library(stringr)
library(patchwork)
library(Matrix)
library(writexl)


rawdata <- read.delim("C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NAFLD_option2/RAW_DATA/arranged_before_surgery_P53.csv", sep=",")

rownames(rawdata) = make.names(rawdata[,1], unique=TRUE)
rawdata <- rawdata[,-1]

SingleCellSeq <- CreateSeuratObject(counts = rawdata,project = "SingleCellSeq", assay = "RNA", min.cells = 3, min.features = 200)

#Pre-Processing
SingleCellSeq[["percent.mt"]] <- PercentageFeatureSet(object = SingleCellSeq, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(SingleCellSeq, features = c("nCount_RNA"), ncol = 3)

plot1 <- FeatureScatter(SingleCellSeq, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot2 <- FeatureScatter(SingleCellSeq, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
plot1 + plot2

FinalData <- NormalizeData(SingleCellSeq, normalization.method = "LogNormalize", scale.factor = 10000)

FinalData <- FindVariableFeatures(FinalData, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(FinalData), 10)
print(top10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(FinalData)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot3 + plot4

# Scaling the data
options(ggrepel.max.overlaps = Inf)
all.genes <- rownames(FinalData)
FinalData <- ScaleData(FinalData, features = all.genes)

# Perform linear dimensional reduction
FinalData <- RunPCA(FinalData,assay = NULL,features = NULL,npcs = 50,
                    rev.pca = FALSE,weight.by.var = TRUE,verbose = TRUE,
                    ndims.print = 1:5,nfeatures.print = 30,reduction.name = "pca",
                    reduction.key = "PC_",seed.use = 42)

# Examine and visualize PCA results a few different ways
print(FinalData[["pca"]], dims = 1:7, nfeatures = 5)
VizDimLoadings(FinalData, dims = 1:7, reduction = "pca")
DimPlot(FinalData, reduction = "pca")
DimHeatmap(FinalData, dims = 1:7, cells = 54, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
FinalData <- JackStraw(FinalData, num.replicate = 100)
FinalData <- ScoreJackStraw(FinalData, dims = 1:7)
JackStrawPlot(FinalData, dims = 1:7)
ElbowPlot(FinalData)

# Cluster the cells
FinalData <- FindNeighbors(FinalData,reduction = "pca",dims = 1:7,k.param = 20,
                           prune.SNN = 1/15,nn.method = "annoy",n.trees = 50,
                           annoy.metric = "euclidean",nn.eps = 0,verbose = TRUE)
FinalData <- FindClusters(FinalData,modularity.fxn = 1,resolution = 1.0,
                          method = "matrix",algorithm = 1,n.start = 10,
                          n.iter = 10,random.seed = 0,group.singletons = TRUE,
                          verbose = TRUE)

head(Idents(FinalData), 54)
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
                     reduction.name = "umap",reduction.key = "UMAP_", dims=1:3)

DimPlot(FinalData, reduction = "umap", label = TRUE, pt.size = 2.5)
umap <- DimPlot(FinalData, reduction = "umap", label = TRUE, pt.size = 0.5)
LabelPoints(plot = umap, points = TopCells(object = FinalData), repel = TRUE)

fibrosis_vector <- c("GSM1178970" = "0", "GSM1178972" = "0",
                     "GSM1178973" = "0", "GSM1178974" = "0", "GSM1178975" = "1",
                     "GSM1178976" = "1", "GSM1178977" = "0", "GSM1178978" = "0",
                     "GSM1178979" = "0", "GSM1178981" = "0",
                     "GSM1178982" = "0", "GSM1178983" = "0", "GSM1178984" = "1",
                     "GSM1178985" = "0", "GSM1178986" = "0", "GSM1178987" = "0",
                     "GSM1178988" = "0", "GSM1178989" = "1", "GSM1178990" = "0",
                     "GSM1178991" = "0", "GSM1178992" = "0", "GSM1178993" = "0",
                     "GSM1178995" = "0", "GSM1178996" = "1",
                     "GSM1178999" = "1", "GSM1179001" = "3",
                     "GSM1179003" = "3", "GSM1179004" = "0",
                     "GSM1179005" = "0", "GSM1179006" = "1", "GSM1179007" = "0",
                     "GSM1179010" = "1", "GSM1179012" = "0",
                     "GSM1179015" = "1", "GSM1179016" = "0", "GSM1179017" = "1",
                     "GSM1179018" = "1", "GSM1179021" = "1", "GSM1179024" = "0",
                     "GSM1179025" = "0", "GSM1179026" = "1", "GSM1179027" = "0",
                     "GSM1179030" = "0", "GSM1179033" = "1", "GSM1179035" = "1",
                     "GSM1179036" = "1", "GSM1179037" = "0", "GSM1179038" = "0")

nas_vector <- c("GSM1178970" = "0", "GSM1178971" = "1", "GSM1178972" = "0",
                "GSM1178973" = "0", "GSM1178974" = "0", "GSM1178975" = "4",
                "GSM1178976" = "0", "GSM1178977" = "0", "GSM1178978" = "0",
                "GSM1178979" = "0", "GSM1178980" = "5", "GSM1178981" = "0",
                "GSM1178982" = "0", "GSM1178983" = "0", "GSM1178984" = "1",
                "GSM1178985" = "0", "GSM1178986" = "3", "GSM1178987" = "6",
                "GSM1178988" = "1", "GSM1178989" = "1", "GSM1178990" = "2",
                "GSM1178991" = "0", "GSM1178992" = "0", "GSM1178993" = "1",
                "GSM1178994" = "0", "GSM1178995" = "5", "GSM1178996" = "6",
                "GSM1178998" = "0", "GSM1178999" = "1", "GSM1179001" = "5",
                "GSM1179002" = "4", "GSM1179003" = "5", "GSM1179004" = "7",
                "GSM1179005" = "0", "GSM1179006" = "5", "GSM1179007" = "0",
                "GSM1179008" = "5", "GSM1179010" = "0", "GSM1179012" = "0",
                "GSM1179015" = "5", "GSM1179016" = "0", "GSM1179017" = "5",
                "GSM1179018" = "0", "GSM1179021" = "3", "GSM1179024" = "0",
                "GSM1179025" = "1", "GSM1179026" = "5", "GSM1179027" = "3",
                "GSM1179030" = "0", "GSM1179033" = "5", "GSM1179035" = "6",
                "GSM1179036" = "5", "GSM1179037" = "1", "GSM1179038" = "1")

age_vector <- c("GSM1178970" = "50", "GSM1178971" = "50", "GSM1178972" = "70",
                "GSM1178973" = "20", "GSM1178974" = "80", "GSM1178975" = "40",
                "GSM1178976" = "40", "GSM1178977" = "60", "GSM1178978" = "40",
                "GSM1178979" = "40", "GSM1178980" = "30", "GSM1178981" = "30",
                "GSM1178982" = "50", "GSM1178983" = "40", "GSM1178984" = "40",
                "GSM1178985" = "40", "GSM1178986" = "40", "GSM1178987" = "50",
                "GSM1178988" = "20", "GSM1178989" = "30", "GSM1178990" = "50",
                "GSM1178991" = "40", "GSM1178992" = "60", "GSM1178993" = "30",
                "GSM1178994" = "40", "GSM1178995" = "30", "GSM1178996" = "50",
                "GSM1178998" = "20", "GSM1178999" = "40", "GSM1179001" = "40",
                "GSM1179002" = "40", "GSM1179003" = "40", "GSM1179004" = "30",
                "GSM1179005" = "40", "GSM1179006" = "40", "GSM1179007" = "30",
                "GSM1179008" = "40", "GSM1179010" = "40", "GSM1179012" = "50",
                "GSM1179015" = "30", "GSM1179016" = "40", "GSM1179017" = "50",
                "GSM1179018" = "70", "GSM1179021" = "30", "GSM1179024" = "40",
                "GSM1179025" = "30", "GSM1179026" = "50", "GSM1179027" = "40",
                "GSM1179030" = "30", "GSM1179033" = "30", "GSM1179035" = "40",
                "GSM1179036" = "50", "GSM1179037" = "30", "GSM1179038" = "50")


FinalData[["fibrosis_Identity"]]<- fibrosis_vector
FinalData[["nas_Identity"]]<- nas_vector
FinalData[["age_Identity"]]<- age_vector

# Finding differentially expressed features (cluster biomarkers)
FinalData.markers <- FindAllMarkers(FinalData, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
markers <- FinalData.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write_xlsx(markers,"C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NAFLD_option2/PLOTS/genes/P53/ClusterMarkers.xlsx")

# find all markers distinguishing clusters from other clusters
cluster6_11.markers <- FindMarkers(FinalData, ident.1 = c(6, 11), ident.2 = c(0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14), min.pct = 0.25)
markers6_11 <- cluster6_11.markers %>%
  slice_max(n = 20, order_by = avg_log2FC)

markers6_11 <- cbind(" "=rownames(markers6_11), markers6_11)
write_xlsx(markers6_11,"C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NAFLD_option2/PLOTS/genes/P53/ClusterMarkers.xlsx")

##############################################################################
# TP53
VlnPlot(object = FinalData, features = 'TP53', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'TP53', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'TP53', group.by = 'age_Identity')

# TP53
VlnPlot(FinalData, features = c("TP53"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("TP53"), pt.size = 1.0)
FeaturePlot(FinalData, features = c("TP53"), pt.size = 1.0, min.cutoff = 0.5, max.cutoff = 0.8)
