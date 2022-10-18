library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(rlang)
library(stringr)
library(patchwork)
library(Matrix)
library(writexl)

rawdata <- read.delim("C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NASH_option1/RAW_DATA/EGR1_targets/GSE162694_raw_counts_EGR1.csv", sep=",")

rownames(rawdata) = make.names(rawdata[,1], unique=TRUE)
rawdata <- rawdata[,-1]

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
                    ndims.print = 1:7,nfeatures.print = 30, approx=FALSE,
                    reduction.name = "pca",reduction.key = "PC_",seed.use = 42)

# Examine and visualize PCA results a few different ways
print(FinalData[["pca"]], dims = 1:7, nfeatures = 5)
VizDimLoadings(FinalData, dims = 1:7, reduction = "pca")
DimPlot(FinalData, reduction = "pca")
DimHeatmap(FinalData, dims = 1:7, cells = 143, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
FinalData <- JackStraw(FinalData, num.replicate = 100)
FinalData <- ScoreJackStraw(FinalData, dims = 1:7)
JackStrawPlot(FinalData, dims = 1:7)
ElbowPlot(FinalData)

# Cluster the cells
FinalData <- FindNeighbors(FinalData,reduction = "pca",dims = 1:7,k.param = 10,
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
                     reduction.name = "umap",reduction.key = "UMAP_", dims=1:4)

DimPlot(FinalData, reduction = "umap", label = TRUE, pt.size = 2.5)
umap <- DimPlot(FinalData, reduction = "umap", label = TRUE, pt.size = 2.5)
LabelPoints(plot = umap, points = TopCells(object = FinalData), repel = TRUE)

fibrosis_vector <- c("548nash1" = "0", "548nash10" = "1", "548nash100" = "normal",
                     "548nash101" = "3", "548nash102" = "1", "548nash103" = "normal",
                     "548nash104" = "0", "548nash105" = "2", "548nash106" = "1",
                     "548nash108" = "normal", "548nash109" = "normal", "548nash11" = "1",
                     "548nash110" = "2", "548nash111" = "2", "548nash1119" = "4",
                     "548nash112" = "normal", "548nash1120" = "4", "548nash1121" = "4",
                     "548nash1122" = "4", "548nash1123" = "4", "548nash1124" = "4",
                     "548nash1125" = "4", "548nash1126" = "2", "548nash1127" = "2",
                     "548nash1128" = "3", "548nash1129" = "2", "548nash113" = "1",
                     "548nash1130" = "2", "548nash1131" = "4", "548nash1132" = "2",
                     "548nash1133" = "2", "548nash1134" = "2", "548nash1135" = "2",
                     "548nash1136" = "2", "548nash1137" = "3", "548nash1138" = "3",
                     "548nash1139" = "2", "548nash114" = "normal", "548nash115" = "1",
                     "548nash116" = "4", "548nash117" = "1", "548nash12" = "2",
                     "548nash13" = "0", "548nash14" = "0", "548nash143" = "1",
                     "548nash144" = "4", "548nash145" = "1", "548nash146" = "0",
                     "548nash148" = "0", "548nash149" = "1", "548nash15" = "0",
                     "548nash150" = "0", "548nash151" = "0", "548nash153" = "3",
                     "548nash154" = "0", "548nash155" = "1", "548nash156" = "1",
                     "548nash157" = "0", "548nash158" = "0", "548nash159" = "0",
                     "548nash16" = "1", "548nash160" = "0", "548nash17" = "1",
                     "548nash19" = "2", "548nash2" = "2", "548nash21" = "0",
                     "548nash22" = "1", "548nash23" = "0", "548nash24" = "0",
                     "548nash25" = "2", "548nash26" = "normal", "548nash27" = "1",
                     "548nash28" = "1", "548nash29" = "0", "548nash3" = "4",
                     "548nash30" = "2", "548nash31" = "0", "548nash32" = "normal",
                     "548nash33" = "0", "548nash34" = "2", "548nash35" = "0",
                     "548nash36" = "normal", "548nash37" = "3", "548nash38" = "normal",
                     "548nash39" = "0", "548nash40" = "normal", "548nash41" = "3",
                     "548nash44" = "0", "548nash45" = "normal", "548nash46" = "0",
                     "548nash47" = "1", "548nash48" = "1", "548nash49" = "normal",
                     "548nash5" = "0", "548nash51" = "2", "548nash52" = "0",
                     "548nash53" = "1", "548nash54" = "normal", "548nash55" = "1",
                     "548nash56" = "1", "548nash57" = "normal", "548nash58" = "0",
                     "548nash59" = "normal", "548nash6" = "1", "548nash60" = "normal",
                     "548nash61" = "1", "548nash62" = "normal", "548nash63" = "0",
                     "548nash64" = "normal", "548nash65" = "normal", "548nash66" = "1",
                     "548nash67" = "2", "548nash68" = "2", "548nash69" = "normal",
                     "548nash7" = "normal", "548nash70" = "0", "548nash71" = "normal",
                     "548nash72" = "normal", "548nash73" = "normal", "548nash75" = "normal",
                     "548nash76" = "normal", "548nash77" = "0", "548nash78" = "1",
                     "548nash79" = "normal", "548nash80" = "2", "548nash81" = "2",
                     "548nash82" = "1", "548nash83" = "2", "548nash84" = "0",
                     "548nash85" = "0", "548nash86" = "4", "548nash87" = "0",
                     "548nash88" = "normal", "548nash89" = "1", "548nash9" = "1",
                     "548nash90" = "normal", "548nash91" = "2", "548nash92" = "3",
                     "548nash93" = "0", "548nash95" = "2", "548nash96" = "0",
                     "548nash98" = "1", "548nash99" = "normal")

nas_vector <- c("548nash1" = "3", "548nash10" = "3", "548nash100" = "0",
                "548nash101" = "7", "548nash102" = "5", "548nash103" = "0",
                "548nash104" = "4", "548nash105" = "5", "548nash106" = "4",
                "548nash108" = "0", "548nash109" = "0", "548nash11" = "5",
                "548nash110" = "6", "548nash111" = "5", "548nash1119" = "NA",
                "548nash112" = "0", "548nash1120" = "NA", "548nash1121" = "NA",
                "548nash1122" = "NA", "548nash1123" = "NA", "548nash1124" = "NA",
                "548nash1125" = "NA", "548nash1126" = "7", "548nash1127" = "5",
                "548nash1128" = "7", "548nash1129" = "7", "548nash113" = "4",
                "548nash1130" = "5", "548nash1131" = "4", "548nash1132" = "6",
                "548nash1133" = "6", "548nash1134" = "NA", "548nash1135" = "7",
                "548nash1136" = "8", "548nash1137" = "NA", "548nash1138" = "NA",
                "548nash1139" = "7", "548nash114" = "0", "548nash115" = "3",
                "548nash116" = "6", "548nash117" = "6", "548nash12" = "6",
                "548nash13" = "2", "548nash14" = "5", "548nash143" = "NA",
                "548nash144" = "NA", "548nash145" = "NA", "548nash146" = "NA",
                "548nash148" = "NA", "548nash149" = "NA", "548nash15" = "3",
                "548nash150" = "NA", "548nash151" = "NA", "548nash153" = "NA",
                "548nash154" = "NA", "548nash155" = "NA", "548nash156" = "NA",
                "548nash157" = "NA", "548nash158" = "NA", "548nash159" = "NA",
                "548nash16" = "5", "548nash160" = "NA", "548nash17" = "4",
                "548nash19" = "4", "548nash2" = "4", "548nash21" = "3",
                "548nash22" = "5", "548nash23" = "2", "548nash24" = "1",
                "548nash25" = "5", "548nash26" = "0", "548nash27" = "5",
                "548nash28" = "6", "548nash29" = "3", "548nash3" = "0",
                "548nash30" = "7", "548nash31" = "4", "548nash32" = "0",
                "548nash33" = "1", "548nash34" = "7", "548nash35" = "3",
                "548nash36" = "0", "548nash37" = "5", "548nash38" = "0",
                "548nash39" = "1", "548nash40" = "0", "548nash41" = "6",
                "548nash44" = "2", "548nash45" = "0", "548nash46" = "2",
                "548nash47" = "2", "548nash48" = "2", "548nash49" = "0",
                "548nash5" = "3", "548nash51" = "1", "548nash52" = "1",
                "548nash53" = "1", "548nash54" = "0", "548nash55" = "4",
                "548nash56" = "4", "548nash57" = "0", "548nash58" = "3",
                "548nash59" = "0", "548nash6" = "6", "548nash60" = "0",
                "548nash61" = "5", "548nash62" = "0", "548nash63" = "5",
                "548nash64" = "0", "548nash65" = "0", "548nash66" = "1",
                "548nash67" = "3", "548nash68" = "4", "548nash69" = "0",
                "548nash7" = "0", "548nash70" = "2", "548nash71" = "0",
                "548nash72" = "0", "548nash73" = "0", "548nash75" = "0",
                "548nash76" = "0", "548nash77" = "1", "548nash78" = "1",
                "548nash79" = "0", "548nash80" = "3", "548nash81" = "5",
                "548nash82" = "5", "548nash83" = "5", "548nash84" = "2",
                "548nash85" = "1", "548nash86" = "4", "548nash87" = "2",
                "548nash88" = "0", "548nash89" = "6", "548nash9" = "1",
                "548nash90" = "0", "548nash91" = "6", "548nash92" = "4",
                "548nash93" = "1", "548nash95" = "6", "548nash96" = "5",
                "548nash98" = "5", "548nash99" = "0")

age_vector <- c("548nash1" = "30", "548nash10" = "50", "548nash100" = "20",
                "548nash101" = "40", "548nash102" = "40", "548nash103" = "50",
                "548nash104" = "60", "548nash105" = "40", "548nash106" = "20",
                "548nash108" = "50", "548nash109" = "30", "548nash11" = "20",
                "548nash110" = "50", "548nash111" = "60", "548nash1119" = "50",
                "548nash112" = "30", "548nash1120" = "60", "548nash1121" = "60",
                "548nash1122" = "50", "548nash1123" = "50", "548nash1124" = "60",
                "548nash1125" = "40", "548nash1126" = "40", "548nash1127" = "40",
                "548nash1128" = "50", "548nash1129" = "10", "548nash113" = "30",
                "548nash1130" = "50", "548nash1131" = "60", "548nash1132" = "50",
                "548nash1133" = "60", "548nash1134" = "50", "548nash1135" = "30",
                "548nash1136" = "50", "548nash1137" = "60", "548nash1138" = "60",
                "548nash1139" = "50", "548nash114" = "50", "548nash115" = "50",
                "548nash116" = "60", "548nash117" = "20", "548nash12" = "30",
                "548nash13" = "40", "548nash14" = "30", "548nash143" = "30",
                "548nash144" = "50", "548nash145" = "40", "548nash146" = "50",
                "548nash148" = "40", "548nash149" = "20", "548nash15" = "30",
                "548nash150" = "50", "548nash151" = "40", "548nash153" = "50",
                "548nash154" = "50", "548nash155" = "20", "548nash156" = "40",
                "548nash157" = "50", "548nash158" = "20", "548nash159" = "50",
                "548nash16" = "30", "548nash160" = "50", "548nash17" = "60",
                "548nash19" = "40", "548nash2" = "50", "548nash21" = "60",
                "548nash22" = "50", "548nash23" = "50", "548nash24" = "20",
                "548nash25" = "20", "548nash26" = "60", "548nash27" = "40",
                "548nash28" = "60", "548nash29" = "30", "548nash3" = "60",
                "548nash30" = "20", "548nash31" = "50", "548nash32" = "50",
                "548nash33" = "20", "548nash34" = "30", "548nash35" = "70",
                "548nash36" = "40", "548nash37" = "40", "548nash38" = "20",
                "548nash39" = "50", "548nash40" = "20", "548nash41" = "30",
                "548nash44" = "50", "548nash45" = "30", "548nash46" = "40",
                "548nash47" = "60", "548nash48" = "40", "548nash49" = "40",
                "548nash5" = "40", "548nash51" = "30", "548nash52" = "20",
                "548nash53" = "30", "548nash54" = "40", "548nash55" = "30",
                "548nash56" = "60", "548nash57" = "60", "548nash58" = "40",
                "548nash59" = "40", "548nash6" = "40", "548nash60" = "20",
                "548nash61" = "50", "548nash62" = "20", "548nash63" = "40",
                "548nash64" = "40", "548nash65" = "40", "548nash66" = "40",
                "548nash67" = "60", "548nash68" = "40", "548nash69" = "50",
                "548nash7" = "50", "548nash70" = "30", "548nash71" = "30",
                "548nash72" = "50", "548nash73" = "40", "548nash75" = "30",
                "548nash76" = "30", "548nash77" = "50", "548nash78" = "20",
                "548nash79" = "30", "548nash80" = "40", "548nash81" = "30",
                "548nash82" = "60", "548nash83" = "30", "548nash84" = "50",
                "548nash85" = "40", "548nash86" = "60", "548nash87" = "30",
                "548nash88" = "50", "548nash89" = "60", "548nash9" = "50",
                "548nash90" = "50", "548nash91" = "30", "548nash92" = "50",
                "548nash93" = "20", "548nash95" = "40", "548nash96" = "40",
                "548nash98" = "40", "548nash99" = "50")

FinalData[["fibrosis_Identity"]]<- fibrosis_vector
FinalData[["nas_Identity"]]<- nas_vector
FinalData[["age_Identity"]]<- age_vector


LabelPoints(plot = umap, points = cells_vector, repel = TRUE)
DimPlot(FinalData, reduction = "pca", cells = cells_vector, 
        cols = type_colors_vector, pt.size = 2.5)

# Finding differentially expressed features (cluster biomarkers)
FinalData.markers <- FindAllMarkers(FinalData, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
FinalData.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write_xlsx(FinalData.markers,"C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NASH_option1/PLOTS/EGR1/new_merged/ClusterMarkers.xlsx")


# EGR1
VlnPlot(object = FinalData, features = 'EGR1', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'EGR1', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'EGR1', group.by = 'age_Identity')

# RYR2
VlnPlot(object = FinalData, features = 'RYR2', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'RYR2', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'RYR2', group.by = 'age_Identity')

# LRFN5
VlnPlot(object = FinalData, features = 'LRFN5', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'LRFN5', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'LRFN5', group.by = 'age_Identity')

# SIK1
VlnPlot(object = FinalData, features = 'SIK1', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'SIK1', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'SIK1', group.by = 'age_Identity')

# PTPRT
VlnPlot(object = FinalData, features = 'PTPRT', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'PTPRT', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'PTPRT', group.by = 'age_Identity')

# DUSP1
VlnPlot(object = FinalData, features = 'DUSP1', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'DUSP1', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'DUSP1', group.by = 'age_Identity')

# CYP1A1
VlnPlot(object = FinalData, features = 'CYP1A1', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'CYP1A1', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'CYP1A1', group.by = 'age_Identity')

# IGFBP2
VlnPlot(object = FinalData, features = 'IGFBP2', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'IGFBP2', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'IGFBP2', group.by = 'age_Identity')

# KDM5D
VlnPlot(object = FinalData, features = 'KDM5D', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'KDM5D', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'KDM5D', group.by = 'age_Identity')

VlnPlot(FinalData, features = c("EGR1"), slot = "counts", log = TRUE)
FeaturePlot(FinalData, features = c("EGR1"), pt.size = 2.5)
FeaturePlot(FinalData, features = c("EGR1"), pt.size = 2.5, min.cutoff = 0.5, max.cutoff = 0.8)
