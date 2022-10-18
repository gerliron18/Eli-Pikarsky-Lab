library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(rlang)
library(stringr)
library(patchwork)
library(Matrix)

rawdata <- read.delim("C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NAFLD_option2/RAW_DATA/arranged_before_surgery_genes.csv", sep=",")

rownames(rawdata) = make.names(rawdata[,1], unique=TRUE)
rawdata <- rawdata[,-1]

SingleCellSeq <- CreateSeuratObject(counts = rawdata,project = "SingleCellSeq", assay = "RNA", min.cells = 3, min.features = 200)

FinalData <- NormalizeData(SingleCellSeq, normalization.method = "LogNormalize", scale.factor = 10000)

FinalData <- FindVariableFeatures(FinalData, selection.method = "vst", nfeatures = 2000)

# Scaling the data
all.genes <- rownames(FinalData)
FinalData <- ScaleData(FinalData, features = all.genes)

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


# EGR1
VlnPlot(object = FinalData, features = 'ENSG00000120738', group.by = 'fibrosis_Identity')
VlnPlot(object = FinalData, features = 'ENSG00000120738', group.by = 'nas_Identity')
VlnPlot(object = FinalData, features = 'ENSG00000120738', group.by = 'age_Identity')
