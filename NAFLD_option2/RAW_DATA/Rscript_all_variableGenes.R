library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(rlang)
library(stringr)
library(patchwork)
library(Matrix)

rawdata <- read.delim("C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/RAW_DATA/arranged_before_surgery_genes.csv", sep=",")

#rownames(rawdata) <- rawdata[,1]
rownames(rawdata) = make.names(rawdata[,1], unique=TRUE)
rawdata <- rawdata[,-1]

SingleCellSeq <- CreateSeuratObject(counts = rawdata,project = "SingleCellSeq", assay = "RNA", min.cells = 3, min.features = 200)

FinalData <- NormalizeData(SingleCellSeq, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling the data
all.genes <- rownames(FinalData)
FinalData <- ScaleData(FinalData, features = all.genes)

Idents(object = FinalData)

Nash_labels <- c("GSM1178975", "GSM1178980", "GSM1178987", "GSM1178995",
                  "GSM1178996", "GSM1179001", "GSM1179002", "GSM1179003", 
                  "GSM1179004", "GSM1179006", "GSM1179008", "GSM1179015",
                  "GSM1179017", "GSM1179026", "GSM1179033", "GSM1179035", 
                  "GSM1179036")


cells_vector = c('GSM1178970' = 'Control', 'GSM1178971' = 'Control', 'GSM1178972' = 'Control', 
                  'GSM1178973' = 'Control', 'GSM1178974' = 'Control', 'GSM1178975' = 'Nash',
                  'GSM1178976' = 'Healthy Obese', 'GSM1178977' = 'Control', 'GSM1178978' = 'Control',
                  'GSM1178979' = 'Control', 'GSM1178980' = 'Nash', 'GSM1178981' = 'Healthy Obese',
                  'GSM1178982' = 'Healthy Obese', 'GSM1178983' = 'Healthy Obese', 'GSM1178984' = 'Healthy Obese', 
                  'GSM1178985' = 'Healthy Obese', 'GSM1178986' = 'Steatosis', 'GSM1178987' = 'Nash',
                  'GSM1178988' = 'Steatosis', 'GSM1178989' = 'Steatosis', 'GSM1178990' = 'Healthy Obese',
                  'GSM1178991' = 'Healthy Obese', 'GSM1178992' = 'Healthy Obese', 'GSM1178993' = 'Steatosis',
                  'GSM1178994' = 'Healthy Obese', 'GSM1178995' = 'Nash', 'GSM1178996' = 'Nash',
                  'GSM1178998' = 'Control', 'GSM1178999' = 'Steatosis', 'GSM1179001' = 'Nash',
                  'GSM1179002' = 'Nash', 'GSM1179003' = 'Nash', 'GSM1179004' = 'Nash', 
                  'GSM1179005' = 'Healthy Obese', 'GSM1179006' = 'Nash', 'GSM1179007' = 'Healthy Obese',
                  'GSM1179008' = 'Nash', 'GSM1179010' = 'Control', 'GSM1179012' = 'Healthy Obese',
                  'GSM1179015' = 'Nash', 'GSM1179016' = 'Healthy Obese', 'GSM1179017' = 'Nash',
                  'GSM1179018' = 'Control', 'GSM1179021' = 'Steatosis', 'GSM1179024' = 'Control',
                  'GSM1179025' = 'Steatosis', 'GSM1179026' = 'Nash', 'GSM1179027' = 'Steatosis',
                  'GSM1179030' = 'Healthy Obese', 'GSM1179033' = 'Nash', 'GSM1179035' = 'Nash',
                  'GSM1179036' = 'Nash', 'GSM1179037' = 'Steatosis', 'GSM1179038' = 'Healthy Obese')


#Idents(object = FinalData) <- cells_vector
#head(x = Idents(object = FinalData))


FinalData[["Identity"]]<- cells_vector

VlnPlot(object = FinalData, features = 'ENSG00000134184', group.by = 'Identity')

# ALL
VlnPlot(FinalData, features = c("ENSG00000134184", "ENSG00000184674.1", 
                                "ENSG00000170345", "ENSG00000164266", 
                                "ENSG00000167910", "ENSG00000120738",
                                "ENSG00000146678", "ENSG00000188257",
                                "ENSG00000133048", "ENSG00000122121"), 
        slot = "counts", log = TRUE, group.by = 'Identity')


#FOS
FOS_vlnPlot <- VlnPlot(object = FinalData, features = 'ENSG00000170345', group.by = 'Identity')
LabelPoints(plot = FOS_vlnPlot, points = Nash_labels, repel = TRUE)

#EGR1
EGR1_vlnPlot <- VlnPlot(object = FinalData, features = 'ENSG00000120738', group.by = 'Identity')
LabelPoints(plot = EGR1_vlnPlot, points = Nash_labels, repel = TRUE)

