# Header 
# Title: Standart Seurat Pipeline on our Reference Breast Cancer Dataset
# Type of data: scRNA-seq of a breast cqncer atlas
# Brief description: Build Reference Seurat Object from callranger output. Then, set QCs and visualize author's pdata information & stemness index
# Date: 12/19/2022
# Author: Maycon Marção 	

# Packages
library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(filesstrings)
library(RColorBrewer)
library(GEOquery)
library(readr)
library(HelpersMG)
library(ggpubr)
library(viridis)


# Download scRNA-seq data from GEO 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
setwd("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject")
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2) # turning to 
getOption('timeout') # 60
options(timeout=300) 
wget("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE176078&format=file&file=GSE176078%5FWu%5Fetal%5F2021%5FBRCA%5FscRNASeq%2Etar%2Egz")

untar("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/?acc=GSE176078&format=file&file=GSE176078%5FWu%5Fetal%5F2021%5FBRCA%5FscRNASeq%2Etar%2Egz")

file_names = dir("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/Wu_etal_2021_BRCA_scRNASeq")
file.rename(file_names[1],"barcodes.tsv")
file.rename(file_names[2],"genes.tsv")
file.rename(file_names[3],"matrix.mtx")

# Part 1 (start) ---------------------
## Create Seurat object 
  # Load sparse matrix
GSE176078_data = Read10X(data.dir = "/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/Wu_etal_2021_BRCA_scRNASeq",gene.column=1)
dim(GSE176078_data) #29733 100064
range(GSE176078_data) #0 29831
  # Edit strings on rownames and colnames
colnames(GSE176078_data) = gsub("_", "-", colnames(GSE176078_data)) 
rownames(GSE176078_data) = gsub("_", "-", rownames(GSE176078_data)) 
  # Create Seurat Object
Sobj = CreateSeuratObject(counts = GSE176078_data, min.cells = 3, min.features = 200)
dim(Sobj[["RNA"]]@data) # 27719 100064


## Quality control 
Sobj[["percent.mt"]] = PercentageFeatureSet(Sobj, pattern = "^MT-")
Sobj$Patient_ID <- str_split_fixed(as.character(names(Sobj$orig.ident)), "[-]", 2)[,1]
Idents(Sobj) <- "Patient_ID"
VlnPlot(Sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
saveRDS(Sobj, file = "/mnt/plummergrp/maycon/Figures_to_paper/Ref_dset/Sobj_RefDataset_NO_QCs.rds")

Sobj = subset(Sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA > 200)
VlnPlot(Sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
dim(Sobj[["RNA"]]@data) #27719 81376

## Normalize the data
Sobj = NormalizeData(Sobj, normalization.method = "LogNormalize", scale.factor = 10000)
range(Sobj[["RNA"]]@data) # 0.000000 9.111652

##**NO - Now we predict for every dataset together and then mmovie on with separated Sobj**##
# ## Score stemness ##**NO**##
# matrix = as.matrix(GetAssayData(object = Sobj, slot = "data"))
# dim(matrix) # 27719 81376
#   # Loading stemness genes wights vector (- the model)
# load("/mnt/storage1/mmaycon/downloads/model_RNA_MALTA.2018.Rda") # mm = modelo
# w = mm$w
# w[1:5]
# length(w) #12953 (number of genes on the model)
#   # Filter matrix expression by the genes on stemenss model
# length(intersect(rownames(matrix), names(w))) #12380 
# predict.DATA = matrix[rownames(matrix) %in% names(w) ,]
# length(rownames(predict.DATA)) #12380 
# w = w[ rownames(predict.DATA) ]
# length(intersect(names(w),rownames(predict.DATA))) # 12380 
# w[1:5]
# length(names(w)) # 12380
# is.vector(w) #TRUE
# 
#   # Score the Matrix `X` using Spearman correlation
# s = apply( predict.DATA, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
# s[1:5]
# 
#   # Scale the scores to be between 0 and 1
# s = s - min(s)
# s = s / max(s)
# s[1:5]
# 
#   # stemness index histogram
#   ggplot(s, aes(x=s)) + 
#   geom_histogram(bins = 50) +
#   xlab("Stemness index") + ylab("n Cells") + theme(axis.title.x = element_text(colour="black", size = 12), axis.title.y = element_text(colour="black", size = 15)) + ggtitle("") +theme(axis.text.x = element_text(angle = 45, size = 15 ,vjust=1,hjust=1))
# 
#   # Merge stemness index to Seurat metadata 
# s = as.data.frame(s)
# names(s) = "stemness"
# Sobj = AddMetaData(
#   object = Sobj,
#   metadata = s
# )



## Download pdata 
# Single Cell Portal - need to login before download any data
# https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers#study-download


library(readr)
Whole_miniatlas_meta = read_csv("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/Whole_miniatlas_meta.csv")
head(Whole_miniatlas_meta)
Whole_miniatlas_meta = as.data.frame(Whole_miniatlas_meta)
Whole_miniatlas_meta$NAME = gsub("_", "-", Whole_miniatlas_meta$NAME) 
rownames(Whole_miniatlas_meta) = Whole_miniatlas_meta$NAME

  # Merge pdata to Seurat object
Sobj = AddMetaData(
  object = Sobj,
  metadata = Whole_miniatlas_meta
)


## Standard Seurat pipeline
# NormalizeData() has been used before
Sobj = FindVariableFeatures(Sobj, selection.method = 'vst', nfeatures = 3000, verbose = FALSE)
Sobj = ScaleData(Sobj, verbose = FALSE) #scales only most varible features
Sobj = RunPCA(Sobj, verbose = FALSE, features = VariableFeatures(object = Sobj)) 
ElbowPlot(Sobj) # 20 PCs
Sobj = FindNeighbors(Sobj, reduction = 'pca', dims = 1:20, verbose = FALSE)
Sobj = FindClusters(Sobj, resolution = 0.5, verbose = FALSE)
Sobj = RunUMAP(Sobj, reduction = 'pca', dims = 1:20, verbose = FALSE)


## Seurat visualization
  # UMAP seurat clusters
n_colors <- length(table(Sobj@meta.data$seurat_clusters))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj) <- "seurat_clusters"
DimPlot(Sobj) + scale_color_manual(values = pal)

  # UMAP major celltypes
n_colors <- length(table(Sobj@meta.data$celltype_major))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj) <- "celltype_major"
DimPlot(Sobj) + scale_color_manual(values = pal)

  # UMAP patients
Sobj$Patient_ID <- str_split_fixed(as.character(names(Sobj$orig.ident)), "[-]", 2)[,1]
n_colors <- length(table(Sobj@meta.data$Patient_ID))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj) <- "Patient_ID"
DimPlot(Sobj) + scale_color_manual(values = pal)

  # UMAP stemness
Idents(Sobj) <- "stemness"
FeaturePlot(object = Sobj, features = c("stemness")) & 
scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

  # Ridge plots - stemness by major celltype
Idents(Sobj) <- "celltype_major"
RidgePlot(Sobj, features = "stemness", ncol = 1)

  # Vln - stemness by major celltype
Idents(Sobj) <- "celltype_major"
VlnPlot(Sobj, features = "stemness")

  # UMAP of clinical subtypes
    # Recovering clinical subtype information from our merged Ref_Valid_iPS Sobj
# Sobj_ref_valid_iPS_0.9_transferedlabel <- readRDS("/mnt/plummergrp/maycon/Sobj_ref_valid_iPS_0.9_transferedlabel.rds")

pdata = as.data.frame(Sobj_ref_valid_iPS_0.9_transferedlabel@meta.data)
pdata$cell_barcode = str_split_fixed(as.character(pdata$cell_barcode), "[_]", 2)[,2]

pdata = pdata[,c("cell_barcode","clinical_subtype")]
rownames(pdata) = pdata$cell_barcode
Sobj = Sobj[,colnames(Sobj) %in% pdata$cell_barcode]

pdata <- pdata[rownames(pdata) %in% colnames(Sobj),]

dim(Sobj) #27719 80682
dim(pdata) #80682     2

    # Merge clinical subtype info. to Seurat object
Sobj = AddMetaData(
  object = Sobj,
  metadata = pdata
)

dim(Sobj) # 27719 80682
Sobj@meta.data$clinical_subtype[1:5]

n_colors <- length(table(Sobj@meta.data$clinical_subtype))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj) <- "clinical_subtype"
DimPlot(Sobj) + scale_color_manual(values = pal)

saveRDS(Sobj, file = "/mnt/plummergrp/maycon/Reference_dataset_standart_seut_pipeline.rds")

meta_data = as.data.frame(Sobj@meta.data)
save(meta_data, file = "/mnt/plummergrp/maycon/Reference_dataset_standart_seut_pipeline_MetaData.rda")

n_colors <- length(table(Sobj@meta.data$celltype_minor))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj) <- "celltype_minor"
DimPlot(Sobj) + scale_color_manual(values = pal)

DimPlot(Sobj) + scale_color_manual(values = c(rep("gray",12),"gray","red","gray","gray","gray",rep("gray",6),"black",rep("gray",2),"orange",rep("gray",2)))


Idents(Sobj) <- "celltype_minor"
RidgePlot(Sobj, features = "stemness", ncol = 1)

Idents(Sobj) <- "clinical_subtype"
VlnPlot(Sobj, features = "stemness", split.by = "Patient_ID")

Idents(Sobj) <- "clinical_subtype"
VlnPlot(Sobj, features = "stemness", split.by = "Patient_ID")

Idents(Sobj) <- "clinical_subtype"
VlnPlot(Sobj, features = "stemness", split.by = "Patient_ID")

Sobj_SCSubtype = subset(Sobj, subset = celltype_minor %in% c("Cancer Basal SC","Cancer LumA SC","Cancer LumB SC","Cancer Her2 SC"))
Idents(Sobj_SCSubtype) <- "celltype_minor"
VlnPlot(Sobj_SCSubtype, features = "stemness", split.by = "Patient_ID")


table(Sobj@meta.data$clinical_subtype, Sobj@meta.data$Patient_ID, )



