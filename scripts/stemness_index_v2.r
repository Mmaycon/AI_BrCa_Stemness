######******            ******###### 


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


### Reference expression matrix  ###

# Load Seurat Object
Reference_dataset_standart_seut_pipeline <- readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Ref_dset/Reference_dataset_standart_seut_pipeline.rds")
Sobj_Ref = Reference_dataset_standart_seut_pipeline
dim(Sobj_Ref) #21843 276371
head(Sobj_Ref@meta.data)
names(Sobj_Ref@meta.data)
VlnPlot(Sobj_Ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
range(Sobj_Ref[["RNA"]]@) # 0.000000 9.063545

### Validation expression matrix  ###
Seurat_object = readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Valid_Merged_No_lynph_samples_yes_clusterization_CELLTYPE.rds")
Sobj_Val = Seurat_object
dim(Sobj_Val) #21843 276371
head(Sobj_Val@meta.data)
names(Sobj_Val@meta.data)
VlnPlot(Sobj_Val, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
range(Sobj_Val[["RNA"]]@data) #0.000000 9.020428

# Load Sobj_Val clinical data
load("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/GSE161529_pData_sub_temp_from_GEO_all_tumors_except_lynph.rda")
dim(pData_sub)

# Add clinical info into metadata  
# First merge, then cbind
meta_data = Sobj_Val@meta.data
meta_data_toadd = merge(meta_data, pData_sub) 

identical(meta_data_toadd$Barcode, meta_data$Barcode) 
meta_data_toadd_final = cbind(meta_data, meta_data_toadd)

identical(Sobj_Val@meta.data$Barcode,
          meta_data_toadd_final$Barcode) #TRUE

Sobj_Val@meta.data = meta_data_toadd_final
head(Sobj_Val@meta.data) # okay

Idents(Sobj_Val) = "cancer_type_ch1"
Sobj_Val = subset(x = Sobj_Val, idents = c("Normal"), invert = TRUE)
dim(Sobj_Val) #21843 180315

### iPS expression matrix ###

# Creating seurat object
  # iPS experiment 1
iPS_smp1 <- fread("/mnt/storage1/mmaycon/downloads/sc_iPS/IPSCMedia_scRNA_Sample1_V1_Count.txt") 
iPS_smp1 <- iPS_smp1 %>% as.data.frame() 
rownames(iPS_smp1) <- iPS_smp1$V1
iPS_smp1[1:4,1:4]
iPS_smp1$V1 <- NULL
iPS_smp1 <- iPS_smp1 %>% as.matrix() 
dim(iPS_smp1) # 32838  8740
colnames(iPS_smp1) <- gsub("_", "-", colnames(iPS_smp1)) 
rownames(iPS_smp1) <- gsub("_", "-", rownames(iPS_smp1)) 
Sobj_1 <- CreateSeuratObject(counts = iPS_smp1, min.cells = 3, min.features = 200)

  # iPS experiment 2
iPS_smp2 <- fread("/mnt/storage1/mmaycon/downloads/sc_iPS/IPSCMedia_scRNA_Sample2_V1_Count.txt") 
iPS_smp2 <- iPS_smp2 %>% as.data.frame() 
rownames(iPS_smp2) <- iPS_smp2$V1
iPS_smp2[1:4,1:4]
iPS_smp2$V1 <- NULL
iPS_smp2 <- iPS_smp2 %>% as.matrix() 
dim(iPS_smp2) # 32838 12859
colnames(iPS_smp2) <- gsub("_", "-", colnames(iPS_smp2)) 
rownames(iPS_smp2) <- gsub("_", "-", rownames(iPS_smp2)) 
Sobj_2 <- CreateSeuratObject(counts = iPS_smp2, min.cells = 3, min.features = 200)

Sobj_iPS <- merge(x = Sobj_1, 
              y = Sobj_2, add.cell.ids = c("iPS_1","iPS_2"), 
              project = "sc_iPS_dataset")

Sobj_iPS[["percent.mt"]] <- PercentageFeatureSet(Sobj_iPS, pattern = "^MT-")
VlnPlot(Sobj_iPS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)  
range(Sobj_iPS[["RNA"]]@data) #0 482 

### Combine them all in one matrix ###

# Combine seurat objects
dim(Sobj_Ref) #27719 80682
dim(Sobj_Val) #21843 276371
dim(Sobj_iPS) #20202 21599
Sobj_all <- merge(x = Sobj_Ref, 
                  y = c(Sobj_Val, 
                        Sobj_iPS),
                  add.cell.ids = c("REF","VAL","IPS"))
                  #project = "")
dim(Sobj_all) #31404 378652
range(Sobj_all[["RNA"]]@data)  #0 482 (iPS must have the highst gene count that's why this range is the same on Sobj_iPS only)
colnames(Sobj_all)



# LogNormalize method  ---------------------------------------
Sobj_all <- NormalizeData(Sobj_all, normalization.method = "LogNormalize", scale.factor = 10000)
range(Sobj_all[["RNA"]]@data) #

# Extract matrix expression to apply Stemness model (Malta et. al. 2018)
matrix_all <- as.matrix(GetAssayData(object = Sobj_all, slot = "data"))
dim(matrix_all) #31404 282596


###Predict stemness index

# Loading stemness genes wights vector (- the model)
load("/mnt/storage1/mmaycon/downloads/model_RNA_MALTA.2018.Rda") # mm = modelo
w = mm$w
w[1:5]
length(w) #12953 (number of genes on the model)
# Filter matrix expression by the genes on stemenss model
matrix = matrix_all
length(intersect(rownames(matrix), names(w))) #12922 
predict.DATA = matrix[rownames(matrix) %in% names(w) ,]
length(rownames(predict.DATA)) # 12922
w = w[ rownames(predict.DATA) ]
length(intersect(names(w),rownames(predict.DATA))) # 12922 
w[1:5]
length(names(w)) # 12380
is.vector(w) #TRUE

# Score the Matrix `X` using Spearman correlation
s = apply( predict.DATA, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
s[1:5]

# Scale the scores to be between 0 and 1
s = s - min(s)
s = s / max(s)
s[1:5]
s = as.data.frame(s)
names(s) = "stemness"
library(stringr)
rownames(s) = str_split_fixed(as.character(rownames(s)), "[_]", 2)[,2]

stemnes_index = s
save(stemnes_index, file = "/mnt/plummergrp/maycon/Figures_to_paper/stemness_index_prediction/2023_03_09_Stemness_index/reference_validation_and_iPS_predicted_together_except_Normalsamples.rda")

# stemness index histogram
ggplot(s, aes(x=stemness)) + 
  geom_histogram(bins = 50) +
  xlab("Stemness index") + ylab("n Cells") + theme(axis.title.x = element_text(colour="black", size = 12), axis.title.y = element_text(colour="black", size = 15)) + ggtitle("") +theme(axis.text.x = element_text(angle = 45, size = 15 ,vjust=1,hjust=1))


### Merge each subset meta_data into its respectively seurat object
# Load stemness index predicted 
# load("/mnt/plummergrp/maycon/Figures_to_paper/stemness_index_prediction/2023_03_09_Stemness_index/reference_validation_and_iPS_predicted_together.rda")
# head(stemnes_index)

# # Reference dataset - adding stemness index
# dim(Sobj_Ref) #27719 80682
# Sobj_Ref = AddMetaData(
#   object = Sobj_Ref,
#   metadata = s
# )
# dim(Sobj_Ref) #27719 81376
# head(Sobj_Ref@meta.data)
# 
# # Validation dataset - adding stemness index
# dim(Sobj_Val) #21843 276371
# Sobj_Val = AddMetaData(
#   object = Sobj_Val,
#   metadata = s
# )
# dim(Sobj_Val) #21843 276371
# head(Sobj_Val@meta.data)
# 
# 
# # iPS single cell dataset - adding stemness index
# dim(Sobj_iPS) #20202 21599
# Sobj_iPS = AddMetaData(
#   object = Sobj_iPS,
#   metadata = s
# )
# dim(Sobj_iPS) #20202 21599
# head(Sobj_iPS@meta.data)

# # Saving the newest Seurat object version 
# saveRDS(Sobj_Ref, "/mnt/plummergrp/maycon/Figures_to_paper/Ref_dset/Reference_dataset_standart_seut_pipeline_STEMNESScomputed.rds")
# 
# saveRDS(Sobj_Val, "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Valid_Merged_No_lynph_samples_yes_clusterization_CELLTYPE.STEMNESScomputed.rds")
# 
# saveRDS(Sobj_iPS, "/mnt/plummergrp/maycon/Figures_to_paper/iPS/Sobj_iPSsinglecell_no_clusterization.STEMNESScomputed.rds")

