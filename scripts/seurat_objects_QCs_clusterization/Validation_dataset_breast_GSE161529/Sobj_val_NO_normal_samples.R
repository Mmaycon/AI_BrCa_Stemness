### TAKE NORMAL SAMPLES OUT OF SEURAT OBJECT ###
### RUN SEURAT PIPELINE OVER SOBJ_VAL ###
### PLOT ALL THAT AGAIN ###

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(filesstrings)
library(RColorBrewer)
library(GEOquery)
library(readr)

# 0. load Sobj_Val
# 1. add clinical metadata 
# 2. subset Sobj_Val without Normal cells
# 3. run seurat pipeline
# 4. save it
# 5. come back to "/mnt/plummergrp/maycon/Figures_to_paper/Compare_stemness/Compare_stemness.R" to plot all those visualizations again

### Load Sobj_Val
Sobj_Val = readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Valid_Merged_No_lynph_samples_yes_clusterization_CELLTYPE.STEMNESScomputed.rds")
dim(Sobj_Val) #21843 276371
head(Sobj_Val@meta.data)
names(Sobj_Val@meta.data)

# Fix Patient ID label
names(Sobj_Val@meta.data)[names(Sobj_Val@meta.data) == "PatientID"] <- "Patient_ID"


### Add clinical metadata 
load("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/GSE161529_pData_sub_temp_from_GEO_all_tumors_except_lynph.rda")
dim(pData_sub)

# Fix Patient ID label
names(pData_sub)[names(pData_sub) == "PatientID"] <- "Patient_ID"
# Check $Patient_ID in both metadata
all(pData_sub$Patient_ID %in% Sobj_Val@meta.data$Patient_ID) # TRUE

# First merge, then cbind, and over write it
meta_data = Sobj_Val@meta.data
meta_data_toadd = merge(meta_data, pData_sub) # merge 
identical(meta_data_toadd$Barcode, meta_data$Barcode) 
meta_data_toadd_final = cbind(meta_data, meta_data_toadd) #cbind

identical(Sobj_Val@meta.data$Barcode,
          meta_data_toadd_final$Barcode) #TRUE

Sobj_Val@meta.data = meta_data_toadd_final #over write it
head(Sobj_Val@meta.data) # okay



### Subset Sobj_Val
Idents(Sobj_Val) = "cancer_type_ch1"
Sobj_Val_tumors = subset(x = Sobj_Val, idents = c("Normal"), invert = TRUE)
dim(Sobj_Val_tumors) #21843 180315

table(Sobj_Val_tumors@meta.data$cancer_type_ch1) #no "Normal" samples, okay


### Run seurat pipeline
Sobj_Val_tumors = NormalizeData(Sobj_Val_tumors, normalization.method = "LogNormalize", scale.factor = 10000)
#range(Sobj_Val_tumors[["RNA"]]@data) # 0.000000 9.111652
Sobj_Val_tumors = FindVariableFeatures(Sobj_Val_tumors, selection.method = 'vst', nfeatures = 3000, verbose = FALSE)
Sobj_Val_tumors = ScaleData(Sobj_Val_tumors, verbose = FALSE) #scales only most varible features
Sobj_Val_tumors = RunPCA(Sobj_Val_tumors, verbose = FALSE, features = VariableFeatures(object = Sobj_Val_tumors)) 
#ElbowPlot(Sobj_Val_tumors) # at least 20 PCs
Sobj_Val_tumors = FindNeighbors(Sobj_Val_tumors, reduction = 'pca', dims = 1:30, verbose = FALSE)
Sobj_Val_tumors = FindClusters(Sobj_Val_tumors, resolution = 0.5, verbose = FALSE)
Sobj_Val_tumors = RunUMAP(Sobj_Val_tumors, reduction = 'pca', dims = 1:30, verbose = FALSE)

DimPlot(Sobj_Val_tumors)

bk_Sobj_Val_tumors = Sobj_Val_tumors
### Add  new stemness index (predicted withou Normal samples)
load("/mnt/plummergrp/maycon/Figures_to_paper/stemness_index_prediction/2023_03_09_Stemness_index/reference_validation_and_iPS_predicted_together_except_Normalsamples.rda")
dim(stemnes_index)
stemnes_index$Barcode = rownames(stemnes_index)
stemnes_index$cell_barcode = NULL
names(stemnes_index)[names(stemnes_index) == "stemness"] <- "noNormal_stemness"
names(Sobj_Val_tumors@meta.data)[names(Sobj_Val_tumors@meta.data) == "stemness"] <- "all_stemness"

##**NOOO IT IS TAKING TO MUCH MEMORY**## - GO TO ADDMETADATA OPTION 2
# Merge, cbind, and over writte (merge option 1)
meta_data = Sobj_Val_tumors@meta.data
meta_data_to_add = merge(meta_data, stemnes_index)

identical(meta_data_to_add$cell_barcode, meta_data$cell_barcode) #FALSE
# Put these dfs in order
rownames(meta_data_to_add) = meta_data_to_add$cell_barcode
meta_data_to_add = meta_data_to_add[meta_data$cell_barcode, ]
identical(meta_data_to_add$cell_barcode, meta_data$cell_barcode) #TRUE, okay
meta_data_to_add_indeed = cbind(meta_data_to_add, meta_data)
Sobj_Val_tumors@meta.data = meta_data_to_add_indeed
##**NOOO IT IS TAKING TO MUCH MEMORY**##


# AddMetaData (merge option 2)
dim(Sobj_Val_tumors) #21843 180315
Sobj_Val_tumors = AddMetaData(
  object = Sobj_Val_tumors,
  metadata = stemnes_index
)
dim(Sobj_Val_tumors) #21843 180315
head(Sobj_Val_tumors@meta.data)
Sobj_Val_tumors$noNormal_stemness


# Desnsity plot to compare stemenss versions
dens_plot_1 = data.frame(stemness = Sobj_Val_tumors@meta.data$all_stemness,
                         stemness_version = "all_stemness")

dens_plot_2 = data.frame(stemness = Sobj_Val_tumors@meta.data$noNormal_stemness,
                         stemness_version = "noNormal_stemness")

dens_plot = plyr::rbind.fill(dens_plot_1, dens_plot_2)

ggplot(dens_plot, aes(x = stemness, fill = stemness_version)) + 
  geom_density(alpha = 0.5) + 
  theme_classic() 





### Save it
saveRDS(Sobj_Val_tumors, file = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Val_onlyd_tumors_No_lynph_samples_yes_clusterization_yes_stemness_without_normalsmp.rds")


### Come back to "/mnt/plummergrp/maycon/Figures_to_paper/Compare_stemness/Compare_stemness.R" to plot all those visualizations again


