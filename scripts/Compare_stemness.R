### Compare stemnessin both datasets Reference and Validation 
### along with clinical data (clinical subtype etc) 
#OBS: sciPS will not be part of this analysis 

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(filesstrings)
library(RColorBrewer)
library(GEOquery)
library(readr)


# Load Ref Sobj
# Check meta.data info - If necessary, add it
# Sobj_Ref = readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Ref_dset/Reference_dataset_standart_seut_pipeline_noNormal_stemnesscomputed.rds") - stemnesspredicted to normal samples together with tumors
Sobj_Ref = readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Ref_dset/Reference_dataset_standart_seut_pipeline_NEW_noNormal_stemnesscomputed.rds")
dim(Sobj_Ref) # 27719 80682
head(Sobj_Ref@meta.data)
names(Sobj_Ref@meta.data)


# Load Valid Sobj 
# Check meta.data info - If necessary, add it
# Sobj_Val = readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Valid_Merged_No_lynph_samples_yes_clusterization_CELLTYPE.noNormal_stemnesscomputed.rds") - without clinical info; there are normal samples; stemnesspredicted to normal samples together with tumors
##**Only tumors**##
Sobj_Val = readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Val_onlyd_tumors_No_lynph_samples_yes_clusterization_yes_noNormal_stemness_without_normalsmp.rds")
dim(Sobj_Val) #21843 180315
head(Sobj_Val@meta.data)
names(Sobj_Val@meta.data)



# Plot noNormal_stemness, clinical info. and stack barplot for cell proportion - see Felipe package suggestion


### Stack barplot - cell proportion by patient
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dittoSeq")

library(dittoSeq)
dittoBarPlot(Sobj_Ref, "celltype_major", group.by = "Patient_ID")
dittoBarPlot(Sobj_Val, "Major_celltypes_sctyped", group.by = "Patient_ID")

table(Sobj_Ref@meta.data$Patient_ID) # to many patients to visualize in a table
table(Sobj_Val@meta.data$Patient_ID) # to many patients to visualize in a table

### Rank plot of patient cells - instead of tables

# Sobj_Ref
rankplot = data.frame(table(Sobj_Ref@meta.data$Patient_ID))
names(rankplot) = c("Patient_ID", "cells")
rankplot <- rankplot[order(rankplot$cells),] #ordena

library(viridis)
n_colors <- length(table(rankplot$Patient_ID))
pal <- viridis(n = n_colors, option = "G", direction = -1)
library(ggplot2)
ggplot(rankplot, aes(reorder(Patient_ID, +cells, sum), cells)) +
  geom_bar(stat="identity", aes(fill=factor(Patient_ID))) +
  scale_fill_manual(values=c(pal), name="Legend") +
  ggtitle("Cells by patient") +
  ylab("Number of cells") + 
  xlab("Patient ID") + 
  theme(axis.title.x = element_text(colour="black", size = 18), 
      axis.title.y = element_text(colour="black", size = 18),
      legend.text = element_text(colour="black", size = 16),
      legend.title = element_text(colour="black", size = 18)) + 
  #geom_jitter(alpha=0.4) +
  theme(axis.text.x = element_text(angle = 45, size = 10 ,vjust=1,hjust=1),
        axis.text.y = element_text(angle = 0, size = 16)) +
  guides(fill = FALSE)  


# Sobj_Val
rankplot = data.frame(table(Sobj_Val@meta.data$Patient_ID))
names(rankplot) = c("Patient_ID", "cells")
rankplot <- rankplot[order(rankplot$cells),] #ordena

library(viridis)
n_colors <- length(table(rankplot$Patient_ID))
pal <- viridis(n = n_colors, option = "G", direction = -1)
library(ggplot2)
ggplot(rankplot, aes(reorder(Patient_ID, +cells, sum), cells)) +
  geom_bar(stat="identity", aes(fill=factor(Patient_ID))) +
  scale_fill_manual(values=c(pal), name="Legend") +
  ggtitle("Cells by patient") +
  ylab("Number of cells") + 
  xlab("Patient ID") + 
  theme(axis.title.x = element_text(colour="black", size = 18), 
        axis.title.y = element_text(colour="black", size = 18),
        legend.text = element_text(colour="black", size = 16),
        legend.title = element_text(colour="black", size = 18)) + 
  #geom_jitter(alpha=0.4) +
  theme(axis.text.x = element_text(angle = 45, size = 10 ,vjust=1,hjust=1),
        axis.text.y = element_text(angle = 0, size = 16)) +
  guides(fill = FALSE) 

# Cluster proportion
library(dittoSeq)
dittoBarPlot(Sobj_Ref, "seurat_clusters", group.by = "Patient_ID")
dittoBarPlot(Sobj_Val, "seurat_clusters", group.by = "Patient_ID")

# stemnessview
FeaturePlot(object = Sobj_Ref, features = c("noNormal_stemness")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle("noNormal_stemness")

FeaturePlot(object = Sobj_Val, features = c("noNormal_stemness")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle("noNormal_stemness")


# QC vln 
Idents(Sobj_Ref) = "Patient_ID"
VlnPlot(Sobj_Ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

Idents(Sobj_Val) = "Patient_ID"
VlnPlot(Sobj_Val, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

# Seurat clusters
library(viridis)
n_colors <- length(table(Sobj_Ref@meta.data$seurat_clusters))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj_Ref) = "seurat_clusters"
DimPlot(Sobj_Ref) + scale_color_manual(values = pal)

n_colors <- length(table(Sobj_Val@meta.data$seurat_clusters))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj_Val) = "seurat_clusters"
DimPlot(Sobj_Val) + scale_color_manual(values = pal)

# Patients
library(viridis)
n_colors <- length(table(Sobj_Ref@meta.data$Patient_ID))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj_Ref) = "Patient_ID"
DimPlot(Sobj_Ref, label = FALSE) + scale_color_manual(values = pal)

n_colors <- length(table(Sobj_Val@meta.data$Patient_ID))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj_Val) = "Patient_ID"
DimPlot(Sobj_Val, label = FALSE) + scale_color_manual(values = pal)

# stemness~ clusters (RidgePlot)
Idents(Sobj_Ref) <- "seurat_clusters"
RidgePlot(Sobj_Ref, features = "noNormal_stemness", ncol = 1)

Idents(Sobj_Val) <- "seurat_clusters"
RidgePlot(Sobj_Val, features = "noNormal_stemness", ncol = 1)


# Clinical subtypes (vln, noNormal_stemness, by patient,)
Idents(Sobj_Ref) <- "clinical_subtype"
VlnPlot(Sobj_Ref, features = "noNormal_stemness", split.by = "Patient_ID") + 
  theme(axis.title.x = element_text(colour="black", size = 22))

# IT IS NOT PLOTTING. IDK WHY
bk = Sobj_Val
Idents(Sobj_Val) <- "cancer_type_ch1"
VlnPlot(Sobj_Val, features = "noNormal_stemness", split.by = "Patient_ID") +
  theme(axis.title.x = element_text(colour="black", size = 22))


# Clinical subtypes (UMAP)
library(viridis)
n_colors <- length(table(Sobj_Ref@meta.data$clinical_subtype))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj_Ref) = "clinical_subtype"
DimPlot(Sobj_Ref) + scale_color_manual(values = pal)


n_colors <- length(table(Sobj_Val@meta.data$cancer_type_ch1))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj_Val) = "cancer_type_ch1"
DimPlot(Sobj_Val, label = FALSE) + scale_color_manual(values = pal)


### Fixing Sob_Val metadata variables (start)
x = Sobj_Val@meta.data$cancer_type_ch1
Sobj_Val@meta.data$cancer_type_ch1 = NA
Sobj_Val@meta.data$cancer_type_ch1 = x
# Create an other column with the same levels
Sobj_Val@meta.data$clinical_subtype =  NA
Sobj_Val@meta.data$clinical_subtype =  Sobj_Val@meta.data$cancer_type_ch1

n_colors <- length(table(Sobj_Val@meta.data$clinical_subtype))
pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj_Val) = "clinical_subtype"
DimPlot(Sobj_Val, label = FALSE) + scale_color_manual(values = pal)

### Fixing Sob_Val metadata variables (end)


table(Sobj_Val@meta.data$cancer_type_ch1)
head(Sobj_Val@meta.data)
str(Sobj_Val@meta.data$cancer_type_ch1)
dim(Sobj_Val)
dim(Sobj_Val@meta.data)
table(is.na(Sobj_Val@meta.data$cancer_type_ch1))
summary(Sobj_Val@meta.data$cancer_type_ch1)


class(Sobj_Val@meta.data$cancer_type_ch1)
Sobj_Val@meta.data$

  


# Minor celltypes (UMAP)

# Minor celltypes (UMAP) + Major celltypes (UMAP)

# RidgePlot for Minor celltypes


# REMINDER: WE MUST SEE CONCORNDACE ON HIGH stemnessMINOR CELL TYPES IN BOTH DATASETS; PAY ATTENTION HOW THOSE 3 CELLTYPES ARE CLUSTERED; WHAT ABOUT NORMAL CELLTYPES ON VALIDATIOND DATASET??? 


# Give colors only to Cancer Cycling, Cancer Basal, Cycling Tcells - the high stemnesscells ones
DimPlot(Sobj) + scale_color_manual(values = c(rep("gray",12),"gray","red","gray","gray","gray",rep("gray",6),"black",rep("gray",2),"orange",rep("gray",2)))



