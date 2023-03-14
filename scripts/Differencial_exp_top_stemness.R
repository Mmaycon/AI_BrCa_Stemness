
### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
### Differential Expression - Top vs Bottom stemness cells within 
### the following cell types: Cancer Cycling, Cancer Basal, Cycling Tcells
### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 


# Load Sobj_Ref (reference) and Sobj_Val (validation) - follow all further analysis with these two datasets in paralel
# Plot Ridge plot by stemness ~ minor cell types - make sure those 3 cell types have indeed the highst stemness index 

# Subset Sobj by the three cell types
# eg: Cycling_tcells_Valid = subset() etc

# Generate stemness percentile label (top, bottom, in between)  

# FindMarkers group 1 = top stemness; group 2 = bottom stemenss
# Filter significant markers 
# upset plot: i) among cell types within same dataset; ii) among celltypes within different datasets (Reference and Validation on the same upsetplot)

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(filesstrings)
library(RColorBrewer)
library(GEOquery)
library(readr)

# Loading Seurat Objects
Sobj_Ref = readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Ref_dset/Reference_dataset_standart_seut_pipeline_NEW_noNormal_stemnesscomputed.rds")
dim(Sobj_Ref) # 27719 80682
head(Sobj_Ref@meta.data)
names(Sobj_Ref@meta.data)
bk_Sobj_Ref = Sobj_Ref

Sobj_Val = readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Val_onlyd_tumors_No_lynph_samples_yes_clusterization_yes_noNormal_stemness_without_normalsmp.rds")
dim(Sobj_Val) #21843 180315
head(Sobj_Val@meta.data)
names(Sobj_Val@meta.data)
bk_Sobj_Val = Sobj_Val
Sobj_Ref$celltype_minor

# Ridge Plot: temness~ clusters 
Idents(Sobj_Ref) <- "celltype_minor"
RidgePlot(Sobj_Ref, features = "noNormal_stemness", ncol = 1)

Idents(Sobj_Val) <- "Minor_celltypes_sctyped_only_3_celltype_ofInterest_10_markers"
RidgePlot(Sobj_Val, features = "noNormal_stemness", ncol = 1)


### Subset cell types
## Reference data "Cancer Basal SC"

# Subset
Idents(Sobj_Ref) = "celltype_minor"
Sobj_Ref_CB = subset(x = Sobj_Ref, idents = "Cancer Basal SC")

# Quantile 
# 0.4838806 stemness >= represent top 75 percentile
# 0.4003201 stemness =<  represent bottom 25 percentile
quantile(Sobj_Ref_CB@meta.data[Sobj_Ref_CB@meta.data$celltype_minor %in% "Cancer Basal SC",]$noNormal_stemness,na.rm = T,probs = c(0.25,0.50,0.75,1.0))

# Add label top and bottom stemness based on quantiles
Sobj_Ref_CB@meta.data$top_bottom_stemness = NA
Sobj_Ref_CB@meta.data[Sobj_Ref_CB@meta.data$noNormal_stemness >= 0.4838806,]$top_bottom_stemness = "top"
Sobj_Ref_CB@meta.data[Sobj_Ref_CB@meta.data$noNormal_stemness <= 0.4003201,]$top_bottom_stemness = "bottom"
Sobj_Ref_CB@meta.data[!Sobj_Ref_CB@meta.data$top_bottom_stemness %in% c("top", "bottom"), ]$top_bottom_stemness = "middle"

# FindMarkers (Differential Expression)
Idents(Sobj_Ref_CB) <- "top_bottom_stemness"
high_stemness_markers <- FindMarkers(Sobj_Ref_CB, ident.1 = "top", ident.2 = "bottom", min.pct = 0.25)

head(high_stemness_markers)
dim(high_stemness_markers) #1725    5
findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,] 
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #685   5
hist(findmarkers_sig$avg_log2FC)

findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"

save(findmarkers_sig, file = "/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Reference_Cancer_Basal_top_stemness_DE_genes.rda")



## Reference data "Cancer Cycling"

# Subset
Idents(Sobj_Ref) = "celltype_minor"
Sobj_Ref_CC = subset(x = Sobj_Ref, idents = "Cancer Cycling")

# Quantile 
# 0.5337746 stemness >= represent top 75 percentile
# 0.4442086 stemness =<  represent bottom 25 percentile
quantile(Sobj_Ref_CC@meta.data[Sobj_Ref_CC@meta.data$celltype_minor %in% "Cancer Cycling",]$noNormal_stemness,na.rm = T,probs = c(0.25,0.50,0.75,1.0))

# Add label top and bottom stemness based on quantiles
Sobj_Ref_CC@meta.data$top_bottom_stemness = NA
Sobj_Ref_CC@meta.data[Sobj_Ref_CC@meta.data$noNormal_stemness >= 0.5337746,]$top_bottom_stemness = "top"
Sobj_Ref_CC@meta.data[Sobj_Ref_CC@meta.data$noNormal_stemness <= 0.4442086,]$top_bottom_stemness = "bottom"
Sobj_Ref_CC@meta.data[!Sobj_Ref_CC@meta.data$top_bottom_stemness %in% c("top", "bottom"), ]$top_bottom_stemness = "middle"

# FindMarkers (Differential Expression)
Idents(Sobj_Ref_CC) <- "top_bottom_stemness"
high_stemness_markers <- FindMarkers(Sobj_Ref_CC, ident.1 = "top", ident.2 = "bottom", min.pct = 0.25)

head(high_stemness_markers)
dim(high_stemness_markers) # 1395    5
findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,] 
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #393   5
hist(findmarkers_sig$avg_log2FC)

findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"

save(findmarkers_sig, file = "/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Reference_Cancer_Cycling_top_stemness_DE_genes.rda")



## Reference data "Cycling T-cells"

# Subset
Idents(Sobj_Ref) = "celltype_minor"
Sobj_Ref_CT = subset(x = Sobj_Ref, idents = "Cycling T-cells")

# Quantile 
# 0.5447295 stemness >= represent top 75 percentile
# 0.4729583 stemness =<  represent bottom 25 percentile
quantile(Sobj_Ref_CT@meta.data[Sobj_Ref_CT@meta.data$celltype_minor %in% "Cycling T-cells",]$noNormal_stemness,na.rm = T,probs = c(0.25,0.50,0.75,1.0))

# Add label top and bottom stemness based on quantiles
Sobj_Ref_CT@meta.data$top_bottom_stemness = NA
Sobj_Ref_CT@meta.data[Sobj_Ref_CT@meta.data$noNormal_stemness >= 0.5447295,]$top_bottom_stemness = "top"
Sobj_Ref_CT@meta.data[Sobj_Ref_CT@meta.data$noNormal_stemness <= 0.4729583,]$top_bottom_stemness = "bottom"
Sobj_Ref_CT@meta.data[!Sobj_Ref_CT@meta.data$top_bottom_stemness %in% c("top", "bottom"), ]$top_bottom_stemness = "middle"

# FindMarkers (Differential Expression)
Idents(Sobj_Ref_CT) <- "top_bottom_stemness"
high_stemness_markers <- FindMarkers(Sobj_Ref_CT, ident.1 = "top", ident.2 = "bottom", min.pct = 0.25)

head(high_stemness_markers)
dim(high_stemness_markers) # 1329    5
findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,] 
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #230   5
hist(findmarkers_sig$avg_log2FC)

findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"

save(findmarkers_sig, file = "/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Reference_Cycling_Tcells_top_stemness_DE_genes.rda")




## Validation data "Cancer Basal SC"

# Subset
Idents(Sobj_Val) = "Minor_celltypes_sctyped_only_3_celltype_ofInterest_10_markers"
Sobj_Val_CB = subset(x = Sobj_Val, idents = "Cancer Basal SC")

# Quantile 
# 0.4580021 stemness >= represent top 75 percentile
# 0.3839707 stemness =<  represent bottom 25 percentile
quantile(Sobj_Val_CB@meta.data[Sobj_Val_CB@meta.data$Minor_celltypes_sctyped_only_3_celltype_ofInterest_10_markers %in% "Cancer Basal SC",]$noNormal_stemness,na.rm = T,probs = c(0.25,0.50,0.75,1.0))

# Add label top and bottom stemness based on quantiles
Sobj_Val_CB@meta.data$top_bottom_stemness = NA
Sobj_Val_CB@meta.data[Sobj_Val_CB@meta.data$noNormal_stemness >= 0.4580021,]$top_bottom_stemness = "top"
Sobj_Val_CB@meta.data[Sobj_Val_CB@meta.data$noNormal_stemness <= 0.3839707,]$top_bottom_stemness = "bottom"
Sobj_Val_CB@meta.data[!Sobj_Val_CB@meta.data$top_bottom_stemness %in% c("top", "bottom"), ]$top_bottom_stemness = "middle"

# Checking number os cells
table(Sobj_Val_CB@meta.data$top_bottom_stemness %in% "top") # 5039 top 
table(Sobj_Val_CB@meta.data$top_bottom_stemness %in% "bottom") # 5039 bottom 

# FindMarkers (Differential Expression)
Idents(Sobj_Val_CB) <- "top_bottom_stemness"
high_stemness_markers <- FindMarkers(Sobj_Val_CB, ident.1 = "top", ident.2 = "bottom", min.pct = 0.25)

head(high_stemness_markers)
dim(high_stemness_markers) # 1418    5
findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,] 
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #525   5
hist(findmarkers_sig$avg_log2FC)

findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"

save(findmarkers_sig, file = "/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Validation_Cancer_Basal_top_stemness_DE_genes.rda")






## Validation data "Cancer Cycling"

# Subset
Idents(Sobj_Val) = "Minor_celltypes_sctyped_only_3_celltype_ofInterest_10_markers"
Sobj_Val_CC = subset(x = Sobj_Val, idents = "Cancer Cycling")

# Quantile 
# 0.4519422 stemness >= represent top 75 percentile
# 0.3856006 stemness =<  represent bottom 25 percentile
quantile(Sobj_Val_CC@meta.data[Sobj_Val_CC@meta.data$Minor_celltypes_sctyped_only_3_celltype_ofInterest_10_markers %in% "Cancer Cycling",]$noNormal_stemness,na.rm = T,probs = c(0.25,0.50,0.75,1.0))

# Add label top and bottom stemness based on quantiles
Sobj_Val_CC@meta.data$top_bottom_stemness = NA
Sobj_Val_CC@meta.data[Sobj_Val_CC@meta.data$noNormal_stemness >= 0.4519422,]$top_bottom_stemness = "top"
Sobj_Val_CC@meta.data[Sobj_Val_CC@meta.data$noNormal_stemness <= 0.3856006,]$top_bottom_stemness = "bottom"
Sobj_Val_CC@meta.data[!Sobj_Val_CC@meta.data$top_bottom_stemness %in% c("top", "bottom"), ]$top_bottom_stemness = "middle"

# Checking number os cells
table(Sobj_Val_CC@meta.data$top_bottom_stemness %in% "top") # 3032 top 
table(Sobj_Val_CC@meta.data$top_bottom_stemness %in% "bottom") # 3032 bottom 

# FindMarkers (Differential Expression)
Idents(Sobj_Val_CC) <- "top_bottom_stemness"
high_stemness_markers <- FindMarkers(Sobj_Val_CC, ident.1 = "top", ident.2 = "bottom", min.pct = 0.25)

head(high_stemness_markers)
dim(high_stemness_markers) # 1616    5
findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,] 
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #528   5
hist(findmarkers_sig$avg_log2FC)

findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"

save(findmarkers_sig, file = "/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Validation_Cancer_Cycling_top_stemness_DE_genes.rda")







## Validation data "Cycling T-cells"

# Subset
Idents(Sobj_Val) = "Minor_celltypes_sctyped_only_3_celltype_ofInterest_10_markers"
Sobj_Val_CT = subset(x = Sobj_Val, idents = "Cycling T-cells")

# Quantile 
# 0.4917912 stemness >= represent top 75 percentile
# 0.4296064 stemness =<  represent bottom 25 percentile
quantile(Sobj_Val_CT@meta.data[Sobj_Val_CT@meta.data$Minor_celltypes_sctyped_only_3_celltype_ofInterest_10_markers %in% "Cycling T-cells",]$noNormal_stemness,na.rm = T,probs = c(0.25,0.50,0.75,1.0))

# Add label top and bottom stemness based on quantiles
Sobj_Val_CT@meta.data$top_bottom_stemness = NA
Sobj_Val_CT@meta.data[Sobj_Val_CT@meta.data$noNormal_stemness >= 0.4917912,]$top_bottom_stemness = "top"
Sobj_Val_CT@meta.data[Sobj_Val_CT@meta.data$noNormal_stemness <= 0.4296064,]$top_bottom_stemness = "bottom"
Sobj_Val_CT@meta.data[!Sobj_Val_CT@meta.data$top_bottom_stemness %in% c("top", "bottom"), ]$top_bottom_stemness = "middle"

# Checking number os cells
table(Sobj_Val_CT@meta.data$top_bottom_stemness %in% "top") # 3440 top 
table(Sobj_Val_CT@meta.data$top_bottom_stemness %in% "bottom") # 3440 bottom 

# FindMarkers (Differential Expression)
Idents(Sobj_Val_CT) <- "top_bottom_stemness"
high_stemness_markers <- FindMarkers(Sobj_Val_CT, ident.1 = "top", ident.2 = "bottom", min.pct = 0.25)

head(high_stemness_markers)
dim(high_stemness_markers) # 1634    5
findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,] 
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #632   5
hist(findmarkers_sig$avg_log2FC)

findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"

save(findmarkers_sig, file = "/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Validation_Cycling_Tcells_top_stemness_DE_genes.rda")



### Upsetplot 

# Gene intersected within indiviually datasets 
# Reference dataset
load("/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Reference_Cycling_Tcells_top_stemness_DE_genes.rda")
Reference_Cycling_Tcells = findmarkers_sig
Reference_Cycling_Tcells_UP = Reference_Cycling_Tcells[Reference_Cycling_Tcells$up_down_exp %in% "up", ]$genes
Reference_Cycling_Tcells_DOWN = Reference_Cycling_Tcells[Reference_Cycling_Tcells$up_down_exp %in% "down", ]$genes


load("/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Reference_Cancer_Cycling_top_stemness_DE_genes.rda")
Reference_Cancer_Cycling = findmarkers_sig
Reference_Cancer_Cycling_UP = Reference_Cancer_Cycling[Reference_Cancer_Cycling$up_down_exp %in% "up", ]$genes
Reference_Cancer_Cycling_DOWN = Reference_Cancer_Cycling[Reference_Cancer_Cycling$up_down_exp %in% "down", ]$genes


load("/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Reference_Cancer_Basal_top_stemness_DE_genes.rda")
Reference_Cancer_Basal = findmarkers_sig
Reference_Cancer_Basal_UP = Reference_Cancer_Basal[Reference_Cancer_Basal$up_down_exp %in% "up", ]$genes
Reference_Cancer_Basal_DOWN = Reference_Cancer_Basal[Reference_Cancer_Basal$up_down_exp %in% "down", ]$genes


library(UpSetR)
listInput <- list(Reference_Cycling_Tcells_UP = Reference_Cycling_Tcells_UP,
                  Reference_Cycling_Tcells_DOWN = Reference_Cycling_Tcells_DOWN,
                  Reference_Cancer_Cycling_UP = Reference_Cancer_Cycling_UP, 
                  Reference_Cancer_Cycling_DOWN = Reference_Cancer_Cycling_DOWN,
                  Reference_Cancer_Basal_UP = Reference_Cancer_Basal_UP, 
                  Reference_Cancer_Basal_DOWN = Reference_Cancer_Basal_DOWN)

upset(nsets = 6, fromList(listInput), order.by = "freq")

# Validation dataset
load("/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Validation_Cycling_Tcells_top_stemness_DE_genes.rda")
Validation_Cycling_Tcells = findmarkers_sig
Validation_Cycling_Tcells_UP = Validation_Cycling_Tcells[Validation_Cycling_Tcells$up_down_exp %in% "up", ]$genes
Validation_Cycling_Tcells_DOWN = Validation_Cycling_Tcells[Validation_Cycling_Tcells$up_down_exp %in% "down", ]$genes


load("/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Validation_Cancer_Cycling_top_stemness_DE_genes.rda")
Validation_Cancer_Cycling = findmarkers_sig
Validation_Cancer_Cycling_UP = Validation_Cancer_Cycling[Validation_Cancer_Cycling$up_down_exp %in% "up", ]$genes
Validation_Cancer_Cycling_DOWN = Validation_Cancer_Cycling[Validation_Cancer_Cycling$up_down_exp %in% "down", ]$genes


load("/mnt/plummergrp/maycon/Figures_to_paper/DE_top_bottom_stemness/Validation_Cancer_Basal_top_stemness_DE_genes.rda")
Validation_Cancer_Basal = findmarkers_sig
Validation_Cancer_Basal_UP = Validation_Cancer_Basal[Validation_Cancer_Basal$up_down_exp %in% "up", ]$genes
Validation_Cancer_Basal_DOWN = Validation_Cancer_Basal[Validation_Cancer_Basal$up_down_exp %in% "down", ]$genes

library(UpSetR)
listInput <- list(Validation_Cycling_Tcells_UP = Validation_Cycling_Tcells_UP,
                  Validation_Cycling_Tcells_DOWN = Validation_Cycling_Tcells_DOWN,
                  Validation_Cancer_Cycling_UP = Validation_Cancer_Cycling_UP, 
                  Validation_Cancer_Cycling_DOWN = Validation_Cancer_Cycling_DOWN,
                  Validation_Cancer_Basal_UP = Validation_Cancer_Basal_UP, 
                  Validation_Cancer_Basal_DOWN = Validation_Cancer_Basal_DOWN)

upset(nsets = 6, fromList(listInput), order.by = "freq")



# Reference dataset & Validation dataset
library(UpSetR)
listInput <- list(Reference_Cycling_Tcells_UP = Reference_Cycling_Tcells_UP,
                  Reference_Cycling_Tcells_DOWN = Reference_Cycling_Tcells_DOWN,
                  Reference_Cancer_Cycling_UP = Reference_Cancer_Cycling_UP, 
                  Reference_Cancer_Cycling_DOWN = Reference_Cancer_Cycling_DOWN,
                  Reference_Cancer_Basal_UP = Reference_Cancer_Basal_UP, 
                  Reference_Cancer_Basal_DOWN = Reference_Cancer_Basal_DOWN,
                  Validation_Cycling_Tcells_UP = Validation_Cycling_Tcells_UP,
                  Validation_Cycling_Tcells_DOWN = Validation_Cycling_Tcells_DOWN,
                  Validation_Cancer_Cycling_UP = Validation_Cancer_Cycling_UP, 
                  Validation_Cancer_Cycling_DOWN = Validation_Cancer_Cycling_DOWN,
                  Validation_Cancer_Basal_UP = Validation_Cancer_Basal_UP, 
                  Validation_Cancer_Basal_DOWN = Validation_Cancer_Basal_DOWN)

upset(nsets = 12, fromList(listInput), order.by = "freq")

