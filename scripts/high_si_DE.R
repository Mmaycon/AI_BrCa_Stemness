

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Differential Expression on MAJOR cell types within stemness > 0.7055146 ###

Sobj = readRDS("/mnt/plummergrp/maycon/Reference_dataset_standart_seut_pipeline.rds")
Idents(Sobj) = "celltype_major"
VlnPlot(Sobj, features = "stemness")
Idents(Sobj) = "celltype_minor"
VlnPlot(Sobj, features = "stemness")

# Let's take the iPS bottom 25 percentile as a cutoff for high stemness index
# In this case, out cutoff would be > 0.7055146 of stemness


  # High stemness cell >= 0.7055146 stemness
  load("/mnt/storage1/mmaycon/After_iPS_dataset/meta.data_CNV_CellCycleScore.rda")
  quantile(pdata[pdata$celltype_major %in% "iPS",]$stemness,na.rm = T,probs = c(0.25,0.50,0.75,1.0))

  # There're 2 celltypes with enough cells above 0.7055146 stemness 
  # i) Cancer epithelial with 708 cells
  # ii) T-cells with 365 cells
  table(Sobj@meta.data[Sobj@meta.data$stemness >= 0.7055146,]$celltype_major)

  # Bottom 25 percentile stemness cells < 0.5109824 (T-cells)
  quantile(Sobj@meta.data[Sobj@meta.data$celltype_major %in% "T-cells",]$stemness, na.rm = T,probs = c(0.25,0.50,0.75,1.0))
  # Bottom 25 percentile stemness cells < 0.4647798 (Cancer Epithelial)
  quantile(Sobj@meta.data[Sobj@meta.data$celltype_major %in% "Cancer     Epithelial",]$stemness, na.rm = T,probs = c(0.25,0.50,0.75,1.0))

  # Label high and low stemness based on iPS cutoff 0.7055146 
  # - it will be for group indicator on Differential Expression
  Sobj@meta.data$high_low_stemness = NA
  Sobj@meta.data[Sobj@meta.data$stemness >= 0.7055146,]$high_low_stemness = "high"
  Sobj@meta.data[Sobj@meta.data$stemness < 0.7055146,]$high_low_stemness = "low"


### T-cells subset for Differential Expression (DE) on high stemness cells 

  # Checking identical cellbarcodes and colnames on seurat object
  dim(Sobj) #27719 80682
  length(intersect(Sobj@meta.data$cell_barcode, colnames(Sobj))) #80682 (ok)

  # Filtering T-cellss to DEG of high stemness vs low stemness cells 
  # high stemness T-cells (>= 0.7055146) and low bottom 25% stemness T-cells (< 0.5109824)
Sobje_metadata = as.data.frame(Sobj@meta.data)
Sobje_metadata = Sobje_metadata[Sobje_metadata$celltype_major %in% "T-cells",]
Sobje_metadata_hisi = Sobje_metadata[Sobje_metadata$stemness >= 0.7055146,]
hist(Sobje_metadata_hisi$stemness)

Sobje_metadata_lowsi = Sobje_metadata[Sobje_metadata$stemness < 0.5109824,]
hist(Sobje_metadata_lowsi$stemness)

  # Downsampling low stemness cells
  # 365 HIGH stemness Tcells
  dim(Sobje_metadata_hisi) #365  25
  # Downsample on low stemness Tcells around 365 cells
  lowsi_cells = Sobje_metadata_lowsi$cell_barcode
  Idents(Sobj) <- "Patient_ID"
  Sobje_lowsi = Sobj[,colnames(Sobj) %in% lowsi_cells]
  Idents(Sobje_lowsi) <- "Patient_ID" # donwsampling over patients
  Sobje_lowsi_downspm <- subset(x = Sobje_lowsi, downsample = 14)
  print(dim(Sobje_lowsi_downspm)) # 364 cells (ok)
  # Subset Sobj for high stemness samples 
  highsi_cells = Sobje_metadata_hisi$cell_barcode
  Idents(Sobj) <- "Patient_ID"
  Sobje_highsi = Sobj[,colnames(Sobj) %in% highsi_cells]
  print(dim(Sobje_highsi)) 
  # Merge high stemness cells with bottom stemness dowmsampled cells
  Sobj_Tcells_stemness = merge(Sobje_highsi,Sobje_lowsi_downspm)
    # " Warning message: Some cell names are duplicated across objects provided. Renaming to enforce unique cell names "
    # Is that a problem?
  
  # Differential expression on high stemness cells
  Idents(Sobj_Tcells_stemness) <- "high_low_stemness"
high_stemness_markers <- FindMarkers(Sobj_Tcells_stemness, ident.1 = "high", ident.2 = "low", min.pct = 0.25)
  # Filtering Differential Expression output
  head(high_stemness_markers)
  dim(high_stemness_markers) #2415   5
  findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,] 
  findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
  dim(findmarkers_sig) #518   5
  hist(findmarkers_sig$avg_log2FC)
  
  findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
  findmarkers_sig$genes <- rownames(findmarkers_sig)
  findmarkers_sig$up_down_exp = NA
  findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
  findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"
  
  save(findmarkers_sig, file = "/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/tcells_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")

  

### Cancer Epithelial subset for Differential Expression (DE) 
### on high stemness cells
  
  # Checking identical cellbarcodes and colnames on seurat object
  dim(Sobj) #27719 80682
  length(intersect(Sobj@meta.data$cell_barcode, colnames(Sobj))) #80682 (ok)

  # Filtering Cancer Epithelial to DEG of high stemness vs low stemness cells 
  # high stemness Cancer Epithelial (>= 0.7055146) and low bottom 25% stemness Cancer Epithelial (< 0.4647798)
  Sobje_metadata = as.data.frame(Sobj@meta.data)
  Sobje_metadata = Sobje_metadata[Sobje_metadata$celltype_major %in% "Cancer Epithelial",]
  Sobje_metadata_hisi = Sobje_metadata[Sobje_metadata$stemness >= 0.7055146,]
  hist(Sobje_metadata_hisi$stemness)
  
  Sobje_metadata_lowsi = Sobje_metadata[Sobje_metadata$stemness < 0.4647798,]
  hist(Sobje_metadata_lowsi$stemness)

  # Downsampling low stemness cells
  # 708 HIGH stemness Cancer Epithelial
  dim(Sobje_metadata[Sobje_metadata$stemness >= 0.7055146,]) 
  # Downsample low stemness Cancer Epithelial around 708 cells
  lowsi_cells = Sobje_metadata_lowsi$cell_barcode
  Idents(Sobj) <- "Patient_ID"
  Sobje_lowsi = Sobj[,colnames(Sobj) %in% lowsi_cells]
  Idents(Sobje_lowsi) <- "Patient_ID" # donwsampling over patients
  Sobje_lowsi_downspm <- subset(x = Sobje_lowsi, downsample = 65)
  print(dim(Sobje_lowsi_downspm)) # 706 cells (ok)
  # Subset Sobj for high stemness samples 
  highsi_cells = Sobje_metadata_hisi$cell_barcode
  Idents(Sobj) <- "Patient_ID"
  Sobje_highsi = Sobj[,colnames(Sobj) %in% highsi_cells]
  print(dim(Sobje_highsi)) 
  # Merge high stemness cells with bottom stemness dowmsampled cells
  Sobj_Cancer_Epithelial_stemness = merge(Sobje_highsi,Sobje_lowsi_downspm)
  # " Warning message: Some cell names are duplicated across objects provided. Renaming to enforce unique cell names "
  # Is that a problem?

  table(Sobj_Cancer_Epithelial_stemness$high_low_stemness)
  Idents(Sobj_Cancer_Epithelial_stemness) = "celltype_major"
  VlnPlot(Sobj_Cancer_Epithelial_stemness, features = "stemness")
  
  # Differential expression on high stemness cells
Idents(Sobj_Cancer_Epithelial_stemness) <- "high_low_stemness"
high_stemness_markers <- FindMarkers(Sobj_Cancer_Epithelial_stemness, ident.1 = "high", ident.2 = "low", min.pct = 0.25)

head(high_stemness_markers)
dim(high_stemness_markers) #1256   5
findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,] 
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #256   5
hist(findmarkers_sig$avg_log2FC)

findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"

save(findmarkers_sig, file = "/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Epithelial_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")


### Find the purest gene markers among top 10 DE genes

  ### T-cells 
  # DE output
  load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/tcells_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")
  # 10 top DE genes on high vs low stemness T-cells
  test_genes = findmarkers_sig$genes[1:10] 

# Attention: scatter plot need to be done with Sobj_Tcells_stemness because it has only selected high and low stemness cells
  f1 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[1]); f1 # +
  f2 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[2]); f2 # + -
  f3 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[3]); f3 # + +
  f4 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[4]); f4 # + +
  f5 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[5]); f5 # -
  f6 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[6]); f6 # -
  f7 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[7]); f7 # -
  f8 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[8]); f8 # + +
  f9 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[9]); f9 # +
  f10 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[10]); f10 # +

f1+ f2+ f3+ f4+ f5+ f6+ f7+ f8+ f9+ f10
  
# Attention: feature plot need to be done with Sobj because it has all cells. We wish to view those markers on the entire umap 
FeaturePlot(object = Sobj, features = test_genes)

Idents(Sobj) = "celltype_major"
Sobj <- ScaleData(object = Sobj, features = rownames(Sobj))
test_genes = findmarkers_sig$genes[1:10] 
DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)
test_genes = findmarkers_sig$genes[1:20] 
DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)
  # top 10 DE up genes, top 10 DE down genes
test_genes = findmarkers_sig$genes[c(1:10,508:518)] 
DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)


# Top DOWN regulated genes 
test_genes = findmarkers_sig$genes[c(508:518)] 
f1 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[1]); f1 # +
f2 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[2]); f2 # + -
f3 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[3]); f3 # + +
f4 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[4]); f4 # + +
f5 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[5]); f5 # -
f6 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[6]); f6 # -
f7 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[7]); f7 # -
f8 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[8]); f8 # + +
f9 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[9]); f9 # +
f10 = FeatureScatter(object = Sobj_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[10]); f10 # +

f1+ f2+ f3+ f4+ f5+ f6+ f7+ f8+ f9+ f10

test_genes = findmarkers_sig$genes[c(508:518)] 
Idents(Sobj) = "celltype_major"
FeaturePlot(object = Sobj, features = test_genes)
DotPlot(object = Sobj, features = test_genes)


  ### Cancer Epithelial 
  # DE output
  load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Epithelial_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")
  # 10 top DE genes on high vs low stemness Cancer Epithelial cells
  test_genes = findmarkers_sig$genes[1:10] 

  # Attention: scatter plot need to be done with Sobj_Cancer_Epithelial_stemness because it has only selected high and low stemness cells
  f1 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[1]); f1 # +
  f2 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[2]); f2 # + -
  f3 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[3]); f3 # + +
  f4 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[4]); f4 # + +
  f5 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[5]); f5 # -
  f6 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[6]); f6 # -
  f7 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[7]); f7 # -
  f8 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[8]); f8 # + +
  f9 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[9]); f9 # +
  f10 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[10]); f10 # +
  
f1+ f2+ f3+ f4+ f5+ f6+ f7+ f8+ f9+ f10
  
  # Attention: feature plot need to be done with Sobj because it has all cells. We wish to view those markers on the entire umap 
  FeaturePlot(object = Sobj, features = test_genes)
  DotPlot(object = Sobj, features = test_genes)
  
  Idents(Sobj) = "celltype_major"
  #Sobj <- ScaleData(object = Sobj, features = rownames(Sobj))
  test_genes = findmarkers_sig$genes[1:10] 
  DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)
  test_genes = findmarkers_sig$genes[1:20] 
  DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)
  # top 10 DE up genes, top 10 DE down genes
  test_genes = findmarkers_sig$genes[c(1:10,246:256)] 
  DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)
  
  
  # Top DOWN regulated genes 
  test_genes = findmarkers_sig$genes[c(246:256)] 
  f1 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[1]); f1 # 
  f2 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[2]); f2 # 
  f3 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[3]); f3 # 
  f4 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[4]); f4 # 
  f5 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[5]); f5 # 
  f6 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[6]); f6 # 
  f7 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[7]); f7 # 
  f8 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[8]); f8 # 
  f9 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[9]); f9 # 
  f10 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[10]); f10 # 
  
f1+ f2+ f3+ f4+ f5+ f6+ f7+ f8+ f9+ f10


test_genes = findmarkers_sig$genes[c(246:256)] 
Idents(Sobj) = "celltype_major"
FeaturePlot(object = Sobj, features = test_genes)
DotPlot(object = Sobj, features = test_genes)


  
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Differential Expression on MINOR cell types within stemness > 0.7055146 ###

  # Stemness cutoff == 0.7055146 (based on sciPS stemness index)
  # minor cell types within reasonable quantity of cells > 0.7055146 stemness 
    # 3 minor cell types: "Cancer Cycling" (519 cells), "Cycling T-cells" (347 cells), 
    # and "Cancer Basal SC" (173 cells)
table(Sobj@meta.data[Sobj@meta.data$high_low_stemness %in% "high",]$celltype_minor)

  # Bottom stemness cutoff for the 3 minor cell types 
  # (low stemness index cells to be compared on DE)
    # Bottom stemness cutoff on Cancer Cycling == 0.5734144 stemness
quantile(Sobj@meta.data[Sobj@meta.data$celltype_minor %in% "Cancer Cycling",]$stemness, na.rm = T,probs = c(0.25,0.50,0.75,1.0))
    # Bottom stemness cutoff on Cycling T-cells == 0.6098451 stemness
quantile(Sobj@meta.data[Sobj@meta.data$celltype_minor %in% "Cycling T-cells",]$stemness, na.rm = T,probs = c(0.25,0.50,0.75,1.0))
    # Bottom stemness cutoff on Cancer Basal SC == 0.5167555 stemness
quantile(Sobj@meta.data[Sobj@meta.data$celltype_minor %in% "Cancer Basal SC",]$stemness, na.rm = T,probs = c(0.25,0.50,0.75,1.0))


## Cancer Cycling 
  # Filtering Cancer Cycling to DEG of high stemness vs low stemness cells 
Sobje_metadata = as.data.frame(Sobj@meta.data)
Sobje_metadata = Sobje_metadata[Sobje_metadata$celltype_minor %in% "Cancer Cycling",]
Sobje_metadata_hisi = Sobje_metadata[Sobje_metadata$stemness >= 0.7055146,]
hist(Sobje_metadata_hisi$stemness)

Sobje_metadata_lowsi = Sobje_metadata[Sobje_metadata$stemness < 0.5734144,]
hist(Sobje_metadata_lowsi$stemness)

# Downsampling low stemness cells
# 519 high stemness Cancer Cycling
dim(Sobje_metadata[Sobje_metadata$stemness >= 0.7055146,]) 
lowsi_cells = Sobje_metadata_lowsi$cell_barcode
Idents(Sobj) <- "Patient_ID"
Sobje_lowsi = Sobj[,colnames(Sobj) %in% lowsi_cells]
Idents(Sobje_lowsi) <- "Patient_ID" # donwsampling over patients
Sobje_lowsi_downspm <- subset(x = Sobje_lowsi, downsample = 90)
print(dim(Sobje_lowsi_downspm)) # 520 cells (ok)
# Subset Sobj for high stemness samples 
highsi_cells = Sobje_metadata_hisi$cell_barcode
Idents(Sobj) <- "Patient_ID"
Sobje_highsi = Sobj[,colnames(Sobj) %in% highsi_cells]
print(dim(Sobje_highsi)) 
# Merge high stemness cells with bottom stemness dowmsampled cells
Sobj_Cancer_Cycling_stemness = merge(Sobje_highsi,Sobje_lowsi_downspm)

table(Sobj_Cancer_Cycling_stemness$high_low_stemness)
Idents(Sobj_Cancer_Cycling_stemness) = "celltype_minor"
VlnPlot(Sobj_Cancer_Cycling_stemness, features = "stemness")

  # Differential expression on high stemness cells
  # Label "high_low_stemness" was prepared on the beginning of this script
Idents(Sobj_Cancer_Cycling_stemness) <- "high_low_stemness"
high_stemness_markers <- FindMarkers(Sobj_Cancer_Cycling_stemness, ident.1 = "high", ident.2 = "low", min.pct = 0.25)

head(high_stemness_markers)
dim(high_stemness_markers) #693   5
findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,] 
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #69   5
hist(findmarkers_sig$avg_log2FC)

findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"

save(findmarkers_sig, file = "/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Cycling_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")


## Cycling Tcells 
# Filtering Cycling Tcells to DEG of high stemness vs low stemness cells 
Sobje_metadata = as.data.frame(Sobj@meta.data)
Sobje_metadata = Sobje_metadata[Sobje_metadata$celltype_minor %in% "Cycling T-cells",]
Sobje_metadata_hisi = Sobje_metadata[Sobje_metadata$stemness >= 0.7055146,]
hist(Sobje_metadata_hisi$stemness)

Sobje_metadata_lowsi = Sobje_metadata[Sobje_metadata$stemness < 0.6098451,]
hist(Sobje_metadata_lowsi$stemness)

# Downsampling low stemness cells
# 347 high stemness Cycling Tcells
dim(Sobje_metadata[Sobje_metadata$stemness >= 0.7055146,]) 
lowsi_cells = Sobje_metadata_lowsi$cell_barcode
Idents(Sobj) <- "Patient_ID"
Sobje_lowsi = Sobj[,colnames(Sobj) %in% lowsi_cells]
Idents(Sobje_lowsi) <- "Patient_ID" # donwsampling over patients
Sobje_lowsi_downspm <- subset(x = Sobje_lowsi, downsample = 50)
print(dim(Sobje_lowsi_downspm)) # 343 cells (ok)
# Subset Sobj for high stemness samples 
highsi_cells = Sobje_metadata_hisi$cell_barcode
Idents(Sobj) <- "Patient_ID"
Sobje_highsi = Sobj[,colnames(Sobj) %in% highsi_cells]
print(dim(Sobje_highsi)) 
# Merge high stemness cells with bottom stemness dowmsampled cells
Sobj_Cycling_Tcells_stemness = merge(Sobje_highsi,Sobje_lowsi_downspm)

table(Sobj_Cycling_Tcells_stemness$high_low_stemness)
Idents(Sobj_Cycling_Tcells_stemness) = "celltype_minor"
VlnPlot(Sobj_Cycling_Tcells_stemness, features = "stemness")


# Differential expression on high stemness cells
# Label "high_low_stemness" was prepared on the beginning of this script
Idents(Sobj_Cycling_Tcells_stemness) <- "high_low_stemness"
high_stemness_markers <- FindMarkers(Sobj_Cycling_Tcells_stemness, ident.1 = "high", ident.2 = "low", min.pct = 0.25)

head(high_stemness_markers)
dim(high_stemness_markers) #873   5
findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,]
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #68   5
hist(findmarkers_sig$avg_log2FC)

findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"
dim(findmarkers_sig) #68   7

save(findmarkers_sig, file = "/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cycling_Tcells_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")


## Cancer_Basal_SC 
# Filtering Cancer_Basal_SC to DEG of high stemness vs low stemness cells 
Sobje_metadata = as.data.frame(Sobj@meta.data)
Sobje_metadata = Sobje_metadata[Sobje_metadata$celltype_minor %in% "Cancer Basal SC",]
Sobje_metadata_hisi = Sobje_metadata[Sobje_metadata$stemness >= 0.7055146,]
hist(Sobje_metadata_hisi$stemness)

Sobje_metadata_lowsi = Sobje_metadata[Sobje_metadata$stemness < 0.5167555,]
hist(Sobje_metadata_lowsi$stemness)

# Downsampling low stemness cells
# 173 high stemness Cancer_Basal_SC
dim(Sobje_metadata[Sobje_metadata$stemness >= 0.7055146,]) 
lowsi_cells = Sobje_metadata_lowsi$cell_barcode
Idents(Sobj) <- "Patient_ID"
Sobje_lowsi = Sobj[,colnames(Sobj) %in% lowsi_cells]
Idents(Sobje_lowsi) <- "Patient_ID" # donwsampling over patients
Sobje_lowsi_downspm <- subset(x = Sobje_lowsi, downsample = 20)
print(dim(Sobje_lowsi_downspm)) # 168 cells (ok)
# Subset Sobj for high stemness samples 
highsi_cells = Sobje_metadata_hisi$cell_barcode
Idents(Sobj) <- "Patient_ID"
Sobje_highsi = Sobj[,colnames(Sobj) %in% highsi_cells]
print(dim(Sobje_highsi)) 
# Merge high stemness cells with bottom stemness dowmsampled cells
Sobj_Cancer_Basal_SC_stemness = merge(Sobje_highsi,Sobje_lowsi_downspm)

table(Sobj_Cancer_Basal_SC_stemness$high_low_stemness)
Idents(Sobj_Cancer_Basal_SC_stemness) = "celltype_minor"
VlnPlot(Sobj_Cancer_Basal_SC_stemness, features = "stemness")

# Differential expression on high stemness cells
# Label "high_low_stemness" was prepared on the beginning of this script
Idents(Sobj_Cancer_Basal_SC_stemness) <- "high_low_stemness"
high_stemness_markers <- FindMarkers(Sobj_Cancer_Basal_SC_stemness, ident.1 = "high", ident.2 = "low", min.pct = 0.25)

head(high_stemness_markers)
dim(high_stemness_markers) #1039   5
findmarkers_sig <- high_stemness_markers[high_stemness_markers$p_val_adj < 0.05,]
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #173   5
hist(findmarkers_sig$avg_log2FC)

findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"
dim(findmarkers_sig) #173   7

save(findmarkers_sig, file = "/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Basal_SC_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")




### Gene makers visualization 

  ## Cancer Cycling 
load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Cycling_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")
# 10 top DE genes on high vs low stemness Cancer_Cycling cells
test_genes = findmarkers_sig$genes[1:10] 

f1 = FeatureScatter(object = Sobj_Cancer_Cycling_stemness, feature1 = "stemness", feature2 = test_genes[1]); f1 # 
f2 = FeatureScatter(object = Sobj_Cancer_Cycling_stemness, feature1 = "stemness", feature2 = test_genes[2]); f2 # 
f3 = FeatureScatter(object = Sobj_Cancer_Cycling_stemness, feature1 = "stemness", feature2 = test_genes[3]); f3 # 
f4 = FeatureScatter(object = Sobj_Cancer_Cycling_stemness, feature1 = "stemness", feature2 = test_genes[4]); f4 #
f5 = FeatureScatter(object = Sobj_Cancer_Cycling_stemness, feature1 = "stemness", feature2 = test_genes[5]); f5 # 
f6 = FeatureScatter(object = Sobj_Cancer_Cycling_stemness, feature1 = "stemness", feature2 = test_genes[6]); f6 # 
f7 = FeatureScatter(object = Sobj_Cancer_Cycling_stemness, feature1 = "stemness", feature2 = test_genes[7]); f7 # 
f8 = FeatureScatter(object = Sobj_Cancer_Cycling_stemness, feature1 = "stemness", feature2 = test_genes[8]); f8 # 
f9 = FeatureScatter(object = Sobj_Cancer_Cycling_stemness, feature1 = "stemness", feature2 = test_genes[9]); f9 # 
f10 = FeatureScatter(object = Sobj_Cancer_Cycling_stemness, feature1 = "stemness", feature2 = test_genes[10]); f10 # 

f1+ f2+ f3+ f4+ f5+ f6+ f7+ f8+ f9+ f10

# Attention: feature plot need to be done with Sobj because it has all cells. We wish to view those markers on the entire umap 
Idents(Sobj) = "celltype_minor"
FeaturePlot(object = Sobj, features = test_genes)
DotPlot(object = Sobj, features = test_genes)

Idents(Sobj) = "celltype_minor"
#Sobj <- ScaleData(object = Sobj, features = rownames(Sobj))
test_genes = findmarkers_sig$genes[1:10] 
DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)
test_genes = findmarkers_sig$genes[1:20] 
DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)
# top 10 DE up genes, top 10 DE down genes
test_genes = findmarkers_sig$genes[c(1:10,59:69)] 
DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)


# Top DOWN regulated genes 
test_genes = findmarkers_sig$genes[c(59:69)] 
f1 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[1]); f1 # 
f2 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[2]); f2 # 
f3 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[3]); f3 # 
f4 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[4]); f4 # 
f5 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[5]); f5 # 
f6 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[6]); f6 # 
f7 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[7]); f7 # 
f8 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[8]); f8 # 
f9 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[9]); f9 # 
f10 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[10]); f10 # 

f1+ f2+ f3+ f4+ f5+ f6+ f7+ f8+ f9+ f10


test_genes = findmarkers_sig$genes[c(59:69)] 
Idents(Sobj) = "celltype_minor"
FeaturePlot(object = Sobj, features = test_genes)
DotPlot(object = Sobj, features = test_genes)


## Cycling T-cells
load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cycling_Tcells_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")

# 10 top DE genes on high vs low stemness Cycling_Tcells cells
test_genes = findmarkers_sig$genes[1:10] 

f1 = FeatureScatter(object = Sobj_Cycling_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[1]); f1 # 
f2 = FeatureScatter(object = Sobj_Cycling_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[2]); f2 # 
f3 = FeatureScatter(object = Sobj_Cycling_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[3]); f3 # 
f4 = FeatureScatter(object = Sobj_Cycling_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[4]); f4 #
f5 = FeatureScatter(object = Sobj_Cycling_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[5]); f5 # 
f6 = FeatureScatter(object = Sobj_Cycling_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[6]); f6 # 
f7 = FeatureScatter(object = Sobj_Cycling_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[7]); f7 # 
f8 = FeatureScatter(object = Sobj_Cycling_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[8]); f8 # 
f9 = FeatureScatter(object = Sobj_Cycling_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[9]); f9 # 
f10 = FeatureScatter(object = Sobj_Cycling_Tcells_stemness, feature1 = "stemness", feature2 = test_genes[10]); f10 # 

f1+ f2+ f3+ f4+ f5+ f6+ f7+ f8+ f9+ f10

# Attention: feature plot need to be done with Sobj because it has all cells. We wish to view those markers on the entire umap 
Idents(Sobj) = "celltype_minor"
FeaturePlot(object = Sobj, features = test_genes)
DotPlot(object = Sobj, features = test_genes)

Idents(Sobj) = "celltype_minor"
#Sobj <- ScaleData(object = Sobj, features = rownames(Sobj))
test_genes = findmarkers_sig$genes[1:10] 
DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)
test_genes = findmarkers_sig$genes[1:20] 
DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)
# top 10 DE up genes, top 10 DE down genes
test_genes = findmarkers_sig$genes[c(1:10,58:68)] 
DoHeatmap(subset(Sobj, downsample = 100), features = test_genes, size = 3)


# Top DOWN regulated genes 
test_genes = findmarkers_sig$genes[c(58:68)] 
f1 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[1]); f1 # 
f2 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[2]); f2 # 
f3 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[3]); f3 # 
f4 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[4]); f4 # 
f5 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[5]); f5 # 
f6 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[6]); f6 # 
f7 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[7]); f7 # 
f8 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[8]); f8 # 
f9 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[9]); f9 # 
f10 = FeatureScatter(object = Sobj_Cancer_Epithelial_stemness, feature1 = "stemness", feature2 = test_genes[10]); f10 # 

f1+ f2+ f3+ f4+ f5+ f6+ f7+ f8+ f9+ f10


test_genes = findmarkers_sig$genes[c(58:68)] 
Idents(Sobj) = "celltype_minor"
FeaturePlot(object = Sobj, features = test_genes)
DotPlot(object = Sobj, features = test_genes)



### Intersect DE genes across Celltypes ###
load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/tcells_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")
tcells = findmarkers_sig$genes
tcells_up = findmarkers_sig[findmarkers_sig$up_down_exp %in% "up",]$genes
tcells_down = findmarkers_sig[findmarkers_sig$up_down_exp %in% "down",]$genes

load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Epithelial_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")
cancer_epi = findmarkers_sig$genes
cancer_epi_up = findmarkers_sig[findmarkers_sig$up_down_exp %in% "up",]$genes
cancer_epi_down = findmarkers_sig[findmarkers_sig$up_down_exp %in% "down",]$genes

load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Cycling_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")
cancer_cyc = findmarkers_sig$genes
cancer_cyc_up = findmarkers_sig[findmarkers_sig$up_down_exp %in% "up",]$genes
cancer_cyc_down = findmarkers_sig[findmarkers_sig$up_down_exp %in% "down",]$genes

load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cycling_Tcells_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")
cyc_tcells = findmarkers_sig$genes
cyc_tcells_up = findmarkers_sig[findmarkers_sig$up_down_exp %in% "up",]$genes
cyc_tcells_down = findmarkers_sig[findmarkers_sig$up_down_exp %in% "down",]$genes

load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Basal_SC_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")
cancer_bas = findmarkers_sig$genes
cancer_bas_up = findmarkers_sig[findmarkers_sig$up_down_exp %in% "up",]$genes
cancer_bas_down = findmarkers_sig[findmarkers_sig$up_down_exp %in% "down",]$genes


library(UpSetR)
listInput <- list(tcells = tcells, cancer_epi = cancer_epi, cancer_cyc = cancer_cyc, cyc_tcells = cyc_tcells, cancer_bas = cancer_bas)
upset(fromList(listInput), order.by = "freq")

library(UpSetR)
listInput <- list(tcells = tcells_up, cancer_epi = cancer_epi_up, cancer_cyc = cancer_cyc_up, cyc_tcells = cyc_tcells_up, cancer_bas = cancer_bas_up)
upset(fromList(listInput), order.by = "freq")

library(UpSetR)
listInput <- list(tcells = tcells_down, cancer_epi = cancer_epi_down, cancer_cyc = cancer_cyc_down, cyc_tcells = cyc_tcells_down, cancer_bas = cancer_bas_down)
upset(fromList(listInput), order.by = "freq")


