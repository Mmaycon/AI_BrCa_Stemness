# Validation dataset GSE161529
# Contains breast cancer and other conditions (e.g normal breast)
library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(filesstrings)
library(RColorBrewer)
library(GEOquery)
library(readr)

### Donwload dataset from GEO GSE161529
### Prepare directory 
## The following code was used to store the three files (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz) in a folder for each sample - maybe not be working anymore 


## Put each sample files (1 barcode, 1 feature, 1 matrix) into one folder by sample

# List files
all_files <- list.files("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529", pattern = ".gz$", full = TRUE, recursive = TRUE) # list all the files
matrix_files = all_files[grepl(c("-matrix.mtx.gz"),all_files)]
barcode_files = all_files[grepl(c("-barcodes.tsv.gz"),all_files)]
feature_files = all_files[grepl(c("-features.tsv.gz"),all_files)]
files_paths = c(matrix_files, barcode_files, feature_files)

# Keep only file names
files_splited <- str_split_fixed (as.character(files_paths), "/mnt/storage1/mmaycon/downloads/Valid_dset_GSE161529/", 2)[,2] 
# Get sample name of it 
sample_name = substring(files_splited, 1, 10)
sample_name = sample_name[!duplicated(sample_name)]
length(sample_name) # 69, okay!

# Create directory to store one sample per folder
dir.create("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/all_samples")

# Loop to create 69 sub directories (69 = number of samples)
for(i in 1:length(sample_name)) { 
  dir.create(paste0("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/all_samples/", sample_name[i]))
} 

# Loop to copy GSE161529_features.tsv.gz to all subdirectories - there is only one feature.tsv because it is the same for every sample
# Coping the GSE161529_features.tsv.gz into all directories
all_dirs <- list.dirs("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/all_samples")
all_dirs <- all_dirs[-1]
for(i in 1:length(sample_name)) { 
  file.copy(from = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/GSE161529_features.tsv.gz",
            to   = all_dirs[i])
}


# Coping the respectively barcodes and matrix to each sample folder
for(i in 1:length(sample_name)) { 
  x <- files_splited[grepl(paste0(sample_name[i]),files_splited)]
  file.copy(from = 
              # matrix
              c(paste0("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/",x[1]),
              # barcodes
                paste0("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/",x[2])),
              # to the sample folder
            to   = paste0("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/all_samples/",sample_name[i]))
  
}


## After set each sample files to a single folder
## Fix their names to be possible to run Read10x()

# Sign variables for loop
folders <- list.dirs("/mnt/storage1/mmaycon/downloads/Valid_dset_GSE161529/all_samples")
folders <- folders[-1] #removing the current directory 
samples <- str_split_fixed (as.character(folders), "[e][s][/]", 2)[,2] 
samples = samples[!samples %in% "all_sample"]

# Fix file names for all samples
for(i in 1:length(sample_name)) {
  setwd(paste0("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/all_samples/",sample_name[i]))
  from <- list.files(paste0("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/all_samples/",sample_name[i]))
  
  file.rename(from[grepl("matrix",from)], "matrix.mtx.gz")
  file.rename(from[grepl("barcodes",from)], "barcodes.tsv.gz")
  file.rename(from[grepl("features",from)], "features.tsv.gz")
  
}


## Change gene ID from ENSB to gene symbol
features_genes = read.table(gzfile("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/GSE161529_features.tsv.gz"),header = F)


library(readr)
geneID_converted <- read_csv("/mnt/storage1/mmaycon/downloads/Valid_dset_GSE161529/gProfiler_hsapiens_11-7-2022_1-13-48 PM.csv") # gene names from features
geneID_converted <- as.data.frame(geneID_converted)
geneID_converted <- geneID_converted[!duplicated(geneID_converted$name),]
rownames(geneID_converted) <- geneID_converted$initial_alias
# geneID_converted is the dataframe with ENSB and gene_symbol info
# .csv download upon on of the next "mtx" objects gene ensembl via gProfiler


## Make a list of seurat obj (1 obj/sample) 
Sobj <- list()
for(i in 1:length(sample_name)) {
  mtx <- Read10X(data.dir = paste0("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Valid_dset_GSE161529/all_samples/",sample_name[i]),gene.column=1) # read expression matrix
  mtx <- mtx[geneID_converted$initial_alias,]
  mtx_genes_ordered <- rownames(mtx)
  gene_symbol_to_mtx <- geneID_converted[mtx_genes_ordered,]$name
  rownames(mtx) <- gene_symbol_to_mtx
  colnames(mtx) <- gsub("_", "-", colnames(mtx)) 
  rownames(mtx) <- gsub("_", "-", rownames(mtx)) 
  
  Sobj[[i]] <- CreateSeuratObject(counts = mtx, min.cells = 3, min.features = 200) # create seurat objects in a list
  Sobj[[i]]@meta.data$PatientID =  as.character(sample_name[i])
} 


## Merge seurat objects with reduce() function 
## it is the same of merge 2 by 2 objects until all the objects get merged
Merge_Seurat_List <- function(
    list_seurat,
    add.cell.ids = NULL,
    merge.data = TRUE
    #project = "SeuratProject"
) {
  merged_object <- reduce(list_seurat, function(x, y) {
    merge(x = x, y = y, add.cell.ids = add.cell.ids, merge.data = merge.data)
  })
}

#lapply(Sobj, Reduce, f = merge)
#reduce(Sobj)

Sobj_Valid_Merged <- Merge_Seurat_List(Sobj)
saveRDS(Sobj_Valid_Merged, file = ".rds")
Sobj_Valid_Merged = readRDS(".rds")


# SUBSET ONLY INTERSTING SAMPLES
gset <- getGEO("GSE161529", GSEMatrix =TRUE, destdir = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset")
dim(pData(gset[[1]])) # 69 54
head(pData(gset [[1]])[, 1:3])
pData.GSE161529 <- pData(gset [[1]])
head(pData.GSE161529)
table(pData.GSE161529$`menopause:ch1`)

pData_sub <- pData.GSE161529 %>% 
  filter(`gender:ch1` %in% "Female", 
         `cell population:ch1` %in% c("Epithelial","Total"
         ))

pData_sub$PatientID = rownames(pData_sub)


# Keep only interstant columns 
pData_sub = pData_sub[, c(48:55)]
colnames(pData_sub) = gsub(":","_",colnames(pData_sub))
colnames(pData_sub) = gsub(" ","_",colnames(pData_sub))

save(pData_sub, file = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/GSE161529_pData_sub_temp_from_GEO_all_tumors_except_lynph.rda")
load("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/GSE161529_pData_sub_temp_from_GEO_all_tumors_except_lynph.rda")

library(stringr)
Sobj_Valid_Merged@meta.data$Barcode = colnames(Sobj_Valid_Merged)
meta_data = Sobj_Valid_Merged@meta.data
pData_sub = merge(pData_sub, meta_data) #merging by Patient ID
#identical(Sobj_Valid_Merged@meta.data$Barcode_cut,
          #pData_sub$Barcode_cut) #TRUE - just in case you ant to over write it on Seurat object meta.data like this: Sobj_Valid_Merged@meta.data = pData_sub

Sobj_Valid_Merged = Sobj_Valid_Merged[, colnames(Sobj_Valid_Merged) %in% pData_sub$Barcode]
dim(Sobj_Valid_Merged) #21843 373924

saveRDS(Sobj_Valid_Merged, file = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Valid_Merged_No_lynph_samples_no_clusterization.rds")

Sobj_Valid_Merged = readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Valid_Merged_No_lynph_samples_no_clusterization.rds")
Sobj = Sobj_Valid_Merged
dim(Sobj) #21843 373924

# QUALITY CONTROL
# Check mit % and vln plots
Sobj[["percent.mt"]] = PercentageFeatureSet(Sobj, pattern = "^MT-")
Idents(Sobj) = "PatientID"
No_qcs_vln = VlnPlot(Sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(No_qcs_vln, file = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Figures/Vlnplot_No_QCs.pdf", width = 14, height = 8, dpi = 100)


Sobj = subset(Sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA > 200)
Qcs_vln = VlnPlot(Sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
ggsave(Qcs_vln, file = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Figures/Vlnplot_QCs.pdf", width = 14, height = 8, dpi = 100)
dim(Sobj[["RNA"]]@data) #27719 81376
# Filter patients by pData_sub 


## SEURAT PIPELINE - CLUSTER, UMAP
Sobj = NormalizeData(Sobj, normalization.method = "LogNormalize", scale.factor = 10000)
#range(Sobj[["RNA"]]@data) # 0.000000 9.111652
Sobj = FindVariableFeatures(Sobj, selection.method = 'vst', nfeatures = 3000, verbose = FALSE)
Sobj = ScaleData(Sobj, verbose = FALSE) #scales only most varible features
Sobj = RunPCA(Sobj, verbose = FALSE, features = VariableFeatures(object = Sobj)) 
#ElbowPlot(Sobj) # at least 20 PCs
Sobj = FindNeighbors(Sobj, reduction = 'pca', dims = 1:30, verbose = FALSE)
Sobj = FindClusters(Sobj, resolution = 0.5, verbose = FALSE)
Sobj = RunUMAP(Sobj, reduction = 'pca', dims = 1:30, verbose = FALSE)

saveRDS(Sobj, file = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Valid_Merged_No_lynph_samples_yes_clusterization.rds")


### Seurat visualization
# UMAP seurat clusters
library(viridis)
#n_colors <- length(table(Sobj@meta.data$seurat_clusters))
#pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj) <- "seurat_clusters"
DimPlot(Sobj) #+ scale_color_manual(values = pal)

# UMAP patients
#n_colors <- length(table(Sobj@meta.data$PatientID))
#pal <- viridis(n = n_colors, option = "G", direction = -1)
Idents(Sobj) <- "PatientID"
DimPlot(Sobj) #+ scale_color_manual(values = pal)

# UMAP seurat clusters splited by patients
Idents(Sobj) <- "seurat_clusters"
DimPlot(Sobj, split.by = "PatientID")
table(Sobj@meta.data$PatientID, Sobj@meta.data$seurat_clusters)


# # UMAP major celltypes
# n_colors <- length(table(Sobj@meta.data$celltype_major))
# pal <- viridis(n = n_colors, option = "G", direction = -1)
# Idents(Sobj) <- "celltype_major"
# DimPlot(Sobj) + scale_color_manual(values = pal)
# 
# # UMAP stemness
# Idents(Sobj) <- "stemness"
# FeaturePlot(object = Sobj, features = c("stemness")) &
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# 
# # Ridge plots - stemness by major celltype
# Idents(Sobj) <- "celltype_major"
# RidgePlot(Sobj, features = "stemness", ncol = 1)
# 
# # Vln - stemness by major celltype
# Idents(Sobj) <- "celltype_major"
# VlnPlot(Sobj, features = "stemness")


#**SC-TYPE - MAP CELL TYPES FROM REFERENCE PAPER MARKERS**#
# PLOT STEMNESS AND FIND THE SAME POPULATIONS AS IN THE REF DATASET


# DONT PREDICT STEMNESS !! MERGE STEMNESS FROM A .rda WHICH THERE ARE STEMNESS PREDICTED FOR REF+VALID+IPS



