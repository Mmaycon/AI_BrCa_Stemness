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


### Reference expression matrix

# Creating Seurat object
GSE176078_data = Read10X(data.dir = "/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/Wu_etal_2021_BRCA_scRNASeq",gene.column=1)
dim(GSE176078_data) #29733 100064
range(GSE176078_data) #0 29831
  # Edit strings on rownames and colnames
colnames(GSE176078_data) = gsub("_", "-", colnames(GSE176078_data)) 
rownames(GSE176078_data) = gsub("_", "-", rownames(GSE176078_data)) 
  # Create Seurat Object
Sobj = CreateSeuratObject(counts = GSE176078_data, min.cells = 3, min.features = 200)
dim(Sobj[["RNA"]]@data) # 27719 100064


# QC filter
Sobj[["percent.mt"]] = PercentageFeatureSet(Sobj, pattern = "^MT-")
VlnPlot(Sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
Sobj = subset(Sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA > 200)
VlnPlot(Sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
dim(Sobj[["RNA"]]@data) #27719 81376

Sobj_refer = Sobj






### Validation expression matrix 
# Load seurat object
Sobj_valid = readRDS("/mnt/storage1/mmaycon/downloads/Valid_dset_GSE161529/Sobj_GSE161529_onlytumors_NOclusterization_all_tumors_except_lynph.rds")

# QC filter
Sobj_valid[["percent.mt"]] <- PercentageFeatureSet(Sobj_valid, pattern = "^MT-")
Idents(Sobj_valid) <- "patientID"
VlnPlot(Sobj_valid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
Sobj_valid_QC <- subset(Sobj_valid, 
                        subset = nFeature_RNA > 200 & 
                          nFeature_RNA < 5000 & 
                          nCount_RNA > 200 &
                          percent.mt < 10) # I didn't do it base on papers parameters. I used as the dataset came
VlnPlot(Sobj_valid_QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

dim(Sobj_valid_QC) #21843 152250

Sobj_valid = Sobj_valid_QC 



### iPS expression matrix

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


### Combine them all in one matrix

# Combine seurat objects
dim(Sobj_refer) #27719 81376
dim(Sobj_valid) #21843 152250
dim(Sobj_iPS) #20202 21599
Sobj_all <- merge(x = Sobj_refer, 
                  y = c(Sobj_valid, 
                        Sobj_iPS),
                  add.cell.ids = c("REF","VAL","IPS"))
                  #project = "")
dim(Sobj_all) #31404 255225




# TO BE CONTINUED - last edit 2023/01/18


# LogNormalize method  ---------------------------------------
Sobj_all <- NormalizeData(Sobj_all, normalization.method = "LogNormalize", scale.factor = 10000)
range(Sobj_all[["RNA"]]@data) # 

# Extract matrix expression to apply Stemness model (Malta et. al. 2018)
matrix_all <- as.matrix(GetAssayData(object = Sobj_all, slot = "data"))
dim(matrix_all) 



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
save(s, file = "/mnt/plummergrp/maycon/Figures_to_paper/stemness_index_prediction/reference_validation_and_iPS_predicted_together.rda")

# stemness index histogram
ggplot(s, aes(x=stemness)) + 
  geom_histogram(bins = 50) +
  xlab("Stemness index") + ylab("n Cells") + theme(axis.title.x = element_text(colour="black", size = 12), axis.title.y = element_text(colour="black", size = 15)) + ggtitle("") +theme(axis.text.x = element_text(angle = 45, size = 15 ,vjust=1,hjust=1))


### Merge each subset meta_data into its respectively seurat object
# Load stemness index predicted 
load("/mnt/plummergrp/maycon/Figures_to_paper/stemness_index_prediction/reference_validation_and_iPS_predicted_together.rda")

dim(Sobj_refer) #27719 81376
Sobj_refer = AddMetaData(
  object = Sobj_refer,
  metadata = s
)
dim(Sobj_refer) #27719 81376


dim(Sobj_valid) #21843 152250
Sobj_valid = AddMetaData(
  object = Sobj_valid,
  metadata = s
)
dim(Sobj_valid) #21843 152250



dim(Sobj_iPS) #20202 21599
Sobj_iPS = AddMetaData(
  object = Sobj_iPS,
  metadata = s
)
dim(Sobj_iPS) #20202 21599



### plot stemness on each UMAP 
# (I hope it would not mess with what we alredy got)
Sobj_valid_from_all_si@meta.data$
  val_md
  
# Standard Seurat pipeline
# NormalizeData() has been used before
  # Put it in a function
Sobj_list = list()
Sobj_list = c(Sobj_refer, Sobj_valid, Sobj_iPS)
Sobj_list_UMAPapped = list()
for (i in 1:length(Sobj_list)) {
  Sobj = NormalizeData(Sobj_list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  Sobj = FindVariableFeatures(Sobj, selection.method = 'vst', nfeatures = 3000, verbose = FALSE)
  Sobj = ScaleData(Sobj, verbose = FALSE) #scales only most varible features
  Sobj = RunPCA(Sobj, verbose = FALSE, features = VariableFeatures(object = Sobj)) 
  #ElbowPlot(Sobj) # 20 PCs
  Sobj = FindNeighbors(Sobj, reduction = 'pca', dims = 1:20, verbose = FALSE)
  Sobj = FindClusters(Sobj, resolution = 0.5, verbose = FALSE)
  Sobj = RunUMAP(Sobj, reduction = 'pca', dims = 1:20, verbose = FALSE)
  Sobj_list_UMAPapped [[i]] = Sobj
}

saveRDS(Sobj_list_UMAPapped, file = "/mnt/plummergrp/maycon/Figures_to_paper/stemness_index_prediction/Ref_Valid_iPS_each_sobj_with_stemness_predicted_alltogether.rds")




Sobj_refer_test = Sobj_list_UMAPapped[[1]]
Idents(Sobj_refer_test) <- "stemness"
FeaturePlot(object = Sobj_refer_test, features = c("stemness")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle("Sobj_reference_seuratobject")

Sobj_valid_test = Sobj_list_UMAPapped[[2]]
Idents(Sobj_valid_test) <- "stemness"
FeaturePlot(object = Sobj_valid_test, features = c("stemness")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle("Sobj_validation_seuratobject")


Sobj_iPS_test = Sobj_list_UMAPapped[[3]]
Idents(Sobj_iPS_test) <- "stemness"
FeaturePlot(object = Sobj_iPS_test, features = c("stemness")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle("Sobj_iPS_seuratobject")


### Plot multiple density plot with the 3 stemness distributions

# iPS
df = data.frame(Sobj_iPS_test@meta.data$stemness)
colnames(df) = "Stemness index"
df$dset = "iPS"

# Ref. dataset
df_2 = data.frame(Sobj_refer_test@meta.data$stemness)
colnames(df_2) = "Stemness index"
df_2$dset = "reference"

# Valid. dataset
df_3 = data.frame(Sobj_valid_test@meta.data$stemness)
colnames(df_3) = "Stemness index"
df_3$dset = "validation"

df_dens = rbind.fill(df,df_2,df_3)

ggplot(df_dens, aes(x =`Stemness index`, fill = dset)) + 
  geom_density(alpha = 0.5) + 
  theme_classic() +
  geom_vline(xintercept=0.7055146, lwd=2, colour="black", linetype="dotted")


dim(Sobj_refer_test@meta.data[Sobj_refer_test@meta.data$stemness >= 0.7055146,])
# 18  7
dim(Sobj_valid_test@meta.data[Sobj_valid_test@meta.data$stemness >= 0.7055146,])
#  7 8

quantile(Sobj_iPS_test@meta.data$stemness,na.rm = T,probs = c(0.25,0.50,0.75,1.0)) # 0.6805127 is the 25 percentile

dim(Sobj_refer_test@meta.data[Sobj_refer_test@meta.data$stemness >= 0.6805127,])
# 35  7
dim(Sobj_valid_test@meta.data[Sobj_valid_test@meta.data$stemness >= 0.6805127,])
#  12  8


ggplot(df_dens, aes(x =`Stemness index`, fill = dset)) + 
  geom_density(alpha = 0.5) + 
  theme_classic() +
  geom_vline(xintercept=0.56, lwd=2, colour="black", linetype="dotted")


dim(Sobj_refer_test@meta.data[Sobj_refer_test@meta.data$stemness >= 0.56,])
# 829   7
dim(Sobj_valid_test@meta.data[Sobj_valid_test@meta.data$stemness >= 0.56,])
# 883   8


# TO BE CONTINUED - last edit 2023/01/19

### Add celltype information ###
# Who are those cells ?

library(readr)
Whole_miniatlas_meta = read_csv("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/Whole_miniatlas_meta.csv")
head(Whole_miniatlas_meta)
Whole_miniatlas_meta = as.data.frame(Whole_miniatlas_meta)
Whole_miniatlas_meta$NAME = gsub("_", "-", Whole_miniatlas_meta$NAME) 
rownames(Whole_miniatlas_meta) = Whole_miniatlas_meta$NAME

# Merge pdata to Seurat object
Sobj_refer_test = AddMetaData(
  object = Sobj_refer_test,
  metadata = Whole_miniatlas_meta
)

table(Sobj_refer_test@meta.data[Sobj_refer_test@meta.data$stemness >= 0.56,]$celltype_minor) # thanks God

  # iPS lowest stemness cells 
quantile(Sobj_iPS_test@meta.data$stemness,na.rm = T,probs = c(0.08)) # stemness == 0.618299 

  # Reference dataset highst stemenss 99 percentile
quantile(Sobj_refer_test@meta.data$stemness,na.rm = T,probs = c(0.99)) # stemness == 0.5609754 

  # Validation dataset highst stemenss 99 percentile
quantile(Sobj_valid_test@meta.data$stemness,na.rm = T,probs = c(0.99)) # stemness == 0.5444483 

  # Reference dataset highst stemenss 95 percentile
quantile(Sobj_refer_test@meta.data$stemness,na.rm = T,probs = c(0.95)) # stemness == 0.4832981 

  # Validation dataset highst stemenss 95 percentile
quantile(Sobj_valid_test@meta.data$stemness,na.rm = T,probs = c(0.95)) # stemness == 0.49034 


  # Reference cell types high stemness: above 99 percentile
table(Sobj_refer_test@meta.data[Sobj_refer_test@meta.data$stemness >= 0.5609754,]$celltype_minor) 
  # Reference cell types high stemness: above 95 percentile
table(Sobj_refer_test@meta.data[Sobj_refer_test@meta.data$stemness >= 0.4832981,]$celltype_minor) 


# Do 2 round of DE with high stemness cells: i) 99percentile and ii) 95percentile 
# Compare the same 3 minor cell types - high stemness cells vs bottom stemness cells within the same cell type



saveRDS(Sobj_refer_test, file = "/mnt/plummergrp/maycon/Figures_to_paper/stemness_index_prediction/Ref_sobj_stemness_predicted_alltogether.rds")

saveRDS(Sobj_valid_test, file = "/mnt/plummergrp/maycon/Figures_to_paper/stemness_index_prediction/Valid_sobj_stemness_predicted_alltogether.rds")

saveRDS(Sobj_iPS_test, file = "/mnt/plummergrp/maycon/Figures_to_paper/stemness_index_prediction/iPS_sobj_stemness_predicted_alltogether.rds")



