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
setwd("/mnt/storage1/mmaycon/downloads")
# iPS scRNA-seq download from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6524

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

Sobj <- merge(x = Sobj_1, 
      y = Sobj_2, add.cell.ids = c("iPS_1","iPS_2"), 
            project = "sc_iPS_dataset")

Sobj_iPS <- Sobj

Sobj[["percent.mt"]] <- PercentageFeatureSet(Sobj, pattern = "^MT-")
VlnPlot(Sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

Sobj_QC <- Sobj

# LogNormalize method  ---------------------------------------
options(future.globals.maxSize= 891289600)
Sobj_QC <- NormalizeData(Sobj_QC, normalization.method = "LogNormalize", scale.factor = 10000)
range(Sobj_QC[["RNA"]]@data) # 

# Extract matrix expression to apply Stemness model (Malta et. al. 2018)
matrix_exp <- as.matrix(GetAssayData(object = Sobj_QC, slot = "data"))
dim(matrix_exp) # 
matrix_exp[1:3,1:3]

# Stemness prediction
load("/mnt/storage1/mmaycon/downloads/model_RNA_MALTA.2018.Rda") # mm = modelo
w <- mm$w
w[1:5]
length(w) #12953 (number of genes in the model)

# Filtering commun genes
length(intersect(rownames(matrix_exp), names(w)))  
predict.DATA <- matrix_exp[rownames(matrix_exp)%in%  names(w) ,]
length(rownames(predict.DATA))  
w <- w[ rownames(predict.DATA) ]
length(intersect(names(w),rownames(predict.DATA)))  
w[1:5]
length(names(w)) 
is.vector(w) #TRUE

## Score the Matrix `X` using Spearman correlation
s <- apply( predict.DATA, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
s[1:5]

## Scale the scores to be between 0 and 1
s <- s - min(s)
s <- s / max(s)
s[1:5]

## Merge stemenss index to Seurat metadata 
names(s)
head(Sobj_QC@meta.data, 5) 

length(intersect(colnames(x = Sobj_QC),names(s))) # 36424
identical(colnames(x = Sobj_QC),names(s)) # TRUE
s <- as.data.frame(s)
#save(s, file = "/mnt/storage1/mmaycon/downloads/Valid_dset_GSE161529/stemness_index_breast_GSE16152_scRNAseq.rda")


# stemness index histogram
ggplot(s, aes(x=s)) + 
  geom_histogram(bins = 50) +
  xlab("Stemness index") + ylab("n Cells") + theme(axis.title.x = element_text(colour="black", size = 12), axis.title.y = element_text(colour="black", size = 15)) + ggtitle("") +theme(axis.text.x = element_text(angle = 45, size = 15 ,vjust=1,hjust=1))


# add stemnees index into seurat meta.data
Sobj_QC[["stemness"]] <- s$s


Sobj_QC <- FindVariableFeatures(Sobj_QC, selection.method = "vst", nfeatures = 2000) 
all.genes <- rownames(Sobj_QC)
Sobj_QC <- ScaleData(Sobj_QC, features = all.genes)
Sobj_QC <- RunPCA(Sobj_QC, features = VariableFeatures(object = Sobj_QC)) 
Sobj_QC <- FindNeighbors(Sobj_QC, dims = 1:10) # defaut
Sobj_QC <- FindClusters(Sobj_QC, resolution = 0.5)
Sobj_QC <- RunUMAP(Sobj_QC, dims = 1:10)
DimPlot(Sobj_QC, reduction = "umap")

FeaturePlot(object = Sobj_QC, features = c("stemness")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  xlab("UMAP_1") + 
  ylab("UMAP_2") + 
  labs(colour = "Stemness Index") +
  ggtitle("") + 
  theme(axis.title.x = element_text(colour="black", size = 18), 
        axis.title.y = element_text(colour="black", size = 18),
        legend.text = element_text(colour="black", size = 16),
        legend.title = element_text(colour="black", size = 18)) + 
  #geom_jitter(alpha=0.4) +
  theme(axis.text.x = element_text(angle = 45, size = 10 ,vjust=1,hjust=1)) 
  ggsave(file="/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/Figures_to_paper/iPS/only_iPS_singlecell_stemness_UMAP.pdf",width = 6, height = 4, dpi = 300) 


save(Sobj_QC, file = "/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/Figures_to_paper/iPS/only_iPS_singlecell_stemness_seuratobject.rds")
meta_data = Sobj_QC@meta.data
save(meta_data, file = "/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/Figures_to_paper/iPS/only_iPS_singlecell_stemness_metadata.rda")

VlnPlot(Sobj_QC, features = "stemness")

