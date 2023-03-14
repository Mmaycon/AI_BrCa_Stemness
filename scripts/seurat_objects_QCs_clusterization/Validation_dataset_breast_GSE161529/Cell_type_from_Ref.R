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
library(readxl)
### Sc-type template - just to take a look

sctype_gene_markers <- read_excel("/mnt/plummergrp/maycon/Miami_Ovarian_scrnaseq/sctype_gene_markers.xlsx")
sctype_gene_markers = sctype_gene_markers %>%
  as.data.frame()

head(sctype_gene_markers)
names(sctype_gene_markers)
# This should be the structure of sc-type input on gene_sets_prepare function
sctype_gene_markers[1:5,] 



### Retriving cell type markers from Reference paper (PMID34493872) based on sctype_gene_markers template
# Loading markers
celltype_markers = read_excel("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/major_and_minor_celltypes_from_ref_paper_PMID34493872.xlsx")
celltype_markers = celltype_markers %>%
  as.data.frame()

# Fixing colnames
names(celltype_markers)[names(celltype_markers) == "classification tier"] <- "classification_tier"
names(celltype_markers)[names(celltype_markers) == "avg_logFC"] <- "avg_log2FC"

# There are no negative markers 
celltype_markers[celltype_markers$avg_log2FC < 0, ] #none negative
celltype_markers[celltype_markers$avg_log2FC > 0, ] #only positive


###########################
###** Major Cell Type **###
###########################
major_celltypes = celltype_markers %>%
  dplyr::filter(classification_tier %in% "major")

## Filter top 10 markers WITHIN each celltype
# Cell types
celltype = as.vector(names(table(major_celltypes$cluster)))
# Loop to subset top and bottom markers
df_list = list()
#i = 1
for(i in 1:length(celltype)) {
  # Odering df by Fold Change
  ordered_df = major_celltypes[order(major_celltypes$avg_log2FC, decreasing = TRUE), ]
  # Subset cell type
  celltype_select = ordered_df[ordered_df$cluster %in% celltype[i], ]
  # Taking top cell markers
  top_markers = celltype_select[c(1:30), ]
  # Store it in a list
  df_list[[i]] = top_markers
  
}

# Bind all df from the list
markers_sub = do.call(bind_rows, df_list)
length(intersect(markers_sub[markers_sub$cluster %in% "Endothelial", ]$gene, major_celltypes[major_celltypes$cluster %in% "Endothelial", ]$gene)) 

names(table(markers_sub$cluster)) == celltype # All TRUE, okay 
table(markers_sub$cluster) 
high_avg_pos_celltype = markers_sub[markers_sub$avg_log2FC > 0,]


# Load seurat object
Seurat_object <- readRDS("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Valid_Merged_No_lynph_samples_yes_clusterization.rds")
table(Seurat_object@meta.data$seurat_clusters) 
dim(Seurat_object) #21843 276371 (Sobj after QCs)

# Differential expression on high stemness cells - it takes some time
# Idents(Seurat_object) <- "seurat_clusters"
# all_markers <- FindAllMarkers(object = Seurat_object, min.pct = 0.25)
# save(all_markers, file = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/all_markers.rda")
# Filtering Differential Expression output
load("/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/all_markers.rda")
head(all_markers)
dim(all_markers) #
findmarkers_sig <- all_markers[all_markers$p_val_adj < 0.05,] 
findmarkers_sig <- findmarkers_sig[findmarkers_sig$avg_log2FC < -0.5 | findmarkers_sig$avg_log2FC > 0.5,]
dim(findmarkers_sig) #
hist(findmarkers_sig$avg_log2FC)


findmarkers_sig <- findmarkers_sig[order(findmarkers_sig$avg_log2FC, decreasing = TRUE),]
findmarkers_sig$genes <- rownames(findmarkers_sig)
findmarkers_sig$up_down_exp = NA
findmarkers_sig[findmarkers_sig$avg_log2FC > 0,]$up_down_exp = "up"
findmarkers_sig[findmarkers_sig$avg_log2FC < 0,]$up_down_exp = "down"
table(findmarkers_sig$up_down_exp) 

# Intersect genes from cell type markers (paper) and cluster markers (that Sobj Validation we are working on)
length(intersect(findmarkers_sig[findmarkers_sig$avg_log2FC > 0.5, ]$gene, high_avg_pos_celltype$gene)) # n genes intersected
high_intersect = intersect(findmarkers_sig[findmarkers_sig$avg_log2FC > 0.5, ]$gene, high_avg_pos_celltype$gene)
# Percentage of high markers out of length(high_intersect) on clusters
hist(table(findmarkers_sig[findmarkers_sig$gene %in% high_intersect, ]$cluster) / length(high_intersect))


## Try to call celltypes with both those postive and negative markers OR only with positive markers

## Variables dictionary
# high_avg_pos_celltype contains cell type markers from Paper reference (and genes etc)
# high_intersect contains genes exp > 0.5 in Validation data clusters and among top genes we had choosen across each cell type from Paper Reference
# celltype is the number of cell types from Paper Reference
# Loop to create sc-type markers input based on our chosen markers
df_list = list()
for(i in 1:length(celltype)) {
  # Subset cell type specific positive markers
  gen_marker_pos = high_avg_pos_celltype[high_avg_pos_celltype$cluster %in% celltype[i], ]$gene
  # ###***Keep positive markers that match ***###
  gen_marker_pos = gen_marker_pos[gen_marker_pos %in% high_intersect]
  # remove genes with "orf" on their names - otherwise the function gene_sets_prepare() would not work
  gen_marker_pos = gen_marker_pos[!grepl(c("orf"),gen_marker_pos)]
  # Create a vector with all 10 markers separeted by comma with no space in between 
  gen_marker_pos_c = do.call(paste, c(as.list(gen_marker_pos), sep = ","))
  
  df = data.frame(tissueType = "Breast", cellName = celltype[i], 
                  geneSymbolmore1 = gen_marker_pos_c, geneSymbolmore2 = "NA",
                  shortName = "NA")
  df_list[[i]] = df
}

# Bind ao lists in one dataframe
df = do.call(bind_rows, df_list)


# Save it. It is one of inputs of sc-type
library(writexl)
write_xlsx(df,"/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/high_intersect_top30PaperMarkers_MoreThan0.5FCqueryData_without_orf_gene_all_MAJOR_celltypes.xlsx")



### Sc-type cell type mapping
# Packages
install.packages("HGNChelper")
install.packages("openxlsx")
library(openxlsx)
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Load gene markers we previous prepared 
db_ = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/high_intersect_top30PaperMarkers_MoreThan0.5FCqueryData_without_orf_gene_all_MAJOR_celltypes.xlsx"
# String with the tissue type we chose 
tissue = "Breast"



# Prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

##**The following steps are exactly how it is on sc-type githbub**##
#sctype_function <- function(Seurat_object) {
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = Seurat_object[["RNA"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = NULL) # gs2 = NULL is really important whenever there is no negative marker in the table
  
  
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(Seurat_object@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(Seurat_object@meta.data[Seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Seurat_object@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  
  # set low-confident (low ScType score) clusters to "unknown"
  ## gs2 = NULL when there is no negative marker
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  
  Seurat_object@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    Seurat_object@meta.data$customclassif[Seurat_object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
DimPlot(Seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') + ggtitle("")
  
#}

# Changing colnames
names(Seurat_object@meta.data)[names(Seurat_object@meta.data) == "customclassif"] <- "Major_celltypes_sctyped"






###########################
###** Minor Cell Type **###
###########################

minor_celltypes = celltype_markers %>%
  dplyr::filter(classification_tier %in% "minor", 
                # cluster %in% c("Cancer Basal SC",
                #                "Cancer Cycling",
                #                "Cancer LumA SC",
                #                "Cancer LumB SC",
                #                "Cancer Her2 SC",
                #                "Cycling T-cells",
                #                "T cells CD4+",
                #                "T cells CD8+",
                #                "NK cells",
                #                "NKT cells"))
                # cluster %in% c("Cancer Basal SC",
                #                "Cancer Cycling",
                #                "Cycling T-cells"))
                  cluster %in% c("Cancer Her2 SC",
                                 "Cancer Cycling"))
                               

## Filter top markers WITHIN each celltype
# Cell types
celltype = as.vector(names(table(minor_celltypes$cluster)))
# Loop to subset top and bottom markers
df_list = list()
#i = 1
for(i in 1:length(celltype)) {
  # Odering df by Fold Change
  ordered_df = minor_celltypes[order(minor_celltypes$avg_log2FC, decreasing = TRUE), ]
  # Subset cell type
  celltype_select = ordered_df[ordered_df$cluster %in% celltype[i], ]
  # Taking first 10 markers and bottom 10 markers
  top_markers = celltype_select[c(1:10), ]
  # Store it in a list
  df_list[[i]] = top_markers
  
}

# Bind all df from the list
markers_sub = do.call(bind_rows, df_list)
length(intersect(markers_sub[markers_sub$cluster %in% "Cancer Cycling", ]$gene, minor_celltypes[minor_celltypes$cluster %in% "Cancer Cycling", ]$gene))

names(table(markers_sub$cluster)) == celltype # All TRUE, okay 
table(markers_sub$cluster) # n markers from each celltype
high_avg_pos_celltype = markers_sub[markers_sub$avg_log2FC > 0,]



# Intersect genes from cell type markers (paper) and cluster markers (that Sobj Validation we are working on)
length(intersect(findmarkers_sig[findmarkers_sig$avg_log2FC > 0.5, ]$gene, high_avg_pos_celltype$gene)) # 
high_intersect = intersect(findmarkers_sig[findmarkers_sig$avg_log2FC > 0.5, ]$gene, high_avg_pos_celltype$gene)
# Percentage of high markers out of length(high_intersect) on clusters
hist(table(findmarkers_sig[findmarkers_sig$gene %in% high_intersect, ]$cluster) / length(high_intersect))




## Try to call celltypes with both those postive and negative markers OR only with positive markers

## Variables dictionary
# high_avg_pos_celltype contains cell type markers from Paper reference (and genes etc)
# high_intersect contains genes exp > 0.5 in Validation data clusters and among top genes we had choosen across each cell type from Paper Reference
# celltype is the number of cell types from Paper Reference
# Loop to create sc-type markers input based on our chosen markers

df_list = list()
for(i in 1:length(celltype)) {
  # Subset cell type specific positive markers
  gen_marker_pos = high_avg_pos_celltype[high_avg_pos_celltype$cluster %in% celltype[i], ]$gene
  # ###***Keep positive markers that match ***###
  gen_marker_pos = gen_marker_pos[gen_marker_pos %in% high_intersect]
  # remove genes with "orf" on their names
  gen_marker_pos = gen_marker_pos[!grepl(c("orf"),gen_marker_pos)]
  # Create a vector with  markers separeted by comma with no space in between 
  gen_marker_pos_c = do.call(paste, c(as.list(gen_marker_pos), sep = ","))
  
  df = data.frame(tissueType = "Breast", cellName = celltype[i], 
                  geneSymbolmore1 = gen_marker_pos_c, geneSymbolmore2 = "NA",
                  shortName = "NA")
  df_list[[i]] = df
}

df = do.call(bind_rows, df_list)



library(writexl)
write_xlsx(df,"/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/high_intersect_top10PaperMarkers_MoreThan0.5FCqueryDataUsed_without_orf_gene_2_MINORcelltypes.xlsx")



### Sc-type cell type mapping
# Packages
install.packages("HGNChelper")
install.packages("openxlsx")
library(openxlsx)
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Load gene markers we previous prepared 
db_ = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/high_intersect_top10PaperMarkers_MoreThan0.5FCqueryDataUsed_without_orf_gene_2_MINORcelltypes.xlsx"
# String with the tissue type we chose 
tissue = "Breast"



# Prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# The following steps are exactly how it is on sc-type githbub 
#sctype_function <- function(Seurat_object) {

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = Seurat_object[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = NULL) # gs2 = NULL is really important whenever there is no negative marker in the table


# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(Seurat_object@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(Seurat_object@meta.data[Seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Seurat_object@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
## gs2 = NULL when there is no negative marker
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

Seurat_object@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  Seurat_object@meta.data$customclassif[Seurat_object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(Seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') + ggtitle("")

#}

# Changing colnames
names(Seurat_object@meta.data)[names(Seurat_object@meta.data) == "customclassif"] <- "Minor_celltypes_sctyped_only_2_celltypes_10markers"

DimPlot(Seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Minor_celltypes_sctyped_only_10_celltypes_80markers') + ggtitle("")

DimPlot(Seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Minor_celltypes_sctyped_only_2_celltypes_80markers') + ggtitle("")

DimPlot(Seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Minor_celltypes_sctyped_only_10_celltypes_10_markers') + ggtitle("")

DimPlot(Seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Minor_celltypes_sctyped_only_3_celltype_ofInterest_10_markers') + ggtitle("")

DimPlot(Seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Minor_celltypes_sctyped_only_2_celltypes_10markers') + ggtitle("")

DimPlot(Seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Major_celltypes_sctyped') + ggtitle("")


DimPlot(Seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Celltype_to_use') + ggtitle("")

DimPlot(Seurat_object, reduction = "umap", label = FALSE, repel = TRUE, group.by = 'PatientID') + ggtitle("")
table(Seurat_object$PatientID)

# Rank plot of patient cells
rankplot = data.frame(table(Seurat_object$PatientID))
names(rankplot) = c("patientID", "cells")
rankplot <- rankplot[order(rankplot$cells),] #ordena

library(viridis)
n_colors <- length(table(rankplot$patientID))
pal <- viridis(n = n_colors, option = "G", direction = -1)
library(ggplot2)
ggplot(rankplot, aes(reorder(patientID, +cells, sum), cells)) +
  geom_bar(stat="identity", aes(fill=factor(patientID))) +
  scale_fill_manual(values=c(pal), name="Legend") +
  ggtitle("Cells by patient") +
  ylab("Number of cells") + 
  xlab("Patient ID") + 
  theme(axis.text.x = element_blank()) +
  theme_classic()


saveRDS(Seurat_object, file = "/mnt/plummergrp/maycon/Figures_to_paper/Validation_dset/Sobj_Valid_Merged_No_lynph_samples_yes_clusterization_CELLTYPE.rds")

  