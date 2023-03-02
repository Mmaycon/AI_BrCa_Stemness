### Pathway analysis for high stemness DE 
  # Cancer Eptithelial
    # Cancer Cycling
    # Cancer Basal
  # T-cells
    # Cyling T-cells

# Total of 5 DEs 
# Thus, 5 genes lists 
# Obs: run Gene. Ontology with up and down genes separated 

library(clusterProfiler)
library(org.Hs.eg.db)


# Cancer Epithelil high stemness DEG
  # Load DE output
load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Epithelial_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")

  # up_regulated_genes
dim(findmarkers_sig)
up_genes = findmarkers_sig[findmarkers_sig$up_down_exp %in% "up",]$genes

ego2 <- enrichGO(gene = up_genes,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = 'SYMBOL',
                                # ont	One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                                ont           = "ALL",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05)

ego2@result[order(ego2@result$Count, decreasing = TRUE),][1:10,]
barplot(ego2, drop=TRUE, showCategory=12, main = "Cancer_epi up genes high stemness DE")
clusterProfiler::dotplot(ego2) + ggtitle("Cancer_epi up genes high stemness DE")

  # down_regulated_genes
dim(findmarkers_sig)
down_genes = findmarkers_sig[findmarkers_sig$up_down_exp %in% "down",]$genes

ego2 <- enrichGO(gene = down_genes,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 # ont	One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego2@result[order(ego2@result$Count, decreasing = TRUE),][1:10,]
barplot(ego2, drop=TRUE, showCategory=12, main = "Cancer_epi down genes high stemness DE")
clusterProfiler::dotplot(ego2) + ggtitle("Cancer_epi down genes high stemness DE")


  # Cancerr Cycling high stemness DEG
    # Load DE output
Pathways_up = function(findmarkers_sig) {
  
  # up_regulated_genes
  dim(findmarkers_sig)
  up_genes = findmarkers_sig[findmarkers_sig$up_down_exp %in% "up",]$genes
  
  ego2 <- enrichGO(gene = up_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   # ont	One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
  
  ego2@result[order(ego2@result$Count, decreasing = TRUE),][1:10,]
  barplot(ego2, drop=TRUE, showCategory=12, main = "up genes high stemness DE")
  clusterProfiler::dotplot(ego2) + ggtitle("up genes high stemness DE")
  
}


Pathways_down = function(findmarkers_sig) {
# down_regulated_genes
dim(findmarkers_sig)
down_genes = findmarkers_sig[findmarkers_sig$up_down_exp %in% "down",]$genes

ego2 <- enrichGO(gene = down_genes,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 # ont	One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego2@result[order(ego2@result$Count, decreasing = TRUE),][1:10,]
barplot(ego2, drop=TRUE, showCategory=12, main = "down genes high stemness DE")
clusterProfiler::dotplot(ego2) + ggtitle("down genes high stemness DE")

}

load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Cycling_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")
Pathways_up(findmarkers_sig)
Pathways_down(findmarkers_sig)


  # Cancerr Basal high stemness DEG
    # Load DE output
load("/mnt/storage1/mmaycon/Back_to_beginning_BCscProject/high_stemness_DE/Cancer_Basal_SC_high_si_Donwsampled_BOTTOM25%_si_DE_genes.rda")





