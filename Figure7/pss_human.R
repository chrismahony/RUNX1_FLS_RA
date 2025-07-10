



Idents(pss_obj) <- 'gross_celltype'
cols_umap <- ArchR::paletteDiscrete(pss_obj@meta.data[, "gross_celltype"])
DimPlot(pss_obj, cols=cols_umap)

#corr with RUNX1

fibs_pss <- subset(pss_obj, idents="Fibroblast")

Idents(fibs_pss) <- 'subcluster'
cols_umap <- ArchR::paletteDiscrete(fibs_pss@meta.data[, "subcluster"])


ncol(fibs_pss)
cols_df <- colours$cluster %>% as.data.frame()
DimPlot(fibs_pss, cols=c("#208A42", "#F47D2B", "#89288F","#D51F26", "#272E6A" ))

DotPlot(fibs_pss, features=c("RUNX1"))


