

DimPlot(pss_obj)+NoLegend()
DotPlot(pss_obj, features = c("RUNX1", "CTHRC1", "IGF1", "MMP14"), idents = levels(pss_obj)[grep("Fibroblast", levels(pss_obj))])
DotPlot(pss_obj, features = c("RUNX1", "CXCL11"), idents = levels(pss_obj)[grep("Fibroblast", levels(pss_obj))])

library(babelgene)

inhibitor_DEGs_f <- inhibitor_DEGs %>% filter(padj < 0.05 & log2FoldChange > 0.5)

mouse_inhibitor <- orthologs(genes=inhibitor_DEGs_f$gene, species = "mouse")

DEGs_human_RA <- subres_human_RA %>% filter(comparison == "EV_vs_R1C" & log2FoldChange < -2 &padj < 0.05)
mouse_RA <- orthologs(genes = DEGs_human_RA$gene, species = "mouse")


pss_obj <- AddModuleScore(pss_obj, features = all_markers_medfibros %>% dplyr::filter(cluster== "SPARC+COL3A1+ C4") %>% head(40) %>% rownames() %>% list(), name="endo_interacting_new")

pss_obj <- AddModuleScore(pss_obj, features = all_markers_medfibros %>% dplyr::filter(cluster== "CXCL10+CCL19+ C11") %>% head(40) %>% rownames() %>% list(), name="Tcell_interacting_new")


pss_obj <- AddModuleScore(pss_obj, features = list(mouse_inhibitor$human_symbol), name = "mouse_inhibitor")
pss_obj <- AddModuleScore(pss_obj, features = list(mouse_RA$human_symbol), name = "human_RA_fibs")



dotplot<-DotPlot(pss_obj, features=c( "Tcell_interacting_new1","endo_interacting_new1", "RUNX1"),idents=levels(pss_obj)[grep("Fibroblast", levels(pss_obj))])

dotplot<-dotplot$data
library(tidyr)
dotplot<-dotplot %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)
library(splitstackshape)
anno <- colnames(dotplot) %>% as.data.frame()
colnames(anno) <- c("cluster")

colours <- list('cluster' = ArchR::paletteDiscrete(anno$cluster))
col_ann <- HeatmapAnnotation(df = anno, col=colours)

Heatmap(dotplot, top_annotation = col_ann)





Idents(pss_obj) <- 'gross_celltype'
cols_umap <- ArchR::paletteDiscrete(pss_obj@meta.data[, "gross_celltype"])
DimPlot(pss_obj, cols=cols_umap)

fibs_pss <- subset(pss_obj, idents="Fibroblast")

Idents(fibs_pss) <- 'subcluster'
cols_umap <- ArchR::paletteDiscrete(fibs_pss@meta.data[, "subcluster"])


ncol(fibs_pss)
cols_df <- colours$cluster %>% as.data.frame()
DimPlot(fibs_pss, cols=c("#208A42", "#F47D2B", "#89288F","#D51F26", "#272E6A" ))



