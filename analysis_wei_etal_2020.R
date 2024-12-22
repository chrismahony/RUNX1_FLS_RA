```{r}
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
#library(EnsDb.Mmusculus.v75)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(cowplot)
library(patchwork) #latest version is required!
library(TFBSTools)
library(JASPAR2020)
library(gsfisher)
library(EnhancedVolcano)
#library(monocle)
setwd("/rds/projects/c/croftap-stia-atac/CM_multiome/STIA_andATAC/")
options(bitmapType='cairo')
library(CellChat)
#load()
```
```{r}
mouse <- readRDS("/rds/projects/c/croftap-stia-atac/CM_multiome/data_from_Ilya/mouse.rds")
counts<-mouse[["exprs_norm"]]
metadata<-mouse[["meta_data"]]
row.names(metadata) <- metadata$cell_id
metadata <- subset(metadata, select = -cell_id)
umap<-mouse[["umap_harmony"]]

mouse_seurat<-CreateSeuratObject(counts = counts, meta.data = metadata)
mouse_seurat <- FindVariableFeatures(mouse_seurat, selection.method = "vst", nfeatures = 2000)
mouse_seurat <- ScaleData(mouse_seurat, verbose = FALSE)
mouse_seurat <- RunPCA(mouse_seurat, npcs = 30, verbose = FALSE)
mouse_seurat[['UMAP']] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", global = T, assay = "RNA")
mouse_seurat <- RunUMAP(mouse_seurat, reduction = "pca", dims = 1:30)
DimPlot(mouse_seurat, reduction = "umap", group.by = "cell_type")
```
```{r}

library(harmony)
mouse_seurat<-RunHarmony(mouse_seurat, dims.use = 1:30,group.by.vars="donor")
DimPlot(mouse_seurat, reduction = "umap", group.by = "donor")


```
```{r}
Idents(mouse_seurat)<-"cell_type"
mouse_seurat$cell_type_label <- paste(Idents(mouse_seurat), mouse_seurat$label, sep = "_")
Idents(mouse_seurat)<-"cell_type_label"

sublining_N3vsISO<-FindMarkers(mouse_seurat, ident.1 = "sublining_N3", ident.2 = "sublining_ISO", logfc.threshold=-Inf)
sublining_KOvsWT<-FindMarkers(mouse_seurat, ident.1 = "sublining_KO", ident.2 = "sublining_WT", logfc.threshold=-Inf)

sublining_N3vsISO$gene<-rownames(sublining_N3vsISO)
sublining_KOvsWT$gene<-rownames(sublining_KOvsWT)

sublining_N3vsISO<-sublining_N3vsISO[order(rownames(sublining_N3vsISO)),] 
sublining_KOvsWT<-sublining_KOvsWT[order(rownames(sublining_KOvsWT)),] 


sublining_N3vsISO_common <- sublining_N3vsISO[sublining_N3vsISO$gene %in% sublining_KOvsWT$gene,]
sublining_KOvsWT_common <- sublining_KOvsWT[sublining_KOvsWT$gene %in% sublining_N3vsISO$gene,]

df_for_sctterplot<-subset(sublining_N3vsISO_common, select = -c(1, 3, 4, 5, 6))
df_for_sctterplotKO_vs_WT<-subset(sublining_KOvsWT_common, select = -c(1, 3, 4, 5, 6))
df_for_sctterplot$avg_log2FC_KOWT<-df_for_sctterplotKO_vs_WT$avg_log2FC

genes.to.label = c("Runx1")
p1 <- ggplot(data=df_for_sctterplot, aes(x=avg_log2FC_KOWT,y=avg_log2FC)) + geom_point() + ggtitle("X")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1

lining_N3vsISO<-FindMarkers(mouse_seurat, ident.1 = "lining_N3", ident.2 = "lining_ISO")
lining_KOvsWT<-FindMarkers(mouse_seurat, ident.1 = "lining_KO", ident.2 = "lining_WT")
intermediate_N3vsISO<-FindMarkers(mouse_seurat, ident.1 = "intermediate_N3", ident.2 = "intermediate_ISO")
intermediate_KOvsWT<-FindMarkers(mouse_seurat, ident.1 = "intermediate_KO", ident.2 = "intermediate_WT")

lining_N3vsISO$gene<-rownames(lining_N3vsISO)
lining_KOvsWT$gene<-rownames(lining_KOvsWT)

intermediate_N3vsISO$gene<-rownames(intermediate_N3vsISO)
intermediate_KOvsWT$gene<-rownames(intermediate_KOvsWT)
```


```{r}

Idents(mouse_seurat)<-"cell_type_label"

DotPlot(mouse_seurat, features = c("Il33"), idents = c("lining_WT", "lining_KO","lining_ISO", "lining_N3"))
DotPlot(mouse_seurat, features = c("Il33"), idents = c("sublining_WT", "sublining_KO","sublining_ISO", "sublining_N3"))

DotPlot(mouse_seurat, features = c("Runx1", "Nfatc2", "Nfic"), idents = c("sublining_ISO", "sublining_N3"))

DotPlot(mouse_seurat, features = c("Cbfb"), idents = c("sublining_ISO", "sublining_N3"))



median.stat <- function(x){
    out <- quantile(x, probs = c(0.5))
    names(out) <- c("ymed")
    return(out) 
}

levels(mouse_seurat)
VlnPlot(mouse_seurat, features =c("Runx1"), idents = c("sublining_WT", "sublining_KO"), pt.size = 0)+
    stat_summary(fun.y = median.stat, geom='point', size = 2, colour = "black") 

VlnPlot(mouse_seurat, features =c("Runx1"), idents = c("sublining_ISO", "sublining_N3"), pt.size = 0)+
    stat_summary(fun.y = median.stat, geom='point', size = 2, colour = "black") 



markers_cell_type_label<-FindAllMarkers(mouse_seurat, only.pos = T)
markers_cell_type_label_f<-markers_cell_type_label[markers_cell_type_label$p_val_adj  <0.05,]

library(splitstackshape)
markers_cell_type_label_f<-cSplit(markers_cell_type_label_f, splitCols = "cluster", sep="_")
markers_cell_type_label_f_SL<-markers_cell_type_label_f[markers_cell_type_label_f$cluster_1=="sublining",]
markers_cell_type_label_f_SL_WT<-markers_cell_type_label_f_SL[markers_cell_type_label_f_SL$cluster_2 == "WT",]



Idents(mouse_seurat)<-"cell_type_label"
levels(mouse_seurat)


dotplot<-DotPlot(mouse_seurat, features = unique(markers_cell_type_label_f_SL_WT$gene), idents = c("sublining_WT", "sublining_KO", "sublining_ISO", "sublining_N3"))

dotplot<-dotplot$data

dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

Heatmap(dotplot, cluster_columns = F)


dotplot["Runx1",]
dotplot["Igf1",]

```
```{r}
VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('sublining_KO', 'sublining_WT'), stack = T)

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('sublining_N3', 'sublining_ISO'), stack = T)+geom_jitter(size=0.01)

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('lining_KO', 'lining_WT'), stack = T)+geom_jitter(size=0.01)

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('lining_N3', 'lining_ISO'), stack = T)+geom_jitter(size=0.01)

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('intermediate_KO', 'intermediate_WT'), stack = T)+geom_jitter(size=0.01)

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('intermediate_N3', 'intermediate_ISO'), stack = T)+geom_jitter(size=0.01)
```


```{r}









```



```{r}

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('sublining_KO', 'sublining_WT'), stack = T)

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('sublining_N3', 'sublining_ISO'), stack = T)

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('lining_KO', 'lining_WT'), stack = T)

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('lining_N3', 'lining_ISO'), stack = T)

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('intermediate_KO', 'intermediate_WT'), stack = T)

VlnPlot(mouse_seurat, features = c("Runx1", "Cbfb"), idents = c('intermediate_N3', 'intermediate_ISO'), stack = T)
```

```{r}
mouse_seurat<-AddModuleScore(mouse_seurat, features = "Runx1", name = "Runx1_score")
mouse_seurat_metadata<-mouse_seurat@meta.data
mouse_seurat_metadata <- subset(mouse_seurat_metadata, select = c(score_notch, Runx1_score1))
ggplot(mouse_seurat_metadata, aes(x = Runx1_score1, y = score_notch)) +
        geom_point() +
        stat_smooth(method = "lm",
        col = "#C42126", se = FALSE, size = 1
)
```

```{r}
mouse_seurat_SL_ISO <- mouse_seurat[,grepl("sublining_ISO", mouse_seurat$cell_type_label, ignore.case=TRUE)]
mouse_seurat_SL_KO <- mouse_seurat[,grepl("sublining_KO", mouse_seurat$cell_type_label, ignore.case=TRUE)]
mouse_seurat_SL_N3 <- mouse_seurat[,grepl("sublining_N3", mouse_seurat$cell_type_label, ignore.case=TRUE)]
mouse_seurat_SL_WT <- mouse_seurat[,grepl("sublining_WT", mouse_seurat$cell_type_label, ignore.case=TRUE)]


mouse_seurat_SL_ISO<-mouse_seurat_SL_ISO@meta.data
mouse_seurat_SL_ISO <- subset(mouse_seurat_SL_ISO, select = c(score_notch, Runx1_score1))
ggplot(mouse_seurat_SL_ISO, aes(x = Runx1_score1, y = score_notch)) +
        geom_point() +
        stat_smooth(method = "lm",
        col = "#C42126", se = FALSE, size = 1
)

```
```{r}
mouse_seurat_SL_N3<-mouse_seurat_SL_N3@meta.data
mouse_seurat_SL_N3 <- subset(mouse_seurat_SL_N3, select = c(score_notch, Runx1_score1))
ggplot(mouse_seurat_SL_N3, aes(x = Runx1_score1, y = score_notch)) +
        geom_point() +
        stat_smooth(method = "lm",
        col = "#C42126", se = FALSE, size = 1
)
```

```{r}
mouse_seurat <- RunUMAP(mouse_seurat, dims = 1:30)

DimPlot(mouse_seurat, group.by = "cell_type")




```

