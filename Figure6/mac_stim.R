
```{r}
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(gsfisher)
library(EnhancedVolcano)
options(bitmapType='cairo')
```




```{r}
DimPlot(aggr, group.by = "integrated_snn_res.0.01")
```

```{r}
DefaultAssay(aggr)<-'RNA'
FeaturePlot(aggr, features = "COL1A1")
FeaturePlot(aggr, features = "PDGFRA")
FeaturePlot(aggr, features = "CD248")

FeaturePlot(aggr, features = "CD14")
FeaturePlot(aggr, features = "CD68")

```
```{r}
DimPlot(aggr, group.by = "orig.ident")

```
```{r}
Idents(aggr)<-"orig.ident"
VlnPlot(aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


```
```{r}
Idents(aggr)<-'orig.ident'
DimPlot(aggr, split.by = "orig.ident")

```

```{r}


levels(aggr)
current.sample.ids <- c("AH1", "AH2", "AH3", "AH4", "AH5")
new.sample.ids <- c("Fb", "DexM", "LPSM", "Fb_DexM", "Fb_LPSM")

aggr@meta.data[["orig.ident"]] <- plyr::mapvalues(x = aggr@meta.data[["orig.ident"]], from = current.sample.ids, to = new.sample.ids)

DimPlot(aggr,split.by = "orig.ident")

```
```{r}


library(clustree)
clustree(aggr, prefix = "integrated_snn_res.")

DimPlot(aggr, group.by = "integrated_snn_res.0.3")
aggr<-FindClusters(aggr, resolution = c(0.4, 0.5), graph.name = "integrated_snn")
DimPlot(aggr, group.by = "integrated_snn_res.0.4")
DimPlot(aggr, group.by = "integrated_snn_res.0.5", label=T)
DimPlot(aggr, group.by = "orig.ident")

```

```{r}

Idents(aggr)<-'integrated_snn_res.0.5'
almarkers_0.5<-FindAllMarkers(aggr, only.pos = T, logfc.threshold = 1)


aggr$cm_clusters<-aggr$integrated_snn_res.0.5
levels(aggr)

current.sample.ids <- c("0",  "1",  "2" , "3"  ,"4" , "5" , "6",  "7" , "8",  "9" , "10", "11", "12", "13")
new.sample.ids <- c("fb_PRG4_STMN2",  "M_F13A1",  "M_SLC16A10" , "M_TPRG1_FABP5"  ,"M_APOC1" , "fb_CXCL12" , "fb_CEMIP",  "M_WTAP" , "fb_SLIT3",  "fb_COMP" , "M_TNFAIP6", "fb_MKI67", "NotSure", "M_CTSK")

aggr@meta.data[["cm_clusters"]] <- plyr::mapvalues(x = aggr@meta.data[["cm_clusters"]], from = current.sample.ids, to = new.sample.ids)
DimPlot(aggr, group.by = "cm_clusters", label = T)


Idents(aggr)<-'cm_clusters'

DotPlot(aggr, features = "RUNX1")

```

```{r}

DefaultAssay(aggr)<-'RNA'
#M1 markers
FeaturePlot(aggr, features = c("IL1A", "IL1B", "IL6", "NOS2", "TLR2", "TLR4", "CD80", "CD86"))
FeaturePlot(aggr, features = c("CD38", "GPR18", "FPR2"))

#M2 markers
FeaturePlot(aggr, features = c("CD115", "CD206", "PPARG", "ARG1", "CD163", "CD301", "DECTIN1", "PDL2", "FIZZ1"))
FeaturePlot(aggr, features = c("CD206",  "IL10", "YM1", "ARG1", "AMAC1", "MGL1"))

FeaturePlot(aggr, features = c("MERTK"))


```
```{r}
FeaturePlot(aggr, features = c("PI16", "CD34", "RUNX1", "RUNX2"))

```
```{r}
Idents(aggr)<-'orig.ident'

All_markers_samples<-FindAllMarkers(aggr, logfc.threshold = 0.5, only.pos = T)

DefaultAssay(aggr)<-'RNA'
aggr<-ScaleData(aggr)

All_markers_samples_top10<-All_markers_samples %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
DotPlot(aggr, features = unique(All_markers_samples_top10$gene)) + NoLegend()+RotatedAxis()



```
```{r}

```


```{r}
aggr$global<-aggr$integrated_snn_res.0.01
levels(aggr$global)

current.sample.ids <- c("0",  "1")
new.sample.ids <- c("Mac",  "Fb")

aggr@meta.data[["global"]] <- plyr::mapvalues(x = aggr@meta.data[["global"]], from = current.sample.ids, to = new.sample.ids)

aggr$orig.ident_condition<-paste(aggr$global, aggr$orig.ident, sep="_")
DimPlot(aggr, split.by = "orig.ident_condition")

```


```{r}
Idents(aggr)<-'orig.ident_condition'
levels(aggr)
Fb_LPSMvsFb<-FindMarkers(aggr, ident.1 = "Fb_Fb_LPSM", ident.2 = "Fb_Fb")
Fb_DexvsFb<-FindMarkers(aggr, ident.1 = "Fb_Fb_DexM", ident.2 = "Fb_Fb")
Fb_LPSMvsFb$gene<-rownames(Fb_LPSMvsFb)
Fb_DexvsFb$gene<-rownames(Fb_DexvsFb)


EnhancedVolcano(Fb_LPSMvsFb,
    lab = rownames(Fb_LPSMvsFb),
    selectLab = c("CD74", "APOE", "HLA-DRA", "CD68", "EBF1", "COL6A3", "COL1A2"),
    x = 'avg_log2FC',
    y = 'p_val_adj',
           title = 'Fb_LPSMvsFb',
    subtitle = "GEX, red=p_adj<0.05 & FC > 0.5",
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
    max.overlaps = 600,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5
    )


EnhancedVolcano(Fb_DexvsFb,
    lab = rownames(Fb_DexvsFb),
    selectLab = c("FGF7", "EBF1", "COL6A3", "APOE", "HLA-DRA", "CD74"),
    x = 'avg_log2FC',
    y = 'p_val_adj',
           title = 'Fb_LPSMvsFb',
    subtitle = "GEX, red=p_adj<0.05 & FC > 0.5",
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
    max.overlaps = 600,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5
    )
```


```{r}

#unique to fb_dex_vvs_fb compared to lPS_fb_vs_fb
Fb_DexMvsFb_unique<-Fb_DexvsFb[!Fb_DexvsFb$gene %in% Fb_LPSMvsFb$gene,]


#unique to fb_LPS_vvs_fb compared to Dex_fb_vs_fb
Fb_LPSMvsFb_unique<-Fb_LPSMvsFb[!Fb_LPSMvsFb$gene %in% Fb_DexvsFb$gene,]


```

```{r}
Idents(aggr)<-'orig.ident'
fibs<-subset(x = aggr, idents = "Fb")

fibs <- fibs %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims=1:30)
DimPlot(fibs)



fibs <-FindNeighbors(fibs, dims = 1:30)
fibs <-FindClusters(fibs, resolution = c(0.01, 0.05, 0.1, 0.2))

DimPlot(fibs, group.by = "RNA_snn_res.0.1", label=F)

Idents(fibs)<-'RNA_snn_res.0.1'

DefaultAssay(fibs)<-'RNA'
fibs_markers<-FindAllMarkers(fibs, only.pos = T)

VlnPlot(fibs, features = c("RUNX1", "RUNX2", "RUNX3", "CBFB"), pt.size = 0)

Idents(aggr)<-'integrated_snn_res.0.05'
aggr<-CellSelector(p1, object = aggr, ident = "Chondrocytes")
aggr$new_clusters<-aggr@active.ident

FeaturePlot(fibs, features = c("PRG4","FGF10","IGFBP3","COMP", "MKI67", "CD34"))


fibs$named_CM <- fibs@meta.data[["RNA_snn_res.0.1"]]
Idents(fibs) <- 'named_CM'
levels(fibs)
current.sample.ids <- c("0","1","2", "3", "4", "5")
new.sample.ids <- c("LL", "FGF10","IGFBP3","COMP", "MKI67", "CD34")

fibs@meta.data[["named_CM"]] <- plyr::mapvalues(x = fibs@meta.data[["named_CM"]], from = current.sample.ids, to = new.sample.ids)

DimPlot(fibs, group.by = "named_CM", label=F)


VlnPlot(fibs, features = c("RUNX1", "RUNX2", "CBFB"), pt.size = 0, stack=T)
VlnPlot(fibs, features = c("RUNX3"), pt.size = 0, stack=F)+ coord_flip()

VlnPlot(fibs, features = c("PRG4","FGF10","IGFBP3","COMP", "MKI67", "CD34", "RUNX1", "RUNX2", "CBFB"), pt.size = 0, stack=T)

FeaturePlot(fibs, features = c("RUNX1", "RUNX2", "RUNX3", "CBFB"), pt.size = 0)
ncol(fibs)
```

```{r}
Idents(aggr)<-'orig.ident_condition'
levels(aggr)
fibs_all<-subset(x = aggr, idents = c("Fb_Fb", "Fb_Fb_DexM", "Fb_Fb_LPSM"))

fibs_all <- fibs_all %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims=1:30)
DimPlot(fibs_all)

Idents(fibs_all)<-'orig.ident_condition'

DefaultAssay(fibs_all)<-'RNA'
#fibs_markers<-FindAllMarkers(fibs, only.pos = T)

VlnPlot(fibs_all, features = c("RUNX1", "RUNX2", "RUNX3", "CBFB"), pt.size = 0, stack)
```
```{r}

FindMarkers(fibs_all, ident.1 = "")


```




```{r}
Idents(fibroblasts_har)<-'orig.ident'
all.markers<-FindAllMarkers(fibroblasts_har ,only.pos = T)
all.markers<-all.markers[all.markers$p_val_adj < 0.05,]

dotplot<-DotPlot(fibroblasts_har, features = unique(all.markers$gene))

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


DimPlot(fibroblasts_har)

```




