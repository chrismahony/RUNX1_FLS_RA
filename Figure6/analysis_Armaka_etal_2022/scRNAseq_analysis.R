

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
library(monocle)
setwd("/rds/projects/c/croftap-stia-atac/CM_multiome/STIA_andATAC/")
options(bitmapType='cairo')
library(CellChat)
#load()
```

```{r}
wt4_counts <- Read10X(data.dir = "/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/data_from_koliaspaper/scRNAseq/wt4_count/outs/filtered_feature_bc_matrix")
wt4 <- CreateSeuratObject(counts = wt4_counts, assay = "RNA", min.cells = 3, min.features = 0)
rm(wt4_counts)
wt4$sample <- 'wt4'
wt4[["percent.mt"]] <- PercentageFeatureSet(wt4, pattern = "^mt-")
VlnPlot(wt4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
wt4 <- subset(wt4, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)


wt4 <- wt4 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

```

```{r}


tg4_counts <- Read10X(data.dir = "/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/data_from_koliaspaper/scRNAseq/tg4_count/outs/filtered_feature_bc_matrix")
tg4 <- CreateSeuratObject(counts = tg4_counts, assay = "RNA", min.cells = 3, min.features = 0)
rm(tg4_counts)
tg4$sample <- 'tg4'
tg4[["percent.mt"]] <- PercentageFeatureSet(tg4, pattern = "^mt-")
VlnPlot(tg4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
tg4 <- subset(tg4, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)


tg4 <- tg4 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
```

```{r}
tg8_counts <- Read10X(data.dir = "/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/data_from_koliaspaper/scRNAseq/tg8_count/outs/filtered_feature_bc_matrix")
tg8 <- CreateSeuratObject(counts = tg8_counts, assay = "RNA", min.cells = 3, min.features = 0)
rm(tg8_counts)
tg8$sample <- 'tg8'
tg8[["percent.mt"]] <- PercentageFeatureSet(tg8, pattern = "^mt-")
VlnPlot(tg8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
tg8 <- subset(tg8, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)


tg8 <- tg8 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
```
```{r}

reference.list <- c(tg4, wt4, tg8)
anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
aggr <- IntegrateData(anchorset = anchors, dims = 1:30)

aggr <- aggr %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

DimPlot(aggr, group.by = 'sample')
```
```{r}
library(harmony)
aggr_safe<-aggr
my_harmony_embeddings <- HarmonyMatrix(
  data_mat  = as.matrix(aggr@reductions$umap@cell.embeddings),
  meta_data = aggr@meta.data,
  vars_use  = "sample",
  do_pca = FALSE
)
rownames(my_harmony_embeddings) <- rownames(aggr@reductions$umap@cell.embeddings)

#store the harmony reduction as a custom dimensional reduction called 'harmony' in the default assay
aggr[["harmony"]] <- CreateDimReducObject(embeddings = my_harmony_embeddings, key = "harmony_", assay = DefaultAssay(aggr))
aggr <- FindNeighbors(object = aggr, reduction = "harmony", dims = 1:2)
aggr <- FindClusters(object = aggr, verbose = TRUE, algorithm = 1) # Louvain algorithm
aggr <- RunUMAP(object = aggr, reduction = "harmony", dims = 1:2)
DimPlot(aggr_safe, group.by = "orig.ident", reduction = "umap")
DimPlot(aggr_safe, split.by = "orig.ident")
DimPlot(aggr, group.by = "sample", reduction = "harmony")
```
```{r}
#batch correct no good!
rm(aggr)
aggr=aggr_safe

aggr<-RunUMAP(aggr, dims = 1:10)
DimPlot(aggr, group.by = 'sample')


```

```{r}
aggr <- FindNeighbors(object = aggr, dims = 1:20)
aggr <- FindClusters(object = aggr, verbose = TRUE, resolution = c(0.05, 0.1, 0.2, 0.3))
library(clustree)
clustree(aggr)
```
```{r}
aggr <- FindClusters(object = aggr, verbose = TRUE, resolution = c(0.02, 0.04))
clustree(aggr)

DimPlot(aggr, group.by = 'integrated_snn_res.0.1', label=T)



```
```{r}

Idents(aggr)<-'integrated_snn_res.0.1'
aggr$seurat_new <- aggr@active.ident
current.sample.ids <- c("0", "1", "2", "3", "4", "5", "6", "7")
new.sample.ids <- c("0", "1", "2", "3", "4", "5", "0", "4")
aggr@meta.data[["seurat_new"]] <- plyr::mapvalues(x = aggr@meta.data[["seurat_new"]], from = current.sample.ids, to = new.sample.ids)

DimPlot(aggr, group.by = 'seurat_new', label = F)

```
```{r}
Idents(aggr)<-'seurat_new'
all_markers<-FindAllMarkers(aggr, only.pos = T)
```
```{r}
DotPlot(aggr, features=c("Pi16", "Ccl11","Gsn", "Pdpn", "F13a1", "Prg4", "Gas6","Sfrp2", "Comp", "Mgp","C7", "Dkk3", "Myod1", "Myog","Des", "Thbs4", "Fmod", "Tnmd"))+RotatedAxis()+coord_flip()
```
```{r}
DotPlot(aggr, features=c("Pi16", "Ccl11", "Pdpn", "Prg4", "Gas6", "Comp", "Mgp", "Dkk3", "Myod1", "Myog", "Thbs4", "Fmod"))+RotatedAxis()+coord_flip()
```


```{r}

```

```{r}

aggr$bulk <- aggr$seurat_new
current.sample.ids <- c("0", "1", "2", "3", "4", "5")
new.sample.ids <- c("SL", "LL", "SL", "SL", "peri", "SL")
aggr@meta.data[["bulk"]] <- plyr::mapvalues(x = aggr@meta.data[["bulk"]], from = current.sample.ids, to = new.sample.ids)
DimPlot(aggr, group.by = 'bulk', label = F)




```
```{r}
Idents(aggr)<-'bulk'
FeaturePlot(aggr, features="F13a1")
p1<-FeaturePlot(aggr, features="Prg4")
aggr<-CellSelector(aggr, plot = p1, ident="lining")
aggr$bulk<-aggr@active.ident
DimPlot(aggr, group.by = 'bulk', label = F)

Idents(aggr)<-'bulk'
levels(aggr)
current.sample.ids <- c("lining", "SL" ,    "LL"  ,   "peri" )
new.sample.ids <- c("lining", "SL" ,    "SL"  ,   "peri")
aggr@meta.data[["bulk"]] <- plyr::mapvalues(x = aggr@meta.data[["bulk"]], from = current.sample.ids, to = new.sample.ids)
DimPlot(aggr, group.by = 'bulk', label = F)

DimPlot(aggr, group.by = 'bulk', label = F)

```

```{r}
DefaultAssay(aggr)<-'RNA'
aggr$bulk_sample<-paste(aggr$bulk, aggr$sample, sep="_")
Idents(aggr)<-'bulk_sample'

levels(aggr)<-c("SL_wt4","SL_tg4","SL_tg8","lining_wt4","lining_tg4","lining_tg8","peri_wt4","peri_tg4","peri_tg8")

VlnPlot(aggr, features = "Runx1", idents = c("SL_wt4", "SL_tg4", "SL_tg8"), pt.size = 0.01)
VlnPlot(aggr, features = "Runx1", idents = c("lining_wt4", "lining_tg4", "lining_tg8"), pt.size = 0.01)
VlnPlot(aggr, features = "Runx1", idents = c("peri_wt4", "peri_tg4", "peri_tg8"), pt.size = 0.01)


SL_tg4_wt4_markers<-FindMarkers(aggr, ident.1 = "SL_tg4", ident.2 = "SL_wt4")
SL_tg8_wt4_markers<-FindMarkers(aggr, ident.1 = "SL_tg8", ident.2 = "SL_wt4")

LL_tg4_wt4_markers<-FindMarkers(aggr, ident.1 = "LL_tg4", ident.2 = "LL_wt4")
LL_tg8_wt4_markers<-FindMarkers(aggr, ident.1 = "LL_tg8", ident.2 = "LL_wt4")



DotPlot(aggr, features = "Runx1", idents = c("SL_wt4", "SL_tg4", "SL_tg8"))
DotPlot(aggr, features = "Runx1", idents = c("LL_wt4", "LL_tg4", "LL_tg8"))

  
                      

```


```{r}

Idents(aggr)<-'bulk'
all_markers_bulk<-FindAllMarkers(aggr, only.pos = T)

DotPlot(aggr, features=c("F13a1", "Clic5", "Tspan15", "Prg4", "Apod", "Igf1", "Lum", "Cd34", "Des", "Myod1", "Cdh15", "Actn3"))+RotatedAxis()+coord_flip()

```



```{r}
aggr$bulk_sample2 <- aggr$bulk_sample
current.sample.ids <- c("peri_tg4", "SL_tg4"  , "LL_tg4" ,  "SL_wt4" ,  "peri_wt4" ,"LL_wt4" ,  "SL_tg8" ,  "LL_tg8" ,  "peri_tg8")
new.sample.ids <- c("peri_d", "SL_d"  , "LL_d" ,  "SL_wt4",   "peri_wt4", "LL_wt4",   "SL_d"  , "LL_d" ,  "peri_d")
aggr@meta.data[["bulk_sample2"]] <- plyr::mapvalues(x = aggr@meta.data[["bulk_sample2"]], from = current.sample.ids, to = new.sample.ids)

Idents(aggr)<-'bulk_sample2'
DefaultAssay(aggr)<-'RNA'
SL_d_wt4_markers<-FindMarkers(aggr, ident.1 = "SL_d", ident.2 = "SL_wt4")
LL_d_wt4_markers<-FindMarkers(aggr, ident.1 = "LL_d", ident.2 = "LL_wt4")
SL_d_wt4_markers$gene<-rownames(SL_d_wt4_markers)
LL_d_wt4_markers$gene<-rownames(LL_d_wt4_markers)

levels(aggr)<-c("SL_wt4" ,  "peri_wt4" ,"LL_wt4","peri_d",   "SL_d"  ,   "LL_d"     )

VlnPlot(aggr, features="Runx1", idents = c("SL_d", "SL_wt4"))
VlnPlot(aggr, features="Runx1", idents = c("LL_d", "LL_wt4"))


```

