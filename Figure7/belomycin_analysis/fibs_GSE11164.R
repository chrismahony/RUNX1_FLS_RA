
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
control_1_1 <- CreateSeuratObject(
  counts = GSM3036808_Control_1_1_Mouse_lung_digital_gene_expression_400.dge,
  assay = "RNA"
)

control_1_1 <- control_1_1 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

control_1_2 <- CreateSeuratObject(
  counts = control_1_2,
  assay = "RNA"
)

control_1_2 <- control_1_2 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)


Control_2 <- CreateSeuratObject(
  counts = Control_2,
  assay = "RNA"
)

Control_2 <- Control_2 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

Control_3 <- CreateSeuratObject(
  counts = Control_3,
  assay = "RNA"
)

Control_3 <- Control_3 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

Control_4 <- CreateSeuratObject(
  counts = Control_4,
  assay = "RNA"
)

Control_4 <- Control_4 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

Control_5 <- CreateSeuratObject(
  counts = Control_5,
  assay = "RNA"
)

Control_5 <- Control_5 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

Control_6 <- CreateSeuratObject(
  counts = Control_6,
  assay = "RNA"
)

Control_6 <- Control_6 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)



Bleo_1 <- CreateSeuratObject(
  counts = Bleo_1,
  assay = "RNA"
)

Bleo_1 <- Bleo_1 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)


Bleo_2 <- CreateSeuratObject(
  counts = Bleo_2,
  assay = "RNA"
)

Bleo_2 <- Bleo_2 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

Bleo_3 <- CreateSeuratObject(
  counts = Bleo_3,
  assay = "RNA"
)

Bleo_3 <- Bleo_3 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)


control_1_1$condition<-'control'
control_1_2$condition<-'control'
Control_2$condition<-'control'
Control_3$condition<-'control'
Control_4$condition<-'control'
Control_5$condition<-'control'
Control_6$condition<-'control'
Bleo_1$condition<-'bleo'
Bleo_2$condition<-'bleo'
Bleo_3$condition<-'bleo'

```
```{r}
list=c(control_1_1, control_1_2, Control_2, Control_3, Control_4, Control_5, Control_6, Bleo_1, Bleo_2, Bleo_3)
anchors <- FindIntegrationAnchors(object.list = list, dims = 1:30)
aggr <- IntegrateData(anchorset = anchors, dims = 1:30)
aggr <- ScaleData(aggr, verbose = FALSE)
aggr <- RunPCA(aggr, verbose = FALSE)
aggr <- RunUMAP(aggr, dims = 1:50)
```

```{r}
DefaultAssay(aggr)<-'RNA'
FeaturePlot(aggr, features = "Col1a1")
FeaturePlot(aggr, features = "Col1a2")

aggr<-FindNeighbors(aggr, dims = 1:30)
aggr<-FindClusters(aggr, resolution = 0.2, graph.name = "integrated_snn")
DimPlot(aggr, label=T)

fibs_GSE11164<- aggr[,grepl("8", aggr$integrated_snn_res.0.2, ignore.case=TRUE)]
fibs_GSE11164<- fibs_GSE11164[,!grepl("18", fibs_GSE11164$integrated_snn_res.0.2, ignore.case=TRUE)]
DimPlot(fibs_GSE11164, label=T)
saveRDS(fibs_GSE11164, "/rds/projects/m/mahonyc-kitwong-runx1/GSE111664_Lung_bleomycin/fibs_GSE11164.rds")

```


