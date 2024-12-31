
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
library(OmnipathR)
library(dorothea)
library(decoupleR)
library(SCpubr)
```


```{r}

setwd("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE132771_lung_fibrosis")
data.10x = list()
dirs <- list.dirs(".", recursive = FALSE)

for (i in 1:length(dirs)) {
    data.10x[[i]] <- Read10X(data.dir = dirs[[i]])
}
names(data.10x) <- sub('./', '', dirs)
  

scrna.list = list()
samples<-sub('./', '', dirs)

for (i in 1:length(data.10x)) {
    scrna.list[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=3, min.features=200, project=samples[i]);
    scrna.list[[i]] =NormalizeData(object = scrna.list[[i]]);
    scrna.list[[i]] =ScaleData(object = scrna.list[[i]]);
    scrna.list[[i]] =FindVariableFeatures(object = scrna.list[[i]]);
    scrna.list[[i]] =RunPCA(object = scrna.list[[i]], verbose = FALSE)
    scrna.list[[i]][["percent.mt"]] = PercentageFeatureSet(object=scrna.list[[i]], pattern = "^mt-")
    }
names(scrna.list) <- sub('./', '', dirs)


anchors <- FindIntegrationAnchors(object.list = scrna.list, reduction = "rpca",   dims = 1:50)
aggr <- IntegrateData(anchorset = anchors, dims = 1:50)
aggr <- ScaleData(aggr, verbose = FALSE)
aggr <- RunPCA(aggr, verbose = FALSE)
aggr <- RunUMAP(aggr, dims = 1:50)

DimPlot(aggr)

FeaturePlot(aggr, features="Col1a1")
FeaturePlot(aggr, features="Pdgfra")
VlnPlot(aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
aggr <- subset(aggr, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)


aggr <- ScaleData(aggr, verbose = FALSE)
aggr <- RunPCA(aggr, verbose = FALSE)
aggr <- RunUMAP(aggr, dims = 1:50)
DimPlot(aggr)
FeaturePlot(aggr, features="Col1a1")
FeaturePlot(aggr, features="Pdgfra")

aggr<-FindNeighbors(aggr, dims = 1:30)
aggr<-FindClusters(aggr, resolution = 0.2)
DimPlot(aggr, label=T)

fibs_GSE132771<- aggr[,grepl("1|0|3", aggr$integrated_snn_res.0.2, ignore.case=TRUE)]
fibs_GSE132771<- fibs_GSE132771[,!grepl("11|10|13|12|14", fibs_GSE132771$integrated_snn_res.0.2, ignore.case=TRUE)]
DimPlot(fibs_GSE132771, label=T)
saveRDS(fibs_GSE132771, "/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE132771_lung_fibrosis/fibs_GSE132771.rds")


```



