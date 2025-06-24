


```{r}
library(ArchR)#needs r4.1
library(Nebulosa)#needs r4.1
library(scMEGA)#needs to be r4.1
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
#setwd("/rds/projects/c/croftap-stia-atac/CM_multiome/STIA_andATAC/")
options(bitmapType='cairo')
library(CellChat)
library(chromVAR)
library(pheatmap)
library(viridis)
```

```{r}
metadata.tsv <- read.delim("/rds/projects/m/mahonyc-cesar-data/all_cells/metadata.dir/metadata.dir/metadata.tsv.gz")

metadata.tsv <- metadata.tsv %>% remove_rownames %>% column_to_rownames(var="barcode_id.1")
metadata.tsv$barcode_id<-NULL

umap.0.3.tsv <- read.delim("/rds/projects/m/mahonyc-cesar-data/all_cells/out.30.comp.dir/out.30.comp.dir/umap.dir/umap.0.3.tsv.gz", row.names = 3)

library(SeuratDisk)
Convert("/rds/projects/m/mahonyc-cesar-data/all_cells/out.30.comp.dir/out.30.comp.dir/cellxgene_nMD.h5ad", "/rds/projects/m/mahonyc-cesar-data/all_cells/out.30.comp.dir/out.30.comp.dir/cellxgene.h5seurat")
seuratObject <- LoadH5Seurat("/rds/projects/m/mahonyc-cesar-data/all_cells/out.30.comp.dir/out.30.comp.dir/cellxgene.h5seurat")


library(reticulate)
sc <- import("scanpy")

atlas.data <- sc$read_h5ad("/rds/projects/m/mahonyc-cesar-data/all_cells/out.30.comp.dir/out.30.comp.dir/cellxgene.h5ad")

counts <- t(atlas.data$X)
colnames(counts) <-  atlas.data$obs_names$to_list()
rownames(counts) <-  atlas.data$var_names$to_list()
counts <- Matrix::Matrix(as.matrix(counts), sparse = T)

head(colnames(counts))
head(rownames(counts))

 meta_data_all_cells <- read.csv("/rds/projects/m/mahonyc-cesar-data/all_cells/out.30.comp.dir/out.30.comp.dir/meta_data_all_cells.csv", row.names=1)

library(loomR)  
lfile <- connect(filename = "/rds/projects/m/mahonyc-cesar-data/all_cells/loom.dir/loom.dir/X.loom", mode = "r+", skip.validate = T)
lfile

#data.subset <- lfile[["matrix"]]

all_cells <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0, meta.data = meta_data_all_cells)

umap.0.3.tsv<-as.matrix(umap.0.3.tsv)

all_cells[['umap']] <- CreateDimReducObject(embeddings = umap.0.3.tsv, key = "UMAP_", global = T, assay = "RNA")

all_cells <- all_cells %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) 

DimPlot(all_cells, group.by = "cluster_0.01")

table(all_cells$cluster_0.01)


umap_fibros <- read.delim("/rds/projects/m/mahonyc-cesar-data/umap_fibros.gz", row.names=3)

umap_fibros$ID<-'fibro'
umap_fibros$UMAP_1<-NULL
umap_fibros$UMAP_2<-NULL

all_cells_meta<-all_cells@meta.data
all_cells_meta<-select(all_cells_meta, c("cluster_0.01", "cluster_0.005"))

all_cells_meta_sub<-all_cells_meta[!rownames(all_cells_meta) %in% rownames(umap_fibros),]
colnames(all_cells_meta_sub)<-c("cluster_0.01_fib", "cluster_0.005_fib")

umap_fibros$cluster_0.01_fib<-'fib'
umap_fibros$cluster_0.005_fib<-'fib'
umap_fibros$ID<-NULL

meta_new<-rbind(all_cells_meta_sub, umap_fibros)

all_cells<-AddMetaData(all_cells, meta_new)

DimPlot(all_cells, group.by = "cluster_0.01_fib")

cols <- ArchR::paletteDiscrete(all_cells@meta.data[, "cluster_0.01_fib"])

DimPlot(all_cells, group.by = "cluster_0.01_fib", cols=cols)


cols <- ArchR::paletteDiscrete(all_cells@meta.data[, "condition"])

DimPlot(all_cells, group.by = "condition", cols=cols)

```

```{r}
counts_all_genes<-Read10X("/rds/projects/m/mahonyc-cesar-data/all_cells/GEX.mtx.full.dir-20230711T211614Z-001/GEX.mtx.full.dir")

nrow(counts_all_genes)




all_allgenes<-CreateSeuratObject(counts_all_genes, min.cells = 0, min.features = 0, meta.data = all_cells@meta.data)

all_allgenes[['umap']] <- CreateDimReducObject(embeddings = umap.0.3.tsv, key = "UMAP_", global = T, assay = "RNA")

all_allgenes <- all_allgenes %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) 

DimPlot(all_allgenes, group.by = "cluster_0.01")
max(all_allgenes$nFeature_RNA)
min(all_allgenes$nFeature_RNA)

all_allgenes_check<-CreateSeuratObject(counts_all_genes, min.cells = 0, min.features = 0)
max(all_allgenes_check$nFeature_RNA)
min(all_allgenes_check$nFeature_RNA)

all_allgenes$nFeature_RNA<-all_allgenes_check$nFeature_RNA
min(all_allgenes$pct_mitochondrial)

max(all_allgenes$nCount_RNA)
max(all_allgenes_check$nCount_RNA)

all_allgenes$nCount_RNA<-all_allgenes_check$nCount_RNA
rm(all_allgenes_check)

cols <- ArchR::paletteDiscrete(all_cells@meta.data[, "cluster_0.01_fib"])

DimPlot(all_allgenes, group.by = "cluster_0.01_fib", cols=cols)


```



```{r}

FeaturePlot(all_cells, features = c("Chad", "Chadl", "Sox9", "Cilp"))

FeaturePlot(all_cells, features = c("Lyve1", "Acta2", "Sox9", "Cilp") )

DotPlot(all_cells, features = c("Myoz1", "Actn3", "Aldoa",   "Cd248", "Thy1", "Fap",  "Chadl", "Chad", "Sox9", "Cdh5", "Pecam1", "Flt1", "S100a9", "Lcp1", "Coro1a", "Mpz", "Kcna1", "Mbp", "F2r", "Erfe", "H2-Q4"))+RotatedAxis()


DotPlot(all_cells, features = c("Myoz1", "Actn3", "Aldoa", "Pdgfra",   "Cd248", "Thy1", "Fap",  "Chadl", "Chad", "Sox9", "Cdh5", "Pecam1", "Flt1", "S100a9", "Lcp1", "Coro1a", "Mpz", "Kcna1", "Mbp", "F2r", "Erfe", "H2-Q4")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") + scale_size(range = c(2, 8)) +RotatedAxis()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))

Idents(all_cells)<-'cluster_0.01_fib'

marker<-FindAllMarkers(all_cells, only.pos = T, logfc.threshold = 1)


FeaturePlot(all_cells, features = c("Notch3", "Acta2", "Rgs5") )


FeaturePlot(all_cells, features = c("Clic5", "Col22a1", "Tspan15") )

Idents(all_allgenes)<-'cluster_0.01_fib'
DotPlot(all_allgenes, features = c("Myoz1", "Actn3", "Aldoa","Col1a1", "Pdgfra",   "Cd248", "Thy1", "Fap",  "Chadl", "Chad", "Sox9", "Cdh5", "Pecam1", "Flt1", "S100a9", "Lcp1", "Coro1a", "Mpz", "Kcna1", "Mbp", "F2r", "Erfe", "H2-Q4")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") + scale_size(range = c(2, 8)) +RotatedAxis()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))


```


```{r}

metadata.tsv <- read.delim("/rds/projects/m/mahonyc-cesar-data/all_cells/all_cells_full_data/metadata.tsv.gz")


library(tidyverse)
metadata.tsv <- metadata.tsv %>% remove_rownames %>% column_to_rownames(var="barcode_id.1")
metadata.tsv$barcode_id<-NULL

umap.0.1.tsv <- read.delim("/rds/projects/m/mahonyc-cesar-data/all_cells/all_cells_full_data/umap.0.1.tsv.gz", row.names = 3)

umap.0.3.tsv <- read.delim("/rds/projects/m/mahonyc-cesar-data/all_cells/all_cells_full_data/umap.0.3.tsv.gz", row.names = 3)

umap.0.5.tsv <- read.delim("/rds/projects/m/mahonyc-cesar-data/all_cells/all_cells_full_data/umap.0.5.tsv.gz", row.names = 3)

umap.0.7.tsv <- read.delim("/rds/projects/m/mahonyc-cesar-data/all_cells/all_cells_full_data/umap.0.7.tsv.gz", row.names = 3)

counts<-Read10X(data.dir="/rds/projects/m/mahonyc-cesar-data/all_cells/all_cells_full_data/GEX.mtx.full.dir-20230630T084154Z-001/GEX.mtx.full.dir")


ncol(counts)
nrow(counts)


all_cells<-CreateSeuratObject(counts = counts, meta.data = metadata.tsv)

umap.0.3.tsv<-as.matrix(umap.0.3.tsv)

all_cells[['umap']] <- CreateDimReducObject(embeddings = umap.0.3.tsv, key = "UMAP_", global = T, assay = "RNA")

DimPlot(all_cells)

all_cells <- all_cells %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) 

all_cells<-FindNeighbors(all_cells, dims = 1:30)
all_cells<-FindClusters(all_cells, resolution = c(0.05, 0.1, 0.2))
DimPlot(all_cells, group.by = "RNA_snn_res.0.2")
all_cells_f=all_cells

```

```{r}
load("/rds/projects/m/mahonyc-cesar-data/all_cells/analysis_CM.RData")
ncol(all_cells)
meta_old<-all_cells@meta.data
meta_old<-subset(meta_old, select=c(cluster_0.01_fib))

all_cells_f<-AddMetaData(all_cells_f, meta_old)

cols <- ArchR::paletteDiscrete(all_cells_f@meta.data[, "cluster_0.01_fib"])
cols <- ArchR::paletteDiscrete(all_cells_f@meta.data[, "condition"])

DimPlot(all_cells_f, group.by = "condition", cols=cols)+NoAxes()

ncol(all_cells_f)

```

```{r}

Idents(all_cells_f)<-'cluster_0.01_fib'
DotPlot(all_cells_f, features = c( "Pdgfra", "Cd248", "Col1a1", "Acta2", "Chadl", "Chad", "Sox9","Runx2", "Omd", "Bglap",  "Cdh5", "Pecam1", "Flt1",  "Myoz1", "Actn3", "Aldoa")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") + scale_size(range = c(2, 8)) +RotatedAxis()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))



```
```{r}

FeaturePlot(all_cells_f, features = c("Pdgfra", "Cd248", "Col1a1", "Acta2"))

FeaturePlot(all_cells_f, features = c( "Chadl", "Chad", "Sox9","Runx2") )


FeaturePlot(all_cells_f, features = c("Omd", "Bglap",  "Cdh5", "Pecam1") )


FeaturePlot(all_cells_f, features = c("Flt1",  "Myoz1", "Actn3", "Aldoa") )




p <- FeaturePlot(all_cells_f, c("Pdgfra", "Cd248", "Col1a1", "Acta2"), combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + NoAxes()
}


p <- FeaturePlot(all_cells_f, c( "Chadl", "Chad", "Sox9","Runx2"), combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + NoAxes()
}
cowplot::plot_grid(plotlist = p)
p <- FeaturePlot(all_cells_f, c( "Omd", "Bglap",  "Cdh5", "Pecam1"), combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + NoAxes()
}
cowplot::plot_grid(plotlist = p)
p <- FeaturePlot(all_cells_f, c( "Flt1",  "Myoz1", "Actn3", "Aldoa"), combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + NoAxes()
}
cowplot::plot_grid(plotlist = p)

```


