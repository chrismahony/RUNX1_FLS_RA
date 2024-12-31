
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

#https://rnabio.org/module-08-scrna/0008/02/01/scRNA/


#setwd first!
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

pct <- aggr[["pca"]]@stdev / sum(aggr[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2



DimPlot(aggr, group.by = "orig.ident")
FeaturePlot(aggr, features = c("Col1a1", "Pdgfra"))


```
```{r}
Idents(aggr)<-'1'
p1<-FeaturePlot(aggr, features="Col1a1")
aggr<-CellSelector(p1, aggr, ident="fibs")
aggr$fibs_selected<-aggr@active.ident
Idents(aggr)<-'orig.ident'
levels(aggr)
aggr$condition<-aggr$orig.ident
current.sample.ids <- c( "bleo_N1", "bleo_N2", "bleo_N3", "bleo1" ,  "bleo2" ,  "bleo3"  , "bleo4" ,  "saline1", "saline2", "saline3", "saline4")
new.sample.ids <- c( "bleo_N", "bleo_N", "bleo_N", "bleo" ,  "bleo" ,  "bleo"  , "bleo" ,  "saline", "saline", "saline", "saline")
aggr@meta.data[["condition"]] <- plyr::mapvalues(x = aggr@meta.data[["condition"]], from = current.sample.ids, to = new.sample.ids)

DimPlot(aggr, group.by = "fibs_selected")

aggr$fibs_conditon<-paste(aggr$fibs_selected, aggr$condition, sep="_")
table(aggr$fibs_conditon)
Idents(aggr)<-'fibs_conditon'
DefaultAssay(aggr)<-'RNA'
fibs_belo_saline<-FindMarkers(aggr, ident.1 = "fibs_bleo", ident.2 = "fibs_saline")
fibs_belo_saline$gene<-rownames(fibs_belo_saline)
Idents(aggr)<-'fibs_conditon'
levels(aggr)<-c("1_bleo_N"  ,  "fibs_bleo_N", "1_bleo"   ,     "1_saline" ,   "fibs_saline", "fibs_bleo")
VlnPlot(aggr, features = "Runx1", idents = c("fibs_saline", "fibs_bleo"))



DotPlot(aggr, features = c("Runx1", "Mmp14", "Igf1", "Cthrc1"), idents = c("fibs_saline", "fibs_bleo"))




```




```{r}
network <- decoupleR::get_dorothea(organism = "human",
                                   levels = c("A", "B", "C"))

fibs <- aggr[,grepl("fibs", aggr$fibs_selected, ignore.case=TRUE)]

fibs <- fibs %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

activities <- decoupleR::run_wmean(mat = as.matrix(fibs@assays[['RNA']]@data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "mor",
                                   times = 100,
                                   minsize = 1)
Idents(fibs)<-'condition'
out <- SCpubr::do_TFActivityPlot(sample = fibs,
                                 activities = activities)
p <- out$heatmaps$average_scores
p

saveRDS(fibs, "/rds/projects/m/mahonyc-kitwong-runx1/GSE129605_lung_bleo/fibs_GSE129605.rds")
```

