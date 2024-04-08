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
```

```{r}
counts <- Read10X_h5("/rds/projects/c/croftap-stia-atac/2021_Timecourse_wATACseq/ALL_aggr_MACS3peaks/outs/filtered_feature_bc_matrix.h5")

samples_ID <- read.csv(file.path("/rds/projects/c/croftap-stia-atac/2021_Timecourse_wATACseq/ALL_aggr_MACS3peaks/outs/", "LibraryID.csv"))


fragpath <- "/rds/projects/c/croftap-stia-atac/2021_Timecourse_wATACseq/ALL_aggr_MACS3peaks/outs/atac_fragments.tsv.gz"

#annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#annotation taken from load("/rds/projects/c/croftap-stia-atac/CM_multiome/MainAnalysis/old analysis/CM_analysis.RData")
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"



chrom_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

aggr <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC")

aggr$sampleID=(samples_ID[match(rownames(aggr@meta.data),samples_ID$Barcode),2])

head(aggr$sampleID)
table(aggr$sampleID)
table(samples_ID[,2])
aggr$orig.ident <- aggr$sampleID
table(aggr@meta.data[["orig.ident"]])
aggr$sampleID <- NULL

```

```{r}

aggr <- NucleosomeSignal(object = aggr)
aggr <- TSSEnrichment(object = aggr, fast = FALSE)
VlnPlot(aggr, c("TSS.enrichment", "nCount_ATAC", "nucleosome_signal"), pt.size = 0, ncol=4)

```
```{r}

aggr <- subset(
  x = aggr,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    TSS.enrichment > 3 & 
    nucleosome_signal < 4
)
VlnPlot(aggr, c("TSS.enrichment", "nCount_ATAC", "nucleosome_signal"), pt.size = 0, ncol = 4)

aggr$all<-'1'
Idents(aggr)<-'all'
VlnPlot(aggr, c("TSS.enrichment", "nCount_ATAC", "nucleosome_signal"), pt.size = 0.0, ncol = 3)


```

```{r}
aggr$blacklist_fraction <- FractionCountsInRegion(
  object = aggr, 
  assay = 'ATAC',
  regions = blacklist_mm10
)

Idents(aggr)<-'all'
VlnPlot(
  object = aggr,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 4
)


VlnPlot(
  object = aggr,
  features = c('blacklist_fraction', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 2
)

VlnPlot(
  object = aggr,
  features = c('nCount_ATAC', 'TSS.enrichment'),
  pt.size = 0.1,
  ncol = 2
)
```



```{r}
aggr <- RunTFIDF(aggr)
aggr <- FindTopFeatures(aggr, min.cutoff = 20)
aggr <- RunSVD(aggr)

```

```{r}
aggr <- FindNeighbors(object = aggr, reduction = 'lsi', dims = 2:50)
aggr <- RunUMAP(aggr, dims = 2:50, reduction = 'lsi')
aggr <- RunTSNE(aggr, dims = 2:50, reduction = 'lsi')

p1 <- DimPlot(aggr, group.by = "orig.ident")
p2 <- DimPlot(aggr, group.by = "orig.ident", reduction = "tsne")
plot_grid(p1, p2)
```

```{r}
library(harmony)

aggr <- RunHarmony(aggr, group.by.vars = "orig.ident", 
                         reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE)

aggr <- FindNeighbors(aggr, dims = 2:50, reduction = "harmony", verbose = FALSE)
                                         
aggr <- RunUMAP(aggr, dims = 2:50, reduction = "harmony")
aggr <- RunTSNE(aggr, dims = 2:50, reduction = "harmony")


p3 <- DimPlot(aggr, group.by = "orig.ident")
p4 <- DimPlot(aggr, group.by = "orig.ident", reduction = "tsne")
plot_grid(p3, p4)

#findclusters
aggr@graphs

for (res in c(0.05, 0.1)) {
   aggr <- FindClusters(aggr, resolution = res, algorithm = 3, graph.name = "ATAC_snn")
                          
}

library(clustree)
clustree(aggr, prefix ="ATAC_snn_res.")


for (res in c(0.1, 0.2, 0.3, 0.4)) {
   aggr <- FindClusters(aggr, resolution = res, algorithm = 3, graph.name = "ATAC_nn")
                          
}

clustree(aggr, prefix ="ATAC_nn_res.")
```
```{r}
aggr <- FindClusters(aggr, resolution = 0.3, algorithm = 3, graph.name = "ATAC_snn")
p5 <- DimPlot(aggr)
p6 <- DimPlot(aggr, reduction = "tsne")
plot_grid(p5, p6)           
```
```{r}
aggr <- FindClusters(aggr, resolution = 0.2, algorithm = 3, graph.name = "ATAC_nn")
p7 <- DimPlot(aggr)
p8 <- DimPlot(aggr, reduction = "tsne")
plot_grid(p7, p8)    
```
```{r}

library(cicero)
library(SeuratWrappers)
DefaultAssay(aggr) <- "ATAC"
aggr.cds <- as.cell_data_set(x = aggr)
umap_coords <- reducedDims(aggr.cds)$UMAP
cicero_cds <- make_cicero_cds(aggr.cds, reduced_coordinates = umap_coords)
genome <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns <- run_cicero(cicero_cds, genomic_coords = genome.df, sample_num = 100)
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(aggr) <- links


gene.activities <- GeneActivity(aggr)
aggr[['ATACexpr']] <- CreateAssayObject(counts = gene.activities)
aggr <- NormalizeData(
  object = aggr,
  assay = 'ATACexpr',
  normalization.method = 'LogNormalize',
  scale.factor = median(aggr$nCount_ATACexpr)
)

DefaultAssay(aggr) <- 'ATAC'
Idents(aggr)<-'ATAC_nn_res.0.2'
all.markers<-FindAllMarkers(aggr)
all.markers_repeat <- FindAllMarkers(aggr, min.pct = 0.05, test.use = 'LR',  latent.vars = 'nCount_ATAC', only.pos = T)

Idents(aggr)<-'cm_clusters'
all.markers_cm_clusters <- FindAllMarkers(aggr, min.pct = 0.05, test.use = 'LR',  latent.vars = 'nCount_ATAC', only.pos = T)

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14 Mono",
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

allmarkerspeaks <- rownames(all.markers)
closest_allmarkerspeaks <- ClosestFeature(aggr, regions = allmarkerspeaks)
all.markers$gene_name <- closest_allmarkerspeaks$gene_name 

DefaultAssay(aggr) <- 'ATACexpr'
all.markers_ATACexpr <- FindAllMarkers(aggr, logfc.threshold = 0.15)
```

```{r}
DefaultAssay(aggr) <- 'ATACexpr'
VlnPlot(aggr, features = c("Cd248", "Cd34", "Cdh11", "Col1a1", "Thy1", "Pdgfra", "Pdpn", "Fap"), stack = T)
VlnPlot(aggr, features=c("Acta2", "Des", "Flt1", "Mcam", "Notch3", "Pdgfrb", "Rgs5"), stack = T)
VlnPlot(aggr, features=c("Actn3", "Aldoa", "Tnnt3"), stack = T)
VlnPlot(aggr, features = c("Alpl", "Bglap", "Bglap2", "Omd", "Runx2", "Sp7", "Ostn"), stack = T, pt.size = 0.01)
VlnPlot(aggr, features = c("Sox6", "Cd14", "Chad", "Chadl", "Cilp", "Clu", "Sox9", "Matn3"), stack = T)
VlnPlot(aggr, features = c("Cd55", "Clic5", "Col22a1", "Hbegf", "Htra4", "Tspan15"), stack = T)
VlnPlot(aggr, features = c("Kcna1", "Plp1"), stack = T)
VlnPlot(aggr, features = c("Hbb-bs", "Alas2"), stack = T)
```


```{r}
#1-fibroblasts
StackedVlnPlot(aggr, features = c("Col14a1", "Gas7"))

#2=linning
StackedVlnPlot(aggr, features = c("Col22a1", "F13a1"))
#3-osteoblast
StackedVlnPlot(aggr, features = c("Alpl", "Runx3"))
#4-chondrocytes?check Chd9 peak

#5-smooth muscle
StackedVlnPlot(aggr, features = c("Myh11", "Notch3", "Rbpms"))


#9 vascular
StackedVlnPlot(aggr, features = c("Flt1"))
 
 
```

```{r}
#linning layer
VlnPlot(aggr, features = c("Cd55", "Clic5", "Col22a1", "Hbegf", "Htra4", "Tspan15"))
StackedVlnPlot(aggr, features = c("Cd55", "Clic5", "Col22a1", "Hbegf", "Htra4", "Tspan15"))
```



```{r}
#osteoblasts
VlnPlot(aggr, features = c("Alpl", "Bglap", "Bglap2", "Omd", "Runx2", "Sp7", "Ostn"))
StackedVlnPlot(aggr, features = c("Alpl", "Bglap", "Bglap2", "Omd", "Runx2", "Sp7", "Ostn"))



```


```{r}
#chondorcytes
VlnPlot(aggr, features = c("Sox6", "Cd14", "Chad", "Chadl", "Cilp", "Clu", "Sox9", "Matn3"))
FeaturePlot(aggr, features = c("Sox6", "Cd14", "Chad", "Chadl", "Cilp", "Clu", "Sox9", "Matn3"))
StackedVlnPlot(aggr, features = c("Cd14", "Chad", "Chadl", "Cilp", "Clu", "Sox9", "Matn3"))

```
```{r}
#vascular
StackedVlnPlot(aggr, features = c("Cdh5", "Emcn", "Pecam1"))
```

```{r}
#pericytes
VlnPlot(aggr, features=c("Acta2", "Des", "Flt1", "Mcam", "Notch3", "Pdgfrb", "Rgs5"))
```


```{r}
VlnPlot(aggr, features=c("Actn3", "Aldoa", "Tnnt3"))
```

```{r}
#cell cylel
VlnPlot(aggr, features = c("Cdk1", "Cenpa", "Top2a", "Mki67"))
```

```{r}
VlnPlot(aggr, features = c("H2-Aa", "H2-Ab1"))


```

```{r}

aggr$cm_clusters <- aggr@active.ident
aggr$cm_clusters <- aggr$ATAC_nn_res.0.2
Idents(aggr) <- 'cm_clusters'
current.sample.ids <- c("0","1","2", "3", "4", "5","6","7", "8", "9", "10")
new.sample.ids <- c("fibro","fibro","Linning", "osteoblast", "chondrocyte", "pericyte", "contam1","contam2", "contam3", "vascular", "muscle")

aggr@meta.data[["cm_clusters"]] <- plyr::mapvalues(x = aggr@meta.data[["cm_clusters"]], from = current.sample.ids, to = new.sample.ids)

Idents(aggr) <- "cm_clusters"
levels(aggr) <- c("osteoblast", "fibro", "vascular", "Linning", "pericyte","chondrocyte", "muscle", "contam1", "contam2", "contam3")

DimPlot(aggr, reduction="umap", pt.size=0.2)

```
```{r}
#re-run
DefaultAssay(aggr) <- 'ATAC'
all.markers_cmclusters <- FindAllMarkers(aggr)
allmarkerspeaks <- rownames(all.markers_cmclusters)
closest_allmarkerspeaks <- ClosestFeature(aggr, regions = allmarkerspeaks)
all.markers_cmclusters$gene_name <- closest_allmarkerspeaks$gene_name 

DefaultAssay(aggr) <- 'ATACexpr'
all.markers_ATACexpr_cmclusters <- FindAllMarkers(aggr, logfc.threshold = 0.15)

DefaultAssay(aggr) <- 'ATAC'
all.markers_cmclusters_fc0.15 <- FindAllMarkers(aggr, logfc.threshold = 0.15, only.pos = T)
allmarkerspeaks <- rownames(all.markers_cmclusters_fc0.15)
closest_allmarkerspeaks <- ClosestFeature(aggr, regions = allmarkerspeaks)
all.markers_cmclusters_fc0.15$gene_name <- closest_allmarkerspeaks$gene_name

Idents(aggr)<-'cm_clusters'
all.markers_cmclusters_named <- FindAllMarkers(aggr, logfc.threshold = 0.15, only.pos = T)


DotPlot(aggr, features = c ("Gas7",   "Col14a1", "Pdgfra", "Col1a1", "Pdpn", "Col22a1", "F13a1", "Tspan15", "Alpl", "Runx2", "Runx3", "Anxa8", "Ank", "Chd9",  "Notch3", "Acta2", "Rgs5", "Flt1", "Cdh5", "Pecam1"), idents=c("fibro", "Linning", "osteoblast", "chondrocyte", "pericyte", "vascular")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


```

```{r}

Idents(aggr) <- "cm_clusters"
DefaultAssay(aggr) <- "ATAC"
CoveragePlot(
  object = aggr,
  region = "Col1a1",
    links = T
  
)


```
```{r}
#re-run
DefaultAssay(aggr) <- 'ATACexpr'

DotPlot(aggr, features=c("Col14a1", "F13a1", "Alpl", "Runx3", "Myh11", "Notch3", "Rbpms", "Flt1"))

```
```{r}

#re-run
library(viridis)
Idents(aggr) <- "cm_clusters"

DotPlot(aggr, features = c("Col14a1", "Tshz2", "Col22a1", "F13a1", "Alpl", "Runx3", "Anxa8", "Ank", "Notch3", "Rbpms", "Flt1", "Cdh5"), idents=c("fibro", "Linning", "osteoblast", "chondrocyte", "pericyte", "vascular")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```
```{r}
#re-run

#rna comparison
heatmapgns=c("Sox6", "Cd14", "Chad", "Chadl", "Cilp", "Clu", "Sox9", "Matn3", "Alpl", "Bglap", "Bglap2", "Omd", "Runx2", "Sp7", "Ostn", "Cdh5", "Emcn", "Pecam1", "Myh11", "Cd248", "Cd34", "Cdh11", "Col1a1", "Thy1", "Pdgfra", "Pdpn", "Fap", "Cd55", "Clic5", "Col22a1", "Hbegf", "Htra4", "Tspan15")

DotPlot(aggr, features = heatmapgns) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma")  +RotatedAxis()+coord_flip()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))

DotPlot(aggr, features = heatmapgns) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


```
```{r}
#re-run
fibro_linning_peri <- aggr[,grepl("fibro|linning|pericyte", aggr$cm_clusters, ignore.case=TRUE)]

fibro_linning_peri <- RunTFIDF(fibro_linning_peri)
fibro_linning_peri <- FindTopFeatures(fibro_linning_peri, min.cutoff = 20)
fibro_linning_peri <- RunSVD(fibro_linning_peri)

```

```{r}
fibro_linning_peri <- FindNeighbors(object = fibro_linning_peri, reduction = 'lsi', dims = 2:50)
fibro_linning_peri <- RunUMAP(fibro_linning_peri, dims = 2:50, reduction = 'lsi')
fibro_linning_peri <- RunTSNE(fibro_linning_peri, dims = 2:50, reduction = 'lsi')

p9 <- DimPlot(fibro_linning_peri, group.by = "orig.ident")
p10 <- DimPlot(fibro_linning_peri, group.by = "orig.ident", reduction = "tsne")
plot_grid(p9, p10)

```


```{r}
library(harmony)

fibro_linning_peri_safe <- fibro_linning_peri

fibro_linning_peri <- RunHarmony(fibro_linning_peri, group.by.vars = "orig.ident", 
                         reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE)

fibro_linning_peri <- FindNeighbors(fibro_linning_peri, dims = 2:50, reduction = "harmony", verbose = FALSE)
                                         
fibro_linning_peri <- RunUMAP(fibro_linning_peri, dims = 2:50, reduction = "harmony")
fibro_linning_peri <- RunTSNE(fibro_linning_peri, dims = 2:50, reduction = "harmony")


p11 <- DimPlot(fibro_linning_peri, group.by = "orig.ident")
p12 <- DimPlot(fibro_linning_peri, group.by = "orig.ident", reduction = "tsne")
plot_grid(p11, p12)

fibro_linning_peri$ATAC_snn_res.0.05 <- NULL
fibro_linning_peri$ATAC_snn_res.0.1 <- NULL
fibro_linning_peri$ATAC_snn_res.0.2 <- NULL
fibro_linning_peri$ATAC_snn_res.0.3 <- NULL
fibro_linning_peri$ATAC_snn_res.0.4 <- NULL
fibro_linning_peri$ATAC_snn_res.0.5 <- NULL
fibro_linning_peri$ATAC_snn_res.0.6 <- NULL
fibro_linning_peri$ATAC_nn_res.0.1 <- NULL
fibro_linning_peri$ATAC_nn_res.0.2 <- NULL
fibro_linning_peri$ATAC_nn_res.0.3 <- NULL
fibro_linning_peri$ATAC_nn_res.0.4 <- NULL

#findclusters
aggr@graphs

for (res in c(0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2)) {
   fibro_linning_peri <- FindClusters(fibro_linning_peri, resolution = res, algorithm = 3, graph.name = "ATAC_snn")
                          
}

clustree(fibro_linning_peri, prefix ="ATAC_snn_res.")

fibro_linning_peri <- FindClusters(fibro_linning_peri, resolution = 0.06, algorithm = 3, graph.name = "ATAC_snn")

```
