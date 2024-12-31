
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
library(harmony)
library(clustree)
library(plyr)
library(dplyr)
library(EnhancedVolcano)
library(monocle)
library(SoupX)
library(DoubletFinder)
options(bitmapType='cairo')
```

```{r}
ra1_counts<-Read10X("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/SDY2213/ra1/")
ra1_RNA<-ra1_counts$`Gene Expression`
ra1_ATAC<-ra1_counts$Peaks

ncol(ra1_ATAC)
ncol(ra1_RNA)

head(colnames(ra1_ATAC))
head(colnames(ra1_RNA))

library(EnsDb.Hsapiens.v86)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))


ra1 <- CreateSeuratObject(
  counts = ra1_counts$`Gene Expression`,
  assay = "RNA"
)

ncol(ra1)

# create ATAC assay and add it to the object
atac <- CreateChromatinAssay(
  counts = ra1_counts$Peaks,
  sep = c(":", "-"),
  fragments = NULL,
  annotation = annotation
)

atac_s<-CreateSeuratObject(ra1_ATAC, assay="ATAC")

ra1[["ATAC"]]<-atac_s[["ATAC"]]

DefaultAssay(ra1)<-'RNA'

ra1 <- ra1 %>%
    SCTransform(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)


DefaultAssay(ra1)<-'ATAC'

obj.atac <- obj.atac %>%
    RunTFIDF() %>%
    FindTopFeatures() %>%
    RunSVD() %>%
    RunUMAP(reduction = 'lsi', dims = 2:30, verbose = FALSE)

```





```{r}
ra4_counts<-Read10X("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/SDY2213/ra4/")
ra4_counts$`Gene Expression`

ra4a_counts<-Read10X("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/SDY2213/ra4a/")
ra4a_counts$`Gene Expression`

ra4b_counts<-Read10X("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/SDY2213/ra4b/")
ra4b_counts$`Gene Expression`

ra4c_counts<-Read10X("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/SDY2213/ra4c/")
ra4c_counts$`Gene Expression`


library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)


data.10x = list()
dirs <- dir(".",pattern ="ra")

scrna.list = list()
samples<-sub('./', '', dirs)
dirs <- paste0( "./", dirs) 

for (i in 1:length(dirs)) {
    data.10x[[i]] <- Read10X(data.dir = dirs[[i]])
}

library(EnsDb.Hsapiens.v86)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))


for (i in 1:length(data.10x)) {
    scrna.list[[i]]<-data.10x[[i]]$`Gene Expression`;
    scrna.list[[i]] = CreateSeuratObject(counts = scrna.list[[i]], min.cells=3, min.features=200, project=samples[i]);
    scrna.list[[i]][["ATAC"]] <- CreateChromatinAssay(counts = data.10x[[i]]$Peaks, sep = c(":", "-"),fragments = NULL, annotation = annotation);
    scrna.list[[i]][["percent.mt"]] = PercentageFeatureSet(object=scrna.list[[i]], pattern = "^MT-");
    DefaultAssay(scrna.list[[i]]) <- "ATAC";
    scrna.list[[i]] <- NucleosomeSignal(scrna.list[[i]]);
    scrna.list[[i]] <- TSSEnrichment(scrna.list[[i]]);
    scrna.list[[i]] <- subset(scrna.list[[i]], subset = nFeature_RNA > 500 & nFeature_RNA > 7000 & percent.mt < 10 &
    nCount_ATAC > 1000 & nCount_ATAC < 100000 & nucleosome_signal < 2 & TSS.enrichment > 1);
    scrna.list[[i]] <- scrna.list[[i]] %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 'q5') %>% RunSVD() %>% RunUMAP(reduction = 'lsi', dims = 2:30, verbose = FALSE);
    DefaultAssay(scrna.list[[i]]) <- "RNA";
    scrna.list[[i]] =NormalizeData(object = scrna.list[[i]]);
    scrna.list[[i]] =ScaleData(object = scrna.list[[i]]);
    scrna.list[[i]] =FindVariableFeatures(object = scrna.list[[i]]);
    scrna.list[[i]] =RunPCA(object = scrna.list[[i]], verbose = FALSE);
}
names(scrna.list) <- sub('./', '', dirs)



```
```{r}
library(SeuratDisk)
Convert("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/SDY2213/rna_no_meta.h5ad", dest = "h5seurat", overwrite = F)
rna <- LoadH5Seurat("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/SDY2213/rna_no_meta.h5seurat")

library(tidyverse)
meta_rna<- meta_rna %>% remove_rownames %>% column_to_rownames(var="...1")
rna<-AddMetaData(rna, meta_rna)
meta_rna_obsm<-as.data.frame(meta_rna_obsm)
rownames(meta_rna_obsm)<-rownames(meta_rna)
meta_rna_obsm$...1<-NULL
colnames(meta_rna_obsm)<-c("UMAP_1", "UMAP_2")
meta_rna_obsm_mat <- as(meta_rna_obsm, "matrix")

# Create DimReducObject and add to object
rna[['UMAP']] <- CreateDimReducObject(embeddings = meta_rna_obsm_mat, key = "UMAP_", global = T, assay = "RNA")

cols <- ArchR::paletteDiscrete(rna@meta.data[, "subtype_fine"])

DimPlot(rna, group.by = "subtype_fine", cols=cols)
Idents(rna)<-'subtype_fine'
DotPlot(rna, features = "RUNX1")

Idents(rna)<-'subtype_fine'
DotPlot(rna, features = c("RUNX1", "MMP14", "IGF1", "CTHRC1"))


Idents(rna)<-'cluster_harmony_k100'
DimPlot(rna)+NoAxes()
VlnPlot(rna, features = c("RUNX1", "RUNX2", "RUNX3", "CBFB"), stack = T)

Idents(rna)<-'subtype_fine'
DimPlot(rna)+NoAxes()
levels(rna)<-c("Resting Sublining", "Inflamed Sublining"    ,               "Resting Lining"         ,              "Inflamed Lining", "Intermediate Lining/Sublining"  ,      "Intermediate Inflamed/Resting Lining"  )
DotPlot(rna, features = c("RUNX1"))+RotatedAxis()+coord_flip()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


```


