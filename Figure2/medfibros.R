


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
counts <- readRDS("/rds/projects/c/croftap-stia-atac/CM_multiome/data_from_Ilya/med_fibroblasts/exprs_raw.rds")
meta_data <- readRDS("/rds/projects/c/croftap-stia-atac/CM_multiome/data_from_Ilya/med_fibroblasts/meta_data.rds")


row.names(meta_data) <- meta_data$CellID
meta_data <- subset(meta_data, select = -CellID)

med_fibros<-CreateSeuratObject(counts = counts, meta.data = meta_data)
med_fibros<-NormalizeData(med_fibros)
med_fibros <- FindVariableFeatures(med_fibros, selection.method = "vst", nfeatures = 2000)
med_fibros <- ScaleData(med_fibros, verbose = FALSE)
med_fibros <- RunPCA(med_fibros, npcs = 30, verbose = FALSE)
med_fibros <- RunUMAP(med_fibros, reduction = "pca", dims = 1:30)
DimPlot(med_fibros, reduction = "umap", group.by = "Tissue")

```
```{r}



Idents(med_fibros)<-'Case'
VlnPlot(med_fibros, features=c("CBFB"))
levels(med_fibros)<-c("GutControl", "GutNonisnflamed", "GutInflamed",    "PSS", "SICCA",    "LungControl" ,   "LungEarlyStageILD" ,  "LungEndStageILD", "Osteoarthritis",      "RheumatoidArthritis")

DotPlot(med_fibros, features = c("RUNX1"))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Idents(med_fibros)<-"Cluster_name"
med_fibros$Cluster_name_case <- paste(Idents(med_fibros), med_fibros$Case, sep = "_")
Idents(med_fibros)<-"Cluster_name_case"


Idents(med_fibros)<-'Tissue'


dotplot_.data_avgacc_markers<-DotPlot(med_fibros, features = c("RUNX1"))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


dotplot_.data_avgacc_markers <- dotplot_.data_avgacc_markers[["data"]]
ggplot(dotplot_.data_avgacc_markers, aes(x=id, y=avg.exp.scaled)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_classic()

VlnPlot(med_fibros, features = c("RUNX1", "RUNX2", "RUNX3", "CBFB"), stack=T)+geom_jitter(size=0.05)


#to do
SCTransform data and check expression again
check DE of RUNX1 in Synovium/lung/gut
```
```{r}

Idents(med_fibros)<-"Cluster_name"
med_fibros$Cluster_name_tissue <- paste(Idents(med_fibros), med_fibros$Tissue, sep = "_")
Idents(med_fibros)<-"Cluster_name_tissue"

Dotplot_heatmap<-DotPlot(med_fibros, features = "RUNX1")
Dotplot_heatmap_data <- Dotplot_heatmap[["data"]]


Dotplot_heatmap_data<-cSplit(Dotplot_heatmap_data, 'id', sep="_", type.convert=FALSE)
Dotplot_heatmap_data <- subset(Dotplot_heatmap_data, select = -c (avg.exp, pct.exp)) 

Dotplot_heatmap_data_Salv = subset(Dotplot_heatmap_data, id_2 == "SalivaryGland")
Dotplot_heatmap_data_Salv<-Dotplot_heatmap_data_Salv[order(Dotplot_heatmap_data_Salv$id_1), ]
Dotplot_heatmap_data_Salv <- Dotplot_heatmap_data_Salv %>% remove_rownames %>% column_to_rownames(var="id_1")
Dotplot_heatmap_data_Salv = subset(Dotplot_heatmap_data_Salv, select = -c(1,3) )
names(Dotplot_heatmap_data_Salv)[names(Dotplot_heatmap_data_Salv) == 'avg.exp.scaled'] <- 'SalivaryGland'

Dotplot_heatmap_data_Gut = subset(Dotplot_heatmap_data, id_2 == "Gut")
Dotplot_heatmap_data_Gut<-Dotplot_heatmap_data_Gut[order(Dotplot_heatmap_data_Gut$id_1), ]
Dotplot_heatmap_data_Gut <- Dotplot_heatmap_data_Gut %>% remove_rownames %>% column_to_rownames(var="id_1")
Dotplot_heatmap_data_Gut = subset(Dotplot_heatmap_data_Gut, select = -c(1,3) )
names(Dotplot_heatmap_data_Gut)[names(Dotplot_heatmap_data_Gut) == 'avg.exp.scaled'] <- 'Gut'

Dotplot_heatmap_data_Lung = subset(Dotplot_heatmap_data, id_2 == "Lung")
Dotplot_heatmap_data_Lung<-Dotplot_heatmap_data_Lung[order(Dotplot_heatmap_data_Lung$id_1), ]
Dotplot_heatmap_data_Lung <- Dotplot_heatmap_data_Lung %>% remove_rownames %>% column_to_rownames(var="id_1")
Dotplot_heatmap_data_Lung = subset(Dotplot_heatmap_data_Lung, select = -c(1,3) )
names(Dotplot_heatmap_data_Lung)[names(Dotplot_heatmap_data_Lung) == 'avg.exp.scaled'] <- 'Lung'


Dotplot_heatmap_data_Synovium = subset(Dotplot_heatmap_data, id_2 == "Synovium")
Dotplot_heatmap_data_Synovium<-Dotplot_heatmap_data_Synovium[order(Dotplot_heatmap_data_Synovium$id_1), ]
Dotplot_heatmap_data_Synovium <- Dotplot_heatmap_data_Synovium %>% remove_rownames %>% column_to_rownames(var="id_1")
Dotplot_heatmap_data_Synovium = subset(Dotplot_heatmap_data_Synovium, select = -c(1,3) )
names(Dotplot_heatmap_data_Synovium)[names(Dotplot_heatmap_data_Synovium) == 'avg.exp.scaled'] <- 'Synovium'

Dotplot_heatmap_data_Salv$Gut<-Dotplot_heatmap_data_Gut$Gut
Dotplot_heatmap_data_Salv$Lung<-Dotplot_heatmap_data_Lung$Lung
Dotplot_heatmap_data_Salv$Synovium<-Dotplot_heatmap_data_Synovium$Synovium

Dotplot_heatmap_data_Salv<-as.matrix(Dotplot_heatmap_data_Salv)

library(heatmaply)
heatmaply(Dotplot_heatmap_data_Salv, row_dend_left = FALSE, show_dendrogram = c(F, F), Rowv=F, Colv=F)

library(ComplexHeatmap)
Heatmap(Dotplot_heatmap_data_Salv[,-3], border=T)


library(pheatmap)

cat_df2=data.frame("cluster"=c("SalivaryGland", "Gut", "Lung", "Synovium"))
rownames(cat_df2)=colnames(Dotplot_heatmap_data_Salv)


col=c("green","black","red")

color=colorRampPalette(col)(50)

pheatmap(mat,scale="row",show_rownames = F,color=color,cellwidth = cellwidth,cellheight = cellheight,border_color = "black")

pheatmap(Dotplot_heatmap_data_Salv, main="Title", cluster_rows=F, cluster_cols = F, annotation_col=cat_df2, cellwidth = 40, show_colnames = F, color = inferno(100), border_color = "black")



pheatmap(Dotplot_heatmap_data_Salv)

Dotplot_heatmap_data_Salvdf<-as.data.frame(Dotplot_heatmap_data_Salv)
Dotplot_heatmap_data_Salvdf$Lung<-NULL
pheatmap(Dotplot_heatmap_data_Salvdf)

Heatmap(Dotplot_heatmap_data_Salv)


```



```{r}
Idents(med_fibros)<-"Case"
med_fibros$Case_tissue <- paste(Idents(med_fibros), med_fibros$Tissue, sep = "_")
Idents(med_fibros)<-"Case_tissue"

Dotplot_heatmap<-DotPlot(med_fibros, features = "FAP")
Dotplot_heatmap_data <- Dotplot_heatmap[["data"]]

library(splitstackshape)
library(dplyr)
library(tidyverse)
Dotplot_heatmap_data<-cSplit(Dotplot_heatmap_data, 'id', sep="_", type.convert=FALSE)
Dotplot_heatmap_data <- subset(Dotplot_heatmap_data, select = -c (avg.exp, pct.exp)) 

Dotplot_heatmap_data_Salv = subset(Dotplot_heatmap_data, id_2 == "SalivaryGland")
Dotplot_heatmap_data_Salv<-Dotplot_heatmap_data_Salv[order(Dotplot_heatmap_data_Salv$id_1), ]
Dotplot_heatmap_data_Salv <- Dotplot_heatmap_data_Salv %>% remove_rownames %>% column_to_rownames(var="id_1")
Dotplot_heatmap_data_Salv = subset(Dotplot_heatmap_data_Salv, select = -c(1,3) )
names(Dotplot_heatmap_data_Salv)[names(Dotplot_heatmap_data_Salv) == 'avg.exp.scaled'] <- 'SalivaryGland'

Dotplot_heatmap_data_Gut = subset(Dotplot_heatmap_data, id_2 == "Gut")
Dotplot_heatmap_data_Gut<-Dotplot_heatmap_data_Gut[order(Dotplot_heatmap_data_Gut$id_1), ]
Dotplot_heatmap_data_Gut <- Dotplot_heatmap_data_Gut %>% remove_rownames %>% column_to_rownames(var="id_1")
Dotplot_heatmap_data_Gut = subset(Dotplot_heatmap_data_Gut, select = -c(1,3) )
names(Dotplot_heatmap_data_Gut)[names(Dotplot_heatmap_data_Gut) == 'avg.exp.scaled'] <- 'Gut'

Dotplot_heatmap_data_Lung = subset(Dotplot_heatmap_data, id_2 == "Lung")
Dotplot_heatmap_data_Lung<-Dotplot_heatmap_data_Lung[order(Dotplot_heatmap_data_Lung$id_1), ]
Dotplot_heatmap_data_Lung <- Dotplot_heatmap_data_Lung %>% remove_rownames %>% column_to_rownames(var="id_1")
Dotplot_heatmap_data_Lung = subset(Dotplot_heatmap_data_Lung, select = -c(1,3) )
names(Dotplot_heatmap_data_Lung)[names(Dotplot_heatmap_data_Lung) == 'avg.exp.scaled'] <- 'Lung'


Dotplot_heatmap_data_Synovium = subset(Dotplot_heatmap_data, id_2 == "Synovium")
Dotplot_heatmap_data_Synovium<-Dotplot_heatmap_data_Synovium[order(Dotplot_heatmap_data_Synovium$id_1), ]
Dotplot_heatmap_data_Synovium <- Dotplot_heatmap_data_Synovium %>% remove_rownames %>% column_to_rownames(var="id_1")
Dotplot_heatmap_data_Synovium = subset(Dotplot_heatmap_data_Synovium, select = -c(1,3) )
names(Dotplot_heatmap_data_Synovium)[names(Dotplot_heatmap_data_Synovium) == 'avg.exp.scaled'] <- 'Synovium'


rownames(Dotplot_heatmap_data_Salv)<-c("Control", "Inflamation")

Dotplot_heatmap_data_Gut$all <- '1'
Dotplot_heatmap_data_Lung$all <- '1'

Dotplot_heatmap_data_Gut <- Dotplot_heatmap_data_Gut[c(1,2),]
Dotplot_heatmap_data_Lung <- Dotplot_heatmap_data_Lung[c(1,3),]

Dotplot_heatmap_data_Salv$Gut<-Dotplot_heatmap_data_Gut$Gut
Dotplot_heatmap_data_Salv$Lung<-Dotplot_heatmap_data_Lung$Lung
Dotplot_heatmap_data_Salv$Synovium<-Dotplot_heatmap_data_Synovium$Synovium

Dotplot_heatmap_data_Salv<-as.matrix(Dotplot_heatmap_data_Salv)

library(heatmaply)
heatmaply(Dotplot_heatmap_data_Salv, row_dend_left = FALSE, show_dendrogram = c(F, F), Rowv=F, Colv=F)

library(ComplexHeatmap)

Heatmap(Dotplot_heatmap_data_Salv, border=T)

```


```{r}
Idents(med_fibros)<-'Case_tissue'

med_fibros$Case_tissue_sample <- paste(med_fibros$Case_tissue, med_fibros$SampleID, sep=".")

cts_fibs<-AggregateExpression(med_fibros, group.by = c("Case_tissue_sample"), assays = "RNA", slot = "counts", return.seurat = F)

cts_fibs<-cts_fibs$RNA
cts_fibs<-as.data.frame(cts_fibs)
meta_data=colnames(cts_fibs)
meta_data<-as.data.frame(meta_data)
library(splitstackshape)
meta_data$to_split<-meta_data$meta_data
meta_data<-cSplit(meta_data, splitCols = "to_split", sep=".")
meta_data<-cSplit(meta_data, splitCols = "to_split_1", sep="_")
colnames(meta_data)<-c("all","sample", "tissue_cond", "tissue")
meta_data$all<-as.factor(meta_data$all)
meta_data$donor<-as.factor(meta_data$sample)
meta_data$treatment<-as.factor(meta_data$tissue_cond)
meta_data$treatment<-as.factor(meta_data$tissue)


dds <- DESeqDataSetFromMatrix(countData = cts_fibs,
                                  colData = meta_data,
                                  design = ~1)
    
dds <- scran::computeSumFactors(dds)
print(dds)
print(quantile(rowSums(counts(dds))))

#mingenecount <- quantile(rowSums(counts(dds)), 0.5)
mingenecount <- 200
maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount & rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

dds@colData[['treatment']] <- as.factor(dds@colData[['treatment']])

design(dds) <- formula(~ treatment)
print(design(dds))
dds <- DESeq(dds, test = "Wald")
```


```{r}
ggplots <- list()
for (i in 1:length(unique(colData(dds)$tissue))){


dds_gut <- dds[, colData(dds)$tissue == unique(colData(dds)$tissue)[[i]]]


dds_norm <- counts(dds_gut, normalized=TRUE)
dds_norm_df <- as.data.frame(dds_norm)
dds_norm_df$Gene <- rownames(dds_norm_df)
dds_long <- melt(dds_norm_df, id.vars = "Gene", variable.name = "Sample", value.name = "NormalizedCount")

library(splitstackshape)
dds_long <- cSplit(dds_long, splitCols="Sample",sep = ".")

gene = "FAP"
ggplots[[i]] <-dds_long %>% filter(Gene == gene) %>% 
 ggplot(aes(x = Sample_1, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_1), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip()+  # Clean theme
  theme(
     # Remove y-axis ticks
    legend.position = "none",axis.title.y = element_blank(),  # Remove legend
    plot.title = element_text(hjust = 0.5)  # Center title
  ) 

}


plot_grid(ggplots[[1]], ggplots[[2]], ggplots[[3]], ggplots[[4]])



```



```{r}



Idents(med_fibros)<-'Cluster_name_case'
C4_Lung_end_vs_ctrl<-FindMarkers(med_fibros, ident.1 = "SPARC+COL3A1+ C4_LungEndStageILD", ident.2 = "SPARC+COL3A1+ C4_LungControl")
C4Lung_early_vs_ctrl<-FindMarkers(med_fibros, ident.1 = "SPARC+COL3A1+ C4_LungEarlyStageILD", ident.2 = "SPARC+COL3A1+ C4_LungControl")
C4_Lung_end_vs_ctrl$gene<-rownames(C4_Lung_end_vs_ctrl)
C4Lung_early_vs_ctrl$gene<-rownames(C4Lung_early_vs_ctrl)

C8_Lung_end_vs_ctrl<-FindMarkers(med_fibros, ident.1 = "PTGS2+SEMA4A+ C8_LungEndStageILD", ident.2 = "PTGS2+SEMA4A+ C8_LungControl")
C8Lung_early_vs_ctrl<-FindMarkers(med_fibros, ident.1 = "PTGS2+SEMA4A+ C8_LungEarlyStageILD", ident.2 = "PTGS2+SEMA4A+ C8_LungControl")
C8_Lung_end_vs_ctrl$gene<-rownames(C8_Lung_end_vs_ctrl)
C8Lung_early_vs_ctrl$gene<-rownames(C8Lung_early_vs_ctrl)



C4_Fibro_RA_vs_OA<-FindMarkers(med_fibros, ident.1 = "SPARC+COL3A1+ C4_RheumatoidArthritis", ident.2 = "SPARC+COL3A1+ C4_Osteoarthritis")
C4_Fibro_RA_vs_OA$gene<-rownames(C4_Fibro_RA_vs_OA)

C8_Fibro_RA_vs_OA<-FindMarkers(med_fibros, ident.1 = "PTGS2+SEMA4A+ C8_RheumatoidArthritis", ident.2 = "PTGS2+SEMA4A+ C8_Osteoarthritis")
C8_Fibro_RA_vs_OA$gene<-rownames(C8_Fibro_RA_vs_OA)


C8_SG_PSS_vs_SICCA<-FindMarkers(med_fibros, ident.1 = "PTGS2+SEMA4A+ C8_PSS", ident.2 = "PTGS2+SEMA4A+ C8_SICCA")
C4_SG_PSS_vs_SICCA<-FindMarkers(med_fibros, ident.1 = "SPARC+COL3A1+ C4_PSS", ident.2 = "SPARC+COL3A1+ C4_SICCA")
C8_SG_PSS_vs_SICCA$gene<-rownames(C8_SG_PSS_vs_SICCA)
C4_SG_PSS_vs_SICCA$gene<-rownames(C4_SG_PSS_vs_SICCA)



```
```{r}
Idents(med_fibros)<-"Tissue"
med_fibros$tissue_case <- paste(Idents(med_fibros), med_fibros$Case, sep = "_")
Idents(med_fibros)<-"tissue_case"

Synovium_RAvsOA<-FindMarkers(med_fibros, ident.1 = "Synovium_RheumatoidArthritis", ident.2 = "Synovium_Osteoarthritis")
Synovium_RAvsOA$gene<-rownames(Synovium_RAvsOA)


Lung_EarlyvsCtrl<-FindMarkers(med_fibros, ident.1 = "Lung_LungEarlyStageILD", ident.2 = "Lung_LungControl")
Lung_EarlyvsCtrl$gene<-rownames(Lung_EarlyvsCtrl)

Lung_EndvsCtrl<-FindMarkers(med_fibros, ident.1 = "Lung_LungEndStageILD", ident.2 = "Lung_LungControl")
Lung_EndvsCtrl$gene<-rownames(Lung_EndvsCtrl)





```

```{r}
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
```


```{r}
#pseudo bulk donors based on this method
#https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

med_fibros_synoviumonly <- med_fibros[,grepl("Osteoarthritis|RheumatoidArthritis", med_fibros@meta.data[["Case"]], ignore.case=TRUE)]


Idents(med_fibros_synoviumonly)<-'DonorID'
levels(med_fibros_synoviumonly)

med_fibros_synoviumonly$DonorID_orig<-med_fibros_synoviumonly$DonorID

current.sample.ids <- c("OA180104_A",     "OA180108_B"  ,   "BWH064_CD45N",   "RA174_CD45N" ,  
 "OA213_CD45N" ,   "RA178_CD45N"   , "OA214L_CD45N" ,  "RA195_CD45N",   
 "300_0488_45N" ,  "BX233_CD45N"  ,  "OA180430_CD45N", "OA180507_CD45N",
"RA_186"     ,    "RA_190"   ,      "BX250_CD45N",    "BX254_CD45N",   
"BX256_CD45N"  ,  "BX230"    ,      "BWH075"  ,       "BWH076CD45n",   
 "BWH078")


new.sample.ids <- c("s1",     "s2"  ,   "s3",   "s4" ,  
 "s5" ,   "s6"   , "s7" ,  "s8",   
 "s9" ,  "s10"  ,  "s11", "s12",
"s13"     ,    "s14"   ,      "s15",    "s16",   
"s17"  ,  "s18"    ,      "s19"  ,       "s20",   
 "s21")


med_fibros_synoviumonly@meta.data[["DonorID"]] <- plyr::mapvalues(x = med_fibros_synoviumonly@meta.data[["DonorID"]], from = current.sample.ids, to = new.sample.ids)



seurat=med_fibros_synoviumonly
Idents(seurat)<-'Cluster_name'
counts <- seurat@assays$RNA@counts 
metadata <- seurat@meta.data
metadata$DonorID<-gsub("_","",as.character(metadata$DonorID))
metadata$Cluster_name <- factor(seurat@active.ident)
metadata$DonorID <- factor(seurat$DonorID)
sce <- SingleCellExperiment(assays = list(counts = counts), 
                           colData = metadata)
groups <- colData(sce)[, c("DonorID", "Cluster_name")]
kids <- purrr::set_names(levels(sce$Cluster_name))
nk <- length(kids)
sids <- purrr::set_names(levels(sce$DonorID))
ns <- length(sids)
n_cells <- as.numeric(table(sce$DonorID))
m <- match(sids, sce$DonorID)

ei <- data.frame(colData(sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
                select(-"Cluster_name")


dim(sce)
#sce$is_outlier <- isOutlier(
#        metric = sce$total_features_by_counts,
#        nmads = 2, type = "both", log = TRUE)
#sce <- sce[, !sce$is_outlier]
#dim(sce)
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)
```

```{r}
groups <- colData(sce)[, c("Cluster_name", "DonorID")]
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)



pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
        lapply(function(u) 
                set_colnames(t(u), 
                             stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))


get_sample_ids <- function(x){
        pb[[x]] %>%
                colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
        unlist()



samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
        rep(names(pb)[x], 
            each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
        unlist()



gg_df <- data.frame(Cluster_name = de_cluster_ids,
                    DonorID = de_samples)


#library(splitstackshape)
#ei<-cSplit(ei, 'DonorID', sep="_", type.convert=FALSE)

gg_df <- left_join(gg_df, ei[, c("DonorID", "Case")]) 


metadata <- gg_df %>%
        dplyr::select(Cluster_name, DonorID, Case) 
        
metadata        
```
```{r}
Idents(med_fibros)<-'Cluster_name'
clusters <- levels(med_fibros)
cluster_metadata <- metadata[which(metadata$Cluster_name == clusters[3]), ]
#cluster_metadata$ID <- seq_along(cluster_metadata[,1])
#cluster_metadata$DonorID_2 <- paste(cluster_metadata$DonorID_1, cluster_metadata$ID, sep="-")
rownames(cluster_metadata) <- cluster_metadata$DonorID
counts <- pb[[clusters[5]]]
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% cluster_metadata$DonorID)])
all(rownames(cluster_metadata) == colnames(cluster_counts))       

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ Case)


rld <- rlog(dds, blind=TRUE)
DESeq2::plotPCA(rld, intgroup = "DonorID")


rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor, annotation = cluster_metadata[, c("Case"), drop=F])
```
```{r}
dds <- DESeq(dds)
plotDispEsts(dds)
```

```{r}
Idents(med_fibros_synoviumonly)<-'Case'
levels(med_fibros_synoviumonly)

contrast <- c("Case", levels(med_fibros_synoviumonly)[2], levels(med_fibros_synoviumonly)[1])


res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

#did not work!
#res <- lfcShrink(dds, 
#                 contrast =  contrast,
#                 res=res)


res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()


p_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, pvalue < p_cutoff) %>%
        dplyr::arrange(pvalue)

normalized_counts <- counts(dds, 
                            normalized = TRUE)

```

```{r}
top20_sig_genes <- sig_res %>%
        dplyr::arrange(padj) %>%
        dplyr::pull(gene) %>%
        head(n=20)


top20_sig_norm <- data.frame(normalized_counts) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% top20_sig_genes)

gathered_top20_sig <- top20_sig_norm %>%
        gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
        
gathered_top20_sig <- inner_join(ei[, c("DonorID", "Case" )], gathered_top20_sig, by = c("DonorID" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig) +
        geom_point(aes(x = gene, 
                       y = normalized_counts, 
                       color = DonorID), 
                   position=position_jitter(w=0.1,h=0)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle("Top 20 Significant DE Genes") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        theme(plot.title = element_text(hjust = 0.5))
```



```{r}
sig_norm <- data.frame(normalized_counts) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% sig_res$gene)
        
# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
    color = heat_colors, 
    cluster_rows = T, 
    show_rownames = F,
    annotation = cluster_metadata[, c("DonorID", "Case")], 
    border_color = NA, 
    fontsize = 10, 
    scale = "row", 
    fontsize_row = 10, 
    height = 20)        
```
```{r}

```



```{r}
se <- aggregateAcrossCells(sce, ids = colData(sce)[, c("DonorID", "Cluster_name")])
y.all <- DGEList(counts(se), samples = colData(se))

dea <- pseudoBulkDGE(se, label = se$Cluster_name, condition = se$DonorID)
#try this! https://support.bioconductor.org/p/132898/
```



```{r}


#med_fibros <- AddModuleScore(med_fibros, features = "RUNX1", name = "RUNX1mod")

#med_fibros_meta <- med_fibros@meta.data %>% select(RUNX1mod1, )




#create meta data splot
med_fibros$sample_tissue_cluster<-paste(med_fibros$SampleID, med_fibros$Tissue, med_fibros$Cluster_name, sep=".")


paste(med_fibros$InflamScore, med_fibros$LibraryID, sep="_") %>% table()

#extract avg. scaled expression for your gene in each cluster and sample
Idents(med_fibros)<-"sample_tissue_cluster"
dotplot<-DotPlot(med_fibros, features = "RUNX1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"


dotplot_data <- dotplot_data %>% cSplit(splitCols = "id", sep=".")

inflam_df <- paste(med_fibros$SampleID, med_fibros$InflamScore, sep=".") %>% unique() %>% as.data.frame() %>% cSplit(splitCols = ".", sep=".")

colnames(inflam_df) <- c("id_1", "inflam", "extra")

inflam_df <- inflam_df %>% replace(is.na(.), 0)

inflam_df$inflam <- paste(inflam_df$inflam, inflam_df$extra, sep=".")
inflam_df$inflam <- as.double(inflam_df$inflam)

final_df <- dotplot_data %>% 
  left_join(inflam_df, by="id_1")


#plot
#dotplot_data_synovium <- dotplot_data[dotplot_data$condition_2 == "Synovium",]
ggplot(final_df, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") +
        facet_wrap(~id_2)


#ml = lm(Runx1~InflamScore, data = dotplot_data)
#summary(ml)$r.squared

library(ggpubr)
ggscatter(final_df, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)+
        facet_wrap(~id_2)



#filter if you want a specific cluster
dotplot_data_C4 <- final_df[final_df$id_3 == "SPARC+COL3A1+ C4",]
dotplot_data_C4 <- dotplot_data_C4[dotplot_data_C4$id_2 == "Synovium",]

ggplot(dotplot_data_C4, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") +
        facet_wrap(~id_2)




ml = lm(Runx1~InflamScore, data = dotplot_data_C4)
summary(ml)$r.squared

library(ggpubr)

dotplot_data_syn <- final_df[final_df$id_2 == "Synovium",]

ggscatter(dotplot_data_C4, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 1)+
        facet_wrap(~id_2)


ggplot(data=dotplot_data_C4, aes(Runx1,inflam)) + geom_point(alpha=0.6, color="grey", size=0.1) + ggtitle("C4 cluster") +theme_minimal() +
    geom_point(data = dotplot_data_C4, color = "darkred",size=1.5)+ theme(axis.title.x = element_blank())+
    stat_smooth(method = "lm",
        col = "black",
        se = T,
        size = 0.5)+theme_ArchR()

```

```{r}


med_fibros$sample_tissue<-paste(med_fibros$SampleID, med_fibros$Tissue, sep=".")


#extract avg. scaled expression for your gene in each cluster and sample
Idents(med_fibros)<-"sample_tissue"
dotplot<-DotPlot(med_fibros, features = "RUNX1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"


dotplot_data <- dotplot_data %>% cSplit(splitCols = "id", sep=".")

inflam_df <- paste(med_fibros$SampleID, med_fibros$InflamScore, sep=".") %>% unique() %>% as.data.frame() %>% cSplit(splitCols = ".", sep=".")

colnames(inflam_df) <- c("id_1", "inflam", "extra")

inflam_df <- inflam_df %>% replace(is.na(.), 0)

inflam_df$inflam <- paste(inflam_df$inflam, inflam_df$extra, sep=".")
inflam_df$inflam <- as.double(inflam_df$inflam)

final_df <- dotplot_data %>% 
  left_join(inflam_df, by="id_1")


#plot
#dotplot_data_synovium <- dotplot_data[dotplot_data$condition_2 == "Synovium",]
ggplot(final_df, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") +
        facet_wrap(~id_2)


ml = lm(Runx1~InflamScore, data = dotplot_data)
summary(ml)$r.squared

library(ggpubr)

final_df %>% filter(id_2 == "Synovium") %>% 
ggscatter( x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)+
        facet_wrap(~id_2)








```






```{r}


cols <- ArchR::paletteDiscrete(med_fibros_synoviumonly@meta.data[,"Cluster_name"]) %>% as.data.frame()

med_fibros_synoviumonly$Cluster_name %>% table() %>% as.data.frame() %>% 
ggplot(aes(y=Freq, x=., fill=.))  + 
        geom_bar(stat = 'identity', position = position_dodge(), aes(group = .),width = 1) + 
                theme(
            axis.text.x = element_text(angle = 45, hjust=1),
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
            # strip.text = element_blank()
        ) + 
        guides(color = 'none', fill = 'none') + 
        labs(y = '# cells')+
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = c(0, 0)))+theme_ArchR()+scale_fill_manual(values =cols$.)+
  theme(strip.background = element_rect(fill="white", size=1, color="white"))+RotatedAxis()
```




```{r}


med_fibros$tissue_case %>% unique()
Idents(med_fibros) <- 'tissue_case'
levels(med_fibros)[c(1,2)]

DotPlot(med_fibros, idents=levels(med_fibros)[c(1,2)], features="RUNX1")

```



```{r}

library(readr)
grn <- read_csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/grn.csv")
grn_RUNX1 <- grn %>% filter(Source == "RUNX1")

med_fibros <- AddModuleScore(med_fibros, features=list(grn_RUNX1$Target), name="RUNX1_grn")

Idents(med_fibros) <- 'Tissue'
med_fibros_no_lung <- subset(med_fibros, idents= levels(med_fibros)[-3])
med_fibros_no_lung <- med_fibros_no_lung %>% ScaleData()

Idents(med_fibros_no_lung)<-"Cluster_name"
Dotplot_heatmap<-DotPlot(med_fibros_no_lung, features = c("RUNX1", "RUNX1_grn1"))
Dotplot_heatmap_data <- Dotplot_heatmap[["data"]]

dotplot<-Dotplot_heatmap_data %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

Heatmap(dotplot, border=T)


```

```{r}

med_fibros_no_lung$SampleID %>% unique()
#extract avg. scaled expression for your gene in each cluster and sample
Idents(med_fibros_no_lung)<-"SampleID"
dotplot<-DotPlot(med_fibros_no_lung, features = "RUNX1_grn1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"


#dotplot_data <- dotplot_data %>% cSplit(splitCols = "id", sep=".")
library(splitstackshape)
inflam_df <- paste(med_fibros$SampleID, med_fibros$InflamScore, sep=".") %>% unique() %>% as.data.frame() %>% cSplit(splitCols = ".", sep=".")

colnames(inflam_df) <- c("id_1", "inflam", "extra")

inflam_df <- inflam_df %>% replace(is.na(.), 0)

inflam_df$inflam <- paste(inflam_df$inflam, inflam_df$extra, sep=".")
inflam_df$inflam <- as.double(inflam_df$inflam)

dotplot_data$id_1 <- dotplot_data$id

final_df <- dotplot_data %>% 
  left_join(inflam_df, by="id_1")


#plot
#dotplot_data_synovium <- dotplot_data[dotplot_data$condition_2 == "Synovium",]
ggplot(final_df, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_1))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") #+
        #facet_wrap(~id_2)


ml = lm(Runx1~inflam, data = final_df)
summary(ml)$r.squared

library(ggpubr)

ggscatter(final_df, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = -2, label.y = 1)



ggscatter(final_df, x = "Runx1", y = "inflam",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "black", fill = NA),
          shape = 21,             # Hollow circles
          color = "red",        # Border color
          fill = "white",         # Inside fill color (optional)
          size = 2 ,
           stroke = 1.6# Adjust size as needed
) +
  stat_cor(method = "pearson", label.x = 1, label.y = 1)


```

```{r}

Idents(med_fibros_no_lung)<-"SampleID"
dotplot<-DotPlot(med_fibros_no_lung, features = "RUNX1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"


#dotplot_data <- dotplot_data %>% cSplit(splitCols = "id", sep=".")

inflam_df <- paste(med_fibros$SampleID, med_fibros$InflamScore, sep=".") %>% unique() %>% as.data.frame() %>% cSplit(splitCols = ".", sep=".")

colnames(inflam_df) <- c("id_1", "inflam", "extra")

inflam_df <- inflam_df %>% replace(is.na(.), 0)

inflam_df$inflam <- paste(inflam_df$inflam, inflam_df$extra, sep=".")
inflam_df$inflam <- as.double(inflam_df$inflam)

dotplot_data$id_1 <- dotplot_data$id

final_df <- dotplot_data %>% 
  left_join(inflam_df, by="id_1")


#plot
#dotplot_data_synovium <- dotplot_data[dotplot_data$condition_2 == "Synovium",]
ggplot(final_df, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_1))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") #+
        #facet_wrap(~id_2)


ml = lm(Runx1~inflam, data = final_df)
summary(ml)$r.squared

library(ggpubr)

ggscatter(final_df, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = -1, label.y = 1.1)



ggscatter(final_df, x = "Runx1", y = "inflam",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "black", fill = NA),
          shape = 21,             # Hollow circles
          color = "red",        # Border color
          fill = "white",         # Inside fill color (optional)
          size = 2 ,
           stroke = 1.6# Adjust size as needed
) +
  stat_cor(method = "pearson", label.x = 1, label.y = 1)


```


