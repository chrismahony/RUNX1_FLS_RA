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

VlnPlot(merging, features = c("RUNX1", "RUNX2", "RUNX3", "CBFB"), stack = T)+geom_jitter(size=0.01)



```
```{r}

merging$cluster_condition <- paste(Idents(merging), merging$condition, sep = "_")
Idents(merging)<-'cluster_condition'

COMP_RAvsOA<-FindMarkers(merging, ident.1 = "COMP+ Fibroblast Niche_RA", ident.2 ="COMP+ Fibroblast Niche_OA" )
lining_RAvsOA<-FindMarkers(merging, ident.1 = "Lining Layer Cells_RA", ident.2 ="Lining Layer Cells_OA" )
APOD_RAvsOA<-FindMarkers(merging, ident.1 = "APOD+ GAS5+ FABP4+_RA", ident.2 ="APOD+ GAS5+ FABP4+_OA" )

lining_RAvsOA$gene<-rownames(lining_RAvsOA)
APOD_RAvsOA$gene<-rownames(APOD_RAvsOA)

```
```{r}
Idents(merging)<-'number'
levels(merging)
current.sample.ids <- c("T Cell Rich Niche",      "B Cell Rich Niche"  ,    "Vascular Niche"    ,     "Erythrocytes"   ,        "Lining Layer Cells",  "COMP+ Fibroblast Niche", "APOD+ GAS5+ FABP4+")
new.sample.ids <- c("SL",      "SL"  ,    "EC"    ,     "SL"   ,        "LL",  "SL", "SL")

merging$pseudo_bulk<-merging$number
merging@meta.data[["pseudo_bulk"]] <- plyr::mapvalues(x = merging@meta.data[["pseudo_bulk"]], from = current.sample.ids, to = new.sample.ids)
table(merging$pseudo_bulk)

DimPlot(merging, group.by = "pseudo_bulk")
```
```{r}
Idents(merging)<-'pseudo_bulk'
merging$pseudo_bulk_condition <- paste(Idents(merging), merging$condition, sep = "_")
Idents(merging)<-'pseudo_bulk_condition'
levels(merging)

pseudo_SL_RA_VS_OA<-FindMarkers(merging, ident.1 = "SL_RA", ident.2 = "SL_OA")
VlnPlot(merging, features="RUNX1", idents=c("SL_RA", "SL_OA"))

```
```{r}
merging_largebits_RA  <- merging_largebits[,grepl("RA", merging_largebits$condition, ignore.case=TRUE)]
merging_largebits_OA  <- merging_largebits[,grepl("OA", merging_largebits$condition, ignore.case=TRUE)]


all_coords_dist_norm_largebits_RA<-rbind(coords_slice1.1_norm, coords_slice1.17_norm,  coords_slice1.21_norm, coords_slice1.22_norm, coords_slice1.25_norm)

all_coords_dist_norm_largebits_OA<-rbind( coords_slice1.11_norm,  coords_slice1.13_norm, coords_slice1.14_norm, coords_slice1.15_norm, coords_slice1.16_norm)


merging_largebits_RA<-AddModuleScore(merging_largebits_RA, "RUNX1", name="RUNX1_module")
merging_largebits_OA<-AddModuleScore(merging_largebits_OA, "RUNX1", name="RUNX1_module")


RUNX1_RA<-merging_largebits_RA@meta.data
RUNX1_OA<-merging_largebits_OA@meta.data


RUNX1_RA <- subset(RUNX1_RA, select = c("RUNX1_module1"))
RUNX1_OA <- subset(RUNX1_OA, select = c("RUNX1_module1"))


RUNX1_RA_df <- data.frame(Dist=all_coords_dist_norm_largebits_RA$dist1, 
				Score=as.numeric(RUNX1_RA$RUNX1_module1) / as.numeric(max(RUNX1_RA$RUNX1_module1)), #scaling cell type signature to max here 
				CellType="RA")

RUNX1_OA_df <- data.frame(Dist=all_coords_dist_norm_largebits_OA$dist1, 
				Score=as.numeric(RUNX1_OA$RUNX1_module1) / as.numeric(max(RUNX1_OA$RUNX1_module1)), #scaling cell type signature to max here 
				CellType="OA")

RUNX1_largebits<-rbind(RUNX1_RA_df, RUNX1_OA_df)


ggplot(RUNX1_largebits, aes(Dist, Score, lty=CellType)) + geom_smooth(alpha=.1) + labs(x="Distance (towards lumen)", y="Cell Type Signal") + theme_light(base_size = 16)+theme_classic()+ scale_fill_manual(values=c("#E69F00", "#56B4E9"))


pval <- t.test(RUNX1_largebits$Score[RUNX1_largebits$CellType=="RA"],
               RUNX1_largebits$Score[RUNX1_largebits$CellType=="OA"])$p.value

```
