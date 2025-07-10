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
grn <- read.csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/df_grn2_new_FINAL.csv")
grn_RUNX1 <- grn %>% filter(tf == "RUNX1")

merging_A<- AddModuleScore(merging_A, features = list(grn_runx1$gene), name = "RUNX1_grn")
SpatialFeaturePlot(merging_A, features = "RUNX1_grn1", images = "slice1.17", min.cutoff = "q5",
  max.cutoff = "q95" )

SpatialFeaturePlot(merging_A, features = "RUNX1_grn1", images = "slice1.11", min.cutoff = "q5",
  max.cutoff = "q95")

SpatialFeaturePlot(merging_A, features = "RUNX1_grn1", images = "slice1.16", min.cutoff = "q5",
  max.cutoff = "q95")

merging_largebits_RA<-AddModuleScore(merging_largebits_RA, features = list(grn_runx1$gene), name = "RUNX1_grn")
merging_largebits_OA<-AddModuleScore(merging_largebits_OA, features = list(grn_runx1$gene), name = "RUNX1_grn")


RUNX1_RA<-merging_largebits_RA@meta.data
RUNX1_OA<-merging_largebits_OA@meta.data


RUNX1_RA <- subset(RUNX1_RA, select = c("RUNX1_grn1"))
RUNX1_OA <- subset(RUNX1_OA, select = c("RUNX1_grn1"))


RUNX1_RA_df <- data.frame(Dist=all_coords_dist_norm_largebits_RA$dist1, 
				Score=as.numeric(RUNX1_RA$RUNX1_grn1) / as.numeric(max(RUNX1_RA$RUNX1_grn1)), #scaling cell type signature to max here 
				CellType="RA")

RUNX1_OA_df <- data.frame(Dist=all_coords_dist_norm_largebits_OA$dist1, 
				Score=as.numeric(RUNX1_OA$RUNX1_grn1) / as.numeric(max(RUNX1_OA$RUNX1_grn1)), #scaling cell type signature to max here 
				CellType="OA")

RUNX1_largebits<-rbind(RUNX1_RA_df, RUNX1_OA_df)


ggplot(RUNX1_largebits, aes(Dist, Score, group=CellType, colour=CellType)) + geom_smooth(alpha=.1) + labs(x="Distance (towards lumen)", y="Cell Type Signal") + theme_light(base_size = 16)+theme_classic()+
  scale_colour_manual(values = c("RA" = "red", "OA" = "grey")) 


RUNX1_largebits$CellType_fac <- as.factor(RUNX1_largebits$CellType)

library(mgcv)
model <- gam(Score ~ s(Dist, by = CellType_fac) + CellType_fac, data = RUNX1_largebits)
summary(model)

```
