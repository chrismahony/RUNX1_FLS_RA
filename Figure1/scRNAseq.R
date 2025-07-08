

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

data <- Read10X(data.dir = "/rds/projects/c/croftap-stia-atac/CM_multiome/STIA_andATAC/STIA_2021_cesear_analysis/GEX.mtx.full.mouse.fibro.stia.dir/GEX.mtx.full.mouse.fibro.stia.dir/")


ncol(data)

#meta<-read_tsv("/rds/projects/c/croftap-stia-atac/CM_multiome/STIA_andATAC/STIA_2021_cesear_analysis/GEX.mtx.full.mouse.fibro.stia#.dir/GEX.mtx.full.mouse.fibro.stia.dir/cell_metadata_curated_20220830.tsv")

#meta <- meta %>% remove_rownames %>% column_to_rownames(var="barcode_id")
#meta <- subset(meta, select = -c(barcode))

stia2021_rna <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)

umap_fibros <- read.delim("/rds/projects/m/mahonyc-cesar-data/umap_fibros.gz", row.names=3)

umap_fibros<-as.matrix(umap_fibros)

stia2021_rna[['umap']] <- CreateDimReducObject(embeddings = umap_fibros, key = "UMAP_", global = T, assay = "RNA")

stia2021_rna <- stia2021_rna %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) 

cols <- ArchR::paletteDiscrete(stia2021_rna@meta.data[, "cluster.name"])

DimPlot(stia2021_rna, group.by = "cluster.name", cols=cols)

cols <- ArchR::paletteDiscrete(stia2021_rna@meta.data[, "condition"])

DimPlot(stia2021_rna, group.by = "condition", cols=cols)



DimPlot(stia2021_rna, group.by = "pseudo.bulk.level")

```

```{r}
Idents(stia2021_rna)<-'cluster.name'
DotPlot(stia2021_rna, features = c("F13a1", "Col22a1", "Runx2",   "Bglap", "Rgs5", "Ccl11",  "Clu", "Pi16", "Ccl7", "Ccl2", "Sfrp1", "Cfb", "Serpina3c", "C3", "C1qtnf3", "Col8a1", "Fmo2", "Chodl", "Col23a1", "Crabp1")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") + scale_size(range = c(2, 8)) +RotatedAxis()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))

DoHeatmap(stia2021_rna, features = c("F13a1", "Col22a1", "Runx2",   "Bglap", "Rgs5", "Ccl11",  "Clu", "Pi16", "Ccl7", "Ccl2", "Sfrp1", "Cfb", "Serpina3c", "C3", "C1qtnf3", "Col8a1", "Fmo2", "Chodl", "Col23a1", "Crabp1"))

```

```{r}
FeaturePlot(stia2021_rna, features = c("Runx2", "Bglap", "Alpl", "Ostn") )

FeaturePlot(stia2021_rna, features = c("Sox9", "Cilp", "Clu", "Ostn") )

FeaturePlot(stia2021_rna, features = c("Cd200") )

```

```{r}
counts <- as.data.frame(table(stia2021_rna$condition))
ggplot(data=counts, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="steelblue",color="black")+
  theme_minimal()+ coord_flip()
```



```{r}
Idents(stia2021_rna)<-'cluster.name'
allMarkers<-FindAllMarkers(stia2021_rna, only.pos = T)

library(gsfisher)

expressed_genes<-allMarkers$gene %>% unique()
annotation_gs <- fetchAnnotation(species="mm", ensembl_version=NULL, ensembl_host=NULL)

rm(list=ls()[! ls() %in% c("stia2021_rna","annotation_gs", "allMarkers")])
gc()


allMarkers$cluster %>% table

index <- match(allMarkers$gene, annotation_gs$gene_name)
allMarkers$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- expressed_genes
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]

seurat_obj.res <- allMarkers
seurat_obj <- stia2021_rna
seurat_obj.res <- seurat_obj.res[!is.na(seurat_obj.res$ensembl),]
ensemblUni <- na.omit(ensemblUni)



go.results <- runGO.all(results=seurat_obj.res,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="p_val_adj", p_threshold=0.05,
                  species = "mm")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=3, -p.val)






levels <- c("fibroblast_sublining_C1qtnf3_Col8a1", "fibroblast_sublining_Sfrp1_Cfb", "fibroblast_sublining_lining_Ccl7_Ccl2", "fibroblast__Chodl", "fibroblast_sublining_Fmo2", "fibroblast_sublining_Ccl11", "fibroblast_sublining_Serpina3c_C3", "fibroblast__Crabp1_Col23a1", "fibroblast_sublining_Pi16", "fibroblast_lining_F13a1_Col22a1")

go.results_f <- go.results[go.results$cluster %in% levels,]

sampleEnrichmentDotplot(go.results_f, selection_col = "description", selected_genesets = c( "cytoskeletal regulatory protein binding","collagen binding","regulation of Wnt signaling pathway", "osteoblast differentiation", "hyaluronan metabolic process", "myofibril","bone remodeling",  "response to interferon-beta", "complement activation", "chemokine activity", "regulation of actin filament organization", "positive regulation of phagocytosis", "metalloendopeptidase activity"), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T, sample_levels= )



dev.off()
rm(seurat_obj)
write.csv(go.results, "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/go_results.csv")


go.results_f <- go.results[go.results$description %in% c( "cytoskeletal regulatory protein binding","collagen binding","regulation of Wnt signaling pathway", "osteoblast differentiation", "hyaluronan metabolic process", "myofibril","bone remodeling",  "response to interferon-beta", "complement activation", "chemokine activity", "regulation of actin filament organization", "positive regulation of phagocytosis", "metalloendopeptidase activity"),]


go.results_f <- go.results_f[c(19:20,12:15,10:11,23,21,22,7,16,17,24,25,9, 1),]


#reorder pathways


pathways <- c( "cytoskeletal regulatory protein binding","collagen binding","regulation of Wnt signaling pathway", "osteoblast differentiation", "hyaluronan metabolic process", "myofibril","bone remodeling",  "response to interferon-beta", "complement activation", "chemokine activity", "regulation of actin filament organization", "positive regulation of phagocytosis", "metalloendopeptidase activity")

sampleEnrichmentDotplot(go.results_f, selection_col = "description", selected_genesets = pathways, sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)







library(readxl)
go_terms_selected <- read_excel("./go_terms_seletedNEW.xlsx")
go_terms_selected<-as.data.frame(go_terms_selected)
rownames(go_terms_selected)=go_terms_selected$...1
go_terms_selected$...1<-NULL


sampleEnrichmentDotplot(go_terms_selected, selection_col = "description", selected_genesets = unique(go_terms_selected$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)
```



```{r}
pt <- table(stia2021_rna$sample_id, stia2021_rna$cluster.name)
pt <- as.data.frame(pt)
library(tidyverse)
pt_fibroblast__Chodl = subset(pt, Var2 == "fibroblast__Chodl")
pt_fibroblast__Clu = subset(pt, Var2 == "fibroblast__Clu")
pt_fibroblast__Crabp1_Col23a1 = subset(pt, Var2 == "fibroblast__Crabp1_Col23a1")
pt_fibroblast__Runx2_Bglap = subset(pt, Var2 == "fibroblast__Runx2_Bglap")
pt_fibroblast_lining_F13a1_Col22a1 = subset(pt, Var2 == "fibroblast_lining_F13a1_Col22a1")
pt_fibroblast_mural_Rgs5 = subset(pt, Var2 == "fibroblast_mural_Rgs5")
pt_fibroblast_sublining_C1qtnf3_Col8a1 = subset(pt, Var2 == "fibroblast_sublining_C1qtnf3_Col8a1")
pt_fibroblast_sublining_Ccl11 = subset(pt, Var2 == "fibroblast_sublining_Ccl11")
pt_fibroblast_sublining_Fmo2 = subset(pt, Var2 == "fibroblast_sublining_Fmo2")
pt_fibroblast_sublining_lining_Ccl7_Ccl2 = subset(pt, Var2 == "fibroblast_sublining_lining_Ccl7_Ccl2")
pt_fibroblast_sublining_Pi16 = subset(pt, Var2 == "fibroblast_sublining_Pi16")
pt_fibroblast_sublining_Serpina3c_C3 = subset(pt, Var2 == "fibroblast_sublining_Serpina3c_C3")
pt_fibroblast_sublining_Sfrp1_Cfb = subset(pt, Var2 == "fibroblast_sublining_Sfrp1_Cfb")

pt_fibroblast_mural_Rgs5 = subset(pt, Var2 == "fibroblast_mural_Rgs5")
fibroblast_mural_Rgs5
```


```{r}
colnames(pt_fibroblast__Chodl) <- c("Var1", "Var2", "Chodl")
colnames(pt_fibroblast__Clu)<- c("Var1", "Var2", "Clu")
colnames(pt_fibroblast__Crabp1_Col23a1)<- c("Var1", "Var2", "Crabp1_Col23a1")
colnames(pt_fibroblast__Runx2_Bglap)<- c("Var1", "Var2", "Runx2_Bglap")
colnames(pt_fibroblast_lining_F13a1_Col22a1)<- c("Var1", "Var2", "F13a1_Col22a1")
colnames(pt_fibroblast_sublining_C1qtnf3_Col8a1)<- c("Var1", "Var2", "C1qtnf3_Col8a1")
colnames(pt_fibroblast_sublining_Ccl11)<- c("Var1", "Var2", "Ccl11")
colnames(pt_fibroblast_sublining_Fmo2)<- c("Var1", "Var2", "Fmo2")
colnames(pt_fibroblast_sublining_lining_Ccl7_Ccl2)<- c("Var1", "Var2", "Ccl7_Ccl2")
colnames(pt_fibroblast_sublining_Pi16)<- c("Var1", "Var2", "Pi16")
colnames(pt_fibroblast_sublining_Serpina3c_C3)<- c("Var1", "Var2", "Serpina3c_C3")
colnames(pt_fibroblast_sublining_Sfrp1_Cfb)<- c("Var1", "Var2", "Sfrp1_Cfb")
colnames(pt_fibroblast_mural_Rgs5)<- c("Var1", "Var2", "Rgs5")




#split up for each cluster
pt_master<-pt_fibroblast__Chodl
pt_master$Rgs5<-pt_fibroblast_mural_Rgs5$Freq
pt_master$Clu<-pt_fibroblast__Clu$Clu
pt_master$Crabp1_Col23a1<-pt_fibroblast__Crabp1_Col23a1$Crabp1_Col23a1
pt_master$F13a1_Col22a1<-pt_fibroblast_lining_F13a1_Col22a1$F13a1_Col22a1
pt_master$C1qtnf3_Col8a1<-pt_fibroblast_sublining_C1qtnf3_Col8a1$C1qtnf3_Col8a1
pt_master$Ccl11<-pt_fibroblast_sublining_Ccl11$Ccl11
pt_master$Fmo2<-pt_fibroblast_sublining_Fmo2$Fmo2
pt_master$Ccl7_Ccl2<-pt_fibroblast_sublining_lining_Ccl7_Ccl2$Ccl7_Ccl2
pt_master$Pi16<-pt_fibroblast_sublining_Pi16$Pi16
pt_master$Serpina3c_C3<-pt_fibroblast_sublining_Serpina3c_C3$Serpina3c_C3
pt_master$Sfrp1_Cfb<-pt_fibroblast_sublining_Sfrp1_Cfb$Sfrp1_Cfb
pt_master$Bglap<-pt_fibroblast__Runx2_Bglap$Runx2_Bglap
pt_master$Fmo2<-pt_fibroblast_sublining_Fmo2$Fmo2


pt_master <- pt_master %>% select(-one_of('Var2'))
library(tidyverse)
pt_master <- pt_master %>% remove_rownames %>% column_to_rownames(var="Var1")

pt_master <- pt_master/rowSums(pt_master)

pt_master <- as.data.frame(pt_master)

pt_master$condition<-rownames(pt_master)
```


```{r}
stia2021_rna$condition_corrected<-stia2021_rna$sample_id
Idents(stia2021_rna)<-'sample_id'
levels(stia2021_rna)
current.sample.ids <- c("control1_stia2017"   ,      "control2_stia2017" ,       
  "control3_stia2017"     ,    "control_cd45n_s1_stia2021",
  "control_cd45n_s2_stia2021", "control_cd45n_s3_stia2021",
  "control_d0_1_stia2018"   ,  "control_d0_2_stia2018"    ,
  "control_d0_3_stia2018"   ,  "day15_cd45n_s1_stia2021"  ,
 "day15_cd45n_s2_stia2021"  , "day15_cd45n_s3_stia2021"  ,
 "day1_cd45n_s1_stia2021"   , "day1_cd45n_s2_stia2021"   ,
 "day1_cd45n_s3_stia2021"   , "day22_cd45n_s1_stia2021"  ,
 "day22_cd45n_s2_stia2021"  , "day22_cd45n_s3_stia2021"  ,
 "day28_cd45n_s1_stia2021"  , "day28_cd45n_s2_stia2021"  ,
 "day28_cd45n_s3_stia2021"  , "day8_cd45n_s1_stia2021"   ,
 "day8_cd45n_s2_stia2021"   , "day8_cd45n_s3_stia2021"   ,
 "inflamedpeak1_stia2017"    ,"inflamedpeak2_stia2017"   ,
 "inflamedpeak3_stia2017"    ,"peak_d9_1_stia2018"       ,
 "peak_d9_2_stia2018"        ,"peak_d9_3_stia2018"       ,
 "resolved_d22_1_stia2018"   ,"resolved_d22_2_stia2018"  ,
 "resolved_d22_3_stia2018"   ,"resolving1_stia2017"      ,
 "resolving2_stia2017"       ,"resolving3_stia2017"      ,
 "resolving_d15_1_stia2018"  ,"resolving_d15_2_stia2018" ,
 "resolving_d15_3_stia2018" )
new.sample.ids <- c("control"   ,      "control" ,       
  "control"     ,    "control",
  "control", "control",
  "control"   ,  "control"    ,
  "control"   ,  "resolving"  ,
 "resolving"  , "resolving"  ,
 "initiation"   , "initiation"   ,
 "initiation"   , "resolved"  ,
 "resolved"  , "resolved"  ,
 "persistant"  , "persistant"  ,
 "persistant"  , "peak"   ,
 "peak"   , "peak"   ,
 "peak"    ,"peak"   ,
 "peak"    ,"peak"       ,
 "peak"        ,"peak"       ,
 "resolved"   ,"resolved"  ,
 "resolved"   ,"resolving"      ,
 "resolving"       ,"resolving"      ,
 "resolving"  ,"resolving" ,
 "resolving" )

stia2021_rna@meta.data[["condition_corrected"]] <- plyr::mapvalues(x = stia2021_rna@meta.data[["condition_corrected"]], from = current.sample.ids, to = new.sample.ids)


stia2021_rna$orig.ident_patho<-paste(stia2021_rna$condition_corrected, stia2021_rna$sample_id, sep="_")
Idents(stia2021_rna)<-'orig.ident_patho'
levels<-as.data.frame(levels(stia2021_rna))
colnames(levels)<-'levels_1'
library(splitstackshape)
levels<-cSplit(levels, splitCols="levels_1", sep="_")
pt_master$condition<-levels$levels_1_1

pt_master$check=rownames(pt_master)
```


```{r}
condition<-pt_master$condition

condition[[10]]<-'initiation'
condition[[11]]<-'initiation'
condition[[12]]<-'initiation'

condition[[13]]<-'resolving'
condition[[14]]<-'resolving'
condition[[15]]<-'resolving'

pt_master$condition<-condition


Chodl <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Chodl)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Chodl")+theme_classic()+RotatedAxis()
Chodl
res.aov_Chodl <- aov(Chodl ~ condition, data = pt_master)
summary(res.aov_Chodl)
stats_Chodl <- TukeyHSD(res.aov_Chodl)
stats_Chodl


Clu <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Clu)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Clu")+theme_classic()+RotatedAxis()
Clu
res.aov_Clu <- aov(Clu ~ condition, data = pt_master)
summary(res.aov_Clu)
stats_Clu <- TukeyHSD(res.aov_Clu)
stats_Clu

Crabp1_Col23a1 <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Crabp1_Col23a1)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Crabp1_Col23a1")+theme_classic()+RotatedAxis()
Crabp1_Col23a1
res.aov_Crabp1_Col23a1 <- aov(Crabp1_Col23a1 ~ condition, data = pt_master)
summary(res.aov_Crabp1_Col23a1)
stats_Crabp1_Col23a1 <- TukeyHSD(res.aov_Crabp1_Col23a1)
stats_Crabp1_Col23a1

#liniing
F13a1_Col22a1 <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=F13a1_Col22a1)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("F13a1_Col22a1")+theme_classic()+RotatedAxis()
F13a1_Col22a1
res.aov_F13a1_Col22a1 <- aov(F13a1_Col22a1 ~ condition, data = pt_master)
summary(res.aov_F13a1_Col22a1)
stats_F13a1_Col22a1 <- TukeyHSD(res.aov_F13a1_Col22a1)
stats_F13a1_Col22a1

#interesting
C1qtnf3_Col8a1 <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=C1qtnf3_Col8a1)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("C1qtnf3_Col8a1")+theme_classic()+RotatedAxis()
C1qtnf3_Col8a1
res.aov_C1qtnf3_Col8a1 <- aov(C1qtnf3_Col8a1 ~ condition, data = pt_master)
summary(res.aov_C1qtnf3_Col8a1)
stats_C1qtnf3_Col8a1 <- TukeyHSD(res.aov_C1qtnf3_Col8a1)
stats_C1qtnf3_Col8a1


#interesting
Pi16 <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Pi16)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Pi16")+theme_classic()+RotatedAxis()
Pi16
res.aov_Pi16 <- aov(Pi16 ~ condition, data = pt_master)
summary(res.aov_Pi16)
stats_Pi16 <- TukeyHSD(res.aov_Pi16)
stats_Pi16

Serpina3c_C3 <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Serpina3c_C3)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Serpina3c_C3")+theme_classic()+RotatedAxis()
Serpina3c_C3
res.aov_Serpina3c_C3 <- aov(Serpina3c_C3 ~ condition, data = pt_master)
summary(res.aov_Serpina3c_C3)
stats_Serpina3c_C3 <- TukeyHSD(res.aov_Serpina3c_C3)
stats_Serpina3c_C3

Sfrp1_Cfb <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Sfrp1_Cfb)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Sfrp1_Cfb")+theme_classic()+RotatedAxis()
Sfrp1_Cfb
res.aov_Sfrp1_Cfb <- aov(Sfrp1_Cfb ~ condition, data = pt_master)
summary(res.aov_Sfrp1_Cfb)
stats_Sfrp1_Cfb <- TukeyHSD(res.aov_Sfrp1_Cfb)
stats_Sfrp1_Cfb

Rgs5 <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Rgs5)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Rgs5")+theme_classic()+RotatedAxis()
Rgs5
res.aov_Rgs5 <- aov(Rgs5 ~ condition, data = pt_master)
summary(res.aov_Rgs5)
stats_Rgs5 <- TukeyHSD(res.aov_Rgs5)
stats_Rgs5

Ccl7_Ccl2 <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Ccl7_Ccl2)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Ccl7_Ccl2")+theme_classic()+RotatedAxis()
Ccl7_Ccl2
res.aov_Ccl7_Ccl2 <- aov(Ccl7_Ccl2 ~ condition, data = pt_master)
summary(res.aov_Ccl7_Ccl2)
stats_Ccl7_Ccl2 <- TukeyHSD(res.aov_Ccl7_Ccl2)
stats_Ccl7_Ccl2

Bglap <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Bglap)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Bglap")+theme_classic()+RotatedAxis()
Bglap
res.aov_Bglap <- aov(Bglap ~ condition, data = pt_master)
summary(res.aov_Bglap)
stats_Bglap <- TukeyHSD(res.aov_Bglap)
stats_Bglap

Ccl11 <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Ccl11)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Ccl11")+theme_classic()+RotatedAxis()
Ccl11
res.aov_Ccl11 <- aov(Ccl11 ~ condition, data = pt_master)
summary(res.aov_Ccl11)
stats_Ccl11 <- TukeyHSD(res.aov_Ccl11)
stats_Ccl11

Fmo2 <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=Fmo2)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Fmo2")+theme_classic()+RotatedAxis()
Fmo2
res.aov_Fmo2 <- aov(Fmo2 ~ condition, data = pt_master)
summary(res.aov_Fmo2)
stats_Fmo2 <- TukeyHSD(res.aov_Fmo2)
stats_Fmo2



```

```{r}
Idents(stia2021_rna)<-'cluster.name'
VlnPlot(stia2021_rna, features = c("Runx1", "Runx2", "Runx3", "Cbfb"), pt.size = 0, stack = T)+RotatedAxis()

```






```{r}
levels(stia2021_rna)

cols=c("fibroblast_lining_F13a1_Col22a1"="grey"  ,    
  "fibroblast__Runx2_Bglap"     ="grey"         ,
  "fibroblast_mural_Rgs5"         ="grey"       ,
  "fibroblast_sublining_Ccl11"    ="grey"       ,
  "fibroblast__Clu"               ="grey"       ,
  "fibroblast_sublining_Pi16"         ="black"   ,
  "fibroblast_sublining_lining_Ccl7_Ccl2"="grey",
  "fibroblast_sublining_Sfrp1_Cfb"      ="grey" ,
  "fibroblast_sublining_Serpina3c_C3"   ="grey" ,
 "fibroblast_sublining_C1qtnf3_Col8a1" ="grey" ,
 "fibroblast_sublining_Fmo2"           ="grey" ,
 "fibroblast__Chodl"                  ="grey"  ,
 "fibroblast__Crabp1_Col23a1"="grey" )

DimPlot(stia2021_rna, cols=cols)+NoLegend()



```

```{r}


contplus <- data.frame(stia2021_rna$umap@cell.embeddings)
umap_1 <- stia2021_rna$umap@cell.embeddings[,1]
umap_2 <- stia2021_rna$umap@cell.embeddings[,2]

# Initialize arrays for data subsets for each condition and for plots 
plotData <- list()
modGalaxyPlot <- list()

# Take each subset of data and generate a plot 

# Use grep to subset the data 
plotData[[1]] <- contplus[grepl(as.numeric(1), rownames(contplus)),]

# Generate galaxy plot for each condition
modGalaxyPlot[[1]] <- ggplot(plotData[[1]], aes(umap_1, umap_2)) +
  stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = F) +
  scale_fill_viridis(option = "YlOrRd") +
  coord_cartesian(expand = FALSE, xlim = c(min(umap_1), max(umap_2)), ylim = c(min(umap_1),max(umap_2))) +
  geom_point(shape = '.', col = 'white') 

modGalaxyPlot[[1]]
```
```{r}
library(TreeCorTreat)
genes<-c("Runx2", "Bglap", "Alpl", "Sox9", "Cilp", "Clu")

count<-stia2021_rna@assays[["RNA"]]@counts
rc <- Matrix::colSums(count)
sub_count <- count[genes,] %>% as.matrix
sub_norm <- log2(t(t(sub_count)/rc*1e6 + 1)) %>% as.matrix

df_marker <- t(sub_norm) %>% data.frame %>% mutate(barcode = rownames(.)) %>% gather(gene,expr,-barcode)
umap<-as.data.frame(stia2021_rna@reductions[["umap"]]@cell.embeddings)

df_marker$UMAP_1<-umap$UMAP_1
df_marker$UMAP_2<-umap$UMAP_2

ggplot(df_marker,aes(x = UMAP_1,y = UMAP_2,col = expr)) +
  geom_point(size = 0.01, shape = ".") +
  scale_colour_viridis_c(option = 'C',direction = 1) +
  facet_wrap(~ gene) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom')

```
```{r}

cell_meta=stia2021_rna@meta.data
cell_meta<-select(cell_meta, c(sample_id))
library(tibble)
cell_meta <- tibble::rownames_to_column(cell_meta, "barcode")
colnames(cell_meta)<-c("barcode", "sample")

cell_meta$UMAP_1<-umap$UMAP_1
cell_meta$UMAP_2<-umap$UMAP_2


stia2021_rna$orig.ident_patho<-paste(stia2021_rna$sample_id, stia2021_rna$condition, sep=".")
sample_meta<-as.data.frame(table(stia2021_rna@meta.data[["orig.ident_patho"]]))
library(splitstackshape)
sample_meta<-cSplit(sample_meta, splitCols = "Var1", sep = ".")
sample_meta$Freq<-NULL
colnames(sample_meta)<-c("sample", "condition")
sample_meta$study<-'ALL'

treecor_celldensityplot(cell_meta,
                        sample_meta,
                        row_variable = 'study',
                        col_variable = 'condition',
                        row_combined = F)

```

```{r}
library(TreeCorTreat)
genes<-c("Runx1", "Runx2")

count<-stia2021_rna@assays[["RNA"]]@counts
rc <- Matrix::colSums(count)
sub_count <- count[genes,] %>% as.matrix
sub_norm <- log2(t(t(sub_count)/rc*1e6 + 1)) %>% as.matrix

df_marker <- t(sub_norm) %>% data.frame %>% mutate(barcode = rownames(.)) %>% gather(gene,expr,-barcode)
umap<-as.data.frame(stia2021_rna@reductions[["umap"]]@cell.embeddings)

df_marker$UMAP_1<-umap$UMAP_1
df_marker$UMAP_2<-umap$UMAP_2

ggplot(df_marker,aes(x = UMAP_1,y = UMAP_2,col = expr)) +
  geom_point(size = 0.01, shape = ".") +
  scale_colour_viridis_c(option = 'C',direction = 1) +
  facet_wrap(~ gene) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom')

```



```{r}

library(ggtree)

Idents(stia2021_rna) <- "cluster.name"

rds <- BuildClusterTree(stia2021_rna, 
                        dims = 1:30)

data.tree <- Tool(object = rds, 
                  slot = "BuildClusterTree")



gg_tr <- ggtree(data.tree, 
                layout = "rectangular",
                alpha = .2,
                size = 3) +
  geom_tiplab(hjust = 1.1, size = 3) +
  geom_treescale() +
  theme(plot.margin = unit(c(1,5,1,1), "mm")) 

plot(gg_tr)

```

```{r}


output_file <- "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/mouse_stia.png"

meta<-stia2021_rna@meta.data
table(meta$name)

gs <- unique(meta$name)
gs <- unlist(strsplit(gs, "_"))
gs <- append("Cd34", gs)
gs <- append("Thy1", gs)
gs <- append("Fap", gs)
gs <- append("Col11a1", gs)
gs <- append("Prg4", gs)
gs <- append("Cxcl12", gs)
gs

ggs <- lapply(gs, function(g) {
  FeaturePlot(stia2021_rna,
              features = g) +
  theme_void() +
  theme(legend.position = "bottom", plot.title = element_text(size = 25, face = "bold")) +
  scale_color_gradientn(colours = c("grey90", "#74a9cf", "#023858"))
})

nc <- ceiling(sqrt(length(ggs)))

#png(filename = gsub(".png", paste0("_cell_lineages_top_mkr.png"), output_file),
#    width = 1700*nc, height = 2000*nc, 
#    res = 900)


plot_grid(plotlist = ggs, ncol = nc)




```

```{r}
library(scProportionTest)
test <- sc_utils(stia2021_rna)
table(stia2021_rna$condition)
prop.test <- permutation_test(test, cluster_identity = "cluster.name", sample_1="control", sample_2="initiation", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)

prop.test_peak <- permutation_test(test, cluster_identity = "cluster.name", sample_1="control", sample_2="peak", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test_peak, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)

prop.test_persis <- permutation_test(test, cluster_identity = "cluster.name", sample_1="control", sample_2="persistent", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test_persis, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)

prop.test_res <- permutation_test(test, cluster_identity = "cluster.name", sample_1="control", sample_2="resolved", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test_res, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)

prop.test_resing <- permutation_test(test, cluster_identity = "cluster.name", sample_1="control", sample_2="resolving", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test_resing, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)







```

```{r}
library(presto)
library(scMiko)
stia2021_rna$seurat_clusters <-stia2021_rna$cluster.name

features_expr <- findNetworkFeatures(object = stia2021_rna, method = "expr", min_pct = 0.5)


features_hvg <- findNetworkFeatures(object = stia2021_rna, method = "hvg", n_features = 2000)


features_dev <- findNetworkFeatures(object = stia2021_rna, method = "deviance", n_features = 2000)


feature.list <- list(expr = features_expr, hvg = features_hvg, deviance = features_dev)
ggVennDiagram::ggVennDiagram(feature.list) + scale_color_manual(values = rep("white", 3))

so.gene <- runSSN(object = stia2021_rna, features = unique(c(features_hvg, features_dev)), scale_free = T,
    robust_pca = F, data_type = "pearson", reprocess_sct = T, slot = c("scale"), batch_feature = NULL,
    pca_var_explained = 0.9, optimize_resolution = T, target_purity = 0.8, step_size = 0.05, n_workers = parallel::detectCores(),
    verbose = F)
cowplot::plot_grid(so.gene@misc$scale_free$optimization.plot, so.gene@misc$scale_free$distribution.plot$`2`,
    labels = "AUTO")
```
```{r}
plt_connectivity <- SSNConnectivity(so.gene, quantile_threshold = 0.85, raster_dpi = 500)

plt_connectivity$plot_edge + labs(title = "Network Connectivity")
```


```{r}
# specify pruning threshold [0,1] (low values = less pruning, high values = more pruning)
prune.threshold <- 0.1

# get feature-specific connectivities (wi)
df.wi <- pruneSSN(object = so.gene, graph = "RNA_snn_power", prune.threshold = prune.threshold,
    return.df = T)

# visualize
plt.prune <- df.wi %>%
    ggplot(aes(x = wi_l2)) + geom_histogram(bins = 30) + geom_vline(xintercept = prune.threshold,
    linetype = "dashed", color = "tomato") + labs(x = "Degree (L2 norm)", y = "Count", title = "Network Pruning",
    subtitle = paste0(signif(100 * sum(df.wi$wi_l2 <= prune.threshold)/nrow(df.wi), 3), "% (", sum(df.wi$wi_l2 <=
        prune.threshold), "/", nrow(df.wi), ") genes pruning")) + theme_miko(grid = T)

print(plt.prune)

mod.list <- pruneSSN(object = so.gene, graph = "RNA_snn_power", prune.threshold = prune.threshold)
str(mod.list)
```

```{r}
plt_connectivity_with_features <- SSNConnectivity(so.gene, gene.list = mod.list, quantile_threshold = 0.85,
    raster_dpi = 500, node.size.max = 6, node.size.min = 2, node.alpha = 0.6, node.weights = T,
    node.color = "grey80")

# generate interactive network plot using plotly
plotly::ggplotly(plt_connectivity_with_features$plot_network)
```

```{r}
plt_connectivity_with_features <- SSNConnectivity(so.gene, gene.list = mod.list, quantile_threshold = 0.85,
    raster_dpi = 500, node.size.max = 6, node.size.min = 2, node.alpha = 0.6, node.weights = T,
    node.color = "grey80")

ssn.summary <- summarizeModules(cell.object = stia2021_rna, gene.object = so.gene, gene.list = mod.list,
    module.type = "ssn", n.workers = parallel::detectCores(), connectivity_plot = plt_connectivity_with_features$plot_edge)


# cluster-level heatmap of module activities
plt.ssn.gene.hm.expr <- heatmaply::heatmaply((ssn.summary$data.heatmap), scale = "column", scale_fill_gradient_fun = scale_fill_miko(),
    xlab = "Module", ylab = "Cluster", main = "SSN Module Activity")

plt.ssn.gene.hm.expr
```


```{r}
# get list of module-level summary plots
plt.ssn.gene <- ssn.summary$plt.summary

# assemble figure panels and visualize
x <- plt.ssn.gene$m7

cowplot::plot_grid(cowplot::plot_grid(NULL, x$cell.umap + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)), x$gene.umap + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)), NULL, nrow = 1, labels = c("", "A", "B", ""),
    rel_widths = c(1, 4, 4, 1)), cowplot::plot_grid(x$bp.enrich, x$mf.enrich, x$cc.enrich, nrow = 1,
    labels = c("C", "D", "E")), ncol = 1)
```

```{r}
"Jun" %in% mod.list[["m8"]]
```





```{r}
Idents(stia2021_rna)<-'cluster.name'
levels(stia2021_rna)


input_string <- '@All(@LL(fibroblast_lining_F13a1_Col22a1),@SL(fibroblast_sublining_lining_Ccl7_Ccl2,fibroblast_sublining_Serpina3c_C3,fibroblast_sublining_Serpina3c_C3,fibroblast_sublining_Fmo2,fibroblast__Crabp1_Col23a1,fibroblast_sublining_Ccl11,fibroblast_sublining_Pi16,fibroblast_sublining_Sfrp1_Cfb,fibroblast_sublining_C1qtnf3_Col8a1,fibroblast__Chodl),@Peri(fibroblast_mural_Rgs5),@Chon(fibroblast__Clu),@Osteo(fibroblast__Runx2_Bglap))'
hierarchy_structure <- extract_hrchy_string(input_string,special_character = '@', plot = T)


cell_meta2=stia2021_rna@meta.data

cell_meta$celltype<-cell_meta2$cluster.name


res_ctprop_full <- treecor_ctprop(hierarchy_structure,
                                  cell_meta,
                                  sample_meta,
                                  response_variable = 'condition',
                                  method = "aggregate",
                                  analysis_type = "pearson",
                                  num_permutations = 100)

res_ctprop <- res_ctprop_full[[1]] %>% mutate(condition.absolute_cor = abs(condition.pearson))
head(res_ctprop)


treecortreatplot_p<-treecortreatplot(hierarchy_structure,
                 annotated_df = res_ctprop,
                 response_variable = 'condition',
                 color_variable = 'direction',
                 size_variable = 'absolute_cor',
                 alpha_variable = 'p.sign',
                 font_size = 12,
                 nonleaf_label_pos = 0.4,
                 nonleaf_point_gap = 0.2)
```

```{r}
res_expr_full <- treecor_expr(count,
                              hierarchy_structure,
                              cell_meta,
                              sample_meta,
                              response_variable = 'condition',
                              method = 'aggregate',
                              analysis_type = 'cancor',
                              num_permutations = 100)
names(res_expr_full)
```



```{r}

stia2021_rna$sample_condition<-paste(stia2021_rna$sample_id, stia2021_rna$condition, sep=".")
Idents(stia2021_rna) <- "experiment"
levels(stia2021_rna)


stia2021_rna_only<-subset(stia2021_rna, idents=c(         
  "stia2021"))

Idents(stia2021_rna_only) <- "pseudo.bulk.level"
levels(stia2021_rna_only)

stia2021_rna_only<-subset(stia2021_rna_only, idents=c("sublining", "lining"))


cts<-AggregateExpression(stia2021_rna_only, group.by = c("sample_condition"), assays = "RNA", slot = "counts", return.seurat = F)

cts<-cts$RNA
cts<-as.data.frame(cts)
meta_data=colnames(cts)
meta_data<-as.data.frame(meta_data)
library(splitstackshape)
meta_data$to_split<-meta_data$meta_data
meta_data<-cSplit(meta_data, splitCols = "to_split", sep=".")
colnames(meta_data)<-c("all", "sample", "condition")
meta_data$all<-as.factor(meta_data$all)
meta_data$sample<-as.factor(meta_data$sample)
meta_data$cluster<-as.factor(meta_data$condition)



```


```{r}

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = meta_data,
                                  design = ~1)

dds <- scran::computeSumFactors(dds)
print(dds)
print(quantile(rowSums(counts(dds))))

mingenecount <- 200
#maxgenecount <- quantile(rowSums(counts(dds)), 0.99)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount #& rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

dds@colData[['condition']] <- factor(dds@colData[['condition']],
                                     levels = c("control",
                                                "initiation",
                                                "peak",
                                                "resolving",
                                                "resolved",
                                                "persistent"))

dds@colData[['condition']] <- as.factor(dds@colData[['condition']])

design(dds) <- as.formula(paste0("~", "condition"))

print(design(dds))

dds <- DESeq(dds, test = "Wald")


print(resultsNames(dds))
targetvar <- "condition"
comps <- data.frame(t(combn(unique(as.character(meta[[targetvar]])), 2)))
head(comps)

ress <- apply(comps, 1, function(cp) {
  print(cp)
  res <- data.frame(results(dds, contrast=c(targetvar, cp[1], cp[2])))
  res[["gene"]] <- rownames(res)
  res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
  res
})



res <- Reduce(rbind, ress)

head(res)
comps

res %>% 
  filter(padj < 0.01) %>%
  mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
  arrange(desc(abs(score))) -> subres

head(subres)
dim(subres)

length(unique(subres$gene))


library(ComplexHeatmap)

      if(length(unique(subres$gene)) > 10) {
      vsd <- tryCatch({
        vst(dds, blind=TRUE)
      }, error=function(e) {
        message(e)
        print(e)
        return(NULL)
      })
      
      if(!is.null(vsd)) {
        print(dim(assay(vsd)))
        print(head(assay(vsd), 3))
        vsd_mat <- assay(vsd)
        
        feats <- unique(subres$gene)
        print(length(feats))
        
        # Sub-set matrix to relevant features
        sub_vsd_mat <- vsd_mat[rownames(vsd_mat) %in% feats, ]
        scale_sub_vsd <- t(scale(t(sub_vsd_mat)))
        head(scale_sub_vsd)
        dim(scale_sub_vsd)
      }
      }

meta_data %>%
          mutate('condition' = ifelse(condition == 'control', 't0-control',
                                      ifelse(condition == 'initiation', 't1-initiation',
                                             ifelse(condition == 'peak', 't3-peak',
                                                    ifelse(condition == 'resolving', 't4-resolving',
                                                           ifelse(condition == 'resolved', 't5-resolved',
                                                                  ifelse(condition == 'persistent', 't6-persistent', 'WOH!!!'))))))) %>%
          arrange(condition) %>%
          data.frame -> sub_meta



ss_sm <- sub_meta[, c("sample", "condition")]
ss_sm <- sub_meta[, c("condition")]

col_ann <- HeatmapAnnotation(df = ss_sm)


topedges <- 0.05

        

    ggnet <- Rfast::cora(t(scale_sub_vsd))
    ggnet_gather <- reshape2::melt(ggnet, id.vars = "V1")
    print(dim(ggnet_gather))
    
    ggnet_gather %>%
      filter(value > 0) %>%
      filter(Var1 != Var2) -> ggnet_gather
    
    print(dim(ggnet_gather))
    
    if(nrow(ggnet_gather) > 0) {
      
      edges <- arrange(ggnet_gather, -value)
      print(head(edges))
      print(tail(edges))
      edges <- edges[seq(1, nrow(edges), by = 2), ]
      print(head(edges))
      print(tail(edges))
      qtop <- quantile(edges$value, 1-topedges)
      print(qtop)
      
      edges %>%
        filter(value > qtop) -> top_edges
      
      print(dim(top_edges))
      
      top_edges %>%
        mutate('idx' = paste(Var1, Var2, sep = "_")) -> top_edges# %>%
      #dplyr::distinct(idx, .keep_all = TRUE) -> top_edges
      
      colnames(top_edges) <- c("Source", "Target", "Weight", "Id")
      #top_edges %>%
      #  select(c(Id, Source, Target, Weight)) -> top_edges
      
      print(head(top_edges))
      print(dim(top_edges))
      
      # - community detection
      # - igraph definition
      g <- igraph::graph_from_data_frame(top_edges[, c("Source", "Target")],
                                         directed = FALSE)
      g <- igraph::set_edge_attr(g, "weight", value = top_edges$Weight)
      g <- igraph::set_edge_attr(g, "name", value = top_edges$Id)  
      # leiden
      leiden_mod <- igraph::cluster_leiden(g, objective_function = "modularity")
      mods <- data.frame(cbind(igraph::V(g)$name, leiden_mod$membership))
      colnames(mods) <- c("Id", "leiden_gene_cluster")
      imods <- names(table(mods$leiden_gene_cluster)[table(mods$leiden_gene_cluster) > 10])
      print(imods)
      
      if(length(imods) > 0) {
        mods <- filter(mods, leiden_gene_cluster %in% imods)
        head(mods)
        
        row_cl <- data.frame('gene' = mods$Id,
                             'gene_cluster' = paste0("K", mods$leiden_gene_cluster))
        
        row_cl %>%
          arrange(gene_cluster) %>%
          data.frame -> row_cl
        
        rownames(row_cl) <- row_cl$gene
        
              }
    }
  
colnames(scale_sub_vsd)



    
    


#scale_sub_vsd <- subset(scale_sub_vsd, select=c("control_cd45n_s1_stia2021_fibroblast" ,"control_cd45n_s2_stia2021_fibroblast" #,"control_cd45n_s3_stia2021_fibroblast",
#"day1_cd45n_s1_stia2021_fibroblast"   , "day1_cd45n_s2_stia2021_fibroblast"   , "day1_cd45n_s3_stia2021_fibroblast"  ,
#"day8_cd45n_s1_stia2021_fibroblast"   , "day8_cd45n_s2_stia2021_fibroblast"   , "day8_cd45n_s3_stia2021_fibroblast",
#"day15_cd45n_s1_stia2021_fibroblast"  , "day15_cd45n_s2_stia2021_fibroblast"  , "day15_cd45n_s3_stia2021_fibroblast",  
#  "day22_cd45n_s1_stia2021_fibroblast",   "day22_cd45n_s2_stia2021_fibroblast",   "day22_cd45n_s3_stia2021_fibroblast",  
# "day28_cd45n_s1_stia2021_fibroblast" ,  "day28_cd45n_s2_stia2021_fibroblast"  , "day28_cd45n_s3_stia2021_fibroblast"  ))


scale_sub_vsd <- subset(scale_sub_vsd, select=c("control_cd45n_s1_stia2021.control",  "control_cd45n_s2_stia2021.control",  "control_cd45n_s3_stia2021.control" ,
"day1_cd45n_s1_stia2021.initiation" , "day1_cd45n_s2_stia2021.initiation" , "day1_cd45n_s3_stia2021.initiation", 
 "day8_cd45n_s1_stia2021.peak"       , "day8_cd45n_s2_stia2021.peak"      ,  "day8_cd45n_s3_stia2021.peak" ,
"day15_cd45n_s1_stia2021.resolving" , "day15_cd45n_s2_stia2021.resolving" , "day15_cd45n_s3_stia2021.resolving",  
"day22_cd45n_s1_stia2021.resolved"  , "day22_cd45n_s2_stia2021.resolved"  , "day22_cd45n_s3_stia2021.resolved"  ,
 "day28_cd45n_s1_stia2021.persistent", "day28_cd45n_s2_stia2021.persistent", "day28_cd45n_s3_stia2021.persistent"))
    
        
if(exists("row_cl")) {
    row_ann <- rowAnnotation(df = row_cl[, -1])    
    
    draw(
      Heatmap(scale_sub_vsd[rownames(row_cl),], 
              top_annotation = col_ann, 
              right_annotation = row_ann,
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = FALSE,
              show_row_names = F,
              show_column_names = F))
}



library(dplyr)
library(readr)
grn <- read_csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/df_grn2_new_FINAL.csv")
grn_Runx1 <- grn %>% filter(tf == "RUNX1")
grn_Runx1$Target <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(grn_Runx1$gene), perl=TRUE)

scale_sub_vsd_f <- scale_sub_vsd[rownames(scale_sub_vsd) %in% grn_Runx1$Target,]
Heatmap(scale_sub_vsd_f, 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F,
              border=T)


library(gsfisher)
#annotation_gs<-fetchAnnotation(species = "mm")


index <- match(res$gene, row_cl$gene)
res$cluster <- row_cl$gene_cluster[index]
res_cleaned<-na.omit(res)

index <- match(res_cleaned$gene, annotation_gs$gene_name)
res_cleaned$ensembl <- annotation_gs$ensembl_id[index]



FilteredGeneID <- unique(res$gene)
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]
ensemblUni <- na.omit(ensemblUni)
res_cleaned<-na.omit(res_cleaned)

go.results <- runGO.all(results=res_cleaned,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="padj", p_threshold=0.05,
                  species = "mm")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=5, -p.val)
sampleEnrichmentDotplot(go.results, selection_col = "description", selected_genesets = c("fat cell differentiation", "lipid storage", "embryonic limb morphogenesis", "bone morphogenesis", "metallopeptidase activity", "extracellular matrix organization", "mitotic nuclear division", "nuclear chromosome segregation", "interleukin-27-mediated signaling pathway", "regulation of type I interferon production", "regulation of JUN kinase activity", "lamellipodium","ribosome biogenesis", "ribosomal large subunit biogenesis"), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)


```


```{r}

pc <- prcomp(t(scale_sub_vsd),
             center = TRUE,
            scale. = TRUE)


m <- pc$x %>% as.matrix() 


library(harmony)
harmony_embeddings <- harmony::HarmonyMatrix(
    m, meta_data, c('sample'), do_pca = F, verbose= TRUE
)

colours <- list('condition' = ArchR::paletteDiscrete(stia2021_rna@meta.data[, "condition"]))
col_ann <- HeatmapAnnotation(df = meta_data[,c("condition")], col=colours)

Heatmap(cor(t(harmony_embeddings)),  show_row_names = F,
        show_column_names = F, column_dend_reorder = TRUE, 
        name="harmony space\ncorrelation",
        raster_quality = 3,
        use_raster = TRUE, top_annotation = col_ann)
```











