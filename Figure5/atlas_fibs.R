

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
DimPlot(aggv2, group.by = "CD45")+NoLegend()
FeaturePlot(aggv2, features="Ptprc")
FeaturePlot(aggv2, features="Col1a1")
FeaturePlot(aggv2, features="Cd68")
FeaturePlot(aggv2, features="Pecam1")
FeaturePlot(aggv2, features="Runx2")
FeaturePlot(aggv2, features="Clu")

DimPlot(aggv2, group.by = "orig.ident")+NoLegend()

```
```{r}

aggv2_h<-RunHarmony.Seurat_CM(aggv2, group.by.vars = "orig.ident")
aggv2_h<-RunUMAP(aggv2_h, reduction="harmony", dims=1:20)
DimPlot(aggv2_h, group.by = "orig.ident")+NoLegend()

FeaturePlot(aggv2_h, features="Ptprc")
FeaturePlot(aggv2_h, features="Col1a1")
FeaturePlot(aggv2_h, features="Cd68")
FeaturePlot(aggv2_h, features="Acta2")
FeaturePlot(aggv2_h, features="Runx2")
FeaturePlot(aggv2_h, features="Pecam1")

```
```{r}
pt<-as.data.frame(table(aggv2$orig.ident))
ggplot(data=pt, aes(x=Freq, y=Var1)) +
  geom_bar(stat="identity")+theme_ArchR()+
  geom_text(aes(label=Freq), vjust=0.6, size=2.5, nudge_x = 300)

DimPlot(aggv2, group.by = "integrated_snn_res.0.2")
DimPlot(aggv2, group.by = "integrated_snn_res.0.1")
DimPlot(aggv2, group.by = "integrated_snn_res.0.05")

FeaturePlot(aggv2, features="Clu")
FeaturePlot(aggv2, features="Prg4")


```


```{r}


DimPlot(aggv2, group.by = "CD45")+NoLegend()

write.csv(res0.2markers, "/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/pipeline_all_samples/res0.2markers.csv")
```


```{r}
DimPlot(aggv2, group.by = "integrated_snn_res.0.2", label = T)

```

```{r}
FeaturePlot(aggv2, features="Cd68")
FeaturePlot(aggv2, features="Cd14")
FeaturePlot(aggv2, features="Cd14")

```
```{r}
aggv2$named <- aggv2@meta.data[["integrated_snn_res.0.2"]]
Idents(aggv2) <- 'named'
levels(aggv2)
current.sample.ids <- c("0" , "1" , "10", "11", "12", "13", "14", "15", "16" ,"17", "18", "19" ,"2" , "20", "21" ,"3" , "4",  "5" , "6" , "7" , "8" , "9" )
new.sample.ids <- c("mono" , "mac" , "pericytes", "mono", "Bcells", "Muscle", "fibs", "chondrocytes", "RBCs" ,"Myeloid_proj", "neuron", "mono" ,"fibs" , "fibs", "mac" ,"fibs" , "osteoblasts",  "mac" , "endo_cells" , "chondrocytes" , "prolif" , "Tcells" )

aggv2@meta.data[["named"]] <- plyr::mapvalues(x = aggv2@meta.data[["named"]], from = current.sample.ids, to = new.sample.ids)

cols <- ArchR::paletteDiscrete(aggv2@meta.data[, "named"])
DimPlot(aggv2, group.by = "named", label = F, cols=cols, raster=FALSE)+NoAxes()

ncol(aggv2)

```
```{r}
Idents(aggv2) <- 'named'
DotPlot(aggv2, features = c("Cd14", "Cd68", "Acta2", "Cd79a", "Tnnt3", "Pdgfra", "Clu", "Hbb-bs", "Gata2", "Mpz", "Bglap", "Pecam1", "Top2a", "Cd3d"))+RotatedAxis()



```
```{r}

Idents(aggv2) <- 'named'
DotPlot(aggv2, features = c("Cd14", "Cd68", "Acta2", "Cd79a", "Tnnt3", "Pdgfra", "Clu", "Hbb-bs", "Gata2", "Mpz", "Bglap", "Pecam1", "Top2a", "Cd3d"))+RotatedAxis()
```
```{r}

DimPlot(aggv2, split.by = "CD45", cols=cols)
DimPlot(aggv2, group.by = "CD45")

pt <- table(aggv2$named, aggv2$CD45)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)



ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion")+  theme(legend.title = element_blank())+RotatedAxis()+theme_ArchR()

```
```{r}
library(ggExtra)
p4 <- ggplot(aggv2@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggMarginal(p4, type = "histogram", fill="lightgrey")
p4
```
```{r}
p5 <- ggplot(aggv2@meta.data, aes(x=log10(nCount_RNA), y=log10(percent.mt))) + geom_point() + geom_smooth(method="lm")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggMarginal(p5, type = "histogram", fill="lightgrey")
```


```{r}

Idents(aggv2)<-'all'
VlnPlot(aggv2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)


```




```{r}
Idents(aggv2)<-'orig.ident'
levels(aggv2)
aggv2$inflammation<-aggv2@active.ident


 current.sample.ids <- c("S3_DAY1posB_r"      , "S7_CONTROLposB_r"   , "S1_Day22posA"       , "S3_Day22posB"    ,    "S5_Day22posC"  ,     
  "S1_Day15posB"  ,      "S3_Day15posC" ,       "S5_Day28posA"  ,      "S7_Day28posB"   ,     "S1_Day28posC"       ,
 "S1_DAY1posA"     ,    "S1_DAY8posA"   ,      "S3_DAY8posB"    ,     "S5_DAY1posC"     ,    "S5_DAY8posC"      ,  
 "S7_ControlposC"   ,   "S7_DAY15posA"  ,      "CD45P"          ,     "AIA_infla1_cd45pos",  "AIA_infla2_cd45pos", 
 "AIA_rest1_cd45pos" ,  "D14_1_pos_renamed" ,  "D14_2_pos_renamed",   "D14_3_pos_renamed" ,  "rest_pos_renamed"   ,
 "Day0pos_cr_V5"      , "Day2_pos3"         ,  "Day4_pos1"        ,   "Day4_pos2"         ,  "Day4_pos3"          ,
 "Day7pos1"       ,     "Day7pos2"          ,  "Day7pos3"         ,   "Con1_CD45pos"      ,  "Con3_CD45pos"       ,
 "Infla1_CD45pos"  ,    "CIA_Control_pos_CR5", "Con2_CD45pos"     ,   "Infla2_CD45pos"    ,  "Infla3_CD45pos"     ,
 "Infla4_CD45pos"   ,   "S8_CONTROLnegB_r"  ,  "S4_Day8negB"       ,  "S2_Day22negA"      ,  "S4_Day22negB"       ,
 "S6_Day22negC"      ,  "S2_Day15negB"      ,  "S4_Day15negC"     ,   "S6_Day28negA"      ,  "S8_Day28negB"       ,
 "S2_Day28negC"       , "S2_DAY1negA_cmAB2" ,  "S2_DAY8negA"      ,   "S4_DAY1negB_cmAB2" ,  "S6_DAY1negC_cmAB2"  ,
 "S6_DAY8negC"         ,"S8_ControlnegC"    ,  "S8_DAY15negA"     ,   "CD45N"             ,  "AIA_infla1_cd45neg" ,
 "AIA_infla2_cd45neg" , "AIA_rest1_cd45neg" ,  "D14_1_neg_renamed" ,  "D14_2_neg_renamed" ,  "D14_3_neg_renamed"  ,
 "rest_neg_renamed"   , "Day0neg_cr_V5"     ,  "Day2_neg3"        ,   "Day4_neg1"         ,  "Day4_neg2"          ,
 "Day4_neg3"          , "Day7_neg1"         ,  "Day7_neg2"        ,   "Day7_neg3"         ,  "Con1_CD45neg"       ,
 "Con2_CD45neg"       , "Con3_CD45neg"      ,  "Infla1_CD45neg"   ,   "Infla2_CD45neg"    ,  "Con4_CD45neg"       ,
 "Infla3_CD45neg"   )

new.sample.ids<-c("initiation"      , "rest"   , "resed"       , "resed"    ,    "resed"  ,     
  "resing"  ,      "resing" ,       "persis"  ,      "persis"   ,     "persis"       ,
 "initiation"     ,    "peak"   ,      "peak"    ,     "initiation"     ,    "peak"      ,  
 "rest"   ,   "resing"  ,      "rest"          ,     "peak",  "peak", 
 "rest" ,  "resed" ,  "resed",   "resed" ,  "rest"   ,
 "rest"      , "peak"         ,  "Eresing"        ,   "Eresing"         ,  "Eresing"          ,
 "resing"       ,     "resing"          ,  "resing"         ,   "rest"      ,  "rest"       ,
 "peak"  ,    "rest", "rest"     ,   "peak"    ,  "peak"     ,
 "peak"   ,   "rest"  ,  "peak"       ,  "resed"      ,  "resed"       ,
 "resed"      ,  "resing"      ,  "resing"     ,   "persis"      ,  "persis"       ,
 "persis"       , "initiation" ,  "peak"      ,   "initiation" ,  "initiation"  ,
 "peak"         ,"rest"    ,  "resing"     ,   "rest"             ,  "peak" ,
 "peak" , "rest" ,  "resed" ,  "resed" ,  "resed"  ,
 "rest"   , "rest"     ,  "peak"        ,   "Eresing"         ,  "Eresing"          ,
 "Eresing"          , "resing"         ,  "resing"        ,   "resing"         ,  "rest"       ,
 "rest"       , "rest"      ,  "peak"   ,   "peak"    ,  "rest"       ,
 "peak"   )

aggv2@meta.data[["inflammation"]] <- plyr::mapvalues(x = aggv2@meta.data[["inflammation"]], from = current.sample.ids, to = new.sample.ids)

STIA<-rep("STIA",36)
AIA<-rep("AIA",30)
CIA<-rep("CIA",15)

model<-c(STIA, AIA, CIA)

Idents(aggv2)<-'orig.ident'
current.sample.ids <- levels(aggv2)

new.sample.ids<-model

aggv2$model<-aggv2$orig.ident
aggv2@meta.data[["model"]] <- plyr::mapvalues(x = aggv2@meta.data[["model"]], from = current.sample.ids, to = new.sample.ids)

Idents(aggv2)<-'model'

DimPlot(aggv2, split.by = "model", label = T, raster=F)


```

```{r}
Idents(aggv2)<-'model'
levels(aggv2)
no_AIA <- subset(aggv2, idents=c("STIA", "CIA"))

aggv2$CD45 %>% unique()

Idents(no_AIA)<-'CD45'
no_AIA <- subset(no_AIA, idents=c("neg"))

Idents(no_AIA)<-'named'
levels(no_AIA)

no_AIA <- subset(no_AIA, idents=levels(no_AIA)[c(3,6,7,11,12)])


no_AIA <- no_AIA %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>% RunUMAP(dims = 1:30)

no_AIA <- no_AIA %>% RunUMAP(1:50)

DimPlot(no_AIA, group.by="named")


DimPlot(no_AIA, group.by="orig.ident")


library(harmony)

no_AIA <- RunHarmony(no_AIA, "orig.ident")
no_AIA <- no_AIA %>% RunUMAP(dims=1:50, reduction="harmony")
DimPlot(no_AIA, group.by="orig.ident")
DimPlot(no_AIA, group.by="named")

DimPlot(no_AIA, split.by="model")



```





```{r}
pt_2<-as.data.frame(table(aggv2$orig.ident))
ggplot(data=pt_2, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity")+theme_ArchR()+RotatedAxis()+
  geom_text(aes(label=Freq), vjust=0, nudge_y = 400)+ coord_flip()
```


```{r}


ribo.genes <- grep(pattern = "(^Rpl|^Rps|^Mrp)", x = rownames(x = aggv2@assays[["RNA"]]), value = TRUE)
percent.ribo <- Matrix::colSums(aggv2@assays[["RNA"]][ribo.genes, ])/Matrix::colSums(aggv2@assays[["RNA"]])
aggv2[["percent.ribro"]] <- percent.ribo


VlnPlot(aggv2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribro"), ncol = 4, pt.size = 0)

FeaturePlot(aggv2, "percent.ribro")
```
```{r}
Idents(aggv2) <- 'named'
levels(aggv2)


stromal<-subset(aggv2, idents=c("pericytes"  , "fibs" ,        "chondrocytes",  "osteoblasts",  "endo_cells"))

Idents(aggv2) <- 'CD45'
stromal<-subset(aggv2, idents=c("neg"))

stromal <- stromal %>%
    ScaleData() %>%
    FindVariableFeatures(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)


DimPlot(stromal, group.by = "model")+NoLegend()





#stromal@meta.data$orig.ident <- as.factor(stromal@meta.data$orig.ident)
stromal_h<-RunHarmony.Seurat_CM(stromal, group.by.vars = "model")
stromal_h<-RunUMAP(stromal_h, reduction="harmony", dims = 1:30)
DimPlot(stromal_h, group.by = "orig.ident")+NoLegend()
DimPlot(stromal_h, group.by = "model")



stromal_h<-FindNeighbors(stromal_h, reduction = "harmony",dims = 1:30)
stromal_h<-FindClusters(stromal_h, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4))
stromal_h<-FindClusters(stromal_h, resolution = c(0.5, 0.6, 0.7))

DimPlot(stromal_h, group.by="RNA_snn_res.0.4", label=T)
```


```{r}
table(stromal_h$RNA_snn_res.0.5)

FeaturePlot(stromal_h, features = "Cd34")

FeaturePlot(stromal_h, features = "C1qtnf3")
```


```{r}
DimPlot(stromal_h, group.by="RNA_snn_res.0.4", label=T)
DimPlot(stromal_h, group.by="RNA_snn_res.0.5", label=T)


```
```{r}


Idents(stromal_h)<-'RNA_snn_res.0.5'
markers0.5<-FindAllMarkers(stromal_h, only.pos = T)

Idents(stromal_h)<-'RNA_snn_res.0.4'
markers0.4<-FindAllMarkers(stromal_h, only.pos = T, logfc.threshold = 1)


Idents(stromal_h)<-'RNA_snn_res.0.4'
levels(stromal_h)

current.sample.ids <- c("0" , "1",  "2" , "3" , "4" , "5",  "6",  "7" , "8" , "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19" ,"20")
new.sample.ids <- c("SL_Cd34_Ccl11" , "SL_Cthrc1_C1qtnf3",  "Endo" , "LL_Prg4_Clic5" , "Chron_Clu_Chad" , "Osteo_Bglap_Bglap2",  "Osteo_omd",  "SL_Cxcl9_Cxcl10" , "LS_Pi16_Cd248" , "Peri_Rgs5_Acta2",  "Muscle", "SL_Cthrc1_C1qtnf3", "RBCs", "Myeloid", "Contam", "neuron", "Myeloid", "Muscle", "Contam", "Lymphatics" ,"Osteo_Bglap_Bglap2")

stromal_h$named<-stromal_h$RNA_snn_res.0.4

stromal_h@meta.data[["named"]] <- plyr::mapvalues(x = stromal_h@meta.data[["named"]], from = current.sample.ids, to = new.sample.ids)


DimPlot(stromal_h, group.by="named", label=T)+NoLegend()

cols <- ArchR::paletteDiscrete(stromal_h@meta.data[, "named"])
DimPlot(stromal_h, group.by="named", label=F, cols=cols)

table(stromal_h$model)


Idents(stromal_h)<-'named'
markers_named<-FindAllMarkers(stromal_h, only.pos = T, logfc.threshold = 1)
write.csv(markers_named,  "/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/pipeline_all_samples/markers_named.csv")

```
```{r}
Idents(stromal_h)<-'named'
levels(stromal_h)
DotPlot(stromal_h, features = c("Cd34", "Ccl11", "Cthrc1", "C1qtnf3", "Pecam1", "Prg4", "Clic5", "Clu", "Chad", "Bglap", "Bglap2", "Omd", "Cxcl9", "Cxcl10", "Pi16", "Cd248", "Rgs5", "Acta2"), idents = c("SL_Cd34_Ccl11" ,     "SL_Cthrc1_C1qtnf3",  "Endo" ,              "LL_Prg4_Clic5"  ,    "Chron_Clu_Chad"   ,  "Osteo_Bglap_Bglap2", "Osteo_omd"      ,    "SL_Cxcl9_Cxcl10",   "LS_Pi16_Cd248",      "Peri_Rgs5_Acta2"  ))+RotatedAxis()



```


```{r}
aggv2_meta<-aggv2@meta.data
aggv2_meta<-subset(aggv2_meta, select=c("orig.ident" ,"inflammation"))
aggv2_meta_stromal<-aggv2_meta[rownames(aggv2_meta) %in% colnames(stromal_h),]
aggv2_meta_stromal$orig.ident<-NULL
stromal_h<-AddMetaData(stromal_h, aggv2_meta_stromal)


Idents(stromal_h)<-'inflammation'
DimPlot(stromal_h, split.by = "inflammation")


pt <- table(stromal_h$inflammation, stromal_h$named)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+RotatedAxis()
```
```{r}

Idents(stromal_h)<-'named'
levels(stromal_h)
DotPlot(stromal_h, features = c("Runx1", "Mmp14", "Igf1"), idents = c("SL_Cd34_Ccl11"   ,   "SL_Cthrc1_C1qtnf3",                 "LL_Prg4_Clic5"   ,    "SL_Cxcl9_Cxcl10" ,   "LS_Pi16_Cd248"  ,    "Peri_Rgs5_Acta2" ))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ scale_size(range = c(2, 8))


```
```{r}
stromal_h$sample_named_inflammation<-paste(stromal_h$orig.ident,stromal_h$named, stromal_h$inflammation, sep=".")

Idents(stromal_h)<-"sample_named_inflammation"
dotplot<-DotPlot(stromal_h, features = c("Runx1"))
dotplot_data<-dotplot[["data"]]
library(splitstackshape)
dotplot_data<-cSplit(dotplot_data, splitCols="id", sep=".")
dotplot_data<-dotplot_data[dotplot_data$id_2=="SL_Cthrc1_C1qtnf3",]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id_1, id_2, id_3))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"

dotplot_MMP14<-DotPlot(stromal_h, features = c("Mmp14"))
dotplotmmp14_data<-dotplot_MMP14[["data"]]
dotplotmmp14_data<-cSplit(dotplotmmp14_data, splitCols="id", sep=".")
dotplotmmp14_data<-dotplotmmp14_data[dotplotmmp14_data$id_2=="SL_Cthrc1_C1qtnf3",]
dotplotmmp14_data <- subset(dotplotmmp14_data, select = c(avg.exp.scaled, id_1, id_2, id_3))
names(dotplotmmp14_data)[names(dotplotmmp14_data)=="avg.exp.scaled"] <- "Mmp14"

dotplot_Cthrc1<-DotPlot(stromal_h, features = c("Cthrc1"))
dotplotCthrc1_data<-dotplot_Cthrc1[["data"]]
dotplotCthrc1_data<-cSplit(dotplotCthrc1_data, splitCols="id", sep=".")
dotplotCthrc1_data<-dotplotCthrc1_data[dotplotCthrc1_data$id_2=="SL_Cthrc1_C1qtnf3",]
dotplotCthrc1_data <- subset(dotplotCthrc1_data, select = c(avg.exp.scaled, id_1, id_2, id_3))
names(dotplotCthrc1_data)[names(dotplotCthrc1_data)=="avg.exp.scaled"] <- "Cthrc1"


dotplot_Igf1<-DotPlot(stromal_h, features = c("Igf1"))
dotplotIgf1_data<-dotplot_Igf1[["data"]]
dotplotIgf1_data<-cSplit(dotplotIgf1_data, splitCols="id", sep=".")
dotplotIgf1_data<-dotplotIgf1_data[dotplotIgf1_data$id_2=="SL_Cthrc1_C1qtnf3",]
dotplotIgf1_data <- subset(dotplotIgf1_data, select = c(avg.exp.scaled, id_1, id_2, id_3))
dotplotIgf1_data <- subset(dotplotIgf1_data, select = c(avg.exp.scaled, id_1, id_2, id_3))
names(dotplotIgf1_data)[names(dotplotIgf1_data)=="avg.exp.scaled"] <- "Igf1"

dotplot_data$Mmp14<-dotplotmmp14_data$Mmp14
dotplot_data$Cthrc1<-dotplotCthrc1_data$Cthrc1
dotplot_data$Igf1<-dotplotIgf1_data$Igf1



dotplot_data$samples<-paste(dotplot_data$id_1, dotplot_data$id_2, dotplot_data$id_3, sep=".")
 


library(reshape2)
dotplot_data_long=melt(dotplot_data, id.vars="samples")
dotplot_data_long$id<-dotplot_data_long$samples
dotplot_data_long<-cSplit(dotplot_data_long, splitCols = "samples", sep=".")


STIA<-rep("STIA",18)
AIA<-rep("AIA",15)
CIA<-rep("CIA",7)

model<-c(STIA, AIA, CIA)
dotplot_data$model<-model

ggplot(dotplot_data, aes(x = Mmp14, y = Runx1)) +
    geom_point(aes(shape = factor(model),color=factor(id_3)), size=3) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()

ggplot(dotplot_data, aes(x = Cthrc1, y = Runx1)) +
    geom_point(aes(shape = factor(model),color=factor(id_3)), size=3) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()


ggplot(dotplot_data, aes(x = Igf1, y = Runx1)) +
    geom_point(aes(shape = factor(model),color=factor(id_3)), size=3) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()


library(ggpubr)
ggscatter(dotplot_data, x = "Mmp14", y = "Runx1",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)

ggscatter(dotplot_data, x = "Cthrc1", y = "Runx1",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)

ggscatter(dotplot_data, x = "Igf1", y = "Runx1",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)
```


```{r}
stromal_h$named_inflammation<-paste(stromal_h$named, stromal_h$inflammation, sep="_")
Idents(stromal_h)<-'named_inflammation'
DotPlot(stromal_h, features = c("Runx1", "Mmp14", "Cthrc1", "Igf1", "Cbfb"), idents = c("SL_Cthrc1_C1qtnf3_rest", "SL_Cthrc1_C1qtnf3_resing", "SL_Cthrc1_C1qtnf3_resed", "SL_Cthrc1_C1qtnf3_persis", "SL_Cthrc1_C1qtnf3_initiation", "SL_Cthrc1_C1qtnf3_Eresing", "SL_Cthrc1_C1qtnf3_peak"))


```

```{r}
stromal_h$named_inflammation_model<-paste(stromal_h$named, stromal_h$inflammation,stromal_h$model, sep="_")
Idents(stromal_h)<-'named_inflammation_model'
dotplot<-DotPlot(stromal_h, features = c("Runx1", "Mmp14",  "Igf1"), idents = c("SL_Cthrc1_C1qtnf3_rest_STIA", "SL_Cthrc1_C1qtnf3_resing_STIA", "SL_Cthrc1_C1qtnf3_resed_STIA", "SL_Cthrc1_C1qtnf3_persis_STIA", "SL_Cthrc1_C1qtnf3_initiation_STIA",  "SL_Cthrc1_C1qtnf3_peak_STIA", "SL_Cthrc1_C1qtnf3_rest_CIA", "SL_Cthrc1_C1qtnf3_peak_CIA","SL_Cthrc1_C1qtnf3_rest_AIA", "SL_Cthrc1_C1qtnf3_resing_AIA", "SL_Cthrc1_C1qtnf3_resed_AIA", "SL_Cthrc1_C1qtnf3_persis_AIA", "SL_Cthrc1_C1qtnf3_initiation_AIA",  "SL_Cthrc1_C1qtnf3_peak_AIA", "SL_Cthrc1_C1qtnf3_Eresing_AIA"))+RotatedAxis()
dotplot<-dotplot$data
library(tidyverse)
dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

dotplot<-as.data.frame(dotplot)
dotplot <- dotplot[c("SL_Cthrc1_C1qtnf3_rest_STIA", "SL_Cthrc1_C1qtnf3_initiation_STIA", "SL_Cthrc1_C1qtnf3_peak_STIA", "SL_Cthrc1_C1qtnf3_resing_STIA", "SL_Cthrc1_C1qtnf3_resed_STIA", "SL_Cthrc1_C1qtnf3_persis_STIA", "SL_Cthrc1_C1qtnf3_rest_AIA", "SL_Cthrc1_C1qtnf3_peak_AIA", "SL_Cthrc1_C1qtnf3_Eresing_AIA", "SL_Cthrc1_C1qtnf3_resing_AIA", "SL_Cthrc1_C1qtnf3_resed_AIA", "SL_Cthrc1_C1qtnf3_rest_CIA", "SL_Cthrc1_C1qtnf3_peak_CIA")]


ncol(dotplot)

library(ComplexHeatmap)
split = c(rep(1, each = 6), rep(2, each = 5), rep(3, each = 2))


ha = HeatmapAnnotation(
    type = c(rep("STIA", 6), rep("AIA", 5), rep("CIA", 2)  ), 
    sample = 1:13,
    col = list(type = c("STIA" = "red", "AIA" = "green", "CIA" = "blue"))
)


Heatmap(dotplot, cluster_rows = F,cluster_columns = F, column_split =split,top_annotation=ha)

```

```{r}
stromal_h$orig.ident_named_inflammation_model<-paste(stromal_h$orig.ident, stromal_h$named, stromal_h$inflammation,stromal_h$model, sep=".")
Idents(stromal_h)<-'orig.ident_named_inflammation_model'
levels_complex<-as.data.frame(levels(stromal_h))
levels_complex$split<-levels_complex$`levels(stromal_h)`
library(splitstackshape)
levels_complex<-cSplit(levels_complex, "split", sep=".")
levels_complex_clus=levels_complex

levels_complex_clus_STIA<-levels_complex_clus[levels_complex_clus$split_4 == "STIA",]
levels_complex_clus_STIA_rest<-levels_complex_clus_STIA[levels_complex_clus_STIA$split_3 == "rest",]
levels_complex_clus_STIA_init<-levels_complex_clus_STIA[levels_complex_clus_STIA$split_3 == "initiation",]
levels_complex_clus_STIA_peak<-levels_complex_clus_STIA[levels_complex_clus_STIA$split_3 == "peak",]
levels_complex_clus_STIA_resing<-levels_complex_clus_STIA[levels_complex_clus_STIA$split_3 == "resing",]
levels_complex_clus_STIA_resed<-levels_complex_clus_STIA[levels_complex_clus_STIA$split_3 == "resed",]
levels_complex_clus_STIA_persis<-levels_complex_clus_STIA[levels_complex_clus_STIA$split_3 == "persis",]

levels_complex_clus_AIA<-levels_complex_clus[levels_complex_clus$split_4 == "AIA",]
levels_complex_clus_AIA_rest<-levels_complex_clus_AIA[levels_complex_clus_AIA$split_3 == "rest",]
levels_complex_clus_AIA_peak<-levels_complex_clus_AIA[levels_complex_clus_AIA$split_3 == "peak",]
levels_complex_clus_AIA_Eresing<-levels_complex_clus_AIA[levels_complex_clus_AIA$split_3 == "Eresing",]
levels_complex_clus_AIA_resing<-levels_complex_clus_AIA[levels_complex_clus_AIA$split_3 == "resing",]
levels_complex_clus_AIA_resed<-levels_complex_clus_AIA[levels_complex_clus_AIA$split_3 == "resed",]

levels_complex_clus_CIA<-levels_complex_clus[levels_complex_clus$split_4 == "CIA",]
levels_complex_clus_CIA_rest<-levels_complex_clus_CIA[levels_complex_clus_CIA$split_3 == "rest",]
levels_complex_clus_CIA_peak<-levels_complex_clus_CIA[levels_complex_clus_CIA$split_3 == "peak",]

rbind_list<-as.character(ls(pattern="levels_complex_clu_*"))
rbind_list[!rbind_list %in% c('levels_complex_clus', 'levels_complex_clus_AIA', 'levels_complex_clus_CIA')]
rbind_list<-as.data.frame(rbind_list)
levels_reorder<-rbind(levels_complex_clus_STIA_rest,levels_complex_clus_STIA_init,levels_complex_clus_STIA_peak,   levels_complex_clus_STIA_resing, levels_complex_clus_STIA_resed, levels_complex_clus_STIA_persis, levels_complex_clus_AIA_rest, levels_complex_clus_AIA_peak, levels_complex_clus_AIA_Eresing, levels_complex_clus_AIA_resing, levels_complex_clus_AIA_resed, levels_complex_clus_CIA_rest, levels_complex_clus_CIA_peak)
  
levels(stromal_h)<-levels_reorder$`levels(stromal_h)`

levels_complex_clus_cthrc1<-levels_complex_clus[levels_complex_clus$split_2 == "SL_Cthrc1_C1qtnf3",]

dotplot<-DotPlot(stromal_h, features = c("Runx1", "Mmp14", "Cthrc1", "Igf1"), idents = levels_complex_clus_cthrc1$`levels(stromal_h)`)+RotatedAxis()

dotplot<-dotplot$data
library(tidyverse)
dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

dotplot<-as.data.frame(dotplot)

colnames(dotplot)

dotplot <- dotplot[c("CD45N.SL_Cthrc1_C1qtnf3.rest.STIA",  "S8_CONTROLnegB_r.SL_Cthrc1_C1qtnf3.rest.STIA"   ,  "S8_ControlnegC.SL_Cthrc1_C1qtnf3.rest.STIA"  ,   
                                "S2_DAY1negA_cmAB2.SL_Cthrc1_C1qtnf3.initiation.STIA",   "S4_DAY1negB_cmAB2.SL_Cthrc1_C1qtnf3.initiation.STIA",   "S6_DAY1negC_cmAB2.SL_Cthrc1_C1qtnf3.initiation.STIA",      "S4_Day8negB.SL_Cthrc1_C1qtnf3.peak.STIA"     , "S6_DAY8negC.SL_Cthrc1_C1qtnf3.peak.STIA", "S2_DAY8negA.SL_Cthrc1_C1qtnf3.peak.STIA"       ,   "S2_Day15negB.SL_Cthrc1_C1qtnf3.resing.STIA" ,  "S4_Day15negC.SL_Cthrc1_C1qtnf3.resing.STIA"        ,"S8_DAY15negA.SL_Cthrc1_C1qtnf3.resing.STIA"    ,
 "S2_Day22negA.SL_Cthrc1_C1qtnf3.resed.STIA"       ,    "S4_Day22negB.SL_Cthrc1_C1qtnf3.resed.STIA" ,         
 "S6_Day22negC.SL_Cthrc1_C1qtnf3.resed.STIA"        ,           
   "S6_Day28negA.SL_Cthrc1_C1qtnf3.persis.STIA"  ,       
 "S8_Day28negB.SL_Cthrc1_C1qtnf3.persis.STIA"         , "S2_Day28negC.SL_Cthrc1_C1qtnf3.persis.STIA"   ,
 "AIA_rest1_cd45neg.SL_Cthrc1_C1qtnf3.rest.AIA"    , "Day0neg_cr_V5.SL_Cthrc1_C1qtnf3.rest.AIA"    ,
 "rest_neg_renamed.SL_Cthrc1_C1qtnf3.rest.AIA" ,"AIA_infla1_cd45neg.SL_Cthrc1_C1qtnf3.peak.AIA"  ,     "AIA_infla2_cd45neg.SL_Cthrc1_C1qtnf3.peak.AIA" ,  "Day2_neg3.SL_Cthrc1_C1qtnf3.peak.AIA"  ,
 "Day4_neg1.SL_Cthrc1_C1qtnf3.Eresing.AIA"     ,       
 "Day4_neg2.SL_Cthrc1_C1qtnf3.Eresing.AIA"    ,         "Day4_neg3.SL_Cthrc1_C1qtnf3.Eresing.AIA"     , "Day7_neg1.SL_Cthrc1_C1qtnf3.resing.AIA"      ,        "Day7_neg2.SL_Cthrc1_C1qtnf3.resing.AIA"      ,       
 "Day7_neg3.SL_Cthrc1_C1qtnf3.resing.AIA"       ,"D14_1_neg_renamed.SL_Cthrc1_C1qtnf3.resed.AIA"  ,    
 "D14_2_neg_renamed.SL_Cthrc1_C1qtnf3.resed.AIA"    ,   "D14_3_neg_renamed.SL_Cthrc1_C1qtnf3.resed.AIA"   ,       "Con1_CD45neg.SL_Cthrc1_C1qtnf3.rest.CIA"     ,       
 "Con2_CD45neg.SL_Cthrc1_C1qtnf3.rest.CIA"       ,      "Con3_CD45neg.SL_Cthrc1_C1qtnf3.rest.CIA"     ,  "Con4_CD45neg.SL_Cthrc1_C1qtnf3.rest.CIA"         ,     
 "Infla1_CD45neg.SL_Cthrc1_C1qtnf3.peak.CIA"      ,     "Infla2_CD45neg.SL_Cthrc1_C1qtnf3.peak.CIA"   ,       
     "Infla3_CD45neg.SL_Cthrc1_C1qtnf3.peak.CIA" )]


ncol(dotplot)

library(ComplexHeatmap)
split = c(rep(1, each = 18), rep(2, each = 15), rep(3, each = 7))


ha = HeatmapAnnotation(
    type = c(rep("STIA", 18), rep("AIA", 15), rep("CIA", 7)  ), 
    condition = c(rep("rest", 3), rep("initiation", 3), rep("peak", 3), rep("resing", 3), rep("resed", 3), rep("persis", 3), rep("rest", 3), rep("peak", 3), rep("Eresing", 3), rep("resing", 3), rep("resed", 3),  rep("rest", 4), rep("peak", 3) ),
    col = list(type = c("STIA" = "red", "AIA" = "green", "CIA" = "blue"))
)


Heatmap(dotplot, cluster_rows = F,cluster_columns = F, column_split =split, top_annotation=ha)
```



```{r}
DefaultAssay(object=stromal) <- "RNA"

#stromal@meta.data$orig.ident <- as.factor(stromal@meta.data$orig.ident)
stromal_h_orig <- RunHarmony.Seurat_CM(stromal, "orig.ident", reduction = "pca", assay.use = "integrated", plot_convergence = TRUE, max.iter.harmony = 30)

#not working, object saved and ran on sams setup
saveRDS(stromal, "/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/pipeline_all_samples/stromal.rds")
stromal_h_orig.ident<-readRDS("/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/pipeline_all_samples/stromal.H.rds")

DimPlot(stromal_h_orig.ident, group.by = "orig.ident", reduction="HarmonyUMAP")
DimPlot(stromal_h_orig.ident, group.by = "model", reduction="HarmonyUMAP")
DefaultAssay(stromal_h_orig.ident)<-'RNA'
FeaturePlot(stromal_h_orig.ident, features=c("Runx1", "Cthrc1", "Igf1", "Mmp14"), reduction="HarmonyUMAP")

```
```{r}
stromal_h_orig.ident<-FindNeighbors(stromal_h_orig.ident, reduction = "harmony",dims = 1:30)
stromal_h_orig.ident<-FindClusters(stromal_h_orig.ident, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4), graph.name = "integrated_snn")
stromal_h_orig.ident<-FindClusters(stromal_h_orig.ident, resolution = c(0.5, 0.6, 0.7), graph.name = "integrated_snn")

DimPlot(stromal_h_orig.ident, group.by="integrated_snn_res.0.5", label=T, reduction = "HarmonyUMAP")
FeaturePlot(stromal_h_orig.ident, features="Cd34", reduction = "HarmonyUMAP")


Idents(stromal_h_orig.ident)<-'integrated_snn_res.0.5'
all_markers<-FindAllMarkers(stromal_h_orig.ident, only.pos = T, logfc.threshold = 1)

```

```{r}

current.sample.ids <- c("0" , "1",  "2" , "3" , "4" , "5",  "6",  "7" , "8" , "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19" ,"20")
new.sample.ids <- c("SL_Cxcl14_Cxcl5" , "SL_Pi16_Cd34",  "Endo" , "SL_C1qtnf3_Cthrc1" , "Chon_Clu_Chad" , "LL_Prg4_Clic5",  "Osteo_Bglap_Bglap2",  "Osteo_Omd" , "SL_Pi16_Cd34" , "SL_Smoc2_Ccl11",  "Peri_Rgs5_Acta2", "Mucscle", "Profil", "Chond_Chodl", "Macs", "RBCs", "Muscle2", "GLial", "Chon_Clu_Chad", "Lymphatics" ,"Endo")

stromal_h_orig.ident$named<-stromal_h_orig.ident$integrated_snn_res.0.5

stromal_h_orig.ident@meta.data[["named"]] <- plyr::mapvalues(x = stromal_h_orig.ident@meta.data[["named"]], from = current.sample.ids, to = new.sample.ids)


DimPlot(stromal_h_orig.ident, group.by="named", label=T)+NoLegend()

cols <- ArchR::paletteDiscrete(stromal_h_orig.ident@meta.data[, "named"])
DimPlot(stromal_h_orig.ident, group.by="named", label=F, cols=cols, reduction = "HarmonyUMAP")
DimPlot(stromal_h_orig.ident, group.by="model", label=F,  reduction = "HarmonyUMAP")
DimPlot(stromal_h_orig.ident, group.by="orig.ident", label=F, reduction = "HarmonyUMAP")+NoLegend()


Idents(stromal_h_orig.ident)<-'named'
levels(stromal_h_orig.ident)
DotPlot(stromal_h_orig.ident, features = c("Runx1", "Mmp14", "Cthrc1", "Igf1"), idents = c("SL_Cxcl14_Cxcl5" ,   "SL_Pi16_Cd34",              "SL_C1qtnf3_Cthrc1"  , "LL_Prg4_Clic5"    ,        "SL_Smoc2_Ccl11" , "Peri_Rgs5_Acta2"    ))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+RotatedAxis()+ scale_size(range = c(2, 8))



Idents(stromal_h_orig.ident)<-'named'
markers_named_new<-FindAllMarkers(stromal_h_orig.ident, only.pos = T, logfc.threshold = 1)
write.csv(markers_named,  "/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/pipeline_all_samples/markers_named_new.csv")




```
```{r}


pt <- table(stromal_h_orig.ident$orig.ident , stromal_h_orig.ident$named)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank()) +RotatedAxis()+NoLegend()


pt <- table(stromal_h_orig.ident$model , stromal_h_orig.ident$named)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank()) +RotatedAxis()


agg_meta<-aggv2@meta.data
agg_meta <- subset(agg_meta, select = c(orig.ident, inflammation))
agg_meta<-agg_meta[rownames(agg_meta) %in% colnames(stromal_h_orig.ident),]

stromal_h_orig.ident$inflammation<-agg_meta$inflammation

pt <- table(stromal_h_orig.ident$inflammation , stromal_h_orig.ident$named)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank()) +RotatedAxis()




```
```{r}

stromal_h_orig.ident$named_inflammation<-paste(stromal_h_orig.ident$named, stromal_h_orig.ident$inflammation, sep=".")

library(splitstackshape)
my_patterns=c("SL_", "LL_")
data <-table(stromal_h_orig.ident$named_inflammation) %>% as.data.frame()%>% cSplit(splitCols = "Var1", sep=".") %>% filter(grepl(paste(my_patterns, collapse='|'), Var1_1))

cols <- ArchR::paletteDiscrete(data$Var1_1) %>% as.data.frame()


data %>% ggplot(aes(y=Freq, x=Var1_1, fill=Var1_1))  + 
        geom_bar(stat = 'identity', position = position_dodge(), aes(group = Var1_1),width = 1) + 
        coord_flip()+ 
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
  theme(strip.background = element_rect(fill="white", size=1, color="white"))+facet_wrap("Var1_2", nrow=2)+RotatedAxis()


stromal_h_orig.ident$named_model<-paste(stromal_h_orig.ident$named, stromal_h_orig.ident$model, sep=".")

library(splitstackshape)
my_patterns=c("SL_", "LL_")
data <-table(stromal_h_orig.ident$named_model) %>% as.data.frame()%>% cSplit(splitCols = "Var1", sep=".") %>% filter(grepl(paste(my_patterns, collapse='|'), Var1_1))

cols <- ArchR::paletteDiscrete(stromal_h_orig.ident@meta.data[, "named"]) %>% as.data.frame()
cols <- cols[rownames(cols) %in% data$Var1_1,]

data %>% ggplot(aes(y=Freq, x=Var1_1, fill=Var1_1))  + 
        geom_bar(stat = 'identity', position = position_dodge(), aes(group = Var1_1),width = 1) + 
        coord_flip()+ 
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
  scale_x_discrete(expand = expansion(add = c(0, 0)))+theme_ArchR()+scale_fill_manual(values =c("#FEE500", "#89288F", "#D51F26", "#272E6A", "#D8A767"))+
  theme(strip.background = element_rect(fill="white", size=1, color="white"))+facet_wrap("Var1_2", nrow=1, ncol=3)+RotatedAxis()


data <-table(stromal_h_orig.ident$model) %>% as.data.frame()

#cols <- ArchR::paletteDiscrete(stromal_h_orig.ident@meta.data[, "named"]) %>% as.data.frame()
#cols <- cols[rownames(cols) %in% data$Var1_1,]

data %>% ggplot(aes(y=Freq, x=Var1, fill=Var1))  + 
        geom_bar(stat = 'identity', position = position_dodge(), aes(group = Var1),width = 1) + 
        coord_flip()+ 
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
  scale_x_discrete(expand = expansion(add = c(0, 0)))+theme_ArchR()+scale_fill_manual(values =c("darkred", "darkgreen", "darkblue"))+
  theme(strip.background = element_rect(fill="white", size=1, color="white"))+RotatedAxis()





```

```{r}
DimPlot(stromal_h_orig.ident, group.by = "orig.ident", reduction = "HarmonyUMAP")
```

```{r}
FeaturePlot(stromal_h_orig.ident, features = c("Col1a1"), reduction="HarmonyUMAP")+scale_color_viridis()
FeaturePlot(stromal_h_orig.ident, features = c("Prg4"), reduction="HarmonyUMAP")+scale_color_viridis()
FeaturePlot(stromal_h_orig.ident, features = c("Pdgfra"), reduction="HarmonyUMAP")+scale_color_viridis()
FeaturePlot(stromal_h_orig.ident, features = c("Thy1"), reduction="HarmonyUMAP")+scale_color_viridis()


FeaturePlotSingle<- function(obj, feature, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- levels(obj@meta.data[, metadata_column])
  
  # the minimal and maximal of the value to make the legend scale the same. 
  minimal<- min(obj[['RNA']]@data[feature, ])
  maximal<- max(obj[['RNA']]@data[feature, ])
  ps<- list()
  for (group in groups) {
    subset_indx<- obj@meta.data[, metadata_column] == group
    subset_cells<- all_cells[subset_indx]
    p<- FeaturePlot(obj, features = feature, cells= subset_cells, ...) +
      scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
      ggtitle(group) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[group]]<- p
  }
  
  
  return(ps)
}


p_list<- FeaturePlotSingle(stromal_h_orig.ident, feature= "Col1a1", metadata_column = "model", pt.size = 0.05, order =TRUE)

layout1<-"
ABC
"

wrap_plots(p_list ,guides = 'collect', design = layout1)

```







```{r}
stromal_h_orig.ident$sample_named_inflammation<-paste(stromal_h_orig.ident$orig.ident,stromal_h_orig.ident$named, stromal_h_orig.ident$inflammation, sep=".")



```


```{r}
stromal_h_orig.ident$named_inflammation_model

Idents(stromal_h_orig.ident) <- 'inflammation'

DotPlot(stromal_h_orig.ident, features=c("Runx1", "Cxcl5", "Ccl20", "Igf1"))


```




```{r}
stromal_h_orig.ident$named_inflammation_model<-paste(stromal_h_orig.ident$named, stromal_h_orig.ident$inflammation,stromal_h_orig.ident$model, sep="_")
Idents(stromal_h_orig.ident)<-'named_inflammation_model'
dotplot<-DotPlot(stromal_h_orig.ident, features = c("Runx1", "Mmp14", "Cthrc1", "Igf1"), idents = c("SL_C1qtnf3_Cthrc1_rest_STIA", "SL_C1qtnf3_Cthrc1_resing_STIA", "SL_C1qtnf3_Cthrc1_resed_STIA", "SL_C1qtnf3_Cthrc1_persis_STIA", "SL_C1qtnf3_Cthrc1_initiation_STIA",  "SL_C1qtnf3_Cthrc1_peak_STIA", "SL_C1qtnf3_Cthrc1_rest_CIA", "SL_C1qtnf3_Cthrc1_peak_CIA","SL_C1qtnf3_Cthrc1_rest_AIA", "SL_C1qtnf3_Cthrc1_resing_AIA", "SL_C1qtnf3_Cthrc1_resed_AIA", "SL_C1qtnf3_Cthrc1_persis_AIA", "SL_C1qtnf3_Cthrc1_initiation_AIA",  "SL_C1qtnf3_Cthrc1_peak_AIA", "SL_C1qtnf3_Cthrc1_Eresing_AIA"))+RotatedAxis()
dotplot<-dotplot$data
library(tidyverse)
dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

dotplot<-as.data.frame(dotplot)
dotplot <- dotplot[c("SL_C1qtnf3_Cthrc1_rest_STIA", "SL_C1qtnf3_Cthrc1_initiation_STIA", "SL_C1qtnf3_Cthrc1_peak_STIA", "SL_C1qtnf3_Cthrc1_resing_STIA", "SL_C1qtnf3_Cthrc1_resed_STIA", "SL_C1qtnf3_Cthrc1_persis_STIA", "SL_C1qtnf3_Cthrc1_rest_AIA", "SL_C1qtnf3_Cthrc1_peak_AIA", "SL_C1qtnf3_Cthrc1_Eresing_AIA", "SL_C1qtnf3_Cthrc1_resing_AIA", "SL_C1qtnf3_Cthrc1_resed_AIA", "SL_C1qtnf3_Cthrc1_rest_CIA", "SL_C1qtnf3_Cthrc1_peak_CIA")]


ncol(dotplot)

library(ComplexHeatmap)
split = c(rep(1, each = 6), rep(2, each = 5), rep(3, each = 2))


ha = HeatmapAnnotation(
    type = c(rep("STIA", 6), rep("AIA", 5), rep("CIA", 2)  ), 
    sample = 1:13,
    col = list(type = c("STIA" = "red", "AIA" = "green", "CIA" = "blue"))
)


Heatmap(dotplot, cluster_rows = F,cluster_columns = F, column_split =split,top_annotation=ha)
```



```{r}

fibs_only <- subset(stromal_h_orig.ident, idents=levels(stromal_h_orig.ident)[grep("SL_|LL_", levels(stromal_h_orig.ident))])
fibs_only <- fibs_only %>% ScaleData()


fibs_only$cluster_inflammation_model_smaple<-paste(fibs_only$named, fibs_only$orig.ident, fibs_only$inflammation,fibs_only$model, sep=".")
Idents(fibs_only)<-'cluster_inflammation_model_smaple'
levels(fibs_only)




```

```{r}
fibs_only$inflammation_model_smaple<-paste( fibs_only$orig.ident, fibs_only$inflammation,fibs_only$model, sep=".")


fibs_only$inflammation_model<-paste( fibs_only$model, fibs_only$inflammation, sep=".")

Idents(fibs_only)<-'inflammation_model'


dotplot <- DotPlot(fibs_only, features= c("Cxcl5", "Runx1", "Igf1", "Mmp14"), idents = levels(fibs_only)[grep("STIA|CIA", levels(fibs_only))])

dotplot<-dotplot$data


dotplot<-dotplot %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)


ha = HeatmapAnnotation(
    type = c(rep("STIA", 6), rep("CIA", 2)  ), 
    col = list(type = c("STIA" = "red",  "CIA" = "blue"))
)


Heatmap(dotplot, cluster_rows = F,cluster_columns = F, column_split =split,top_annotation=ha)

dotplot <- dotplot[,c(1,6,2,3,4,5,7,8)]


dotplot <- dotplot[,c(1,3,4,6,5,2,7,8)]

split = c(rep(1, each = 6), rep(3, each = 2))

Heatmap(dotplot, border=T, cluster_columns = F, column_split =split,top_annotation=ha)




```







