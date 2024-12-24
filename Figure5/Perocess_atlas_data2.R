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

load("/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/pipeline_all_samples/all_analysis.RData")


Idents(aggr)<-'orig.ident'
current.sample.ids <- levels(aggr)
new.sample.ids <- c("pos" ,      "pos" ,   "neg" ,   "neg"  ,       "pos" ,      
  "neg"     ,   "pos"    ,    "neg"  ,      "pos"  ,      "neg",       
 "pos"     ,   "neg"  ,      "pos"    ,    "neg"    ,    "pos"  ,     
 "neg"   ,     "pos"    ,    "neg"    ,    "pos"   ,     "neg" ,      
 "pos"    ,     "pos"   ,      "neg",  "neg"   ,      "pos"   ,     
 "neg"  , "pos"     ,    "pos"  ,       "neg" ,  "neg",        
 "pos"   ,   "pos"   ,     "neg"   ,   "neg"  ,      "neg" ,             
 "pos"       ,        "neg",  "pos" , "neg",  "pos" ,
 "neg" ,  "pos" ,  "neg"  , "pos" ,  "neg" , 
 "pos"  , "neg",   "pos" ,  "neg" ,   "pos",   
 "neg"   ,    "pos",       "neg"  ,         "pos" ,          "neg" ,         
 "neg"       ,    "neg"     ,      "pos"    ,       "pos"   ,        "pos"  ,        
 "neg"     ,      "neg"        ,   "neg"     ,      "pos"      ,      "pos" ,          
 "pos"     ,       "pos"  ,      "pos"      ,  "pos"   ,   "pos",
 "pos"   ,     "pos"  ,    "pos"  ,    "pos"  ,    "neg" ,      
 "neg"     ,   "neg"  ,      "neg"  ,    "neg" ,     "neg"  ,     
 "neg"     )

aggr$CD45<-aggr$orig.ident

aggr@meta.data[["CD45"]] <- plyr::mapvalues(x = aggr@meta.data[["CD45"]], from = current.sample.ids, to = new.sample.ids)

aggr$orig.ident_CD45<-paste(aggr$orig.ident, aggr$CD45, sep="_")
table(aggr$orig.ident_CD45)

Idents(aggr)<-'CD45'
levels(aggr)
DefaultAssay(aggr)<-'RNA'
cd45pos<-subset(x = aggr, idents = "pos")
cd45neg<-subset(x = aggr, idents = "neg")

all_split=c(cd45pos, cd45neg)

features <- SelectIntegrationFeatures(object.list = all_split, nfeatures = 20000)


for (i in 1:length(all_split)) {
all_split[[i]] =NormalizeData(object = all_split[[i]]);
all_split[[i]] =ScaleData(object = all_split[[i]]);
all_split[[i]] =FindVariableFeatures(object = all_split[[i]]);
all_split[[i]] =RunPCA(object = all_split[[i]], verbose = FALSE);
all_split[[i]] =RunUMAP(object = all_split[[i]], dims=1:30)
}

rm(list=ls()[! ls() %in% c("all_split", "features")])
gc()


anchors <- FindIntegrationAnchors(object.list = all_split, anchor.features = features, reduction = "rpca",dims = 1:30)

aggv2 <- IntegrateData(anchorset = anchors, dims = 1:30)
aggv2 <- ScaleData(aggv2, verbose = FALSE)
aggv2 <- RunPCA(aggv2, npcs = 30, verbose = FALSE)
aggv2 <- RunUMAP(aggv2, reduction = "pca", dims = 1:30)
aggv2 <- FindNeighbors(aggv2, reduction = "pca", dims = 1:30)


DefaultAssay(aggv2)<-"RNA"

DimPlot(aggv2, group.by="orig.ident")+NoLegend()
FeaturePlot(aggv2, features="Ptprc")
FeaturePlot(aggv2, features="Col1a1")
FeaturePlot(aggv2, features="Cd68")
FeaturePlot(aggv2, features="Acta2")
FeaturePlot(aggv2, features="Runx2")
FeaturePlot(aggv2, features="Pecam1")
dev.off()

aggv2 <- FindClusters(aggv2, resolution = c(0.01, 0.05, 0.1, 0.2), graph.name = 'integrated_snn')
Idents(aggv2)<-'integrated_snn_res.0.01'  
res0.01markers<-FindAllMarkers(aggv2, only.pos = T, logfc.threshold = 1)

Idents(aggv2)<-'integrated_snn_res.0.05'
res0.05markers<-FindAllMarkers(aggv2, only.pos = T, , logfc.threshold = 1)

Idents(aggv2)<-'integrated_snn_res.0.1'
res0.1markers<-FindAllMarkers(aggv2, only.pos = T, , logfc.threshold = 1)


Idents(aggv2)<-'integrated_snn_res.0.2'
res0.2markers<-FindAllMarkers(aggv2, only.pos = T)


rm(list=ls()[! ls() %in% c("aggv2", "res0.01markers", "res0.05markers", "res0.1markers", "res0.02markers")])
gc()


save.image("/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/pipeline_all_samples/all_analysis_selectedfeatures_20k.RData")
