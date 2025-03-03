


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







