

```{r}




dirs <- list.dirs("/rds/projects/c/croftap-sgcar/single_cell_feb_2024/count", recursive = FALSE)
dirs <- paste0(dirs, "/outs/filtered_feature_bc_matrix")

dirs2 <- list.dirs("/rds/projects/c/croftap-sgcar/count", recursive = FALSE)
dirs2 <- dirs2[4:6]

sample <- c("SGH1", "SGH2", "SGH3")
dirs2 <- paste0(dirs2, "/outs/per_sample_outs/", sample, "/count/filtered_feature_bc_matrix")

dirs_all <- c(dirs, dirs2)

rm(list=ls()[! ls() %in% c("dirs_all")])

data.10x = list()

for (i in 1:length(dirs_all)) {
    data.10x[[i]] <- Read10X(data.dir = dirs_all[[i]])
}

samples<-c("P1", "P2", "P3","SGH1", "SGH2", "SGH3" )

names(data.10x) <- samples
  

scrna.list = list()

for (i in 1:length(data.10x)) {
    scrna.list[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=3, min.features=500, project=samples[i]);
    scrna.list[[i]][["percent.mt"]] = PercentageFeatureSet(object=scrna.list[[i]], pattern = "^mt-");
    scrna.list[[i]] <- subset( scrna.list[[i]], subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 10);
    scrna.list[[i]] =NormalizeData(object = scrna.list[[i]]);
    scrna.list[[i]] =ScaleData(object = scrna.list[[i]]);
    scrna.list[[i]] =FindVariableFeatures(object = scrna.list[[i]]);
    scrna.list[[i]] =RunPCA(object = scrna.list[[i]], verbose = FALSE)
    }


scrna.list[[1]][["percent.mt"]] = PercentageFeatureSet(object=scrna.list[[1]], pattern = "^mt-")

max(scrna.list[[1]]$percent.mt)

names(scrna.list) <- samples

anchors <- FindIntegrationAnchors(object.list = scrna.list, reduction = "rpca",   dims = 1:50)
aggr <- IntegrateData(anchorset = anchors, dims = 1:50)
aggr <- FindVariableFeatures(aggr)
aggr <- ScaleData(aggr, verbose = FALSE)
aggr <- RunPCA(aggr, verbose = FALSE)
aggr <- RunUMAP(aggr, dims = 1:50)

aggr <- FindNeighbors(aggr, dims = 1:20)
aggr <- FindClusters(aggr, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3), graph.name = 'integrated_snn')
Idents(aggr)<-'integrated_snn_res.0.01'  
res0.01markers<-FindAllMarkers(aggr, only.pos = T)


Idents(aggr)<-'integrated_snn_res.0.05'
res0.05markers<-FindAllMarkers(aggr, only.pos = T)


Idents(aggr)<-'integrated_snn_res.0.1'
res0.1markers<-FindAllMarkers(aggr, only.pos = T)


Idents(aggr)<-'integrated_snn_res.0.2'
res0.2markers<-FindAllMarkers(aggr, only.pos = T)



```

```{r}
DimPlot(aggr, group.by = 'integrated_snn_res.0.05')
```





```{r}
DimPlot(aggr, group.by="orig.ident")
DefaultAssay(aggr) <- 'RNA'
FeaturePlot(aggr, features = "Pdgfra")
FeaturePlot(aggr, features = "Thy1")
FeaturePlot(aggr, features = "Col1a1")
FeaturePlot(aggr, features = "Runx1")

```

```{r}

aggr$condition <- aggr@meta.data[["orig.ident"]]
Idents(aggr) <- 'condition'
levels(aggr)
current.sample.ids <- levels(aggr)
new.sample.ids <- c("IF" ,  "IF",  "IF" ,  "HC", "HC", "HC")

aggr@meta.data[["condition"]] <- plyr::mapvalues(x = aggr@meta.data[["condition"]], from = current.sample.ids, to = new.sample.ids)



aggr$clus_condition <- paste(aggr$integrated_snn_res.0.05, aggr$condition, sep="_")

Idents(aggr) <- 'clus_condition'
levels(aggr)

VlnPlot(aggr, idents=c("0_IF", "0_HC"), features="Runx1")


```





```{r}
library(harmony)

aggr <- RunHarmony(aggr, group.by.vars = "orig.ident")
aggr <- RunUMAP(aggr, reduction = "harmony", dims = 1:50)
DimPlot(aggr, group.by="orig.ident")

DefaultAssay(aggr) <- 'RNA'
FeaturePlot(aggr, features = "Pdgfra")
FeaturePlot(aggr, features = "Thy1")
FeaturePlot(aggr, features = "Cd248")
FeaturePlot(aggr, features = "Col1a1")

aggr <- FindNeighbors(aggr, dims = 1:30, reduction = "harmony")
aggr <- FindClusters(aggr, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3), graph.name = 'integrated_snn')
DimPlot(aggr, group.by="integrated_snn_res.0.05")
DimPlot(aggr, group.by="orig.ident")

Idents(aggr)<-'integrated_snn_res.0.01'  
res0.01markers<-FindAllMarkers(aggr, only.pos = T, logfc.threshold = 0.75)


p1 <- FeaturePlot(aggr, features = "Pdgfra")+scale_color_viridis()
p2 <- FeaturePlot(aggr, features = "Col1a1")+scale_color_viridis()
p3 <- FeaturePlot(aggr, features = "Cd248")+scale_color_viridis()

FeaturePlot(aggr, features = "Pdpn")+scale_color_viridis()

p4 <- FeaturePlot(aggr, features = "Thy1")+scale_color_viridis()
FeaturePlot(aggr, features = "Runx1")+scale_color_viridis()

DimPlot(aggr, group.by = "orig.ident")

plot_grid(p1,p2,p3,p4, ncol=2, nrow=2)
```

```{r}

aggr$condition <- aggr@meta.data[["orig.ident"]]
Idents(aggr) <- 'condition'
levels(aggr)
current.sample.ids <- levels(aggr)
new.sample.ids <- c("IF" ,  "IF",  "IF" ,  "HC", "HC", "HC")

aggr@meta.data[["condition"]] <- plyr::mapvalues(x = aggr@meta.data[["condition"]], from = current.sample.ids, to = new.sample.ids)


data <-table(aggr$condition) %>% as.data.frame()

cols <- ArchR::paletteDiscrete(data$Var1) %>% as.data.frame()

data %>% ggplot(aes(y=Freq, x=Var1, fill=Var1))  + 
        geom_bar(stat = 'identity', position = position_dodge(), aes(group = Var1),width = 1) + 
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
  theme(strip.background = element_rect(fill="white", size=1, color="white"))
```

```{r}

aggr$clus_condition <- paste(aggr$integrated_snn_res.0.01, aggr$condition, sep="_")
Idents(aggr) <- 'clus_condition'

#pull GEX from seurat obj (mine is all)

meta <- aggr@meta.data

gex <- aggr@assays[["RNA"]]@data %>% t() %>%  as.data.frame()

meta$Runx1 <- gex$Runx1

#get colors (optional)
cols <- as.data.frame(ArchR::paletteDiscrete(aggr@meta.data[, "clus_condition"]))
colnames(cols) <- "colors"

library(tidyverse)

#plot. No outliers and no whiskers
meta  %>% filter(clus_condition %in% c("0_HC", "0_IF")) %>%
ggplot(aes(x=Runx1, y=clus_condition, fill=clus_condition))+
    geom_boxplot(outlier.shape = NA, coef = 0)+theme_minimal()+theme(axis.ticks = element_blank())+scale_fill_manual(values=cols$colors)+coord_flip()+
  theme_cowplot(12)+RotatedAxis()+ theme(legend.position="none")


meta  %>% filter(clus_condition %in% c("0_HC", "0_IF")) %>%
ggplot(aes(x=Runx1, y=clus_condition, fill=clus_condition))+
    geom_boxplot(outlier.shape = 1, coef = 1)+theme_minimal()+theme(axis.ticks = element_blank(), axis.text.x=element_blank())+scale_fill_manual(values=c("#F47D2B", "#C06CAB"))+coord_flip()+
  theme_cowplot(12)+RotatedAxis()+ theme(legend.position="none")






```


```{r}
grn <- read.csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/df_grn2_new_FINAL.csv")

grnRUNX1 <- grn %>% filter(tf=="RUNX1")


grnRUNX1$Target <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(grnRUNX1$gene), perl=TRUE)

DefaultAssay(aggr) <- "RNA"
aggr <- AddModuleScore(aggr, features = list(grnRUNX1$Target), name="RUNX1grn")


meta <- aggr@meta.data

#get colors (optional)
cols <- as.data.frame(ArchR::paletteDiscrete(aggr@meta.data[, "clus_condition"]))
colnames(cols) <- "colors"

library(tidyverse)

#plot. No outliers and no whiskers
meta  %>% filter(clus_condition %in% c("0_HC", "0_IF")) %>%
ggplot(aes(x=RUNX1grn1, y=clus_condition, fill=clus_condition))+
    geom_boxplot(outlier.shape = NA, coef = 0)+theme_minimal()+theme(axis.ticks = element_blank())+scale_fill_manual(values=cols$colors)+coord_flip()+
  theme_cowplot(12)+RotatedAxis()+ theme(legend.position="none")


meta  %>% filter(clus_condition %in% c("0_HC", "0_IF")) %>%
ggplot(aes(x=RUNX1grn1, y=clus_condition, fill=clus_condition))+
    geom_boxplot(outlier.shape = 1, coef = 1)+theme_minimal()+theme(axis.ticks = element_blank(), axis.text.x=element_blank())+scale_fill_manual(values=c("#F47D2B", "#C06CAB"))+coord_flip()+
  theme_cowplot(12)+RotatedAxis()+ theme(legend.position="none")


```


```{r}

Idents(aggr) <- 'integrated_snn_res.0.01'
fibs <- subset(aggr, idents="0")
fibs <- fibs %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>% RunUMAP(dims=1:30)

DimPlot(fibs, group.by="orig.ident")
FeaturePlot(fibs, features = "Runx1")
FeaturePlot(fibs, features = "Cd248")


library(harmony)

fibs <- RunHarmony(fibs, group.by.vars = "orig.ident")
fibs <- RunUMAP(fibs, reduction = "harmony", dims = 1:50)
DimPlot(fibs, group.by="orig.ident")

Idents(fibs) <- 'condition'
IF_vs_HC <- FindMarkers(fibs, ident.1 = "IF", ident.2 = "HC", logfc.threshold = 0.75)


```

