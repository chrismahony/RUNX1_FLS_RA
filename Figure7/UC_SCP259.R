
```{r}


UC <- Read10X(data.dir="/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/SCP259/data", gene.column = 1)

all.meta2 <- all.meta2[-1,]

rownames(all.meta2) <- all.meta2$NAME
all.meta2$NAME <- NULL

all.meta2_f <- all.meta2[rownames(all.meta2) %in% colnames(UC),]

UC <- UC %>% CreateSeuratObject(meta.data = all.meta2_f) %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims=1:30)

UC$Cluster
DimPlot(UC, group.by = "Cluster")


FeaturePlot(UC, features = "PDGFRA")



Idents(UC) <- 'Cluster'

UC_fibs <- subset(UC, idents=levels(UC)[-c(2,7,8,9,11,12)])

UC_fibs <- UC_fibs %>% ScaleData()

UC_fibs <- UC_fibs %>% RunUMAP(dims=1:30)

DimPlot(UC_fibs, group.by = "Cluster")
FeaturePlot(UC_fibs, features = "THY1")


Idents(UC_fibs) <- 'Health'

VlnPlot(UC_fibs, features = "RUNX1")

```

```{r}
grn_Runx1 <- grn %>% filter(Source=="RUNX1")

UC_fibs <- AddModuleScore(UC_fibs, features = list(grn_Runx1$Target), name = "RUNX1grn")



cols <- ArchR::paletteDiscrete(UC_fibs@meta.data[, "Cluster"])

DimPlot(UC_fibs, group.by = "Cluster", cols=cols)


DotPlot(UC_fibs, features = c("RUNX1", "RUNX1grn1"))+RotatedAxis()+ scale_size(range = c(2, 8))



```

```{r}
meta <- UC_fibs@meta.data

gex <- UC_fibs@assays[["RNA"]]@data %>% t() %>%  as.data.frame()

meta$Runx1 <- gex$RUNX1

#get colors (optional)
cols <- as.data.frame(ArchR::paletteDiscrete(UC_fibs@meta.data[, "Cluster"]))
colnames(cols) <- "colors"


meta  %>% 
ggplot(aes(x=RUNX1grn1, y=Health, fill=Health))+
    geom_boxplot(outlier.shape = 1, coef = 1)+theme_minimal()+theme(axis.ticks = element_blank(), axis.text.x=element_blank())+scale_fill_manual(values=c("#F47D2B", "#C06CAB", "#208A42"))+
  theme_cowplot(12)+RotatedAxis()+ theme(legend.position="none")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank())+ggtitle("Runx1 expression Tcell fibroblasts")
```

