
```{r}

setwd("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE138669_Scleradoma")
dirs <-dir("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE138669_Scleradoma/", pattern = ".h5")

samples <-sapply(strsplit(dirs, split='_', fixed=TRUE), function(x) (x[2]))
samples <- gsub("raw","", samples)


dirs <- paste0("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE138669_Scleradoma/", dirs)

data.10x <- list()
for (i in 1:length(dirs)){
  data.10x[[i]] <- Read10X_h5(filename = dirs[[i]])
  data.10x[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=3, min.features=0, project=samples[i]);
data.10x[[i]][["percent.mt"]] = PercentageFeatureSet(object=data.10x[[i]], pattern = "^mt-");
data.10x[[i]]<-subset(x = data.10x[[i]], subset = nFeature_RNA > 200 & nFeature_RNA <7000 & percent.mt < 10);
data.10x[[i]] =NormalizeData(object = data.10x[[i]]);
data.10x[[i]] =ScaleData(object = data.10x[[i]]);
data.10x[[i]] =FindVariableFeatures(object = data.10x[[i]]);
data.10x[[i]] =RunPCA(object = data.10x[[i]], verbose = FALSE);
data.10x[[i]] =RunUMAP(object = data.10x[[i]], dims=1:30)
}

names(data.10x) <- samples


anchors <- FindIntegrationAnchors(object.list = data.10x, reduction = "rpca",   dims = 1:50)
aggr <- IntegrateData(anchorset = anchors, dims = 1:50)
aggr <- FindVariableFeatures(aggr)
aggr <- ScaleData(aggr, verbose = FALSE)
aggr <- RunPCA(aggr, verbose = FALSE)
aggr <- RunUMAP(aggr, dims = 1:50)
aggr <- FindNeighbors(aggr, dims = 1:20)
aggr <- FindClusters(aggr, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3), graph.name = 'integrated_snn')

DimPlot(aggr)
DefaultAssay(aggr) <- 'RNA'
FeaturePlot(aggr, "PDGFRA")
FeaturePlot(aggr, "COL1A1")
FeaturePlot(aggr, "RUNX1")
Idents(aggr)
```
```{r}


aggr$condition <- aggr@meta.data[["orig.ident"]]
Idents(aggr) <- 'condition'
levels(aggr)
current.sample.ids <- levels(aggr)
new.sample.ids <- c(rep("CTRL",8), rep("SCC",14))

aggr@meta.data[["condition"]] <- plyr::mapvalues(x = aggr@meta.data[["condition"]], from = current.sample.ids, to = new.sample.ids)

Idents(aggr) <- 'integrated_snn_res.0.01'
DimPlot(aggr)

aggr <- FindClusters(aggr, resolution = c(0.005), graph.name = 'integrated_snn')
Idents(aggr) <- 'integrated_snn_res.0.01'
DimPlot(aggr)

DimPlot(aggr, group.by = "orig.ident")

FeaturePlot(aggr, features = "PDGFRA")+scale_color_viridis()
FeaturePlot(aggr, features = "COL1A1")+scale_color_viridis()
FeaturePlot(aggr, features = "CD248")+scale_color_viridis()
FeaturePlot(aggr, features = "THY1")+scale_color_viridis()


```
```{r}

skin_fibs <- subset(aggr, idents="2")
Idents(skin_fibs) <-'condition' 
VlnPlot(skin_fibs, features="RUNX1")

skin_fibs <- skin_fibs %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>% RunUMAP(dims=1:30)

skin_fibs <- FindNeighbors()


Idents(skin_fibs) <- 'integrated_snn_res.0.01'

DimPlot(skin_fibs)
FeaturePlot(skin_fibs, features = "RUNX1")
skin_fibs_sclerdoma <- skin_fibs
rm(skin_fibs)

DimPlot(skin_fibs_sclerdoma)

ncol(skin_fibs_sclerdoma)


```

```{r}
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


grn <- read.csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/df_grn2_new_FINAL.csv")

grn <- grn %>% filter(tf=="RUNX1")

skin_fibs_sclerdoma <- AddModuleScore(skin_fibs_sclerdoma, features= list(grn$gene), name="RUNX1grn")

Idents(skin_fibs_sclerdoma) <- 'condition'
DotPlot(skin_fibs_sclerdoma, features=c("RUNX1", "RUNX1grn1"))

                 
```

