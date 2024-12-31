


```{r}
fibs_GSE129605<-readRDS("/rds/projects/m/mahonyc-kitwong-runx1/GSE129605_lung_bleo/fibs_GSE129605.rds")
fibs_GSE11164<-readRDS("/rds/projects/m/mahonyc-kitwong-runx1/GSE111664_Lung_bleomycin/fibs_GSE11164.rds")
fibs_GSE132771<-readRDS("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE132771_lung_fibrosis/fibs_GSE132771.rds")

fibs_GSE129605 <- fibs_GSE129605 %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

fibs_GSE11164 <- fibs_GSE11164 %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

fibs_GSE132771 <- fibs_GSE132771 %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)


head(colnames(fibs_GSE11164))
tail(colnames(fibs_GSE11164))


library_id<-c("control_1_1", "control_1_2", "Control_2", "Control_3", "Control_4", "Control_5", "Control_6", "Bleo_1", "Bleo_2", "Bleo_3")

samples<-as.data.frame(library_id)


cellcodes <- as.data.frame(fibs_GSE11164@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cellcodes) <- "barcodes"
rownames(cellcodes) <- cellcodes$barcodes
cellcodes$libcodes=cellcodes$barcodes
library(splitstackshape)
cellcodes<-cSplit(cellcodes, splitCols="barcodes", sep="_")
cellcodes$barcodes_2<-as.factor(cellcodes$barcodes_2)
cellcodes$samples <- as.vector(samples$library_id[cellcodes$barcodes_2])
cellcodes<-as.data.frame(cellcodes)
rownames(cellcodes) = cellcodes$libcodes

cellcodes<-select(cellcodes, c(samples))
colnames(cellcodes)<-'orig.ident'
fibs_GSE11164<-AddMetaData(fibs_GSE11164, cellcodes)
```


```{r}
fibs_GSE129605$dataset<-'GSE129605'
fibs_GSE11164$dataset<-'GSE11164'
fibs_GSE132771$dataset<-'GSE132771'

FeaturePlot(fibs_GSE11164, features = "Runx1")
#fibs_GSE11164 as is was not seq'ed in this experiment, counts=0!

list=c(fibs_GSE129605, fibs_GSE132771)
anchors <- FindIntegrationAnchors(object.list = list, dims = 1:30)
aggr <- IntegrateData(anchorset = anchors, dims = 1:30)
aggr <- ScaleData(aggr, verbose = FALSE)
aggr <- RunPCA(aggr, verbose = FALSE)
aggr <- RunUMAP(aggr, dims = 1:50)
DimPlot(aggr, group.by = "orig.ident")

```

```{r}
library(harmony)
RunHarmony.Seurat_CM <- function(
  object,
  group.by.vars,
  reduction = 'pca',
  dims.use = NULL,
  theta = NULL,
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4,
  plot_convergence = FALSE,
  verbose = TRUE,
  reference_values = NULL,
  reduction.save = "harmony",
  assay.use = NULL,
  project.dim = TRUE,
  ...
) {
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Harmony on a Seurat object requires Seurat")
  }
  assay.use <- assay.use %||% Seurat::DefaultAssay(object)
  if (reduction == "pca" && !reduction %in% Seurat::Reductions(object = object)) {
    if (isTRUE(x = verbose)) {
      message("Harmony needs PCA. Trying to run PCA now.")
    }
    object <- tryCatch(
      expr = Seurat::RunPCA(
        object = object,
        assay = assay.use,
        verbose = verbose,
        reduction.name = reduction
      ),
      error = function(...) {
        stop("Harmony needs PCA. Tried to run PCA and failed.")
      }
    )
  }
  if (!reduction %in% Seurat::Reductions(object = object)) {
    stop("Requested dimension reduction is not present in the Seurat object")
  }
  embedding <- Seurat::Embeddings(object, reduction = reduction)
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(
    object,
    group.by.vars,
    cells = Seurat::Cells(x = object[[reduction]])
  )

  harmonyEmbed <- HarmonyMatrix(
    embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    FALSE,
    verbose,
    reference_values
  )

  #reduction.key <- Seurat::Key(reduction.save, quiet = TRUE)
  reduction.key <- "harmony_"
  rownames(harmonyEmbed) <- rownames(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.key, seq_len(ncol(harmonyEmbed)))

  object[[reduction.save]] <- Seurat::CreateDimReducObject(
    embeddings = harmonyEmbed,
    stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
    assay = Seurat::DefaultAssay(object = object[[reduction]]),
    key = reduction.key
  )
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  return(object)
}
```

```{r}
aggr_harmony<-RunHarmony.Seurat_CM(aggr, group.by.vars = c("orig.ident", "dataset"), assay.use = "RNA", reduction = "umap")

ncol(aggr_harmony)

DimPlot(aggr_harmony, reduction = "harmony", group.by = "orig.ident")
```
```{r}

aggr_harmony<-FindNeighbors(aggr_harmony, dims = 1:2, reduction = "harmony")
aggr_harmony<-FindClusters(aggr_harmony, resolution=c(0.01, 0.05, 0.1, 0.2, 0.3), graph.name = "integrated_snn")
library(clustree)
clustree(aggr_harmony)
DimPlot(aggr_harmony, group.by = "dataset", reduction = "harmony")
DefaultAssay(aggr_harmony) <- 'RNA'


DimPlot(aggr_harmony, group.by = "integrated_snn_res.0.1")
FeaturePlot(aggr_harmony, reduction = "harmony", features = "Runx1")
FeaturePlot(aggr_harmony, reduction = "harmony", features = "Igf1")
FeaturePlot(aggr_harmony, reduction = "harmony", features = "Mmp14")



Idents(aggr_harmony) <- 'integrated_snn_res.0.1'
markers_0.1 <- FindAllMarkers(aggr_harmony, only.pos = T)


aggr_harmony$named <- aggr_harmony@meta.data[["integrated_snn_res.0.1"]]
Idents(aggr_harmony) <- 'named'
levels(aggr_harmony)
current.sample.ids <- c("0","1","2", "3", "4", "5", "6", "7")
new.sample.ids <- c("Npnt","1","CD34_PI16", "3")

obj@meta.data[["new_ident"]] <- plyr::mapvalues(x = obj@meta.data[["new_ident"]], from = current.sample.ids, to = new.sample.ids)




```
```{r}


grn <- read_csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/grn.csv")


grn$Target <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(grn$Target), perl=TRUE)

grn_RUNX1 <- grn %>% filter(Source == "RUNX1")


Idents(aggr)<-'orig.ident'
levels(aggr)
table(aggr$orig.ident)
aggr$condition<-aggr$orig.ident
current.sample.ids<-c("bleo_N1" ,"bleo_N2", "bleo_N3", "bleo1" ,  "bleo2" ,  "bleo3",   "bleo4" ,  "saline1", "saline2", "saline3", "saline4", "Bleo1",   "Bleo2"  , "UT1" ,    "UT2"     )



new.sample.ids <- c("bleo_N" ,"bleo_N", "bleo_N", "bleo" ,  "bleo" ,  "bleo",   "bleo" ,  "control", "control", "control", "control", "bleo",   "bleo"  , "control" ,    "control"   )
aggr@meta.data[["condition"]] <- plyr::mapvalues(x = aggr@meta.data[["condition"]], from = current.sample.ids, to = new.sample.ids)

aggr_harmony$condition<-aggr$condition

Idents(aggr_harmony)<-'condition'
DotPlot(aggr_harmony, features=c("Runx1", "Mmp14", "Igf1"), idents=levels(aggr_harmony)[2:3])

fibs_bleo_vs_contorl<-FindMarkers(aggr_harmony, ident.1 = "bleo", ident.2 = "control")
fibs_bleo_vs_contorl$gene<-rownames(fibs_bleo_vs_contorl)

EnhancedVolcano(fibs_bleo_vs_contorl,
    lab = rownames(fibs_bleo_vs_contorl),
    x = 'avg_log2FC',
    y = 'p_val_adj',
        selectLab = "Runx1",
    title = 'title',
    subtitle = "GEX, red=p_adj<0.05 & FC > 0.25",
    pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.0,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5
    )


```


```{r}

library(DESeq2)
cts_t<-AggregateExpression(aggr_harmony, group.by = c("orig.ident"), assays = "RNA", slot = "counts", return.seurat = F)

cts_t<-cts_t$RNA
#cts.t<-t(cts)
cts_t<-as.data.frame(cts_t)
samples_2=colnames(cts_t)
samples_2<-as.data.frame(samples_2)
samples_2$condition<-c("bleo_N" ,    "bleo_N" ,    "bleo_N"  ,   "bleo" ,      "bleo"     , 
  "bleo"   ,    "bleo"  ,     "bleo"  ,     "bleo"  ,     "control", "control", "control",   "control" , 
 "control"  , "control" )


dds_t<-DESeqDataSetFromMatrix(countData = cts_t, colData=samples_2, design = ~ condition )
keep<-rowSums(counts(dds_t))>=10
dds_t<-dds_t[keep,]
#dds_t$samples_2 <- relevel(dds_t$samples_2, ref = "Cxcl5 Tnn_Control")
dds_t<-DESeq(dds_t)
res <- results(dds_t)
resultsNames(dds_t)
resLFC <- as.data.frame(lfcShrink(dds_t, coef="condition_control_vs_bleo", type="apeglm"))
resLFC$gene<-rownames(resLFC)
normalized_counts <- as.data.frame(counts(dds_t, normalized=TRUE))
normalized_counts$gene<-rownames(normalized_counts)

 Runx1<- normalized_counts["Runx1", ]
 Runx1$gene<-NULL
 
 Runx1<-as.data.frame(t(Runx1))
 Runx1$condition<-c("bleo_N" ,    "bleo_N" ,    "bleo_N"  ,   "bleo" ,      "bleo"     , 
  "bleo"   ,    "bleo"  ,     "bleo"  ,     "bleo"  ,     "control", "control", "control",   "control" , 
 "control"  , "control" )

 ggplot(Runx1, aes(x=factor(condition, c("control", "bleo")), y=Runx1)) + 
    geom_violin()+   geom_violin(trim=F,fill='#A4A4A4', color="darkred")+ geom_boxplot(width=0.1)+geom_jitter(shape=16, position=position_jitter(0.2), size=3)+ scale_fill_grey() + theme_classic()
```





```{r}
Idents(aggr_harmony)<-"orig.ident"
dotplot<-DotPlot(aggr_harmony, features = "Runx1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"

dotplot_MMP14<-DotPlot(aggr_harmony, features = "Mmp14")
dotplotmmp14_data<-dotplot_MMP14[["data"]]
dotplotmmp14_data <- subset(dotplotmmp14_data, select = c(avg.exp.scaled, id))
names(dotplotmmp14_data)[names(dotplotmmp14_data)=="avg.exp.scaled"] <- "Mmp14"

dotplot_Cxcl5<-DotPlot(aggr_harmony, features = "Cxcl5")
dotplotCxcl5_data<-dotplot_Cxcl5[["data"]]
dotplotCxcl5_data <- subset(dotplotCxcl5_data, select = c(avg.exp.scaled, id))
names(dotplotCxcl5_data)[names(dotplotCxcl5_data)=="avg.exp.scaled"] <- "Cxcl5"

dotplot_Igf1<-DotPlot(aggr_harmony, features = "Igf1")
dotplotIgf1_data<-dotplot_Igf1[["data"]]
dotplotIgf1_data <- subset(dotplotIgf1_data, select = c(avg.exp.scaled, id))
names(dotplotIgf1_data)[names(dotplotIgf1_data)=="avg.exp.scaled"] <- "Igf1"

dotplot_data$Mmp14<-dotplotmmp14_data$Mmp14
dotplot_data$Cxcl5<-dotplotCxcl5_data$Cxcl5
dotplot_data$Igf1<-dotplotIgf1_data$Igf1


dotplot_data$condition<-c("bleo_N" ,    "bleo_N" ,    "bleo_N"  ,   "bleo" ,      "bleo"     , 
  "bleo"   ,    "bleo"  ,     "bleo"  ,     "bleo"  ,     "control", "control", "control",   "control" , 
 "control"  , "control")

ggplot(dotplot_data, aes(x = Mmp14, y = Runx1)) +
    geom_point(aes(color = factor(condition)), size=5) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none")

ggplot(dotplot_data, aes(x = Cxcl5, y = Runx1)) +
    geom_point(aes(color = factor(condition)), size=5) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none")


ggplot(dotplot_data, aes(x = Igf1, y = Runx1)) +
    geom_point(aes(color = factor(condition)), size=5) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none")


ml = lm(Cthrc1~Runx1, data = dotplot_data)
summary(ml)$r.squared
ml = lm(Mmp14~Runx1, data = dotplot_data)
summary(ml)$r.squared
ml = lm(Igf1~Runx1, data = dotplot_data)
summary(ml)$r.squared

library(ggpubr)
ggscatter(dotplot_data, x = "Mmp14", y = "Runx1",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)

ggscatter(dotplot_data, x = "Igf1", y = "Runx1",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)

ggscatter(dotplot_data, x = "Cthrc1", y = "Runx1",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)

```

