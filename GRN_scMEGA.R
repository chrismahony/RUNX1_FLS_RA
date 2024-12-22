

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
load("~/croftap-stia-atac-path/CM_multiome/STIA_andATAC/STIA_2021_cesear_analysis/analysis_cm.RData")
load("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/analysis.RData")

rm(list=ls()[! ls() %in% c("stia2021_rna", "atac")])

FeaturePlot(stia2021_rna, features = c("Chad", "Chadl", "Sox9", "Cilp"))
FeaturePlot(stia2021_rna, features = c("Runx2", "Omd", "Sp7", "Alpl"))


```

```{r}

stia2021_rna<- stia2021_rna[,grepl("lining|sublining", stia2021_rna$pseudo.bulk.level, ignore.case=TRUE)]

stia2021_rna<- stia2021_rna[,!grepl("initiation|resolving", stia2021_rna$condition, ignore.case=TRUE)]

stia2021_rna  <- stia2021_rna[,!grepl("fibroblast__Runx2_Bglap|fibroblast__Clu", stia2021_rna$cluster.name, ignore.case=TRUE)]

```

```{r}
counts <- GetAssayData(stia2021_rna, assay = "RNA",slot = "counts")
rownames(counts) <- toupper(rownames(counts))
obj.rna <- CreateSeuratObject(counts = counts, meta.data = stia2021_rna@meta.data)
gene.activity <- atac@assays$ACTIVITY@counts
rownames(gene.activity) <- toupper(rownames(gene.activity))
counts <- GetAssayData(atac, assay = "ATAC",slot = "counts")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c("-", "-"),
  genome = 'mm10',
  min.cells = 0,
  min.features = 0
)
obj.atac <- CreateSeuratObject(counts = chrom_assay, meta.data = atac@meta.data, assay = "ATAC")

obj.rna <- obj.rna %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

obj.atac <- obj.atac %>% 
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 'q5') %>%
    RunSVD() %>%
    RunUMAP(reduction = 'lsi', dims = 2:30, verbose = FALSE)

obj.rna$orig.ident<-obj.rna$sample_id
obj.atac$cluster.name<-obj.atac$predicted.id

obj.atac[["ATACexpr"]]<-atac[["ATACexpr"]]


rna_control<- obj.rna[,grepl("control", obj.rna$condition, ignore.case=TRUE)]
rna_peak<- obj.rna[,grepl("peak", obj.rna$condition, ignore.case=TRUE)]
rna_persistent<- obj.rna[,grepl("persistent", obj.rna$condition, ignore.case=TRUE)]
rna_resolved<- obj.rna[,grepl("resolved", obj.rna$condition, ignore.case=TRUE)]

atac_control<- obj.atac[,grepl("rest", obj.atac$condition, ignore.case=TRUE)]
atac_peak<- obj.atac[,grepl("Infla", obj.atac$condition, ignore.case=TRUE)]
atac_persistent<- obj.atac[,grepl("Persis", obj.atac$condition, ignore.case=TRUE)]
atac_resolved<- obj.atac[,grepl("resolved", obj.atac$condition, ignore.case=TRUE)]

rm(atac, stia2021_rna)


```


```{r}
gene.activity_control<-gene.activity[,colnames(gene.activity) %in% colnames(atac_control)]
obj.coembed_control <- CoembedData(
  rna_control,
  atac_control, 
  gene.activity_control, 
  weight.reduction = "lsi", 
  verbose = FALSE
)

gene.activity_peak<-gene.activity[,colnames(gene.activity) %in% colnames(atac_peak)]
obj.coembed_peak <- CoembedData(
  rna_peak,
  atac_peak, 
  gene.activity_peak, 
  weight.reduction = "lsi", 
  verbose = FALSE
)

gene.activity_persistent<-gene.activity[,colnames(gene.activity) %in% colnames(atac_persistent)]
obj.coembed_persistent <- CoembedData(
  rna_persistent,
  atac_persistent, 
  gene.activity_persistent, 
  weight.reduction = "lsi", 
  verbose = FALSE
)

gene.activity_resolved<-gene.activity[,colnames(gene.activity) %in% colnames(atac_resolved)]
obj.coembed_resolved <- CoembedData(
  rna_resolved,
  atac_resolved, 
  gene.activity_resolved, 
  weight.reduction = "lsi", 
  verbose = FALSE
)

df.pair_rest <- PairCells(object = obj.coembed_control, reduction = "pca",
                    pair.by = "tech", ident1 = "ATAC", ident2 = "RNA")

df.pair_peak <- PairCells(object = obj.coembed_peak, reduction = "pca",
                    pair.by = "tech", ident1 = "ATAC", ident2 = "RNA")

df.pair_persis <- PairCells(object = obj.coembed_persistent, reduction = "pca",
                    pair.by = "tech", ident1 = "ATAC", ident2 = "RNA")

df.pair_resolved <- PairCells(object = obj.coembed_resolved, reduction = "pca",
                    pair.by = "tech", ident1 = "ATAC", ident2 = "RNA")

save.image("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/analysis.RData")


sel_cells <- c(df.pair_rest$ATAC, df.pair_rest$RNA)
obj.coembed_control <- obj.coembed_control[, sel_cells]
obj.pair_control <- CreatePairedObject(df.pair = df.pair_rest, 
                               object = obj.coembed_control,
                               use.assay1 = "RNA", 
                               use.assay2 = "ATAC")

sel_cells <- c(df.pair_peak$ATAC, df.pair_peak$RNA)
obj.coembed_peak <- obj.coembed_peak[, sel_cells]
obj.pair_peak <- CreatePairedObject(df.pair = df.pair_peak, 
                               object = obj.coembed_peak,
                               use.assay1 = "RNA", 
                               use.assay2 = "ATAC")

sel_cells <- c(df.pair_persis$ATAC, df.pair_persis$RNA)
obj.coembed_persistent <- obj.coembed_persistent[, sel_cells]
obj.pair_persis <- CreatePairedObject(df.pair = df.pair_persis, 
                               object = obj.coembed_persistent,
                               use.assay1 = "RNA", 
                               use.assay2 = "ATAC")

sel_cells <- c(df.pair_resolved$ATAC, df.pair_resolved$RNA)
obj.coembed_resolved <- obj.coembed_resolved[, sel_cells]
obj.pair_resolved <- CreatePairedObject(df.pair = df.pair_resolved, 
                               object = obj.coembed_resolved,
                               use.assay1 = "RNA", 
                               use.assay2 = "ATAC")

all_coembed_merge_nomural <- merge(obj.pair_control, y = c(obj.pair_peak, obj.pair_persis, obj.pair_resolved))

all_coembed_merge_nomural <- all_coembed_merge_nomural %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:30,  verbose = FALSE) %>%
    RunUMAP(dims = 2:25, verbose = FALSE)

obj.pair_all_coembed_merge_nomural <- AddTrajectory(object = all_coembed_merge_nomural, 
                          trajectory = c("sublining", "lining"),
                          group.by = "pseudo.bulk.level", 
                          reduction = "pca",
                          dims = 1:3, 
                          use.all = FALSE)
                          
# we only plot the cells that are in this trajectory
obj.pair_all_coembed_merge_nomural <- obj.pair_all_coembed_merge_nomural[, !is.na(obj.pair_all_coembed_merge_nomural$Trajectory)]

TrajectoryPlot(object = obj.pair_all_coembed_merge_nomural, 
                    reduction = "umap",
                    continuousSet = "blueYellow",
                    size = 1,
                   addArrow = FALSE)

save.image("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/analysis.RData")

DimPlot(obj.pair_all_coembed_merge_nomural, group.by = "pseudo.bulk.level")
DimPlot(obj.pair_all_coembed_merge_nomural, group.by = "cluster.name")

```
```{r}


DefaultAssay(obj.pair_all_coembed_merge_nomural) <- "RNA"
FeaturePlot(obj.pair_all_coembed_merge_nomural, features=c("ZEB1", "EBF1"))


DefaultAssay(obj.pair_all_coembed_merge_nomural) <- "RNA"
FeaturePlot(obj.pair_all_coembed_merge_nomural, features=c("BATF", "FOSL1", "CLIC5"), ncol=3)
FeaturePlot(obj.pair_all_coembed_merge_nomural, features=c("CLIC5"))


grep("MA0477", rownames(obj.pair_all_coembed_merge_nomural@assays[["chromvar"]]))
rownames(obj.pair_all_coembed_merge_nomural@assays[["chromvar"]])[684]


DefaultAssay(obj.pair_all_coembed_merge_nomural) <- "chromvar"
FeaturePlot(obj.pair_all_coembed_merge_nomural, features=c("MA0477.2", "MA1634.1"), max.cutoff = "q90", min.cutoff="q10")



#grep("MA0154", rownames(obj.pair_all_coembed_merge_nomural@assays[["chromvar"]]))
#rownames(obj.pair_all_coembed_merge_nomural@assays[["chromvar"]])[675]

#grep("MA0103", rownames(obj.pair_all_coembed_merge_nomural@assays[["chromvar"]]))
#rownames(obj.pair_all_coembed_merge_nomural@assays[["chromvar"]])[412]


DefaultAssay(obj.pair_all_coembed_merge_nomural) <- "chromvar"
FeaturePlot(obj.pair_all_coembed_merge_nomural, features=c("MA0154.4","MA0103.3"), max.cutoff = "q90", min.cutoff="q10")

Idents(obj.pair_all_coembed_merge_nomural) <- 'name'
DefaultAssay(obj.pair_all_coembed_merge_nomural) <- "RNA"
DotPlot(obj.pair_all_coembed_merge_nomural, features=c("ZEB1", "EBF1"))
DefaultAssay(obj.pair_all_coembed_merge_nomural) <- "chromvar"
DotPlot(obj.pair_all_coembed_merge_nomural, features=c("MA0154.4","MA0103.3"))

```


```{r}
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
obj.pair_all_coembed_merge_nomural <- AddMotifs(
  object = obj.pair_all_coembed_merge_nomural,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm,
    assay = "ATAC"
)


DefaultAssay(obj.pair_all_coembed_merge_nomural)<-'ATAC'

library(BiocParallel)
register(SerialParam())

obj.pair_all_coembed_merge_nomural <- RunChromVAR(
  object = obj.pair_all_coembed_merge_nomural,
  genome = BSgenome.Mmusculus.UCSC.mm10,
    assay = "ATAC"
)

register(MulticoreParam(40, progressbar = F))

res <- SelectTFs(object = obj.pair_all_coembed_merge_nomural, return.heatmap = TRUE)
df.cor <- res$tfs
ht <- res$heatmap
max.overlaps=10000
draw(ht)
```
```{r}
res <- SelectGenes(object = obj.pair_all_coembed_merge_nomural,
                  labelTop1 = 0,
                  labelTop2 = 0)

df.p2g <- res$p2g
ht2 <- res$heatmap

draw(ht2)
```

```{r}
tf.gene.cor <- GetTFGeneCorrelation(object = obj.pair_all_coembed_merge_nomural, 
                                    tf.use = df.cor$tfs, 
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar", 
                                    gene.assay = "RNA",
                                    trajectory.name = "Trajectory")

ht3 <- GRNHeatmap(tf.gene.cor, 
                 tf.timepoint = df.cor$time_point, column_title_gp = gpar(fontsize = 4, fontface = "bold"))




ht3
```


```{r}

ht3
```



```{r}
motif.matching <- obj.pair_all_coembed_merge_nomural@assays$ATAC@motifs@data
colnames(motif.matching) <- obj.pair_all_coembed_merge_nomural@assays$ATAC@motifs@motif.names
motif.matching <-
    motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]


df.grn <- GetGRN(motif.matching = motif.matching, 
                 df.cor = tf.gene.cor, 
                 df.p2g = df.p2g)

df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs

# plot the graph, here we can highlight some genes
df.grn2 <- df.grn %>%
    subset(correlation > 0.4) %>%
    select(c(tf, gene, correlation)) %>%
    rename(weights = correlation)

library(igraph)
library(ggraph)
p <- GRNPlot(df.grn2, 
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = F,
             seed = 42, 
             plot.importance = T,
            min.importance = 2,
            remove.isolated = FALSE,
            genes.highlight=c("RUNX1", "MMP14", "CTHRC1", "IGF1"))

p
```
```{r}
netobj <- graph_from_data_frame(df.grn2,directed = TRUE)
V(netobj)$type <- ifelse(V(netobj)$name %in% df.grn2$tf,"TF/Gene","Gene")
TopEmbGRN(df.grn=netobj)
TopEmbGRN(df.grn=netobj,axis=c(2,3))
max.overlaps=500
NetCentPlot(netobj,"RUNX1", highlights = c("CTHRC1", "C1QTNF3", "ANGPTL1", "MMP14", "MMP3", "COL14A1", "RUNX1"))
```

```{r}
grn <- df.grn2
head(grn)
dim(grn)


grn<-grn[grn$tf=="RUNX1",]

write.csv(grn, "/rds/projects/c/croftap-runx1data01/Dkk3_analysis/bulk_RNA/grn.csv")

#grn<-grn[grn$tf=="RUNX1",]
grn<-grn[grn$weights > 0.6,]

head(grn)
dim(grn)
colnames(grn)[1:3] <- c("Source", "Target", "Weight")
head(grn)

grn[["Type"]] <- "Directed"
grn %>%
  mutate('Id' = paste0(Source, '_', Target)) -> grn
head(grn)

# - community detection
# - igraph definition
g <- graph_from_data_frame(grn[, c("Source", "Target")],
                           directed = FALSE)
g <- set_edge_attr(g, "weight", value = grn$Weight)
g <- set_edge_attr(g, "name", value = grn$Id)

library(igraph)

# leiden
leiden_mod <- cluster_louvain(g)
mods <- data.frame(cbind(V(g)$name, leiden_mod$membership))

colnames(mods) <- c("Id", "leiden_mod")
head(mods)
table(mods$leiden_mod)

#write.table(grn, "~/stroma_stia/atac_rna.dir/scmega.dir/grn_dis_filtered_gephi_edges.tsv",
            #sep = "\t", quote = FALSE, row.names = FALSE)


nodes <- data.frame('Id' = unique(grn$Source))
nodes[["class"]] <- "TF"

target_nodes <- data.frame('Id' = unique(grn$Target))
target_nodes[["class"]] <- "gene"
target_nodes %>%
  filter(!Id %in% nodes$Id) -> target_nodes


nodes <- rbind(nodes, target_nodes)
nodes[["Label"]] <- nodes$Id

nodes <- merge(nodes, mods, by = "Id", all.x = TRUE, sort = FALSE)
head(nodes)
dim(nodes)



e <- get.edgelist(g,names=FALSE)

l <- qgraph::qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                                 area=10*(vcount(g)^2),
                                                 repulse.rad=(vcount(g)^3.1),
                                                 niter = 1000,
                                                 max.delta = 4, 
                                                 cool.exp = 0.3)


plot(leiden_mod, g,
       layout=l,
       vertex.label.cex=0.01, 
       vertex.label.family="Helvetica",
       vertex.label.font=0.5,
       vertex.shape="square", 
       vertex.size=0, 
       col = rgb(1,1,1,0), 
       edge.color = rgb(0,0,0,0.02),
       edge.width = 0.5,
     vertex.label=get.vertex.attribute(g)$name)

nodes$mouse_name <- nodes$Id
nodes$mouse_name <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(nodes$mouse_name), perl=TRUE)
nodes$class <- NULL

write.table(dplyr::select(nodes, c(Id, Label, leiden_mod, mouse_name)), 
            "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/nodes_grn.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.csv(nodes, "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/nodes_grn.csv")

#colnames(grn)<-c("Source", "Target", "")

write.table(grn, "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/gephi_edges_GRN.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


#write.csv(grn, "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/gephi_edges_GRN.csv")


nodes_2<-as.data.frame(nodes[nodes$leiden_mod==2,])


nodes_2$Id<-gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(nodes_2$Id), perl=TRUE)

stia2021_rna<-AddModuleScore(stia2021_rna, features = list(nodes_2$Id), name = "louvain_mod2")

FeaturePlot(stia2021_rna, features="louvain_mod21",  max.cutoff = "q95", min.cutoff ="q5")+NoAxes()

```
```{r}
grn <- df.grn2

if(grn) {
  grn %>%
    dplyr::select(c(tf, gene, weights)) -> grn
} 

if(grn) {
  grn %>%
    dplyr::select(c(tf, gene, correlation, fdr)) %>%
    dplyr::filter(correlation > 0) -> grn
}

if(filter_grn) {
  quant_thr <- 0.6
  quantile(grn$correlation, quant_thr)
  
  grn %>%
    filter(correlation > quantile(grn$correlation, quant_thr)) -> grn 
  
}

head(grn)
dim(grn)
colnames(grn)[1:3] <- c("Source", "Target", "Weight")
head(grn)

grn[["Type"]] <- "Directed"
grn %>%
  mutate('Id' = paste0(Source, '_', Target)) -> grn
head(grn)

# - community detection
# - igraph definition
g <- graph_from_data_frame(grn[, c("Source", "Target")],
                           directed = FALSE)
g <- set_edge_attr(g, "weight", value = grn$Weight)
g <- set_edge_attr(g, "name", value = grn$Id)

# leiden
leiden_mod <- cluster_leiden(g, objective_function = "modularity")
mods <- data.frame(cbind(V(g)$name, leiden_mod$membership))
colnames(mods) <- c("Id", "leiden_mod")
head(mods)
table(mods$leiden_mod)

write.table(grn, "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/grn_dis_filtered_gephi_edges.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


nodes <- data.frame('Id' = unique(grn$Source))
nodes[["class"]] <- "TF"

target_nodes <- data.frame('Id' = unique(grn$Target))
target_nodes[["class"]] <- "gene"
target_nodes %>%
  filter(!Id %in% nodes$Id) -> target_nodes


nodes <- rbind(nodes, target_nodes)
nodes[["Label"]] <- nodes$Id

nodes <- merge(nodes, mods, by = "Id", all.x = TRUE, sort = FALSE)
head(nodes)
dim(nodes)

# mouse gene names
hg2mm <- fread("~/work/synovium_atlas/data/biomart_hg2mm_20211003.tsv")
hg2mm %>%
  dplyr::select(`Gene name`, `Mouse gene name`) %>%
  filter(!duplicated(`Gene name`)) -> hg2mm
head(hg2mm)
dim(hg2mm)

nds <- merge(nodes, hg2mm, by.x = "Id", by.y = "Gene name",  all.x = T, sort = F)

idx <- which(is.na(nds$`Mouse gene name`))
nds[["Mouse gene name"]][idx] <- nds[["Id"]][idx]

dim(nds)
head(nds)

write.table(dplyr::select(nds, c(Id, Label, class, leiden_mod, `Mouse gene name`)), 
            "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/grn_dis_filtered_gephi_nodes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------ 


ndsf <- fread("~/stroma_stia/atac_rna.dir/scmega.dir/grn_full_nodes_gephi_v3_modularity.csv")
ndsf

ndsf %>%
  mutate('module' = paste0('full_M-', modularity_class, "_")) %>%
  dplyr::select(c(Id, class, module)) -> ndsf

ndsf

ndsfi <- fread("~/stroma_stia/atac_rna.dir/scmega.dir/grn_filtered_nodes_gephi_v3_modularity.csv")
ndsfi

ndsfi %>%
  mutate('module' = paste0('filt_M-', modularity_class, "_")) %>%
  dplyr::select(c(Id, class, module)) -> ndsfi

head(ndsfi)

nds <- rbind(ndsf, ndsfi)
head(nds)
dim(nds)

table(nds$module)

hg2mm <- fread("~/work/synovium_atlas/data/biomart_hg2mm_20211003.tsv")
hg2mm %>%
  dplyr::select(`Gene name`, `Mouse gene name`) %>%
  filter(!duplicated(`Gene name`)) -> hg2mm
head(hg2mm)
dim(hg2mm)

nds <- merge(nds, hg2mm, by.x = "Id", by.y = "Gene name",  all.x = T, sort = F)

idx <- which(is.na(nds$`Mouse gene name`))
nds[["Mouse gene name"]][idx] <- nds[["Id"]][idx]

nds %>%
  dplyr::select(`Mouse gene name`, module, class) -> nds

colnames(nds)[1] <- "Id"

dim(nds)
head(nds)

table(nds$module)

write.table(nds, "~/stroma_stia/atac_rna.dir/scmega.dir/grn_modules_full_and_filtered.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
```


```{r}
library(clusterProfiler)

library(gsfisher)
annotation_gs <- fetchAnnotation(species="mm", ensembl_version=NULL, ensembl_host=NULL)

nodes$Id<-gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(nodes$Id), perl=TRUE)

library("org.Mm.eg.db")
nodes_2<-as.data.frame(nodes[nodes$leiden_mod==2,])

index <- match(nodes_2$Id, annotation_gs$gene_name)
nodes_2$entrez_id <- annotation_gs$entrez_id[index]

index <- match(nodes$Id, annotation_gs$gene_name)
nodes$entrez_id <- annotation_gs$entrez_id[index]

nodes_2_f<-na.omit(nodes_2)
nodes_f<-na.omit(nodes)


ego <- enrichGO(gene          = nodes_2_f$entrez_id,
                universe      = nodes_f$entrez_id,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = T)



nodes_1<-as.data.frame(nodes[nodes$leiden_mod==1,])
nodes_1$Id<-gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(nodes_1$Id), perl=TRUE)

index <- match(nodes_1$Id, annotation_gs$gene_name)
nodes_1$entrez_id <- annotation_gs$entrez_id[index]

index <- match(nodes$Id, annotation_gs$gene_name)
nodes$entrez_id <- annotation_gs$entrez_id[index]

nodes_1_f<-na.omit(nodes_1)
nodes_f<-na.omit(nodes)


ego_1 <- enrichGO(gene          = nodes_1_f$entrez_id,
                universe      = nodes_f$entrez_id,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = T)


nodes_3<-as.data.frame(nodes[nodes$leiden_mod==3,])
nodes_3$Id<-gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(nodes_3$Id), perl=TRUE)

index <- match(nodes_3$Id, annotation_gs$gene_name)
nodes_3$entrez_id <- annotation_gs$entrez_id[index]

index <- match(nodes$Id, annotation_gs$gene_name)
nodes$entrez_id <- annotation_gs$entrez_id[index]

nodes_3_f<-na.omit(nodes_3)
nodes_f<-na.omit(nodes)


ego_3 <- enrichGO(gene          = nodes_3_f$entrez_id,
                universe      = nodes_f$entrez_id,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = T)


nodes_4<-as.data.frame(nodes[nodes$leiden_mod==4,])
nodes_4$Id<-gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(nodes_4$Id), perl=TRUE)

index <- match(nodes_4$Id, annotation_gs$gene_name)
nodes_4$entrez_id <- annotation_gs$entrez_id[index]

index <- match(nodes$Id, annotation_gs$gene_name)
nodes$entrez_id <- annotation_gs$entrez_id[index]

nodes_4_f<-na.omit(nodes_4)
nodes_f<-na.omit(nodes)


ego_4 <- enrichGO(gene          = nodes_4_f$entrez_id,
                universe      = nodes_f$entrez_id,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = T)


nodes_5<-as.data.frame(nodes[nodes$leiden_mod==5,])
nodes_5$Id<-gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(nodes_5$Id), perl=TRUE)

index <- match(nodes_5$Id, annotation_gs$gene_name)
nodes_5$entrez_id <- annotation_gs$entrez_id[index]

index <- match(nodes$Id, annotation_gs$gene_name)
nodes$entrez_id <- annotation_gs$entrez_id[index]

nodes_5_f<-na.omit(nodes_5)
nodes_f<-na.omit(nodes)


ego_5 <- enrichGO(gene          = nodes_5_f$entrez_id,
                universe      = nodes_f$entrez_id,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = T)


goplot(ego)
```
```{r}
library(DOSE)

library(enrichplot)
barplot(ego, showCategory=20) 

mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore")
mutate(ego_1, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore")
mutate(ego_3, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore")
mutate(ego_4, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore")
mutate(ego_5, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore")


edox <- setReadable(ego, 'org.Mm.eg.db')


cnetplot(edox, foldChange=nodes_2_f$entrez_id)
p1

p1  <- heatplot(edox, showCategory=5)
p1
```


```{r}


tf<-colnames(as.data.frame(ht3@matrix))

tf_use<-meta_mdata[meta_mdata$motif %in% tf,]
write.csv(tf_use, "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/tf_use.csv")

tf<-as.data.frame(tf)
tf$index=rownames(tf)

aggr$sampleID=(samples_ID[match(rownames(aggr@meta.data),samples_ID$Barcode),2])

tf<-tf[tf$tf %in% tf_use$motif,]
tf_use$index=rownames(tf_use)

tf_use$index=(tf[match(tf$tf,tf_use$motif),2])

rownames(tf_use)=tf_use$motif

tail(tf$tf)
tf_use["USF2",]
```






```{r}
# mouse gene names
hg2mm <- fread("~/work/synovium_atlas/data/biomart_hg2mm_20211003.tsv")
hg2mm %>%
  dplyr::select(`Gene name`, `Mouse gene name`) %>%
  filter(!duplicated(`Gene name`)) -> hg2mm
head(hg2mm)
dim(hg2mm)

nds <- merge(nodes, hg2mm, by.x = "Id", by.y = "Gene name",  all.x = T, sort = F)

idx <- which(is.na(nds$`Mouse gene name`))
nds[["Mouse gene name"]][idx] <- nds[["Id"]][idx]

dim(nds)
head(nds)

write.table(dplyr::select(nds, c(Id, Label, class, leiden_mod, `Mouse gene name`)), 
            "~/stroma_stia/atac_rna.dir/scmega.dir/grn_dis_filtered_gephi_nodes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------ 


ndsf <- fread("~/stroma_stia/atac_rna.dir/scmega.dir/grn_full_nodes_gephi_v3_modularity.csv")
ndsf

ndsf %>%
  mutate('module' = paste0('full_M-', modularity_class, "_")) %>%
  dplyr::select(c(Id, class, module)) -> ndsf

ndsf

ndsfi <- fread("~/stroma_stia/atac_rna.dir/scmega.dir/grn_filtered_nodes_gephi_v3_modularity.csv")
ndsfi

ndsfi %>%
  mutate('module' = paste0('filt_M-', modularity_class, "_")) %>%
  dplyr::select(c(Id, class, module)) -> ndsfi

head(ndsfi)

nds <- rbind(ndsf, ndsfi)
head(nds)
dim(nds)

table(nds$module)

hg2mm <- fread("~/work/synovium_atlas/data/biomart_hg2mm_20211003.tsv")
hg2mm %>%
  dplyr::select(`Gene name`, `Mouse gene name`) %>%
  filter(!duplicated(`Gene name`)) -> hg2mm
head(hg2mm)
dim(hg2mm)

nds <- merge(nds, hg2mm, by.x = "Id", by.y = "Gene name",  all.x = T, sort = F)

idx <- which(is.na(nds$`Mouse gene name`))
nds[["Mouse gene name"]][idx] <- nds[["Id"]][idx]

nds %>%
  dplyr::select(`Mouse gene name`, module, class) -> nds

colnames(nds)[1] <- "Id"

dim(nds)
head(nds)

table(nds$module)

write.table(nds, "~/stroma_stia/atac_rna.dir/scmega.dir/grn_modules_full_and_filtered.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

```




```{r}


DimPlot(all_coembed_merge_nomural, group.by = "name")
DefaultAssay(all_coembed_merge_nomural) <- 'RNA'
FeaturePlot(all_coembed_merge_nomural, features="RUNX1")
Idents()
DefaultAssay(all_coembed_merge_nomural) <- 'ATAC'


DimPlot(obj.atac)

obj.atac %>% ncol

rownames(all_coembed_merge_nomural)


stia2021_rna$cluster_id
FeaturePlot(stia2021_rna, features="Clu")
DimPlot(stia2021_rna, group.by ="cluster.name" )


Idents(stia2021_rna) <- "cluster.name"
stia2021_rna_filt <- subset(stia2021_rna, idents=levels(stia2021_rna)[-c(2,3,5)])
stia2021_rna_filt <- stia2021_rna_filt %>% ScaleData()
DimPlot(stia2021_rna_filt, group.by ="cluster.name" )
FeaturePlot(stia2021_rna_filt, features="Fap")


DimPlot(stia2021_rna_filt, group.by = "cluster.name") +
  theme_minimal(base_size = 14) +  # Use a minimal theme for a clean appearance
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text = element_blank(),   # Remove axis text
    axis.ticks = element_blank(),  # Remove axis ticks
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_blank(),  # Remove default background
    plot.background = element_rect(color = "black", fill = NA, size = 1)  # Add black border
  ) +
  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "blue"
  )) +  # Use custom colors for clusters
  labs(
    title = "STIA scRNAseq",
    color = "Cluster"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  # Center and style the title
  )




```

