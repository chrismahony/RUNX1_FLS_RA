```{r}


DimPlot(aggr, group.by = "integrated_snn_res.0.05")


aggr <- FindClusters(aggr, resolution = 0.08, graph.name= "integrated_snn")

DimPlot(aggr, group.by = "integrated_snn_res.0.08")
Idents(aggr) <- "integrated_snn_res.0.08"
markers_0.08 <- FindAllMarkers(aggr, only.pos = T)


DefaultAssay(aggr) <- 'RNA'
FeaturePlot(aggr, features= "THY1")
FeaturePlot(aggr, features= "PECAM1")
FeaturePlot(aggr, features= "RUNX1", split.by = "orig.ident")


pt <- table(aggr$integrated_snn_res.0.08, aggr$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank())



aggr$cluster_treat <- paste(aggr$integrated_snn_res.0.05, aggr$orig.ident, sep="_")
Idents(aggr) <- 'cluster_treat'

DEGs <- FindMarkers(aggr, ident.1 = "0_cbfbi", ident.2 = "0_DMSO")
DEGs <- DEGs %>% filter(p_val_adj < 0.05)
DEGs$gene <- DEGs %>% rownames()

table(aggr$cluster_treat)


FeaturePlot(aggr, features= "MMP14", split.by = "orig.ident")
FeaturePlot(aggr, features= "CXCL5", split.by = "orig.ident")
FeaturePlot(aggr, features= "IGF1", split.by = "orig.ident")


VlnPlot(aggr, features= "MMP14", group.by = "orig.ident")
VlnPlot(aggr, features= "CXCL5", group.by = "orig.ident")
VlnPlot(aggr, features= "IGF1", group.by = "orig.ident")

FeaturePlot(aggr, features= "COL3A1", split.by = "orig.ident")
FeaturePlot(aggr, features= "COL1A1", split.by = "orig.ident")
FeaturePlot(aggr, features= "THY1", split.by = "orig.ident")
FeaturePlot(aggr, features= "POSTN", split.by = "orig.ident")

```

```{r}



p1 <- DimPlot(aggr, group.by = "integrated_snn_res.0.05")+ 
  theme_void() + # Remove axes, labels, and background
  theme(
    plot.title = element_blank(), # Remove title
    legend.position = "bottom",   # Position the legend below
    legend.box = "horizontal",  
    legend.margin = margin(t = -50, unit = "pt"),# Horizontal legend layout
    panel.border = element_rect(color = "black", fill = NA, size = 1))+ggtitle("integrated_snn_res.0.05")
p2 <- DimPlot(aggr, group.by = "orig.ident")+ 
  theme_void() + # Remove axes, labels, and background
  theme(
    plot.title = element_blank(), # Remove title
    legend.position = "bottom",   # Position the legend below
    legend.box = "horizontal",
    legend.margin = margin(t = -50, unit = "pt"),# Horizontal legend layout
    panel.border = element_rect(color = "black", fill = NA, size = 1))+ggtitle("orig.ident")
p3 <- FeaturePlot(aggr, features= "POSTN", split.by = "orig.ident")+ 
  theme_void() + # Remove axes, labels, and background
  theme(
    plot.title = element_blank(), # Remove title
    legend.position = "bottom",   # Position the legend below
    legend.box = "horizontal",  
    legend.margin = margin(t = -50, unit = "pt"),# Horizontal legend layout
    panel.border = element_rect(color = "black", fill = NA, size = 1))+ggtitle("integrated_snn_res.0.05")

p4 <- FeaturePlot(aggr, features= "COL3A1", split.by = "orig.ident")
p5 <- FeaturePlot(aggr, features= "COL1A1", split.by = "orig.ident")
p6 <- FeaturePlot(aggr, features= "THY1", split.by = "orig.ident")


plot_grid(
  p1, p2, p3,
  p4, p5, p6,
  labels = c("A", "B", "C", "D", "E", "F"),
  ncol = 3, nrow = 2, align = "hv"
)


```



```{r}

# Deconvoluting donors

clusters_CBFB <- read.delim("/rds/projects/m/mahonyc-cesar-data/inhibitor_organoids/analysis/soupourcell/clusters_CBFB.tsv")
clusters_DMSO <- read.delim("/rds/projects/m/mahonyc-cesar-data/inhibitor_organoids/analysis/soupourcell/clusters_DMSO.tsv")

clusters_CBFB
clusters_DMSO


head(colnames(aggr))
tail(colnames(aggr))


clusters_DMSO$barcode <- paste(clusters_DMSO$barcode, "_1", sep="")
clusters_CBFB$barcode <- paste(clusters_CBFB$barcode, "_2", sep="")
clusters_DMSO$assignment <- paste(clusters_DMSO$assignment, "_DMSO", sep="")
clusters_CBFB$assignment <- paste(clusters_CBFB$assignment, "_CBFBI", sep="")


all_clusters <- rbind(clusters_DMSO, clusters_CBFB)

all_clusters <- all_clusters[all_clusters$barcode %in% colnames(aggr),]

all_clusters <- all_clusters[,c(1,3)]
all_clusters

rownames(all_clusters) <- all_clusters$barcode
all_clusters$barcode <- NULL
colnames(all_clusters) <- 'donor'
all_clusters

aggr <- AddMetaData(aggr, all_clusters)

DimPlot(aggr, group.by="donor")


Idents(aggr) <- 'donor'


aggr_f <- subset(aggr, idents=levels(aggr)[c(1,2,3,10,11,12)])
aggr_f %>% ncol

aggr_f <- aggr_f %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>% RunUMAP(dims=1:30)
DimPlot(aggr_f, group.by="donor")

library(harmony)
aggr_f <- RunHarmony(aggr_f, c("donor", "orig.ident"))
aggr_f <- RunUMAP(aggr_f, dims = 1:40, reduction="harmony")
DimPlot(aggr_f, group.by="donor")
DimPlot(aggr_f, group.by="orig.ident")

FeaturePlot(aggr_f, features="CD248")

aggr_f <- FindNeighbors(aggr_f, reduction = "harmony")
aggr_f <- FindClusters(aggr_f, resolution=c(0.05, 0.1, 0.15, 0.2, 0.25))

DimPlot(aggr_f, group.by = "RNA_snn_res.0.15")

Idents(aggr_f) <- 'RNA_snn_res.0.15'
markers_0.15 <- FindAllMarkers(aggr_f, only.pos = T)

DimPlot(aggr_f, group.by = "orig.ident")


library(DropletUtils)
#folder must not already exist i.e. you are creating a new one
write10xCounts(path="/rds/projects/m/mahonyc-cesar-data/inhibitor_organoids/analysis/processed", x=aggr_f@assays$RNA@counts)
write.table(aggr_f@meta.data %>% as.data.frame(), "/rds/projects/m/mahonyc-cesar-data/inhibitor_organoids/analysis/processed/meta.tsv", sep="\t")

```
```{r}


Idents(aggr_f)<-'donor'

cts_fibs<-AggregateExpression(aggr_f, group.by = c("donor"), assays = "RNA", slot = "counts", return.seurat = F)

cts_fibs<-cts_fibs$RNA
cts_fibs<-as.data.frame(cts_fibs)
meta_data=colnames(cts_fibs)
meta_data<-as.data.frame(meta_data)
library(splitstackshape)
meta_data$to_split<-meta_data$meta_data
meta_data<-cSplit(meta_data, splitCols = "to_split", sep="_")
colnames(meta_data)<-c("all","donor", "treatment")
meta_data$all<-as.factor(meta_data$all)
meta_data$donor<-as.factor(meta_data$donor)
meta_data$treatment<-as.factor(meta_data$treatment)


dds <- DESeqDataSetFromMatrix(countData = cts_fibs,
                                  colData = meta_data,
                                  design = ~1)
    
dds <- scran::computeSumFactors(dds)
print(dds)
print(quantile(rowSums(counts(dds))))

#mingenecount <- quantile(rowSums(counts(dds)), 0.5)
mingenecount <- 200
maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount & rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

dds@colData[['treatment']] <- as.factor(dds@colData[['treatment']])

design(dds) <- formula(~ treatment)
print(design(dds))
dds <- DESeq(dds, test = "Wald")



targetvar <- "treatment"

comps1 <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))
      head(comps1)
      
      ress <- apply(comps1, 1, function(cp) {
        print(cp)
        res <- data.frame(results(dds, contrast=c(targetvar, cp[1], cp[2])))
        res[["gene"]] <- rownames(res)
        res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
        res
      })
      
      res1 <- Reduce(rbind, ress)

res1 %>% 
      filter(padj < 0.05) %>%
      mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
      arrange(desc(abs(score))) -> subres

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

ss_sm <- meta_data[, c("donor", "treatment")]
col_ann <- HeatmapAnnotation(df = ss_sm)


filtered_subres <- subres %>% 
  filter(log2FoldChange < -2 | log2FoldChange > 2)

 
  Heatmap(scale_sub_vsd[filtered_subres$gene,], 
              top_annotation = col_ann, 
               col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = T,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F, border=T)




```

```{r}



subres$cluster <- "NO"
subres$cluster[subres$log2FoldChange > 0.001] <- "UP"
subres$cluster[subres$log2FoldChange < -0.001] <- "DOWN"


library(gsfisher)
annotation_gs<-fetchAnnotation(species = "hs")

index <- match(subres$gene, annotation_gs$gene_name)
subres$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- unique(subres$gene)
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]
ensemblUni <- na.omit(ensemblUni)
subres_f <- na.omit(subres)


go.results <- runGO.all(results=subres_f,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="padj", p_threshold=0.05,
                  species = "hs")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=5, -p.val)

sampleEnrichmentDotplot(go.results, selection_col = "description", selected_genesets = c("extracellular matrix organization", "integrin binding", "collagen binding"), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)




```

```{r}


DimPlot(stia2021_rna)


degs <- subres %>% filter(comparison == "CBFBI_vs_DMSO" & log2FoldChange < -2 & padj < 0.05)
degs_up <- subres %>% filter(comparison == "CBFBI_vs_DMSO" & log2FoldChange > 2 & padj < 0.05)


library(babelgene)
mouse <- orthologs(genes = degs$gene, species = "mouse")
mouse_up <- orthologs(genes = degs_up$gene, species = "mouse")

stia2021_rna <- AddModuleScore(stia2021_rna, features=list(mouse_up$symbol), name="inhibotor_up")

stia2021_rna <- AddModuleScore(stia2021_rna, features=list(mouse$symbol), name="inhibotor")
FeaturePlot(stia2021_rna, features="inhibotor_up1", max.cutoff = "q90", min.cutoff = "q10")

FeaturePlot(stia2021_rna, features="inhibotor1", max.cutoff = "q90", min.cutoff = "q10")
FeaturePlot(stia2021_rna, features="Runx1", max.cutoff = "q90", min.cutoff = "q10")



FeaturePlot(stia2021_rna, features="Prg4", max.cutoff = "q90", min.cutoff = "q10")

Idents(stia2021_rna) <- 'pseudo.bulk.level'
levels(stia2021_rna)
DotPlot(stia2021_rna, features="inhibotor1", idents= c("lining" ,"sublining"))

Idents(stia2021_rna) <- 'seurat_clusters'
DotPlot(stia2021_rna, features=c("inhibotor1"))


Idents(stia2021_rna) <- 'condition'
DotPlot(stia2021_rna, features="inhibotor1")
```

```{r}
stia2021_rna$pseudo.bulk.level %>% unique()

DimPlot(stia2021_rna, group.by= "pseudo.bulk.level")





```






```{r}


dds_norm <- counts(dds, normalized=TRUE)
dds_norm_df <- as.data.frame(dds_norm)
dds_norm_df$Gene <- rownames(dds_norm_df)
dds_long <- melt(dds_norm_df, id.vars = "Gene", variable.name = "Sample", value.name = "NormalizedCount")

library(splitstackshape)
dds_long <- cSplit(dds_long, splitCols="Sample",sep = "_")

gene = "THY1"
dds_long %>% filter(Gene == gene) %>% 
ggplot(aes(x = Sample_2, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_2), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip()



gene = "POSTN"
dds_long %>% filter(Gene == gene) %>% 
ggplot(aes(x = Sample_2, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_2), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip()


gene = "CTHRC1"
dds_long %>% filter(Gene == gene) %>% 
ggplot(aes(x = Sample_2, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_2), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip()



gene = "COL3A1"
dds_long %>% filter(Gene == gene) %>% 
ggplot(aes(x = Sample_2, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_2), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip()


gene = "COL1A2"
dds_long %>% filter(Gene == gene) %>% 
ggplot(aes(x = Sample_2, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_2), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip()



gene = "MMP14"
dds_long %>% filter(Gene == gene) %>% 
ggplot(aes(x = Sample_2, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_2), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip()


gene = "MMP9"
dds_long %>% filter(Gene == gene) %>% 
ggplot(aes(x = Sample_2, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_2), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip()





```

```{r}
gene = "PRG4"
dds_long %>% filter(Gene == gene) %>% 
ggplot(aes(x = Sample_2, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_2), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip()



gene = "CLIC5"
dds_long %>% filter(Gene == gene) %>% 
ggplot(aes(x = Sample_2, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_2), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip()


gene = "COL22A1"
dds_long %>% filter(Gene == gene) %>% 
ggplot(aes(x = Sample_2, y = NormalizedCount)) +
  geom_boxplot(aes(fill = Sample_2), outlier.shape = 16, outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x labels for clarity
  labs(title = paste0("Expression of ", gene), x = "Sample_2", y = "Normalized Count") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_ArchR()+coord_flip() 
```

```{r}

library(symphony)


#buiild ref
ref_exp_full = stia2021_rna@assays[["RNA"]]@data
rownames(ref_exp_full) <- toupper(rownames(ref_exp_full))
ref_metadata = stia2021_rna@meta.data

set.seed(0)
reference = symphony::buildReference(
    ref_exp_full,
    ref_metadata,
    vars = c('sample_id'),         # variables to integrate over
    K = 100,                   # number of Harmony clusters
    verbose = TRUE,            # verbose output
    do_umap = TRUE,            # can set to FALSE if want to run umap separately later
    do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
    vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
    vargenes_groups = 'sample_id', # metadata column specifying groups for variable gene selection 
    topn = 2000,               # number of variable genes to choose per group
    d = 20,                    # number of PCs
    save_uwot_path = "/rds/projects/m/mahonyc-cesar-data/inhibitor_organoids/analysis/umap"
)

#extract UMAP embeding
reference_umap<-stia2021_rna@reductions[["umap"]]@cell.embeddings
umap_labels = cbind(ref_metadata, reference_umap)
ggplot(umap_labels, aes(x=UMAP_1, y=UMAP_2, colour = seurat_clusters)) +   geom_point()

reference$umap <- list()
reference$umap$embedding <- reference_umap

#map query
query = symphony::mapQuery(
  aggr[['RNA']]@counts,
  aggr@meta.data,
  reference, 
  vars = 'orig.ident',  # use column names from your meta_data
  do_normalize = TRUE)

query = knnPredict(query, reference, ref_metadata$seurat_clusters, k = 5)
head(query$meta_data)

#rename col to match to combine outputs
reference$meta_data$cell_type_pred_knn=NA
reference$meta_data$cell_type_pred_knn_prob=NA
reference$meta_data$ref_query = 'reference'
query$meta_data$ref_query = 'query'

reference$meta_data = subset(reference$meta_data, select = c(orig.ident,cell_type_pred_knn,cell_type_pred_knn_prob, ref_query))
query$meta_data= subset(query$meta_data, select = c(orig.ident,cell_type_pred_knn,cell_type_pred_knn_prob, ref_query))

reference$meta_data$cell_type_pred_knn <- ref_metadata$seurat_clusters


meta_data_combined = rbind(query$meta_data, reference$meta_data)
umap_combined = rbind(query$umap, reference$umap$embedding)
umap_combined_labels = cbind(meta_data_combined, umap_combined)



ggplot(umap_combined_labels, aes(x=UMAP1, y=UMAP2, colour = ref_query)) +   geom_point(size=1, alpha=0.8) +theme_classic()

p1 <- umap_combined_labels %>% filter(ref_query == "reference") %>% 
ggplot(aes(x = UMAP1, y = UMAP2, colour = ref_query)) + 
  geom_point(size=1, alpha=1) + 
  theme_classic() +
  scale_color_manual(values = c("reference" = "grey", 
                                "query" = "red"))


p2 <- umap_combined_labels %>% filter(ref_query == "query") %>% 
ggplot(aes(x = UMAP1, y = UMAP2, colour = ref_query)) + 
  geom_point(size=1, alpha=1) + 
  theme_classic() +
  scale_color_manual(values = c("reference" = "grey", 
                                "query" = "red"))


ggplot() +
  geom_point(data = umap_combined_labels %>% filter(ref_query == "reference"), 
             aes(x = UMAP1, y = UMAP2, colour = ref_query), 
             size = 1, alpha = 1) + 
  geom_point(data = umap_combined_labels %>% filter(ref_query == "query"), 
             aes(x = UMAP1, y = UMAP2, colour = orig.ident), 
             size = 1, alpha = 1) + 
  theme_classic() +
  scale_color_manual(values = c("reference" = "grey", 
                                "query" = "red", "DMSO"="green", "cbfbi"="yellow"))


umap_combined_labels$orig.ident %>% unique

```

