

```{r}
library(Seurat)

#devtools::install_github("chrismahony/sprokforlife") #install if required
library(sporkforlife)

samples <- c("C1", "C2", "C3", "R1", "R2", "R3")

#optain a list of dirs to samples and smaple names
data.10x = list()
dirs <- list.dirs("/rds/projects/m/mahonyc-runx1-bulk-seq-data/mouse_KO_single_cell/count", recursive = FALSE)  #path to where cellrnager fisished the count step
dirs <- paste(dirs, "/outs/per_sample_outs/", samples, "/count/sample_filtered_feature_bc_matrix", sep="")

sample_names <- samples

process_scrna_data(dirs, sample_names, target_n_clusters = 5,
                               resolution_range = seq(0.05, 0.3, by = 0.5),
                               min_nFeature_RNA = 500, max_nFeature_RNA = 7000,
                               max_percent_mt = 10, n_dims=50)

DefaultAssay(aggr) <- 'RNA'



data.10x <- list()
    scrna.list <- list()
    
    for (i in 1:length(dirs)) {
        data.10x[[i]] <- Read10X(data.dir = dirs[[i]])
    }
    
    
    min_nFeature_RNA = 500
    max_nFeature_RNA = 7000
    max_percent_mt = 10
    n_dims=50
    
    
    for (i in 1:length(data.10x)) {
        scrna.list[[i]] <- CreateSeuratObject(counts = data.10x[[i]], 
            min.cells = 3, min.features = 0, project = sample_names[i])
        scrna.list[[i]][["percent.mt"]] <- PercentageFeatureSet(object = scrna.list[[i]], 
            pattern = "^mt-")
        scrna.list[[i]] <- subset(scrna.list[[i]], subset = nFeature_RNA > 
            min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & 
            percent.mt < max_percent_mt)
        scrna.list[[i]] <- NormalizeData(object = scrna.list[[i]])
        scrna.list[[i]] <- ScaleData(object = scrna.list[[i]])
        scrna.list[[i]] <- FindVariableFeatures(object = scrna.list[[i]])
        scrna.list[[i]] <- RunPCA(object = scrna.list[[i]], verbose = FALSE)
    }
    names(scrna.list) <- sample_names
    anchors <- FindIntegrationAnchors(object.list = scrna.list, 
        dims = 1:50)
    aggr <- IntegrateData(anchorset = anchors, dims = 1:50)
    aggr <- FindVariableFeatures(aggr)
    aggr <- ScaleData(aggr, verbose = FALSE)
    aggr <- RunPCA(aggr, verbose = FALSE)
    aggr <- RunUMAP(aggr, dims = 1:n_dims)
    aggr <- FindNeighbors(aggr, dims = 1:n_dims)
    
    
    aggr <- RunUMAP(aggr, dims = 1:20)
    ElbowPlot(aggr)
    aggr <- FindNeighbors(aggr, dims = 1:20)

    aggr <- FindClusters(aggr, resolution= c(0.01, 0.05, 0.1))
    aggr <- FindClusters(aggr, resolution= c(0.005))

DimPlot(aggr, group.by = "integrated_snn_res.0.1")

DefaultAssay(aggr) <- 'RNA'
DimPlot(aggr, group.by = "integrated_snn_res.0.1")
FeaturePlot(aggr, features = "Prg4")
FeaturePlot(aggr, features = "Clu")
FeaturePlot(aggr, features = "Bglap")
FeaturePlot(aggr, features = "Pecam1")
FeaturePlot(aggr, features = "Cd68")
FeaturePlot(aggr, features = "Pdgfra")
FeaturePlot(aggr, features = "Mki67")
FeaturePlot(aggr, features = "Lyve1")
FeaturePlot(aggr, features = "C1qa")

Idents(aggr) <- "integrated_snn_res.0.1"

markers_aggr_0.1 <- FindAllMarkers(aggr, only.pos = T)
```


```{r}

# Get list of cells that are pdgfra pos/neg

FeaturePlot(aggr, features="Pdgfra")

EXPR_Pdgfra = GetAssayData(object=aggr,assay="RNA",slot="data")["Pdgfra",]
EXPR_Pdgfra_df=data.frame( positive= EXPR_Pdgfra > 0, negative = EXPR_Pdgfra == 0)
names(EXPR_Pdgfra_df)<-paste0( c("positive_","negative_"),"Pdgfra")
aggr <- AddMetaData(aggr,metadata=EXPR_Pdgfra_df)

DimPlot(aggr, group.by="positive_Pdgfra")

Idents(aggr) <- "positive_Pdgfra"
levels(aggr)

pdgfra_pos <- subset(aggr, idents="TRUE")
pdgfra_neg <- subset(aggr, idents="FALSE")

DimPlot(pdgfra_pos, group.by="positive_Pdgfra")
DimPlot(pdgfra_neg, group.by="positive_Pdgfra")

barcodes <- (paste("CB:Z:", colnames(pdgfra_pos), sep=""))

barcodes <- sub("_[^_]+$", "", barcodes)

writeLines(barcodes, "/rds/projects/m/mahonyc-runx1-bulk-seq-data/mouse_KO_single_cell/analysis/subset_Check/barcodes_pdgfra_pos.txt")


barcodes <- (paste("CB:Z:", colnames(pdgfra_neg), sep=""))

barcodes <- sub("_[^_]+$", "", barcodes)

writeLines(barcodes, "/rds/projects/m/mahonyc-runx1-bulk-seq-data/mouse_KO_single_cell/analysis/subset_Check/barcodes_pdgfra_neg.txt")

```





```{r}
DimPlot(aggr, group.by = "integrated_snn_res.0.1", label=T)
Idents(aggr) <- 'integrated_snn_res.0.1'
stromal_cells <- subset(aggr, idents=c("1", "2", "3", "4"))

stromal_cells <- stromal_cells %>% ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:30)

DimPlot(stromal_cells)

library(harmony)
stromal_cells <- RunHarmony(stromal_cells, "orig.ident")
stromal_cells <- RunUMAP(stromal_cells, reduction = "harmony", dims=1:30)

stromal_cells <- stromal_cells %>% FindNeighbors(dims=1:30, reduction="harmony") %>% FindClusters(resolution = 0.05)
DimPlot(stromal_cells, group.by = "orig.ident")

stromal_cells <- stromal_cells %>% FindClusters(resolution = 0.1)

DimPlot(stromal_cells, group.by = "RNA_snn_res.0.1", label=T)
FeaturePlot(stromal_cells, features = "Pdgfra")
FeaturePlot(stromal_cells, features = "Prg4")
FeaturePlot(stromal_cells, features = "Clu")
FeaturePlot(stromal_cells, features = "Bglap")
FeaturePlot(stromal_cells, features = "Mki67")

```

```{r}

Idents(stromal_cells) <- 'RNA_snn_res.0.1'
fibs <- subset(stromal_cells, idents=c("0", "1"))

fibs <- fibs %>% ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:20)

fibs <- RunHarmony(fibs, "orig.ident")
fibs <- RunUMAP(fibs, reduction = "harmony", dims=1:20)
DimPlot(fibs, group.by="orig.ident")
fibs <- fibs %>% FindNeighbors(dims=1:20, reduction="harmony")

fibs <- fibs %>% FindClusters(resolution = c(0.1, 0.2, 0.3))
DimPlot(fibs, group.by="RNA_snn_res.0.3")

FeaturePlot(fibs, features="Pi16")
FeaturePlot(fibs, features="Cxcl5")
FeaturePlot(fibs, features="Cxcl5")

FeaturePlot(fibs, features="Clu")
DimPlot(fibs, group.by="RNA_snn_res.0.3", label=T)


Idents(fibs) <- 'RNA_snn_res.0.3'
markers_fibs <- FindAllMarkers(fibs, only.pos = T, logfc.threshold = 0.5, min.pct = 0.2)


fibs$named_subclusters <- fibs@meta.data[["RNA_snn_res.0.3"]]

current.sample.ids <- c( "0","1","2", "3", "4", "5","6","7", "8")
new.sample.ids <- c("Cd34/Pi16","Comp","Prg4", "Col15a1", "Cxcl5", "Fmo2","Cd34/Pi16","Postn/Cthrc1", "Comp")

fibs@meta.data[["named_subclusters"]] <- plyr::mapvalues(x = fibs@meta.data[["named_subclusters"]], from = current.sample.ids, to = new.sample.ids)

Idents(fibs) <- 'named_subclusters'
DotPlot(fibs, features=c("Pi16","Comp","Prg4", "Col15a1", "Cxcl5","Fmo2", "Cthrc1"))+RotatedAxis()
```
```{r}


fibs$condition <- fibs@meta.data[["orig.ident"]]
Idents(fibs) <- 'condition'

current.sample.ids <- c( "C1", "C2", "C3", "R1", "R2", "R3")
new.sample.ids <- c("Ctrl", "Ctrl", "Ctrl", "no_Runx1", "no_Runx1", "no_Runx1")

fibs@meta.data[["condition"]] <- plyr::mapvalues(x = fibs@meta.data[["condition"]], from = current.sample.ids, to = new.sample.ids)

cols <- ArchR::paletteDiscrete(fibs@meta.data[, "named_subclusters"])
DimPlot(fibs,group.by="named_subclusters", cols=cols, split.by="condition")
DimPlot(fibs,group.by="orig.ident")


DotPlot(fibs, features="Runx1", group.by="named_subclusters")

Idents(fibs) <- 'named_subclusters'
fibs_markers_named <- FindAllMarkers(fibs, only.pos = T)
write.csv(fibs_markers_named, "/rds/projects/m/mahonyc-runx1-bulk-seq-data/mouse_KO_single_cell/analysis/fibs_markers_named.csv")

pt <- table(fibs$named_subclusters, fibs$condition)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank())




```

```{r}
library(scProportionTest)

#create object for test
test <- sc_utils(fibs_SL)


#run test and plot. "cluster.name" is the name of my clusters (e.g. fibroblast, macrophage etc), sample_1/smaple_2 is the name of the conditions you are comparing (e.g. control treated), smaple identity would be somehting like treatments
prop.test <- permutation_test(test, cluster_identity = "named_subclusters", sample_1="no_Runx1", sample_2="Ctrl", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.25, order_clusters = T)
```





```{r}

fibs$cluster_condition <- paste(fibs$named_subclusters, fibs$condition, sep="_")
Idents(fibs) <- 'cluster_condition'

exp <- paste(unique(fibs$named_subclusters), "no_Runx1", sep="_")
ctrl <- paste(unique(fibs$named_subclusters), "Ctrl", sep="_")

DEGs <- list()

for (i in 1:length(exp)){
  DEGs[[i]] <- FindMarkers(fibs, ident.1 =exp[i], ident.2 = ctrl[i] )
}

names(DEGs) <- unique(fibs$named_subclusters)

Idents(fibs) <- 'cluster_condition'
levels(fibs)
VlnPlot(fibs, idents=levels(fibs)[c(5,12)], features=c("Cxcl5","Tnfaip6", "Angptl4", "Tnn" ), ncol=4)

VlnPlot(fibs, idents=levels(fibs)[c(2,10)], features=c("C3","Col5a1", "Cxcl12", "Tnn" ), ncol=4)


```


```{r}
deg_summary <- data.frame(
  condition = character(),
  direction = character(),
  gene_count = numeric()
)

# Fill summary
for (i in 1:length(DEGs)) {
  df <- DEGs[[i]]
  
  up_genes <- sum(df$avg_log2FC > 0 & df$p_val_adj < 0.05)
  down_genes <- sum(df$avg_log2FC < 0 & df$p_val_adj < 0.05)
  
  cond <- names(DEGs)[i]
  deg_summary <- rbind(
    deg_summary,
    data.frame(condition = cond, direction = "Up", gene_count = up_genes),
    data.frame(condition = cond, direction = "Down", gene_count = down_genes)
  )
}



library(ggplot2)

# Flip downregulated genes to negative for left-side bars
deg_summary$plot_count <- ifelse(deg_summary$direction == "Down",
                                 -deg_summary$gene_count,
                                 deg_summary$gene_count)


deg_summary <- deg_summary %>%
  group_by(condition) %>%
  mutate(total = sum(gene_count)) %>%
  ungroup()

# Reorder factor levels of 'condition' based on total
deg_summary$condition <- factor(deg_summary$condition, 
                                levels = deg_summary %>%
                                           distinct(condition, total) %>%
                                           arrange(desc(total)) %>%
                                           pull(condition))

# Bar plot
deg_summary2 <- deg_summary %>% filter(condition != "Prg4")
ggplot(deg_summary2, aes(x = condition, y = plot_count, fill = direction)) +
  geom_bar(stat = "identity", width = 1, color="black") +
  geom_text(aes(label = abs(gene_count)), 
            vjust = ifelse(deg_summary2$plot_count > 0, 0.5, 0.5),
            hjust= ifelse(deg_summary2$plot_count > 0, 2, -0.5),
            size = 3.5) +
  scale_y_continuous(name = "Number of DEGs") +
  scale_fill_manual(values = c("Up" = "red", "Down" = "grey")) +
  labs(title = "Up- and Down-regulated Genes per Condition") +
  theme_minimal()+coord_flip()+theme_ArchR()+geom_vline(xintercept = 0.0, color = "black", size=1.5)


deg_summary2 <- deg_summary %>% filter(condition != "Prg4")
ggplot(deg_summary2, aes(x = condition, y = plot_count, fill = direction)) +
  geom_bar(stat = "identity", width = 1, color="black")  +
  scale_y_continuous(name = "Number of DEGs") +
  scale_fill_manual(values = c("Up" = "red", "Down" = "grey")) +
  labs(title = "Up- and Down-regulated Genes per Condition") +
  theme_minimal()+coord_flip()+theme_ArchR()+geom_vline(xintercept = 0.0, color = "black", size=1.5)


```




```{r}



genes <-DEGs[["Cxcl5"]] %>% filter(avg_log2FC < -0.5 |avg_log2FC > 0.5  & p_val_adj < 0.05) %>% rownames()

DEGs[["Cxcl5"]]$in_list <- ifelse(rownames(DEGs[["Cxcl5"]]) %in% genes, "yes", "no")

mycolors <- c("darkred", "grey")
names(mycolors) <- c("yes", "no")

ggplot(data=DEGs[["Cxcl5"]], aes(x=avg_log2FC, y=-log10(p_val_adj),  col=in_list)) + geom_point(alpha=0.7, size=2) + theme_minimal()+ geom_vline(xintercept=c(-0.6, 0.6), col="darkred", linetype = "dashed") +
    geom_hline(yintercept=-log10(0.05), col="darkred", linetype = "dashed")+theme_ArchR()+ scale_colour_manual(values = mycolors)

DEGs[["Cxcl5"]]$gene_name <- rownames(DEGs[["Cxcl5"]])

ggplot(data = DEGs[["Cxcl5"]], aes(x = avg_log2FC, y = -log10(p_val_adj), col = in_list)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(data = subset(DEGs[["Cxcl5"]], in_list == "yes"),
                  aes(label = gene_name),
                  size = 3.5, max.overlaps = 50) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "darkred", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dashed") +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  theme_ArchR()


ggplot(data = DEGs[["Cxcl5"]], aes(x = avg_log2FC, y = -log10(p_val_adj), col = in_list)) +
  geom_point(alpha = 0.7, size = 2)  +
  geom_vline(xintercept = c(-0.5, 0.5), col = "darkred", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dashed") +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  theme_ArchR()


for (i in 1:length(DEGs)){
  DEGs[[i]]$cluster <- names(DEGs)[i]
  DEGs[[i]]$gene_name <- rownames(DEGs[[i]])
}
all_DEGs <- rbindlist(DEGs, fill=TRUE)
all_DEGs$in_list <- NULL
write.csv(all_DEGs, "/rds/projects/m/mahonyc-runx1-bulk-seq-data/mouse_KO_single_cell/analysis/DEGs.csv")

```

```{r}
genes <-DEGs[["Postn/Cthrc1"]] %>% filter(avg_log2FC < -0.8 |avg_log2FC > 0.8  & p_val_adj < 0.05) %>% rownames()

DEGs[["Postn/Cthrc1"]]$in_list <- ifelse(rownames(DEGs[["Postn/Cthrc1"]]) %in% genes, "yes", "no")

mycolors <- c("darkred", "grey")
names(mycolors) <- c("yes", "no")

ggplot(data=DEGs[["Postn/Cthrc1"]], aes(x=avg_log2FC, y=-log10(p_val_adj),  col=in_list)) + geom_point(alpha=0.7, size=2) + theme_minimal()+ geom_vline(xintercept=c(-0.6, 0.6), col="darkred", linetype = "dashed") +
    geom_hline(yintercept=-log10(0.05), col="darkred", linetype = "dashed")+theme_ArchR()+ scale_colour_manual(values = mycolors)

DEGs[["Postn/Cthrc1"]]$gene_name <- rownames(DEGs[["Postn/Cthrc1"]])

ggplot(data = DEGs[["Postn/Cthrc1"]], aes(x = avg_log2FC, y = -log10(p_val_adj), col = in_list)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(data = subset(DEGs[["Postn/Cthrc1"]], in_list == "yes"),
                  aes(label = gene_name),
                  size = 3.5, max.overlaps = 50) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "darkred", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dashed") +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  theme_ArchR()
```

```{r}
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", 
     "GO_Biological_Process_2023")

enriched_down <- enrichr(DEGs[["Cxcl5"]] %>% filter(avg_log2FC < -0.01) %>% rownames()
, dbs)

plot_down <- enriched_down[["GO_Biological_Process_2023"]]

plotEnrich(plot_down[c(3,44),], showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")

plot_small <- plot[c(3,44),]

goterms <- plot_small$Term

enriched_up <- enrichr(DEGs[["Cxcl5"]] %>% filter(avg_log2FC > 0.01) %>% rownames()
, dbs)

plot_up <- enriched_up[["GO_Biological_Process_2023"]]

plotEnrich(plot_up, showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")


plot_up_new <- plot_up %>% filter(Term %in% goterms)


plotEnrich(plot_up_new, showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")


# Go Term barplot

df <- data.frame(Gene_count=c(23,14,8,2),condition=c("down1", "down2", "up1", "up2"))


df$group <- gsub("[^0-9]", "", df$condition)
df$direction <- ifelse(grepl("down", df$condition), "down", "up")

# Set signed counts: down = negative, up = positive
df$directional_count <- ifelse(df$direction == "up", -df$Gene_count, df$Gene_count)

# Plot
ggplot(df, aes(x = group, y = directional_count, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("down" = "red", "up" = "grey")) +  # Set custom colors
  labs(y = "Gene Count", x = "Group") +
  coord_flip() +
  theme_minimal() +
  theme_ArchR()

get_genes <- plot[c(3,44),]
get_genes$Genes

genes1 <- strsplit(get_genes$Genes[1], ";")[[1]]
genes2 <- strsplit(get_genes$Genes[2], ";")[[1]]

combined_genes <- unique(c(genes1, genes2))

combined_genes <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(combined_genes), perl=TRUE)

DEGs[["Cxcl5"]] %>% filter(gene_name %in% combined_genes)


#genes <-DEGs[["Cxcl5"]] %>% filter(avg_log2FC < -0.2 |avg_log2FC > 0.2  & p_val_adj < 0.05) %>% rownames()

DEGs[["Cxcl5"]]$in_list <- ifelse(rownames(DEGs[["Cxcl5"]]) %in% combined_genes[-c(9,10,13,19,20,21,23,25,27,29,31,36,3,4,12)], "yes", "no")

ggplot(data = DEGs[["Cxcl5"]], aes(x = avg_log2FC, y = -log10(p_val_adj), col = in_list)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(
    data = dplyr::filter(DEGs[["Cxcl5"]], gene_name %in% combined_genes[-c(9,10,13,19,20,21,23,25,27,29,31,36,3,4,12)]),
    aes(label = gene_name),
    size = 3.5, max.overlaps = 50, color="black"
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), col = "darkred", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dashed") +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  theme_ArchR()+
  scale_x_reverse()


ggplot(data = DEGs[["Cxcl5"]], aes(x = avg_log2FC, y = -log10(p_val_adj), col = in_list)) +
  geom_point(alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-0.2, 0.2), col = "darkred", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dashed") +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  theme_ArchR()+
  scale_x_reverse()

enriched <- enrichr(DEGs[["Postn/Cthrc1"]] %>% filter(avg_log2FC < -0.01) %>% rownames()
, dbs)

plot <- enriched[["GO_Biological_Process_2023"]]


plotEnrich(plot, showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")

plotEnrich(plot[c(2,3),], showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")



plot_down <- plot %>% filter(Term %in% plot_small$Term)

plotEnrich(plot_down, showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")

goterms <- plot_small$Term




enriched_up <- enrichr(DEGs[["Postn/Cthrc1"]] %>% dplyr::filter(avg_log2FC > 0.01) %>% rownames(), dbs)

plot_up <- enriched_up[["GO_Biological_Process_2023"]]


plot_downsig <- plot %>% filter(Adjusted.P.value < 0.05)

plot_up_and_down <- plot_up[plot_up$Term %in% plot_downsig$Term,]

plot_up_and_down_ns <- plot_up_and_down %>% filter(Adjusted.P.value > 0.05)


go_terms <- c("Negative Regulation Of Cell Migration (GO:0030336)", "Cartilage Development (GO:0051216)")


plot_down <- plot %>% filter(Term %in% go_terms)

plotEnrich(plot_down, showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")


plot_up_new <- plot_up %>% filter(Term %in% go_terms)

plotEnrich(plot_up_new, showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")


library(data.table)


Idents(fibs)<-'named_subclusters'
fibs_SL<-subset(fibs, idents=levels(fibs)[-c(3)])

fibs_SL$sample_condition <- paste(fibs_SL$orig.ident, fibs_SL$condition, sep=".")

cts<-AggregateExpression(fibs_SL, group.by = c("sample_condition"), assays = "RNA", slot = "counts", return.seurat = F)

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

dds@colData[['condition']] <- as.factor(dds@colData[['condition']])

design(dds) <- as.formula(paste0("~", "condition"))

print(design(dds))

dds <- DESeq(dds, test = "Wald")


print(resultsNames(dds))
targetvar <- "condition"
comps <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))
head(comps)

ress <- apply(comps, 1, function(cp) {
  print(cp)
  res <- data.frame(results(dds, contrast=c(targetvar, cp[2], cp[1])))
  res[["gene"]] <- rownames(res)
  res[["comparison"]] <- paste0(cp[2], "_vs_", cp[1])
  res
})



res <- Reduce(rbind, ress)

head(res)
comps

res %>% 
  filter(padj < 0.05) %>%
  mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
  arrange(desc(abs(score))) -> subres

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




ss_sm <- meta_data[, c("condition")]

col_ann <- HeatmapAnnotation(df = ss_sm)


library(colorRamp2)        

   
      Heatmap(scale_sub_vsd, 
              top_annotation = col_ann,
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F,
              border=T)

      
      

enriched <- enrichr(subres %>% filter(log2FoldChange < -0.01) %>% pull(gene)
, dbs)


plotEnrich(enriched$GO_Biological_Process_2023[c(3,4,82),], showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")

go_terms <- enriched$GO_Biological_Process_2023[c(3,4,82),] %>% pull(Term)


enriched <- enrichr(subres %>% filter(log2FoldChange > 0.01) %>% pull(gene)
, dbs)

all_go <- enriched$GO_Biological_Process_2023
all_go <- all_go[all_go$Term %in% go_terms,]

plotEnrich(all_go, showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")



df <- data.frame( Gene_count=c(5.388808,7.402666,4.692210,1.3389147,0.6051804,1.6517325	),condition=c("down1", "down2","down3", "up1", "up2", "up3"))


df$group <- gsub("[^0-9]", "", df$condition)
df$direction <- ifelse(grepl("down", df$condition), "down", "up")

# Set signed counts: down = negative, up = positive
df$directional_count <- ifelse(df$direction == "up", -df$Gene_count, df$Gene_count)

# Plot
ggplot(df, aes(x = group, y = directional_count, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("down" = "red", "up" = "grey")) +  # Set custom colors
  labs(y = "Odd ratio", x = "Group") +
  coord_flip() +
  theme_minimal() +
  theme_ArchR()



cell_migration_gene <- enriched$GO_Biological_Process_2023[c(3,4,82),]

cell_mig <- strsplit(cell_migration_gene$Genes[1], ";")[[1]]

combined_genes2 <- unique(c(combined_genes, cell_mig))

combined_genes2 <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(combined_genes2), perl=TRUE)


DEGs[["Postn/Cthrc1"]]$in_list <- ifelse(rownames(DEGs[["Postn/Cthrc1"]]) %in% combined_genes[-c(3,4,9,10,12,13,19,20,21,23,25,29,35,36,39,40,4,347,50,51,52)], "yes", "no")

ggplot(data = DEGs[["Postn/Cthrc1"]], aes(x = avg_log2FC, y = -log10(p_val_adj), col = in_list)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(
    data = dplyr::filter(DEGs[["Postn/Cthrc1"]], gene_name %in% combined_genes[-c(3,4,9,10,12,13,19,20,21,23,25,29,35,36,39,40,43,47,50,51,52)]),
    aes(label = gene_name),
    size = 3.5, max.overlaps = 50, color="black"
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), col = "darkred", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dashed") +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  theme_ArchR()+
  scale_x_reverse()


ggplot(data = DEGs[["Postn/Cthrc1"]], aes(x = avg_log2FC, y = -log10(p_val_adj), col = in_list)) +
  geom_point(alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-0.2, 0.2), col = "darkred", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "darkred", linetype = "dashed") +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  theme_ArchR()+
  scale_x_reverse()



library(gsfisher)
annotation_gs<-fetchAnnotation(species = "mm")

subres$cluster <- "NO"
subres$cluster[subres$log2FoldChange > 0.001] <- "UP"
subres$cluster[subres$log2FoldChange < -0.001] <- "DOWN"


library(gsfisher)

index <- match(subres$gene, annotation_gs$gene_name)
subres$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- unique(subres$gene)
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]
ensemblUni <- na.omit(ensemblUni)

subref_f <- na.omit(subres)

go.results <- runGO.all(results=subref_f,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="padj", p_threshold=0.05,
                  species = "mm")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=5, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)





DimPlot(aggr)
DimPlot(fibs)

meta_aggr <- aggr@meta.data

meta_aggr$clusters <- meta_aggr$integrated_snn_res.0.1

meta_aggr <- meta_aggr %>% dplyr::select(orig.ident, clusters)

meta_fibs <- fibs@meta.data

meta_aggr_fibs <- meta_aggr[rownames(meta_aggr) %in% rownames(meta_fibs),]
meta_aggr_fibs$clusters <- 'fibs'

meta_aggr <- meta_aggr[!rownames(meta_aggr) %in% rownames(meta_aggr_fibs),]

meta_aggr <- rbind(meta_aggr, meta_aggr_fibs)

aggr <- AddMetaData(aggr, meta_aggr)

DimPlot(aggr, group.by = "clusters", label = T)


Idents(aggr) <- 'clusters'
markers_aggr_0.1 <- FindAllMarkers(aggr, only.pos = T)

FeaturePlot(aggr, features = "Prg4")
FeaturePlot(aggr, features = "Clu")
FeaturePlot(aggr, features = "Bglap")
FeaturePlot(aggr, features = "Pecam1")
FeaturePlot(aggr, features = "Cd68")
FeaturePlot(aggr, features = "Pdgfra")
FeaturePlot(aggr, features = "Mki67")
FeaturePlot(aggr, features = "Lyve1")
FeaturePlot(aggr, features = "C1qa")


DimPlot(aggr, group.by = "clusters", label = T)
FeaturePlot(aggr, features = "Cd11b")


aggr$clusters <- aggr@meta.data[["clusters"]]
Idents(aggr) <- 'clusters'

current.sample.ids <- levels(aggr)

new.sample.ids <- c( "Endothelial cells"  ,  "fibs"  ,  "lymphoid cells" ,  "Monocytes" ,  "Chondrocytes" ,   "Osteoblasts" ,   "fibs" ,   "Muscle" ,   "Muscle" ,   "Pericytes"   , "Endothelial cells" ,   "Contamination" ,   "fibs")

aggr@meta.data[["clusters"]] <- plyr::mapvalues(x = aggr@meta.data[["clusters"]], from = current.sample.ids, to = new.sample.ids)
DimPlot(aggr, group.by = "clusters", label = T)
table(aggr$orig.ident, aggr$clusters)

aggr$condition <- aggr@meta.data[["orig.ident"]]
Idents(aggr) <- 'condition'

current.sample.ids <- c( "C1", "C2", "C3", "R1", "R2", "R3")
new.sample.ids <- c("Ctrl", "Ctrl", "Ctrl", "no_Runx1", "no_Runx1", "no_Runx1")

aggr@meta.data[["condition"]] <- plyr::mapvalues(x = aggr@meta.data[["condition"]], from = current.sample.ids, to = new.sample.ids)

cols <- ArchR::paletteDiscrete(aggr@meta.data[, "clusters"])
DimPlot(aggr,group.by="clusters", cols=cols, split.by = "condition")
table(aggr$condition)



Idents(aggr) <- 'clusters'
levels(aggr)


DotPlot(aggr, features=c("Pecam1", "Vwf", "Pdgfra","Thy1", "Ccl4","Cxcl2", "Cd68","C1qa", "Clu","Sox9", "Bglap","Bglap2","Acta2", "Des"), idents=levels(aggr)[-c(7,9)])+RotatedAxis()



```
```{r}

aggr$cluster_condition <- paste(aggr$clusters, aggr$condition, sep="_")
Idents(aggr) <- 'cluster_condition'
levels(aggr)

myeloid_markers <- FindMarkers(aggr, ident.1 = "Monocytes_no_Runx1", ident.2 = "Monocytes_Ctrl")
#myeloid_markers <- myeloid_markers %>% filter(p_val_adj < 0.05)
myeloid_markers$gene_name <- rownames(myeloid_markers)

# Too few cel, pseudobul and DE
Idents(aggr)<-'clusters'
macs<-subset(aggr, idents=levels(aggr)[4])

macs$sample_condition <- paste(macs$orig.ident, macs$condition, sep=".")

cts_macs<-AggregateExpression(macs, group.by = c("sample_condition"), assays = "RNA", slot = "counts", return.seurat = F)

cts_macs<-cts_macs$RNA
cts_macs<-as.data.frame(cts_macs)
meta_data=colnames(cts_macs)
meta_data<-as.data.frame(meta_data)
library(splitstackshape)
meta_data$to_split<-meta_data$meta_data
meta_data<-cSplit(meta_data, splitCols = "to_split", sep=".")
colnames(meta_data)<-c("all", "sample", "condition")
meta_data$all<-as.factor(meta_data$all)
meta_data$sample<-as.factor(meta_data$sample)
meta_data$cluster<-as.factor(meta_data$condition)


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts_macs,
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

dds@colData[['condition']] <- as.factor(dds@colData[['condition']])

design(dds) <- as.formula(paste0("~", "condition"))

print(design(dds))

dds <- DESeq(dds, test = "Wald")


print(resultsNames(dds))
targetvar <- "condition"
comps <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))
head(comps)

ress <- apply(comps, 1, function(cp) {
  print(cp)
  res <- data.frame(results(dds, contrast=c(targetvar, cp[2], cp[1])))
  res[["gene"]] <- rownames(res)
  res[["comparison"]] <- paste0(cp[2], "_vs_", cp[1])
  res
})



res_macs <- Reduce(rbind, ress)



res_macs %>% 
  filter(padj < 0.05) %>%
  mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
  arrange(desc(abs(score))) -> subres_macs



library(ComplexHeatmap)

      if(length(unique(subres_macs$gene)) > 10) {
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
        
        feats <- unique(subres_macs$gene)
        print(length(feats))
        
        # Sub-set matrix to relevant features
        sub_vsd_mat_macs <- vsd_mat[rownames(vsd_mat) %in% feats, ]
        scale_sub_vsd_macs <- t(scale(t(sub_vsd_mat_macs)))
      }
      }




ss_sm <- meta_data[, c("condition")]

col_ann <- HeatmapAnnotation(df = ss_sm)


library(colorRamp2)        

   
      Heatmap(scale_sub_vsd_macs, 
              top_annotation = col_ann,
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F,
              border=T)

      

library(DropletUtils)
#folder must not already exist i.e. you are creating a new one
write10xCounts(path="/rds/projects/m/mahonyc-runx1-bulk-seq-data/mouse_KO_single_cell/analysis/geo_upload/", x=fibs@assays$RNA@counts)
write.table(fibs@meta.data %>% as.data.frame(), "/rds/projects/m/mahonyc-runx1-bulk-seq-data/mouse_KO_single_cell/analysis/geo_upload/meta.tsv", sep="\t")      
      






fibs <- AddModuleScore(fibs, features = subres %>% filter(padj < 0.05 & log2FoldChange < -1) %>% pull(gene) %>% list, name= "down_reg_genes")

FeaturePlot(fibs, features = "down_reg_genes1", max.cutoff = "q90", min.cutoff = "q10")
DotPlot(fibs, features = "down_reg_genes1")






