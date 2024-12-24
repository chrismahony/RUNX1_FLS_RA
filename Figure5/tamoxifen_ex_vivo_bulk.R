



```{r}

counts2 <- counts_all_final_n3 %>% dplyr::select(-c(Chr, Start,End, Strand, Length))


colnames(counts2) <- c("Geneid", "CTRL1", "CTRL3", "CTRL5", "TAM1", "TAM3", "TAM5")

rownames(counts2) <- counts2$Geneid

counts2$Geneid <- NULL

meta_data=colnames(counts2)
meta_data<-as.data.frame(meta_data)

condition <- c(rep("CTRL", 3), rep("TAM", 3))
meta_data$condition <- condition
colnames(meta_data)<-c("sample", "condition")


```
```{r}
meta_data$sample<-as.factor(meta_data$sample)
meta_data$condition<-as.factor(meta_data$condition)

library(sva)
batch <- c(1,1,2,1,1,2)


adjusted <- ComBat_seq(counts2, batch=batch, group=NULL)

dds <- DESeqDataSetFromMatrix(countData = adjusted,
                                  colData = meta_data,
                                  design = ~1)

print(quantile(rowSums(counts(dds))))

mingenecount <- 200
#maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

dds@colData[['condition']] <- factor(dds@colData[['condition']])


design(dds) <- formula(~ condition)
print(design(dds))
dds <- DESeq(dds, test = "Wald")
plotDispEsts(dds)

#PCA analysis
rld_v3 <- rlog(dds, blind=TRUE)
plotPCA(rld_v3, intgroup="condition")

pcaData <- plotPCA(rld_v3, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(aes(shape=condition, size=3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


rld_mat_v3 <- assay(rld_v3)    
rld_cor_v3 <- cor(rld_mat_v3)
library(pheatmap)
pheatmap(rld_cor_v3)

targetvar <- "condition"

comps1 <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))
      head(comps1)
      
      ress <- apply(comps1, 1, function(cp) {
        print(cp)
        res <- data.frame(DESeq2::results(dds, contrast=c(targetvar, cp[1], cp[2])))
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

col_ann <- HeatmapAnnotation(df = meta_data)
Heatmap(scale_sub_vsd, 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F, border=T)



runx_grn <- grn %>% filter(Source == "RUNX1")

runx_grn$Target <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(runx_grn$Target), perl=TRUE)


scale_sub_vsd_f <- scale_sub_vsd[rownames(scale_sub_vsd) %in% runx_grn$Target,]

Heatmap(scale_sub_vsd_f, 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F)


```



```{r}

annotation_gs <- fetchAnnotation(species="mm", ensembl_version=NULL, ensembl_host=NULL)


# add 'cluster' column to DEGs ###

subres$cluster <- "NO"
subres$cluster[subres$log2FoldChange > 0.001] <- "UP"
subres$cluster[subres$log2FoldChange < -0.001] <- "DOWN"


library(gsfisher)

FilteredGeneID <- unique(res1$gene)
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]
ensemblUni <- na.omit(ensemblUni)


index <- match(subres$gene, annotation_gs$gene_name)
subres$ensembl <- annotation_gs$ensembl_id[index]
subres <- na.omit(subres)


go.results <- runGO.all(results=subres,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="padj", p_threshold=0.05,
                  species = "mm")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=5, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description",selected_genesets = c(), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)

sampleEnrichmentDotplot(go.results, selection_col = "description", selected_genesets = c("cell chemotaxis", "positive regulation of chemotaxis", "extracellular matrix organization", "phagocytosis", "fatty acid transport", "positive regulation of autophagy"), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)

go.results %>% filter(description == "cell chemotaxis")


```


```{r}

genes=str_split("Adamts1,Col11a1,Col14a1,Col15a1,Col2a1,Col3a1,Col5a1,Col5a2,Col1a1,Col1a2,Mmp16,Adamts5,Mmp17,Mmp23,Postn,Col24a1,Adamts9,Adamts6,Adamts2,Adamts12,Adamts4,Adamtsl3,Col8a2,Ccl7,Cxcl5,Cx3cl1,Cxcl14, Cxcl3, Wnt5a", ",")

subres_plot <- subres

subres_plot$log2FoldChange <- subres_plot$log2FoldChange*-1

EnhancedVolcano(subres_plot,
    lab = subres_plot$gene,
    x = 'log2FoldChange',
    y = 'padj',
        selectLab = "",
    title = 'title',
    subtitle = "GEX, red=p_adj<0.05 & FC > 0.25",
    pCutoff = 0.05,
    FCcutoff = 0.5,
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
    borderWidth = 0.5, max.overlaps = 100
    )
```


```{r}


ECMgenes=str_split("Adamts1,Col11a1,Col14a1,Col15a1,Col2a1,Col3a1,Col5a1,Col5a2,Col1a1,Col1a2,Has2,Lum,Mmp16,Adamts5,Mmp17,Mmp23,Postn,Col24a1,Adamtsl1,Adamts9,Adamts6,Adamts2,Adamts12,Adamts4,Adamtsl3,Col8a2,Col27a1", ",") %>% as.data.frame()


CYtokinesgenes=str_split("Ccl9,Ccl7,Ccl20,Cxcl14,Cxcl3,Cxcl5,Cx3cl1,Igf1,Il6st,Il15", ",")%>% as.data.frame()



col_ann <- HeatmapAnnotation(df = meta_data)
Heatmap(scale_sub_vsd[ECMgenes$c..Adamts1....Col11a1....Col14a1....Col15a1....Col2a1....Col3a1...,], 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 6), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = T,
              show_column_names = F, border=T)

Heatmap(scale_sub_vsd[CYtokinesgenes$c..Ccl9....Ccl7....Ccl20....Cxcl14....Cxcl3....Cxcl5....Cx3cl1...,], 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 12), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = T,
              show_column_names = F, border=T)

Heatmap(scale_sub_vsd[c("Ccl9"),], 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 6), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = T,
              show_column_names = F, border=T)


subres %>% write.csv("/rds/projects/m/mahonyc-runx1-bulk-seq-data/lentivrial_GOF/analysis_bact_corrected/mouse_tamox_n3_DEGs.csv")

```



```{r}
data <-table(subres$cluster) %>% as.data.frame()

ggplot(data, aes(fill=Var1, y=Freq, x=Var1, label = Freq)) + 
    geom_bar(position="stack", stat="identity")+theme_ArchR() +
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
    scale_fill_viridis(discrete = T)+coord_flip()+scale_fill_manual(values = c("lightblue","red"))


data %>% filter(condition== "inhib") %>% ggplot(aes(y=Freq, x=direction, fill=direction)) + 
    geom_bar(position="stack", stat="identity")+theme_ArchR() +
      scale_fill_viridis(discrete = T)+coord_flip()+scale_fill_manual(values = c("lightblue","red"))


data %>% ggplot(aes(y=Freq, x=direction, fill=direction))  + 
        geom_bar(stat = 'identity', position = position_dodge(), aes(group = direction),width = 1) + 
        coord_flip() + 
        facet_wrap(~condition, scales = 'free_y', nrow=2, ncol=1)+ 
        theme(
            axis.text.x = element_text(angle = 45, hjust=1),
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
            # strip.text = element_blank()
        ) + 
        guides(color = 'none', fill = 'none') + 
        labs(y = '# genes')+
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = c(0, 0)))+theme_ArchR()+scale_fill_manual(values = c("blue","red"))+
  theme(strip.background = element_rect(fill="white", size=1, color="white"))

 
```

