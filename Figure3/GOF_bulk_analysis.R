

```{r}

counts2 <- counts_all %>% select(-c(Chr, Start,End, Strand, Length))


counts2_F <- counts2 %>% select(c(1:2, 4:5, 7:8, 10:12))
colnames(counts2_F) <- c("R1A1", "R1A2", "R1C1", "R1C2", "EV1", "EV2", "EV4", "R1A4", "R1C4")

counts2_F <- counts2_F[c(5,6,7,1,2,8,3,4,9)]

meta_data=colnames(counts2_F)
meta_data<-as.data.frame(meta_data)

condition <- c(rep("EV", 3), rep("R1A", 3), rep("R1C", 3))
meta_data$condition <- condition
colnames(meta_data)<-c("sample", "condition")


```
```{r}
meta_data$sample<-as.factor(meta_data$sample)
meta_data$disease<-as.factor(meta_data$condition)

library(sva)
batch <- c(1,1,2,1,1,2,1,1,2)


adjusted <- ComBat_seq(counts2_F, batch=batch, group=NULL)

dds <- DESeqDataSetFromMatrix(countData = adjusted,
                                  colData = meta_data,
                                  design = ~1)

print(quantile(rowSums(counts(dds))))

mingenecount <- quantile(rowSums(counts(dds)), 0.5)
maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount & rowSums(counts(dds)) < maxgenecount
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
              show_column_names = F)
```
```{r}

annotation_gs <- fetchAnnotation(species="hs", ensembl_version=NULL, ensembl_host=NULL)


index <- match(subres$gene, annotation_gs$ensembl_id)
subres$gene_name <- annotation_gs$gene_name[index]




```



```{r}

#plot genes in RUNX1 GRN from mouse stia

grn_runx1_human <- grn %>% filter(Source== "RUNX1")


library(babelgene)

human_for_heatmap <- orthologs(genes = row_cl_k3$gene, species = "mouse", human = F)
heatmap_rows <-  orthologs(genes = rownames(vsd_mat), species = "mouse")

vsd_mat_new <- vsd_mat

subres_new <- res1[rownames(res1) %in% human_for_heatmap$human_ensembl, ]

  
sub_vsd_mat2 <- vsd_mat[rownames(vsd_mat) %in% human_for_heatmap$human_ensembl, ]
        sub_vsd_mat2 <- t(scale(t(sub_vsd_mat2)))  
  
  res1 %>% 
      filter(padj < 0.05) %>%
      mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
      arrange(desc(abs(score))) -> subres

    
  
col_ann <- HeatmapAnnotation(df = meta_data)
Heatmap(sub_vsd_mat2, 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F)




```






```{r}

unique(subres$comparison)

max(subres$padj)
subres_f <- subres %>% filter(comparison == "EV_vs_R1C" & log2FoldChange > 2 &padj < 0.01)
write.csv(subres_f, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/runx1_oversplill/lentiviral_bulk/DEGs_human_RA.csv")



df1 <- str_split("IFI6,IFI16,IFI27,IFIT2,IFIT1,IRF1,IRF5,IRF7,HLA-B,HLA-C,HLA-F,CXCL11,CXCL14,IL23A,CCL1,IL41L,BMP2,BMP7,EGF,IGF1,CXCL10", ",") %>% as.data.frame()
colnames(df1) <- 'col1'

subres_df1 <- res1[res1$gene_name %in% df1$col1,]
#subres_df1 <- subres_df1 %>% filter(comparison == "EV_vs_R1C")

df2 <-str_split("MMP9,MMP8,MMP13,ADAMTS4,ADAM21,ADAMTSL2,ADAM9,COL2A1,COL10A1,COL25A1,COL24A1,MMP14,PDGFRA", ",") %>% as.data.frame()

colnames(df2) <- 'col1'

subres_df2 <- res1[res1$gene_name %in% df2$col1,]
subres_df2 <- subres_df2 %>% filter(comparison == "EV_vs_R1C")

scale_sub_vsd_new <- t(scale(t(vsd_mat)))
scale_sub_vsd_new <- scale_sub_vsd_new[rownames(scale_sub_vsd_new) %in% subres_df1$gene,]
scale_sub_vsd_new <- scale_sub_vsd_new %>% as.data.frame()
index <- match(rownames(scale_sub_vsd_new), annotation_gs$ensembl_id)
scale_sub_vsd_new$gene_name <- annotation_gs$gene_name[index]
rownames(scale_sub_vsd_new) <- scale_sub_vsd_new$gene_name
scale_sub_vsd_new$gene_name <- NULL

Heatmap(scale_sub_vsd_new, 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 8), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = T,
              show_column_names = F, border=T)


scale_sub_vsd_new <- t(scale(t(vsd_mat)))
scale_sub_vsd_new <- scale_sub_vsd_new[rownames(scale_sub_vsd_new) %in% subres_df2$gene,]
scale_sub_vsd_new <- scale_sub_vsd_new %>% as.data.frame()
index <- match(rownames(scale_sub_vsd_new), annotation_gs$ensembl_id)
scale_sub_vsd_new$gene_name <- annotation_gs$gene_name[index]
rownames(scale_sub_vsd_new) <- scale_sub_vsd_new$gene_name
scale_sub_vsd_new$gene_name <- NULL


Heatmap(scale_sub_vsd_new, 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 8), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = T,
              show_column_names = F, border=T)









```


