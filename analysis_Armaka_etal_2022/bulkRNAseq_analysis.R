

```{r}
all_counts <- read.delim("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/data_from_koliaspaper/Bulk_RNAseq/alingment_new/all_counts.txt", comment.char="#")
```


```{r}
#SL

rownames(all_counts) <- all_counts$Geneid
all_counts <- all_counts[,-c(1:6)]

#step 1 create meta data (you will have to adpt this depending on how smaples are named
meta_data_new=colnames(all_counts)
meta_data_new<-as.data.frame(meta_data_new)
library(splitstackshape)
meta_data_new$to_split<-meta_data_new$meta_data_new
meta_data_new<-cSplit(meta_data_new, splitCols = "to_split", sep="_") #check how names are seperated, could be . or _ for eg
meta_data_new <- meta_data_new[,c(1,3)]
meta_data_new$condition <- c(rep("H_LL", 3),rep("H_SL", 3), rep("4wk_LL", 3), rep("4wk_SL", 3), rep("8wk_LL", 3), rep("8wk_SL", 3) )

colnames(meta_data_new)<-c("sample", "time", "condition")  #rename as you want and make sure you have 3 name for 3 cols, 2names for 2 cols etc
meta_data_new$sample<-as.factor(meta_data_new$sample)
meta_data_new$condition<-as.factor(meta_data_new$condition)
meta_data_new$time<-as.factor(meta_data_new$time)


meta_data_new$condition2 <- meta_data_new$condition
meta_data_new<-cSplit(meta_data_new, splitCols = "condition2", sep="_") #check how names are seperated, could be . or _ for eg



meta_data_new_SL <- meta_data_new %>% filter(condition2_2 == "SL")

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = all_counts[,meta_data_new_SL$sample],
                                  colData = meta_data_new_SL,
                                  design = ~condition)

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

meta <- meta_data_new_SL
print(resultsNames(dds))
targetvar <- "condition"
comps <- data.frame(t(combn(unique(as.character(meta[[targetvar]])), 2)))
head(comps)

ress <- apply(comps, 1, function(cp) {
  print(cp)
  res <- data.frame(results(dds, contrast=c(targetvar, cp[1], cp[2])))
  res[["gene"]] <- rownames(res)
  res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
  res
})



res <- Reduce(rbind, ress)

head(res)
comps

res %>% 
  filter(padj < 0.05) %>%
  mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
  arrange(desc(abs(score))) -> subres

head(subres)
dim(subres)

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
ss_sm <- meta_data_new_SL[, c("condition")]
col_ann <- HeatmapAnnotation(df = ss_sm)  


 Heatmap(scale_sub_vsd, 
              top_annotation = col_ann, 
                col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = T,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F)



```

```{r}
meta_data_new_LL <- meta_data_new %>% filter(condition2_2 == "LL")

library(DESeq2)
dds_LL <- DESeqDataSetFromMatrix(countData = all_counts[,meta_data_new_LL$sample],
                                  colData = meta_data_new_LL,
                                  design = ~condition)

mingenecount <- 200
#maxgenecount <- quantile(rowSums(counts(dds)), 0.99)
# Subset low-expressed genes
keep <- rowSums(counts(dds_LL)) > mingenecount #& rowSums(counts(dds)) < maxgenecount
dds_LL <- dds_LL[keep, ]
print(quantile(rowSums(counts(dds_LL))))
dim(dds_LL)

dds_LL@colData[['condition']] <- as.factor(dds_LL@colData[['condition']])

design(dds_LL) <- as.formula(paste0("~", "condition"))

print(design(dds_LL))

dds_LL <- DESeq(dds_LL, test = "Wald")

meta <- meta_data_new_LL

targetvar <- "condition"
comps <- data.frame(t(combn(unique(as.character(meta[[targetvar]])), 2)))
head(comps)

ress <- apply(comps, 1, function(cp) {
  print(cp)
  res <- data.frame(results(dds_LL, contrast=c(targetvar, cp[1], cp[2])))
  res[["gene"]] <- rownames(res)
  res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
  res
})



res_LL <- Reduce(rbind, ress)



res_LL %>% 
  filter(padj < 0.05) %>%
  mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
  arrange(desc(abs(score))) -> subres_LL

library(ComplexHeatmap)

      if(length(unique(subres_LL$gene)) > 10) {
      vsd <- tryCatch({
        vst(dds_LL, blind=TRUE)
      }, error=function(e) {
        message(e)
        print(e)
        return(NULL)
      })
      
      if(!is.null(vsd)) {
        print(dim(assay(vsd)))
        print(head(assay(vsd), 3))
        vsd_mat <- assay(vsd)
        
        feats <- unique(subres_LL$gene)
        print(length(feats))
        
        # Sub-set matrix to relevant features
        sub_vsd_mat <- vsd_mat[rownames(vsd_mat) %in% feats, ]
        scale_sub_vsd_LL <- t(scale(t(sub_vsd_mat)))
        head(scale_sub_vsd)
        dim(scale_sub_vsd)
      }
      }
ss_sm <- meta_data_new_LL[, c("condition")]
col_ann <- HeatmapAnnotation(df = ss_sm)  


 Heatmap(scale_sub_vsd_LL, 
              top_annotation = col_ann, 
                col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = T,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F)


```

```{r}



vsd_matLL_scale <- vsd_matLL %>% scale()
vsd_matSL_scale <- vsd_matSL %>% scale()


vsd_matLL_scale["Runx1",] %>% as.data.frame() %>% rownames_to_column(var="sample") %>% filter(!grepl("unique_wk4_LL", sample))
vsd_matSL_scale["Runx1",] %>% as.data.frame()%>% rownames_to_column(var="sample") %>% filter(!grepl("unique_wk4_SL", sample)) %>%
  cSplit(splitCols = "sample", sep="_") %>% 
ggplot(aes(x=., y=factor(sample_2, levels=c("wk4","wk8")))) +  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+ coord_flip()+ geom_boxplot(width=0.1)+geom_jitter(shape=16, position=position_jitter(0.2))+ scale_fill_grey() + theme_classic()+ggtitle("SL")



vsd_matLL_scale["Runx1",] %>% as.data.frame()%>% rownames_to_column(var="sample") %>% filter(!grepl("unique_wk4_LL", sample)) %>%
  cSplit(splitCols = "sample", sep="_") %>% 
ggplot(aes(x=., y=factor(sample_2, levels=c("wk4","wk8")))) +  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+ coord_flip()+ geom_boxplot(width=0.1)+geom_jitter(shape=16, position=position_jitter(0.2))+ scale_fill_grey() + theme_classic()+ theme_classic()+ggtitle("LL")


```

