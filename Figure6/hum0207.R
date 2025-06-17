hg19.v27.samplefiltered <- read.delim("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/Bulk_human_data/hum0207_RNAseq/hum0207.v1.RNA.v1/RNAseq/hg19.v27.samplefiltered.htseq")

hg19.v27.samplefiltered<-hg19.v27.samplefiltered  %>% distinct(gname, .keep_all = TRUE)

meta.data <- read.delim("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/Bulk_human_data/hum0207_RNAseq/hum0207.v1.RNA.v1/RNAseq/meta.data.txt", header=FALSE)


hg19.v27.samplefiltered_all<-hg19.v27.samplefiltered[,colnames(hg19.v27.samplefiltered) %in% meta.data$V1,]

rownames(hg19.v27.samplefiltered_all)=hg19.v27.samplefiltered$gname

hg19.v27.samplefiltered_all$ENSGID<-NULL
hg19.v27.samplefiltered_all$gname<-NULL

library(splitstackshape)
meta_data_dds=colnames(hg19.v27.samplefiltered_all)
meta_data_dds<-as.data.frame(meta_data_dds)
meta_data_dds<-cSplit(meta_data_dds, "meta_data_dds", sep="_", type.convert=FALSE)
meta_data_dds$sample<-colnames(hg19.v27.samplefiltered_all)
meta_data_dds<- meta_data_dds %>% remove_rownames %>% column_to_rownames(var="sample")
meta_data_dds$meta_data_dds_1<-NULL
colnames(meta_data_dds)<-('group')



dds <- DESeqDataSetFromMatrix(countData=hg19.v27.samplefiltered_all, 
                              colData=meta_data_dds, 
                              design=~group, tidy = F)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


dds$group <- relevel(dds$group, ref = "NS")

dds <- DESeq(dds)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)


deseq2Results <- results(dds)
deseq2ResDF <- as.data.frame(deseq2Results)
summary(deseq2Results)
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)
#strict
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .01,])


deseq2VST <- vst(dds)
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

library(ComplexHeatmap)
Heatmap(deseq2VST,cluster_columns =F)

library(DESeq2)
resultsNames(dds)
IL6_vs_NS <- lfcShrink(dds, coef="group_IL6_vs_NS", type="apeglm")
IL6_vs_NS<-as.data.frame(IL6_vs_NS)
IL6_vs_NS$gene<-rownames(IL6_vs_NS)
IL6_vs_NS<-IL6_vs_NS[order(IL6_vs_NS$padj),]
#IL6_vs_NS<-IL6_vs_NS[IL6_vs_NS$log2FoldChange > 0.2,]
#write.csv(IL6_vs_NS, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/IL6_vs_NS.csv")

#write.csv(IL6_vs_NS, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/zone_annotationIL6_vs_NS_sorted.csv")

TNFa_vs_NS <- lfcShrink(dds, coef="group_TNFa_vs_NS", type="apeglm")
TNFa_vs_NS<-as.data.frame(TNFa_vs_NS)
TNFa_vs_NS$gene<-rownames(TNFa_vs_NS)
TNFa_vs_NS<-TNFa_vs_NS[order(TNFa_vs_NS$padj),]
#TNFa_vs_NS<-TNFa_vs_NS[TNFa_vs_NS$log2FoldChange > 0.2,]
#write.csv(TNFa_vs_NS, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/TNFa_vs_NS.csv")


resultsNames(dds)
ALL_vs_NS <- lfcShrink(dds, coef="group_all_vs_NS", type="apeglm")
ALL_vs_NS<-as.data.frame(ALL_vs_NS)


Il1b_vs_NS <- lfcShrink(dds, coef="group_IL1b_vs_NS", type="apeglm")
Il1b_vs_NS<-as.data.frame(Il1b_vs_NS)
Il1b_vs_NS$gene<-rownames(Il1b_vs_NS)


IFNa_vs_NS <- lfcShrink(dds, coef="group_IFNa_vs_NS", type="apeglm")
IFNa_vs_NS<-as.data.frame(IFNa_vs_NS)
IFNa_vs_NS$gene<-rownames(IFNa_vs_NS)
#normalized_counts["TSPAN6",]

#write.csv(IFNa_vs_NS, "/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/zone_annotation/IFNa_vs_NS.csv")


library(EnhancedVolcano)

genes=c( "CTHRC1", "RUNX1", "MMP14", "IGF1")

EnhancedVolcano(Il1b_vs_NS,
    lab = rownames(Il1b_vs_NS),
    x = 'log2FoldChange',
    y = 'padj',
        selectLab = "",
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



genes=c( "CTHRC1", "RUNX1")

EnhancedVolcano(TNFa_vs_NS,
    lab = rownames(TNFa_vs_NS),
    x = 'log2FoldChange',
    y = 'padj',
        selectLab = "",
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


genes=c( "CTHRC1", "RUNX1", "IGF1", "MMP14")

EnhancedVolcano(ALL_vs_NS,
    lab = rownames(ALL_vs_NS),
    x = 'log2FoldChange',
    y = 'padj',
        selectLab = genes,
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



genes=c( "CTHRC1", "RUNX1", "IGF1", "MMP14")

EnhancedVolcano(IL6_vs_NS,
    lab = rownames(IL6_vs_NS),
    x = 'log2FoldChange',
    y = 'padj',
        selectLab = "",
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



normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts<-as.data.frame(normalized_counts)
normalized_countst<-t(normalized_counts)
normalized_countst<-as.data.frame(normalized_countst)
normalized_countst$group<-meta_data_dds$group


ggplot(normalized_countst, aes(x=RUNX1, y=factor(group, levels=c("NS","TNFa", "IFNg","IFNa", "IL1b", "IL17", "IL18", "TGFb", "IL6", "all")))) +   geom_violin(trim=F,)+ coord_flip()+ geom_boxplot(width=0.1)+geom_jitter(shape=16, position=position_jitter(0.2))+ scale_fill_grey() + theme_classic()

subres %>% filter(gene == "RUNX1")




```


```{r}
#mingenecount <- quantile(rowSums(counts(dds)), 0.5)
mingenecount <- 200
maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount & rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

dds@colData[['group']] <- factor(dds@colData[['group']])


library(genefilter)

top_genes<- genefilter::rowVars(normalized_counts) %>% 
  sort(decreasing = TRUE) %>%
  names() %>%
  head(5000)

top_genes[grep("RUNX1", top_genes)]
top_genes[grep("CTHRC1", top_genes)]
top_genes[grep("MMP14", top_genes)]
top_genes[grep("IGF1", top_genes)]

dds_sub <- dds[rownames(dds) %in% top_genes]

design(dds_sub) <- formula(~ group)
print(design(dds_sub))
dds <- DESeq(dds_sub, test = "Wald")


targetvar <- "group"

comps1 <- data.frame(t(combn(unique(as.character(meta_data_dds[[targetvar]])), 2)))
      head(comps1)
      
      ress <- apply(comps1, 1, function(cp) {
        print(cp)
        res <- data.frame(results(dds_sub, contrast=c(targetvar, cp[1], cp[2])))
        res[["gene"]] <- rownames(res)
        res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
        res
      })
      
      res1 <- Reduce(rbind, ress)
```


```{r}
colnames <- colnames(hg19.v27.samplefiltered_all) %>% as.data.frame()
colnames(colnames) <- 'samples'
library(splitstackshape)
colnames <- cSplit(colnames, splitCols="samples", sep="_")
colnames$samples_3 <- colnames$samples_2


        colnames %>%
          mutate('samples_2' = ifelse(samples_2 == 'NS', 'S0-NS',
                                      ifelse(samples_2 == 'IFNa', 't1-IFNa',
                                             ifelse(samples_2 == 'IFNg', 't3-IFNg',
                                                    ifelse(samples_2 == 'IL1b', 't4-IL1b',
                                                           ifelse(samples_2 == 'IL6', 't5-IL6',
                                                                  ifelse(samples_2 == 'TGFb', 't6-TGFb',
          ifelse(samples_2 == 'TNFa', 't6-TNFa', 
                 ifelse(samples_2 == 'IL17', 't7-IL17', 
                        ifelse(samples_2 == 'IL18', 't8-IL18',    
                               ifelse(samples_2 == 'all', 't9-all',    
                                      'WOH!!!'))))))))))) %>%
          arrange(samples_2) %>%
          data.frame -> colnames_new

colnames_new$final <- paste(colnames_new$samples_1, colnames_new$samples_3, sep="_")


res1 %>% 
      filter(padj < 0.05) %>%
      mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
      arrange(desc(abs(score))) -> subres



```


```{r}
library(ComplexHeatmap)

if(length(unique(subres$gene)) > 10) {
      vsd <- tryCatch({
        vst(dds_sub, blind=TRUE)
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


col_ann <- HeatmapAnnotation(df = colnames_new[,"samples_3"])
scale_sub_vsd<-scale_sub_vsd[ , colnames_new$final]



    draw(
      Heatmap(scale_sub_vsd, 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F))
    
    
 




```
```{r}
topedges <- 0.05
    ggnet <- Rfast::cora(t(scale_sub_vsd))
    ggnet_gather <- reshape2::melt(ggnet, id.vars = "V1")
    print(dim(ggnet_gather))
    
    ggnet_gather %>%
      filter(value > 0) %>%
      filter(Var1 != Var2) -> ggnet_gather
    
    print(dim(ggnet_gather))
    
    if(nrow(ggnet_gather) > 0) {
      
      edges <- arrange(ggnet_gather, -value)
      print(head(edges))
      print(tail(edges))
      edges <- edges[seq(1, nrow(edges), by = 2), ]
      print(head(edges))
      print(tail(edges))
      qtop <- quantile(edges$value, 1-topedges)
      print(qtop)
      
      edges %>%
        filter(value > qtop) -> top_edges
      
      print(dim(top_edges))
      
      top_edges %>%
        mutate('idx' = paste(Var1, Var2, sep = "_")) -> top_edges# %>%
      #dplyr::distinct(idx, .keep_all = TRUE) -> top_edges
      
      colnames(top_edges) <- c("Source", "Target", "Weight", "Id")
      #top_edges %>%
      #  select(c(Id, Source, Target, Weight)) -> top_edges
      
      # - community detection
      # - igraph definition
      g <- igraph::graph_from_data_frame(top_edges[, c("Source", "Target")],
                                         directed = FALSE)
      g <- igraph::set_edge_attr(g, "weight", value = top_edges$Weight)
      g <- igraph::set_edge_attr(g, "name", value = top_edges$Id)  
      # leiden
      leiden_mod <- igraph::cluster_leiden(g, objective_function = "modularity")
      mods <- data.frame(cbind(igraph::V(g)$name, leiden_mod$membership))
      colnames(mods) <- c("Id", "leiden_gene_cluster")
      imods <- names(table(mods$leiden_gene_cluster)[table(mods$leiden_gene_cluster) > 10])
      print(imods)
      
      if(length(imods) > 0) {
        mods <- filter(mods, leiden_gene_cluster %in% imods)
        head(mods)
        
        row_cl <- data.frame('gene' = mods$Id,
                             'gene_cluster' = paste0("K", mods$leiden_gene_cluster))
        
        row_cl %>%
          arrange(gene_cluster) %>%
          data.frame -> row_cl
        
        rownames(row_cl) <- row_cl$gene
        
              }
    }
    
     row_ann <- rowAnnotation(df = row_cl[, -1]) 
     
     
    
     Heatmap(scale_sub_vsd[rownames(row_cl),], 
              top_annotation = col_ann, 
              right_annotation = row_ann,
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = F,
              show_row_names = F,
              show_column_names = F)
    
     
row_cl %>% filter(gene == "RUNX1")     
