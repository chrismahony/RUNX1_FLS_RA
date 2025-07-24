rm(list = setdiff(ls(), "stia2021_rna"))

stia2021_rna$sample_condition<-paste(stia2021_rna$sample_id, stia2021_rna$condition, sep=".")

Idents(stia2021_rna)<-'cluster.name'
stia2021_rna_SLfibs<-subset(stia2021_rna, idents=levels(stia2021_rna)[-c(1,2,3,5,12)])

Idents(stia2021_rna_SLfibs)<-'experiment'
stia2021_rna_SLfibs<-subset(stia2021_rna_SLfibs, idents=levels(stia2021_rna_SLfibs)[2])


cts<-AggregateExpression(stia2021_rna_SLfibs, group.by = c("sample_condition"), assays = "RNA", slot = "counts", return.seurat = F)

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

dds@colData[['condition']] <- factor(dds@colData[['condition']],
                                     levels = c("control",
                                                "initiation",
                                                "peak",
                                                "resolving",
                                                "resolved",
                                                "persistent"))

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

meta_data %>%
          mutate('condition' = ifelse(condition == 'control', 't0-control',
                                      ifelse(condition == 'initiation', 't1-initiation',
                                             ifelse(condition == 'peak', 't3-peak',
                                                    ifelse(condition == 'resolving', 't4-resolving',
                                                           ifelse(condition == 'resolved', 't5-resolved',
                                                                  ifelse(condition == 'persistent', 't6-persistent', 'WOH!!!'))))))) %>%
          arrange(condition) %>%
          data.frame -> sub_meta



ss_sm <- sub_meta[, c("sample", "condition")]
ss_sm <- sub_meta[, c("condition")]

col_ann <- HeatmapAnnotation(df = ss_sm)


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
      
      print(head(top_edges))
      print(dim(top_edges))
      
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
  
colnames(scale_sub_vsd)



    
    


#scale_sub_vsd <- subset(scale_sub_vsd, select=c("control_cd45n_s1_stia2021_fibroblast" ,"control_cd45n_s2_stia2021_fibroblast" #,"control_cd45n_s3_stia2021_fibroblast",
#"day1_cd45n_s1_stia2021_fibroblast"   , "day1_cd45n_s2_stia2021_fibroblast"   , "day1_cd45n_s3_stia2021_fibroblast"  ,
#"day8_cd45n_s1_stia2021_fibroblast"   , "day8_cd45n_s2_stia2021_fibroblast"   , "day8_cd45n_s3_stia2021_fibroblast",
#"day15_cd45n_s1_stia2021_fibroblast"  , "day15_cd45n_s2_stia2021_fibroblast"  , "day15_cd45n_s3_stia2021_fibroblast",  
#  "day22_cd45n_s1_stia2021_fibroblast",   "day22_cd45n_s2_stia2021_fibroblast",   "day22_cd45n_s3_stia2021_fibroblast",  
# "day28_cd45n_s1_stia2021_fibroblast" ,  "day28_cd45n_s2_stia2021_fibroblast"  , "day28_cd45n_s3_stia2021_fibroblast"  ))


scale_sub_vsd <- subset(scale_sub_vsd, select=c("control_cd45n_s1_stia2021.control",  "control_cd45n_s2_stia2021.control",  "control_cd45n_s3_stia2021.control" ,
"day1_cd45n_s1_stia2021.initiation" , "day1_cd45n_s2_stia2021.initiation" , "day1_cd45n_s3_stia2021.initiation", 
 "day8_cd45n_s1_stia2021.peak"       , "day8_cd45n_s2_stia2021.peak"      ,  "day8_cd45n_s3_stia2021.peak" ,
"day15_cd45n_s1_stia2021.resolving" , "day15_cd45n_s2_stia2021.resolving" , "day15_cd45n_s3_stia2021.resolving",  
"day22_cd45n_s1_stia2021.resolved"  , "day22_cd45n_s2_stia2021.resolved"  , "day22_cd45n_s3_stia2021.resolved"  ,
 "day28_cd45n_s1_stia2021.persistent", "day28_cd45n_s2_stia2021.persistent", "day28_cd45n_s3_stia2021.persistent"))
    
library(colorRamp2)        
if(exists("row_cl")) {
    row_ann <- rowAnnotation(df = row_cl[, -1])    
    
    draw(
      Heatmap(scale_sub_vsd[rownames(row_cl),], 
              top_annotation = col_ann,
              right_annotation = row_ann,
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = FALSE,
              show_row_names = F,
              show_column_names = F,
              border=T))
}


library(dplyr)
library(readr)
grn <- read_csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/df_grn2_new_FINAL.csv")
grn_Runx1 <- grn %>% filter(tf == "RUNX1")
grn_Runx1$Target <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(grn_Runx1$gene), perl=TRUE)

scale_sub_vsd_f <- scale_sub_vsd[rownames(scale_sub_vsd) %in% grn_Runx1$Target,]
Heatmap(scale_sub_vsd_f, 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F,
              border=T)



library(gsfisher)
#annotation_gs<-fetchAnnotation(species = "mm")

res1 <- res
index <- match(res1$gene, row_cl$gene)
res1$cluster <- row_cl$gene_cluster[index]
res1_cleaned<-na.omit(res1)

index <- match(res1_cleaned$gene, annotation_gs$gene_name)
res1_cleaned$ensembl <- annotation_gs$ensembl_id[index]



FilteredGeneID <- unique(res1$gene)
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]
ensemblUni <- na.omit(ensemblUni)
res1_cleaned<-na.omit(res1_cleaned)

go.results <- runGO.all(results=res1_cleaned,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="padj", p_threshold=0.05,
                  species = "mm")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=5, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)


stia2021_rna_only <- stia2021_rna_SLfibs

stia2021_rna_only$cluster_sample_condition <- paste(stia2021_rna_only$Key.marker.genes, stia2021_rna_only$sample_id, stia2021_rna_only$condition)
Idents(stia2021_rna_only) <- 'cluster_sample_condition'

stia2021_rna_only <- AddModuleScore(stia2021_rna_only, features=list(row_cl %>% filter(gene_cluster == "K1") %>% rownames()), name="K1_mod")

stia2021_rna_only <- stia2021_rna_only %>% ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA()

levels(stia2021_rna_only)

dotplot <- DotPlot(stia2021_rna_only, features= "K1_mod1")

dotplot <- dotplot$data

library(splitstackshape)
dotplot <- cSplit(dotplot, splitCols = "id", sep="_")

dotplot <- cSplit(dotplot, splitCols = "id_1", sep=" ")



dotplot$id_4 <- gsub("stia2021 control", "0", dotplot$id_4)
dotplot$id_4 <- gsub("stia2021 initiation", "3", dotplot$id_4)
dotplot$id_4 <- gsub("stia2021 peak", "7", dotplot$id_4)
dotplot$id_4 <- gsub("stia2021 resolving", "15", dotplot$id_4)
dotplot$id_4 <- gsub("stia2021 resolved", "22", dotplot$id_4)
dotplot$id_4 <- gsub("stia2021 persistent", "28", dotplot$id_4)
dotplot$id_4 <- as.double(dotplot$id_4)


dotplot$id_1_1 <- gsub("C1qtnf3", "SL_Col8a1", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Ccl11", "SL_Ccl11", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Ccl7", "SL_Ccl2", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Col23a1", "SL_Col23a1", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Fmo2", "SL_Fmo2", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Pi16", "SL_Pi16", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Serpina3c", "SL_C3", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Sfrp1", "SL_Cfb", dotplot$id_1_1)

unique(stia2021_rna_only$cluster.name)
dotplot %>% 
ggplot(aes(x=id_4, y=avg.exp.scaled)) + 
  geom_point(shape = 21)+
  geom_smooth(alpha = .6, fill="#ffb09c", color="black")+RotatedAxis()+theme(
            axis.line = element_line(),panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
        geom_vline(xintercept = c(0,3,7,15,22,28), linetype = 2, color = 'grey')+ 
        guides(color = FALSE, fill = FALSE) + 
        scale_x_continuous(breaks = c(0,3,7,15,22,28))+theme_ArchR()+facet_wrap("id_1_1")


dotplot %>% 
  filter(!id_1_1 %in% c("SL_Bglap", "SL_Clu")) %>%  
  ggplot(aes(x = factor(id_4), y = avg.exp.scaled)) + 
  geom_boxplot(fill = "#ffb09c", color = "black", outlier.shape = NA, alpha = 0.8) + 
  RotatedAxis() + 
  theme(
    axis.line = element_line(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  guides(color = FALSE, fill = FALSE) +
  theme_ArchR() +
  facet_wrap(~factor(id_1_1, levels = c("SL_Ccl2", "SL_Col8a1", "LL_Col22a1", "SL_Cfb",
                                        "SL_Chodl", "SL_Col23a1", "SL_Pi16", "SL_C3", 
                                        "SL_Ccl11", "SL_Fmo2")), nrow = 2)


