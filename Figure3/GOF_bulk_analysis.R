


counts_all <- read.delim("/rds/projects/m/mahonyc-runx1-bulk-seq-data/lentivrial_GOF/count2_star/counts_all.txt", comment.char="#", row.names = 1)

counts2_F <- counts_all %>% select(-c(Chr, Start,End, Strand, Length))


colnames(counts2_F) <- c("EV1", "EV2", "EV4","R1A1", "R1A2","R1A4","R1C1", "R1C2",   "R1C4" )

meta_data=colnames(counts2_F)
meta_data<-as.data.frame(meta_data)

condition <- c(rep("EV", 3), rep("R1A", 3), rep("R1C", 3))
meta_data$condition <- condition
colnames(meta_data)<-c("sample", "condition")



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

subres <- subres %>%
  mutate(direction = if_else(log2FoldChange < 0, "up", "down"))

subres$cluster <- paste(subres$comparison, subres$direction, sep="_")


data <-table(subres$cluster) %>% as.data.frame()
data <- data[-c(5,6),]
data$direction <- c("down", "up", "down", "up")


ggplot(data, aes(fill=direction, y=Freq, x=Var1, label = Freq)) + 
    geom_bar(position="stack", stat="identity")+theme_ArchR() +
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
    scale_fill_viridis(discrete = T)+coord_flip()+scale_fill_manual(values = c("lightblue","red"))


ggplot(data, aes(fill=direction, y=Freq, x=Var1, label = Freq)) + 
    geom_bar(position="stack", stat="identity")+theme_ArchR() +
    scale_fill_viridis(discrete = T)+coord_flip()+scale_fill_manual(values = c("blue","red"))

df1 <- str_split("IFI6,IFI16,IFI27,IFIT2,IFIT1,IRF1,IRF5,IRF7,HLA-B,HLA-C,HLA-F,CXCL11,CXCL14,IL23A,CCL1,IL41L,BMP2,BMP7,EGF,IGF1,CXCL10", ",") %>% as.data.frame()
colnames(df1) <- 'col1'

subres_df1 <- res1[res1$gene %in% df1$col1,]
#subres_df1 <- subres_df1 %>% filter(comparison == "EV_vs_R1C")

df2 <-str_split("MMP9,MMP8,MMP13,ADAMTS4,ADAM21,ADAMTSL2,ADAM9,COL2A1,COL10A1,COL25A1,COL24A1,MMP14,PDGFRL", ",") %>% as.data.frame()
colnames(df2) <- 'col1'

subres_df2 <- res1[res1$gene %in% df2$col1,]
subres_df2 <- subres_df2 %>% filter(comparison == "EV_vs_R1C")

scale_sub_vsd_new <- t(scale(t(vsd_mat)))
scale_sub_vsd_new <- scale_sub_vsd_new[rownames(scale_sub_vsd_new) %in% subres_df1$gene,]
scale_sub_vsd_new <- scale_sub_vsd_new %>% as.data.frame()
#index <- match(rownames(scale_sub_vsd_new), annotation_gs$ensembl_id)
#scale_sub_vsd_new$gene_name <- annotation_gs$gene_name[index]
#rownames(scale_sub_vsd_new) <- scale_sub_vsd_new$gene_name
#scale_sub_vsd_new$gene_name <- NULL

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
#index <- match(rownames(scale_sub_vsd_new), annotation_gs$ensembl_id)
#scale_sub_vsd_new$gene_name <- annotation_gs$gene_name[index]
#rownames(scale_sub_vsd_new) <- scale_sub_vsd_new$gene_name
#scale_sub_vsd_new$gene_name <- NULL


Heatmap(scale_sub_vsd_new, 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 8), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = T,
              show_column_names = F, border=T)



write.csv(subres, "/rds/projects/m/mahonyc-runx1-bulk-seq-data/lentivrial_GOF/count2_star/analysis/DEGs_GOF.csv")






subres_up <- subres %>% filter(log2FoldChange < 0.1 & comparison == "EV_vs_R1C" & gene != "RUNX1")

NatMed_AliverniniEtAl_2020_FLS <- AddModuleScore(NatMed_AliverniniEtAl_2020_FLS, features=list(subres_up$gene), name = "fibs_up")


NatMed_AliverniniEtAl_2020_FLS$group_cluster<-paste(NatMed_AliverniniEtAl_2020_FLS$group, NatMed_AliverniniEtAl_2020_FLS$clusters,  sep=".")
Idents(NatMed_AliverniniEtAl_2020_FLS) <- 'group_cluster'

dotplot <- DotPlot(NatMed_AliverniniEtAl_2020_FLS,features=c("fibs_up1", "RUNX1"))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ scale_size(range = c(2, 8))+RotatedAxis()


dotplot<-dotplot$data

dotplot<-dotplot %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

dotplot <- dotplot[c(1,2,4,6,8,10,3,5,7,9,11)]


row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)


levels(NatMed_AliverniniEtAl_2020_FLS) %>% as.data.frame()


df <- levels(NatMed_AliverniniEtAl_2020_FLS) %>% as.data.frame()
#df$condition <- c("A", "R", "A", "R", "A", "R", "A", "R", "A", "R")

colours <- list('sample' = ArchR::paletteDiscrete(NatMed_AliverniniEtAl_2020_FLS@meta.data[, "group_cluster"]))

colnames(df) <- 'group_cluster'
                  
col_ann <- HeatmapAnnotation(df = df, col=colours)
means <- colMeans(dotplot) %>% as.data.frame()
colnames(means) <- 'mean'

column_ha = HeatmapAnnotation(bar1 = anno_barplot(means), height  = unit(3, "cm"))

Heatmap(dotplot, top_annotation = column_ha, cluster_columns = F)






