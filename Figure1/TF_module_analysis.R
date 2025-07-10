grn <- read.csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/df_grn2_new_FINAL.csv")

grn$Target <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(grn$gene), perl=TRUE)

grn$Source <- grn$tf

tfs <- unique(grn$Source)
grn_list <- list()
for (i in 1:length(tfs)){
grn_list[[i]] <- grn[grn$Source == tfs[[i]],]}


name_tfs=paste("TF_", tfs, sep="")

for (i in 1:length(tfs)){
stia2021_rna <- AddModuleScore(stia2021_rna, features = list(grn_list[[i]]$Target), name = name_tfs[[i]] )}

to_plot <- colnames(stia2021_rna@meta.data)[grep("TF_", colnames(stia2021_rna@meta.data))]

Idents(stia2021_rna) <- 'cluster.name'
#levels(stia2021_rna)[-c(2,3,5,12)]
#stia2021_rna_f <- subset(stia2021_rna, idents=levels(stia2021_rna)[-c(2,3,5,12)])



dotplot<-DotPlot(stia2021_rna, features= to_plot, idents = levels(stia2021_rna)[-c(2,3,5,12)])

dotplot<-dotplot$data
library(tidyr)
dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

ann <- colnames(dotplot) %>% as.data.frame()
colnames(ann) <- 'clusters'
ann$condition <- c("lining", rep("sublining",8))


library(ComplexHeatmap)

colours <- list('condition'= ArchR::paletteDiscrete(ann$condition))
col_ann <- HeatmapAnnotation(df = ann[,c("condition") ], col=colours)



Heatmap(dotplot, row_names_gp = gpar(fontsize = 5), top_annotation = col_ann)
Heatmap(dotplot[,c(7,5,4,8,6,2,3,9,1)], row_names_gp = gpar(fontsize = 5), border=T, cluster_columns = F)




dotplto_1 <- dotplot %>% as.data.frame() %>%  filter(fibroblast_sublining_C1qtnf3_Col8a1 > 0.1  & fibroblast_sublining_Sfrp1_Cfb > 0.1 & fibroblast_sublining_lining_Ccl7_Ccl2 > 0.1)


dotplot["TF_RUNX11",]

dotplto_2 <- dotplot %>% as.data.frame() %>%  filter(fibroblast_sublining_Serpina3c_C3 > 0.1 & fibroblast_sublining_Ccl11 > 0.1 & fibroblast_sublining_Pi16 > 0.1 & fibroblast_sublining_Fmo2 > 0.1)

dotplto_3 <- dotplot %>% as.data.frame() %>%  filter(fibroblast__Crabp1_Col23a1 > 0.1)

dotplto_4 <- dotplot %>% as.data.frame() %>%  filter(fibroblast_lining_F13a1_Col22a1 > 0.1)

tf_1 <- rownames(dotplto_1) %>% as.data.frame() %>% cSplit(splitCols=".", sep="_") 
tf_1 <- substr(tf_1$._2, 1, nchar(tf_1$._2)-1)

tf_2 <- rownames(dotplto_2) %>% as.data.frame() %>% cSplit(splitCols=".", sep="_") 
tf_2 <- substr(tf_2$._2, 1, nchar(tf_2$._2)-1)

tf_3 <- rownames(dotplto_3) %>% as.data.frame() %>% cSplit(splitCols=".", sep="_") 
tf_3 <- substr(tf_3$._2, 1, nchar(tf_3$._2)-1)

tf_4 <- rownames(dotplto_4) %>% as.data.frame() %>% cSplit(splitCols=".", sep="_") 
tf_4 <- substr(tf_4$._2, 1, nchar(tf_4$._2)-1)


tf_1_genes <- grn %>% filter(Source == tf_1)
tf_2_genes <- grn %>% filter(Source == tf_2)
tf_3_genes <- grn %>% filter(Source == tf_3)
tf_4_genes <- grn %>% filter(Source == tf_4)

stia2021_rna <- AddModuleScore(stia2021_rna, features = list(tf_1_genes$Target), name="tf1_mod")
stia2021_rna <- AddModuleScore(stia2021_rna, features = list(tf_2_genes$Target), name="tf2_mod")
stia2021_rna <- AddModuleScore(stia2021_rna, features = list(tf_3_genes$Target), name="tf3_mod")
stia2021_rna <- AddModuleScore(stia2021_rna, features = list(tf_4_genes$Target), name="tf4_mod")


to_plot_tfs <- colnames(stia2021_rna@meta.data)[grep("tf", colnames(stia2021_rna@meta.data))]

Idents(stia2021_rna) <- 'condition'

dotplot<-DotPlot(stia2021_rna, features= to_plot_tfs)
dotplot

dotplot<-dotplot$data

dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

ann <- colnames(dotplot) %>% as.data.frame()
colnames(ann) <- 'clusters'

colours <- list('clusters'= ArchR::paletteDiscrete(ann$clusters))
col_ann <- HeatmapAnnotation(df = ann, col=colours)

Heatmap(dotplot[,c(1,3,6,2,4,5)], row_names_gp = gpar(fontsize = 5), cluster_columns = F,top_annotation = col_ann, border=T)


tf_1_genes$group <- 'mod1'
tf_2_genes$group <- 'mod2'
tf_3_genes$group <- 'mod3'
tf_4_genes$group <- 'mod4'

all_tfs <- rbind(tf_1_genes, tf_2_genes, tf_3_genes,tf_4_genes)

all_tfs$X <- NULL



g <- graph_from_data_frame(all_tfs[,c(5,4,3)], directed = TRUE)

degree_df <- data.frame(
  Node = V(g)$name,
  In_Degree = degree(g, mode = "in"),
  Out_Degree = degree(g, mode = "out")
)

index <- match(degree_df$Node, all_tfs$Source)
degree_df$cluster <- all_tfs$group[index]


degree_df_top <- degree_df %>%
  group_by(cluster) %>% 
  slice_head(n = 40) 

library(ggwordcloud)

gg_wc <- ggplot(degree_df_top,
                aes(label = Node, 
                    size = Out_Degree,
                    color = cluster)) +
  geom_text_wordcloud(shape = "square") +
  scale_size_area(max_size = 12) +
  #scale_color_manual(values = mycolors) +
  theme_minimal() +
  facet_wrap(~cluster)

gg_wc

