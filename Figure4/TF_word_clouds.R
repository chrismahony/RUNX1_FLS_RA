row_cl <- read.csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/DE_clusters_wordclouds/row_cl.csv", row.names = 1)

top_edges <- read.csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/DE_clusters_wordclouds/top_edges.csv", row.names = 1)

TF_list <- read.delim("/rds/projects/c/croftap-stia-atac/TF_list/homer_TFs.txt", col.names = F)

index <- match(top_edges$Target, row_cl$gene)
top_edges$cluster <- row_cl$gene_cluster[index]


#nodes <- dplyr::select(nodes, c(label, class, louvain_mod, in_degree, out_degree))
#head(nodes)
#colnames(nodes)[c(1, 3)] <- c("gene", "gene_module")
#nodes$gene_module <- paste0("M", nodes$gene_module)

mods <- names(which(table(top_edges$cluster) > 10))

nodes <- filter(top_edges, cluster %in% mods) %>%
  arrange(cluster, desc(Weight))


nodes_TF <- nodes[nodes$Source %in% TF_list$FALSE.,]

g <- graph_from_data_frame(nodes_TF, directed = TRUE)

# Calculate in-degree and out-degree
degree_df <- data.frame(
  Node = V(g)$name,
  In_Degree = degree(g, mode = "in"),
  Out_Degree = degree(g, mode = "out")
)


index <- match(degree_df$Node, row_cl$gene)
degree_df$cluster <- row_cl$gene_cluster[index]

#degree_df_f <- degree_df %>% filter(Out_Degree > 80)

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


table(degree_df$cluster)

ggplot(degree_df,
                aes(label = Node, 
                    size = Out_Degree,
                    color = cluster)) +
  geom_text_wordcloud(shape = "square") +
  scale_size_area(max_size = 4) +
  #scale_color_manual(values = mycolors) +
  theme_minimal() +
  facet_wrap(~cluster)

gg_wc
