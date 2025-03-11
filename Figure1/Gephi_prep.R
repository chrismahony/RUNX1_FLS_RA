write.csv(grn, "/rds/projects/c/croftap-runx1data01/Dkk3_analysis/bulk_RNA/grn.csv")

#grn<-grn[grn$tf=="RUNX1",]
grn<-grn[grn$weights > 0.6,]

head(grn)
dim(grn)
colnames(grn)[1:3] <- c("Source", "Target", "Weight")
head(grn)

grn[["Type"]] <- "Directed"
grn %>%
  mutate('Id' = paste0(Source, '_', Target)) -> grn
head(grn)

# - community detection
# - igraph definition
g <- graph_from_data_frame(grn[, c("Source", "Target")],
                           directed = FALSE)
g <- set_edge_attr(g, "weight", value = grn$Weight)
g <- set_edge_attr(g, "name", value = grn$Id)

library(igraph)

# leiden
leiden_mod <- cluster_louvain(g)
mods <- data.frame(cbind(V(g)$name, leiden_mod$membership))

colnames(mods) <- c("Id", "leiden_mod")
head(mods)
table(mods$leiden_mod)

#write.table(grn, "~/stroma_stia/atac_rna.dir/scmega.dir/grn_dis_filtered_gephi_edges.tsv",
            #sep = "\t", quote = FALSE, row.names = FALSE)


nodes <- data.frame('Id' = unique(grn$Source))
nodes[["class"]] <- "TF"

target_nodes <- data.frame('Id' = unique(grn$Target))
target_nodes[["class"]] <- "gene"
target_nodes %>%
  filter(!Id %in% nodes$Id) -> target_nodes


nodes <- rbind(nodes, target_nodes)
nodes[["Label"]] <- nodes$Id

nodes <- merge(nodes, mods, by = "Id", all.x = TRUE, sort = FALSE)
head(nodes)
dim(nodes)



e <- get.edgelist(g,names=FALSE)

l <- qgraph::qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                                 area=10*(vcount(g)^2),
                                                 repulse.rad=(vcount(g)^3.1),
                                                 niter = 1000,
                                                 max.delta = 4, 
                                                 cool.exp = 0.3)


plot(leiden_mod, g,
       layout=l,
       vertex.label.cex=0.01, 
       vertex.label.family="Helvetica",
       vertex.label.font=0.5,
       vertex.shape="square", 
       vertex.size=0, 
       col = rgb(1,1,1,0), 
       edge.color = rgb(0,0,0,0.02),
       edge.width = 0.5,
     vertex.label=get.vertex.attribute(g)$name)

nodes$mouse_name <- nodes$Id
nodes$mouse_name <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(nodes$mouse_name), perl=TRUE)
nodes$class <- NULL

write.table(dplyr::select(nodes, c(Id, Label, leiden_mod, mouse_name)), 
            "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/nodes_grn.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(grn, "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/gephi_edges_GRN.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)



###############################
grn_RUNX1 <- grn %>% filter(Source == "RUNX1")

write.table(grn_RUNX1, "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/gephi_edges_grn_RUNX1.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


nodes_RUNX1 <- nodes_grn[nodes_grn$Label %in% grn_RUNX1$Human_target,]

write.table(nodes_RUNX1, "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/nodes_RUNX1.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

