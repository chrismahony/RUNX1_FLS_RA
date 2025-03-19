library(ggtree)

Idents(stia2021_rna) <- "cluster.name"

stia2021_rna_f <- subset(stia2021_rna, idents=levels(stia2021_rna)[-c(2, 3, 5)])

stia2021_rna_f <- stia2021_rna_f %>% ScaleData()

rds <- BuildClusterTree(stia2021_rna_f, 
                        dims = 1:30)

data.tree <- Tool(object = rds, 
                  slot = "BuildClusterTree")



gg_tr <- ggtree(data.tree, 
                layout = "rectangular",
                alpha = .4,
                size = 3) +
  geom_tiplab(hjust = 1, size = 3) +
  geom_treescale() +
  theme(plot.margin = unit(c(1,5,1,1), "mm")) 

plot(gg_tr)
