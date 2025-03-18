library(rcna)
library(ggthemes)
library(ggplot2)
library(dplyr)
library(glue)
library(harmony)
library(patchwork)
library(purrr)
library(Matrix)
library(data.table)

Idents(stia2021_rna) <- 'cluster.name'
levels(stia2021_rna)[-c(2,3,12)]
stia2021_rna_f <- subset(stia2021_rna, idents=levels(stia2021_rna)[-c(2,3,5)])
Idents(stia2021_rna_f) <- 'cluster.name'
levels(stia2021_rna_f) <- levels(stia2021_rna_f)[c(7,5,4,9,8,2,6,10,3,1)]
stia2021_rna_f$cluster.name <- stia2021_rna_f@active.ident

Idents(stia2021_rna) <- 'condition'
condition <- levels(stia2021_rna)[-1]
ggplots <- list()
Featureplots <- list()
Dotplots <- list()



for (i in 1:length(condition)){
  Idents(stia2021_rna_f) <- 'condition'
stia2021_rna_control_peak <- subset(stia2021_rna_f, idents=c("control", condition[[i]]))

stia2021_rna_control_peak@meta.data$status_val <- as.numeric(factor(stia2021_rna_control_peak@meta.data$condition, c('control', condition[[i]])))

stia2021_rna_control_peak <- FindNeighbors(stia2021_rna_control_peak, dims=1:30)

stia2021_rna_control_peak <- association.Seurat(
    seurat_object = stia2021_rna_control_peak, 
    test_var = 'status_val', 
    samplem_key = 'sample_id', 
    graph_use = 'RNA_nn', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    covs = NULL ## no covariates to include 
)


ord = cbind(stia2021_rna_control_peak@reductions[["umap"]]@cell.embeddings,stia2021_rna_control_peak@meta.data) %>%
    dplyr::select(c(UMAP_1,UMAP_2,cna_ncorrs,cluster.name)) %>%
    dplyr::group_by(cluster.name) %>%
    dplyr::summarize(med = median(cna_ncorrs)) %>%
    .$cluster.name


p <- cbind(stia2021_rna_control_peak@reductions[["umap"]]@cell.embeddings,stia2021_rna_control_peak@meta.data) %>%
    dplyr::select(c(UMAP_1,UMAP_2,cna_ncorrs,cluster.name)) %>%
      ggplot(aes(x=cluster.name, y=cna_ncorrs)) + 
    geom_violin(trim=TRUE, scale = "width")


cor_pos = stia2021_rna_control_peak@meta.data$cna_ncorrs_fdr05
cor_pos = cor_pos[cor_pos > 0]
# summary(cor_pos)
cor_neg = stia2021_rna_control_peak@meta.data$cna_ncorrs_fdr05
cor_neg = cor_neg[cor_neg < 0]


mywidth <- .35 # bit of trial and error
# This is all you need for the fill: 
vl_fill <- data.frame(ggplot_build(p)$data) %>%
  mutate(xnew = x- mywidth*violinwidth, xend = x+ mywidth*violinwidth) 

library(purrr)
library(tidyverse)

vl_poly <- 
  vl_fill %>% 
  dplyr::select(xnew, xend, y, group) %>%
  pivot_longer(-c(y, group), names_to = "oldx", values_to = "x") %>% 
  arrange(y) %>%
  split(., .$oldx) %>%
  map(., function(x) {
    if(all(x$oldx == "xnew")) x <- arrange(x, desc(y))
    x
    }) %>%
  bind_rows()



ggplots[[i]] <- ggplot() +
    geom_polygon(data = vl_poly, aes(x, y, group = group), size = 1, fill = NA) +  
    geom_segment(data = vl_fill, aes(x = xnew, xend = xend, y = y, yend = y, color = y)) +
    scale_color_gradient2(midpoint = 0, low = '#3C5488FF', mid = "white", high = '#DC0000FF', space = "Lab") +
    scale_x_continuous(labels=ord, breaks=1:length(unique(ord)),limits=c(min(vl_poly$x),max(vl_poly$x))) +
    geom_hline(yintercept=c(min(cor_pos), max(cor_neg)), linetype='dashed', color='black') +
    theme_classic() +
    theme(strip.text.x=element_text(size=15, color="black", face="bold"),
          strip.text.y=element_text(size=15, color="black", face="bold"),
          legend.position = "none",
          plot.title = element_text(size=15),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size =15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size =15),
          legend.text =  element_text(size = 15),
          legend.key.size = grid::unit(0.5, "lines"),
          legend.title = element_text(size = 0.8, hjust = 0)) +
    labs(title = paste("CNA control vs", condition[[i]], sep=""),
         x = "",
         y = "Neighborhood\nCorrelation")+RotatedAxis()

print(ggplots[[i]])
library(viridis)
Featureplots[[i]] <- FeaturePlot(stia2021_rna_control_peak, features = "cna_ncorrs", cols =c("blue", "orange", "yellow")) + ggtitle(paste0("control vs ", condition[[i]]))+NoAxes()+NoLegend()
Idents(stia2021_rna_control_peak) <- 'cluster.name'
Dotplots[[i]] <- DotPlot(stia2021_rna_control_peak, features = "cna_ncorrs") + ggtitle(paste0("control vs ", condition[[i]]))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+RotatedAxis()
Idents(stia2021_rna_control_peak) <- 'condition'

print(Dotplots[[i]])


print(Featureplots[[i]])

}

#plot_grid(ggplots[[1]], ggplots[[2]], ggplots[[3]], ggplots[[4]], ggplots[[5]])

FeaturePlot(stia2021_rna_control_peak, features = "cna_ncorrs", cols =c("blue", "orange", "yellow"))

for (i in 1:length(ggplots)){
  print(ggplots[[i]])
#print(Featureplots[[i]])

}


plot_grid(Featureplots[[2]], Featureplots[[5]], Featureplots[[1]], Featureplots[[3]], Featureplots[[4]], ncol = 5)


plot_grid(ggplots[[2]]+coord_flip()+ 
  theme(axis.text.y = element_blank()), ggplots[[5]]+coord_flip()+ 
  theme(axis.text.y = element_blank()), ggplots[[1]]+coord_flip()+ 
  theme(axis.text.y = element_blank()), ggplots[[3]]+coord_flip()+ 
  theme(axis.text.y = element_blank()),ggplots[[4]]+coord_flip()+ 
  theme(axis.text.y = element_blank()), ncol=5)


