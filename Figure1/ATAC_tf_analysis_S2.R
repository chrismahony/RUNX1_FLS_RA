DefaultAssay(obj.pair_all_coembed_merge_nomural) <- 'chromvar'
Idents(obj.pair_all_coembed_merge_nomural) <- 'cluster.name'
All_markers_TFs <- FindAllMarkers(obj.pair_all_coembed_merge_nomural, only.pos=T)


write.csv(All_markers_TFs, "/rds/projects/m/mahonyc-cesar-data/fibro_analysis/TF_subcluster/motifs.csv")

rm(list=ls()[! ls() %in% c("All_markers_TFs","obj.pair_all_coembed_merge_nomural")])

All_markers_TFs_f <- All_markers_TFs %>% filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
      slice_head(n = 3)

colnames(jaspar2) <- c('ID','gene')

index <- match(All_markers_TFs_f$gene, jaspar2$ID)
All_markers_TFs_f$gene_name <- jaspar2$gene[index]
rm(index)


dotplot<-DotPlot(obj.pair_all_coembed_merge_nomural, features = unique(All_markers_TFs_f$gene))

dotplot<-dotplot$data


index <- match(dotplot$features.plot, jaspar2$ID)
dotplot$features.plot <- jaspar2$gene[index]
rm(index)

library(tidyr)
library(dplyr)
dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

Heatmap(dotplot[,c(3,8, 4, 7, 6, 1, 2, 9, 10, 5)], border=T, cluster_columns = F, col=viridis(100))

ncol(dotplot)

heatmap(dotplot)
