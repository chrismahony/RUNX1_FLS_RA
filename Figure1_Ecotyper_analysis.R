


```{r}
rm(list=ls()[! ls() %in% c("stia2021_rna", "state_assignment")])

state_assignment <- read.delim("/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/ecotyper/ecotyper-master/DiscoveryOutput_STIArunx2/fib/state_assignment.txt")
state_assignment<-cSplit(state_assignment, splitCols ="ID", sep="." )
state_assignment$ID<-paste(state_assignment$ID_1, state_assignment$ID_2, sep="-")

stia_ecotyper<-stia2021_rna[,colnames(stia2021_rna) %in% state_assignment$ID]

library(tidyverse)
state_assignment<-state_assignment %>% remove_rownames %>% column_to_rownames(var="ID")
state_assignment<-subset(state_assignment, select=c(State))
stia_ecotyper<-AddMetaData(stia_ecotyper, state_assignment)
Idents(stia_ecotyper)<-'State'

state_markers<-FindAllMarkers(stia_ecotyper, only.pos = T)

library(gsfisher)

expressed_genes<-unique(state_markers$gene)
annotation_gs <- fetchAnnotation(species="mm", ensembl_version=NULL, ensembl_host=NULL)

index <- match(state_markers$gene, annotation_gs$gene_name)
state_markers$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- expressed_genes
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]

state_markers <- state_markers[!is.na(state_markers$ensembl),]
ensemblUni <- na.omit(ensemblUni)



go.results <- runGO.all(results=state_markers,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="p_val_adj", p_threshold=0.05,
                  species = "mm")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=3, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)


levels(stia_ecotyper)<-c('S01', 'S02', 'S03', 'S04', 'S05', 'S06', "S07")
DotPlot(stia_ecotyper, features = "Runx1") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") + scale_size(range = c(2, 8)) +RotatedAxis()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))

dotplot<-dotplot$data

ggplot(data=dotplot, aes(y=factor(id, level=c('S01', 'S02', 'S03', 'S04', 'S05', 'S06', "S07")), x=avg.exp)) +
  geom_bar(stat="identity", color="black", fill="grey")+
        theme(panel.background = element_blank())




gene_info <- read.delim("/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/ecotyper/ecotyper-master/DiscoveryOutput_STIArunx2/fib/gene_info.txt")


DimPlot(stia2021_rna)




pt <- table(stia_ecotyper$condition, stia_ecotyper$State)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)


cols <- ArchR::paletteDiscrete(stia_ecotyper@meta.data[, "State"])


p1<-ggplot(pt, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme(axis.text.y= element_blank(), axis.ticks.y = element_blank()) +
        coord_flip() + scale_fill_brewer(palette="Dark2")+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)

cols <- as.data.frame(ArchR::paletteDiscrete(stia_ecotyper@meta.data[, "cluster.name"]))
colnames(cols)<-"colors"
  
        
pt2 <- table(stia_ecotyper$cluster.name, stia_ecotyper$State)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

p2<-ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme(axis.text.y= element_blank(), axis.ticks.y = element_blank()) +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)


p1+p2
```
```{r}
heatmap_data <- read.delim("/rds/projects/c/croftap-visium-manuscript-01/Visium_CManalysis/ecotyper/ecotyper-master/DiscoveryOutput_STIArunx2/fib/heatmap_data.txt", row.names=1)

Heatmap(heatmap_data, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F)
```
```{r}
cols <- as.data.frame(ArchR::paletteDiscrete(stia2021_rna@meta.data[, "cluster.name"]))
colnames(cols)<-"colors"
cols$colors

cols=c( "grey", "grey", "grey", "#89288F", "grey", "#FEE500", "#8A9FD1", "#C06CAB", "#D8A767", "#90D5E4", "#89C75F", "#F47D2B", "#9983BD")
cols=c( "#272E6A", "grey", "grey", "#D51F26", "grey", "#D51F26", "#D51F26", "#D51F26", "#D51F26", "#D51F26", "#D51F26", "#D51F26", "#D51F26")
DimPlot(stia2021_rna, cols=cols)


```

