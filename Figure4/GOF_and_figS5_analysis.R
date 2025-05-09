
```{r}
#remove RUNX1
subres_human_RA=all_DEGs_human[-2,]


DEGs_human_RA_down <- subres_human_RA %>% filter(comparison == "EV_vs_R1C" & log2FoldChange > 0.5 &padj < 0.05)

mouse_RA_down <- orthologs(genes = DEGs_human_RA_down$gene, species = "mouse")


overlap_down <- mouse_RA_down[mouse_RA_down$symbol %in% mouse_inhibitor_down$symbol,]




```




```{r}



#subres_human_kidney <- subres
subres_f <- all_GEX_human_kidney %>% filter(padj < 0.05)

library(babelgene)

mouse <- orthologs(genes = subres_f$gene, species = "mouse")
DEGs_human_RA <- subres_human_RA %>% filter(comparison == "EV_vs_R1C" & log2FoldChange < -2 &padj < 0.05)


mouse_RA <- orthologs(genes = DEGs_human_RA$gene, species = "mouse")

stia2021_rna <- AddModuleScore(stia2021_rna, features = list(mouse$symbol), name = "human_kidney_fibs")

stia2021_rna <- AddModuleScore(stia2021_rna, features = list(mouse_RA$symbol), name = "human_RA_fibs")


FeaturePlot(stia2021_rna, features="human_kidney_fibs1", min.cutoff = "q10", max.cutoff = "q90")




stia2021_rna$cluster_condition <- paste(stia2021_rna$cluster.name, stia2021_rna$condition, sep="_")


Idents(stia2021_rna) <- 'cluster_condition'
dotplot<-DotPlot(stia2021_rna,features=c("human_RA_fibs1", "mouse_inhibitor1"))

dotplot<-dotplot$data

dotplot<-dotplot %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)


row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)


Heatmap(dotplot)

```


```{r}

#subres_human_kidney <- subres
subres_f <- all_GEX_human_kidney %>% filter(log2FoldChange > 2 & padj < 0.05)


library(babelgene)

mouse <- orthologs(genes = subres_f$gene, species = "mouse")
DEGs_human_RA <- subres_human_RA %>% filter(comparison == "EV_vs_R1C" & log2FoldChange < -2 &padj < 0.05)

DEGs_human_RA[DEGs_human_RA$gene_name %in% inhibitor_DEGs_f$gene_name,]


mouse_RA <- orthologs(genes = DEGs_human_RA$gene, species = "mouse")

overlap <- mouse_RA[mouse_RA$symbol %in% mouse_inhibitor$symbol,]


NatMed_AliverniniEtAl_2020_FLS <- AddModuleScore(NatMed_AliverniniEtAl_2020_FLS, features = list(overlap$human_symbol), name = "overlap")

#NatMed_AliverniniEtAl_2020_FLS <- AddModuleScore(NatMed_AliverniniEtAl_2020_FLS, features = list(subres_f$gene), name = "human_kidney_fibs")

NatMed_AliverniniEtAl_2020_FLS <- AddModuleScore(NatMed_AliverniniEtAl_2020_FLS, features = list(mouse_RA$human_symbol), name = "human_RA_fibs")

NatMed_AliverniniEtAl_2020_FLS <- AddModuleScore(NatMed_AliverniniEtAl_2020_FLS, features = list(DEGs_human_RA_down$gene_name), name = "mouse_RA_down")
NatMed_AliverniniEtAl_2020_FLS <- AddModuleScore(NatMed_AliverniniEtAl_2020_FLS, features = list(inhibitor_DEGs_f_down$gene_name), name = "inhibitor_DEGs_f_down")



NatMed_AliverniniEtAl_2020_FLS$group_cluster<-paste(NatMed_AliverniniEtAl_2020_FLS$group, NatMed_AliverniniEtAl_2020_FLS$clusters,  sep=".")
Idents(NatMed_AliverniniEtAl_2020_FLS) <- 'group_cluster'
dotplot <- DotPlot(NatMed_AliverniniEtAl_2020_FLS,features=c( "human_RA_fibs1", "RUNX1"))+RotatedAxis()


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



```


```{r}
grn <- read_csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/grn.csv")

grn_RUNX1 <- grn %>% filter(Source=="RUNX1")


med_fibros <- AddModuleScore(med_fibros, features = list(grn_RUNX1$Target), name = "RUNX1_regulons")

#extract avg. scaled expression for your gene in each cluster and sample
Idents(med_fibros)<-"sample_tissue_cluster"
dotplot<-DotPlot(med_fibros, features = "RUNX1_regulons1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"


dotplot_data <- dotplot_data %>% cSplit(splitCols = "id", sep=".")

inflam_df <- paste(med_fibros$SampleID, med_fibros$InflamScore, sep=".") %>% unique() %>% as.data.frame() %>% cSplit(splitCols = ".", sep=".")

colnames(inflam_df) <- c("id_1", "inflam", "extra")

inflam_df <- inflam_df %>% replace(is.na(.), 0)

inflam_df$inflam <- paste(inflam_df$inflam, inflam_df$extra, sep=".")
inflam_df$inflam <- as.double(inflam_df$inflam)

final_df <- dotplot_data %>% 
  left_join(inflam_df, by="id_1")


#plot
#dotplot_data_synovium <- dotplot_data[dotplot_data$condition_2 == "Synovium",]
ggplot(final_df, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") +
        facet_wrap(~id_2)


#ml = lm(Runx1~InflamScore, data = dotplot_data)
#summary(ml)$r.squared

library(ggpubr)
ggscatter(final_df, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)+
        facet_wrap(~id_2)



#filter if you want a specific cluster
dotplot_data_C4 <- final_df[final_df$id_3 == "SPARC+COL3A1+ C4",]
dotplot_data_C4 <- dotplot_data_C4[dotplot_data_C4$id_2 == "Synovium",]

ggplot(dotplot_data_C4, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") +
        facet_wrap(~id_2)




ml = lm(Runx1~InflamScore, data = dotplot_data_C4)
summary(ml)$r.squared

library(ggpubr)

dotplot_data_syn <- final_df[final_df$id_2 == "Synovium",]

ggscatter(dotplot_data_C4, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 1)+
        facet_wrap(~id_2)


ggplot(data=dotplot_data_C4, aes(Runx1,inflam)) + geom_point(alpha=0.6, color="grey", size=0.1) + ggtitle("C4 cluster") +theme_minimal() +
    geom_point(data = dotplot_data_C4, color = "darkred",size=1.5)+ theme(axis.title.x = element_blank())+
    stat_smooth(method = "lm",
        col = "black",
        se = T,
        size = 0.5)+theme_ArchR()


```











```{r}
inhibitor_DEGs$cluster <- "NO"
inhibitor_DEGs$cluster[inhibitor_DEGs$log2FoldChange > 0.001] <- "DOWN_inhib"
inhibitor_DEGs$cluster[inhibitor_DEGs$log2FoldChange < -0.001] <- "UP_inhib"

subres_human_RA$cluster <- "NO"
subres_human_RA$cluster[subres_human_RA$log2FoldChange > 0.001] <- "DOWN_lenti"
subres_human_RA$cluster[subres_human_RA$log2FoldChange < -0.001] <- "UP_lenti"

all_GEX_human_kidney$score <- all_GEX_human_kidney$log2FoldChange
all_GEX_human_kidney$gene_name <- all_GEX_human_kidney$gene

index <- match(all_GEX_human_kidney$gene_name, annotation_gs$gene_name)
all_GEX_human_kidney$gene <- annotation_gs$ensembl_id[index]
all_GEX_human_kidney <- na.omit(all_GEX_human_kidney)


all_GEX_human_kidney2 <- all_GEX_human_kidney[,colnames(subres_human_RA)]

all_GEX_human_kidney2$cluster <- "NO"
all_GEX_human_kidney2$cluster[all_GEX_human_kidney2$log2FoldChange > 0.001] <- "UP_kid"
all_GEX_human_kidney2$cluster[all_GEX_human_kidney2$log2FoldChange < -0.001] <- "DOWN_kid"

write.csv(all_GEX_human_kidney, "/rds/projects/m/mahonyc-runx1-bulk-seq-data/lentivrial_GOF/analysis_bact_corrected/all_GEX_human_kidney.csv")



library(gsfisher)
#annotation_gs<-fetchAnnotation(species = "hs")

all_DEGs <- rbind(inhibitor_DEGs, subres_human_RA)

index <- match(all_DEGs$gene_name, annotation_gs$gene_name)
all_DEGs$ensembl <- annotation_gs$ensembl_id[index]

ensemblUni <- unique(all_DEGs$gene)
ensemblUni <- na.omit(ensemblUni)
all_DEGs <- na.omit(all_DEGs)
all_DEGs_f <- all_DEGs %>% filter(padj < 0.05)


go.results <- runGO.all(results=all_DEGs_f,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="padj", p_threshold=0.05,
                  species = "hs")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=20, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(c(go.results.top$description[1:3], "DNA replication", "nuclear chromosome segregation", "mitotic nuclear division", "positive regulation of T cell mediated immunity", "T cell mediated immunity", "regulation of adaptive immune response", "regulation of actin filament organization", "intermediate filament organization", "basal plasma membrane")), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T, p_col = "p.val", sample_levels = c("UP_inhib",   "DOWN_inhib", "UP_lenti" ,  "DOWN_lenti"))

write.csv(go.results.top, "/rds/projects/m/mahonyc-runx1-bulk-seq-data/lentivrial_GOF/analysis_bact_corrected/go.results.top_inhib_lenti.csv")
write.csv(all_DEGs_f, "/rds/projects/m/mahonyc-runx1-bulk-seq-data/lentivrial_GOF/analysis_bact_corrected/all_DEGs_f_inhib_lenti.csv")

```


```{r}

ensemblUni_kid <- unique(all_GEX_human_kidney2$gene)
ensemblUni_kid <- na.omit(ensemblUni_kid)
all_GEX_human_kidney2_f <- all_GEX_human_kidney2 %>% filter(padj < 0.05)


go.results_kid <- runGO.all(results=all_GEX_human_kidney2_f,
                  background_ids = ensemblUni, gene_id_col="gene", gene_id_type="ensembl", sample_col="cluster", p_col="padj", p_threshold=0.05,
                  species = "hs")
go.results_kid <- filterGenesets(go.results_kid)
go.results_kid.top <- go.results_kid %>% group_by(cluster) %>% top_n(n=3, -p.val)
sampleEnrichmentDotplot(go.results_kid.top, selection_col = "description", selected_genesets = unique(go.results_kid.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T, p_col = "p.val", sample_levels = c("UP_kid",   "DOWN_kid"))

```



```{r}

all_DEGs <- rbind(inhibitor_DEGs, subres_human_RA)
all_DEGs <- na.omit(all_DEGs)
all_DEGs_f <- all_DEGs %>% filter(padj < 0.05)


data <-table(all_DEGs_f$cluster) %>% as.data.frame()

data$condition <- c("inhib", "lenti", "inhib", "lenti")
data$direction <- c("down", "down", "up", "up")


ggplot(data, aes(fill=direction, y=Freq, x=condition, label = Freq)) + 
    geom_bar(position="stack", stat="identity")+theme_ArchR() +
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
    scale_fill_viridis(discrete = T)+coord_flip()+scale_fill_manual(values = c("lightblue","red"))


data %>% filter(condition== "inhib") %>% ggplot(aes(y=Freq, x=direction, fill=direction)) + 
    geom_bar(position="stack", stat="identity")+theme_ArchR() +
      scale_fill_viridis(discrete = T)+coord_flip()+scale_fill_manual(values = c("lightblue","red"))


data %>% ggplot(aes(y=Freq, x=direction, fill=direction))  + 
        geom_bar(stat = 'identity', position = position_dodge(), aes(group = direction),width = 1) + 
        coord_flip() + 
        facet_wrap(~condition, scales = 'free_y', nrow=2, ncol=1)+ 
        theme(
            axis.text.x = element_text(angle = 45, hjust=1),
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
            # strip.text = element_blank()
        ) + 
        guides(color = 'none', fill = 'none') + 
        labs(y = '# genes')+
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = c(0, 0)))+theme_ArchR()+scale_fill_manual(values = c("blue","red"))+
  theme(strip.background = element_rect(fill="white", size=1, color="white"))

 


```















```{r}

res1_mouse_kid <- read_csv("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE139950_humankidney_fibs/res1.csv")

write.csv(res1_mouse_kid%>% filter(padj < 0.05),  "/rds/projects/m/mahonyc-runx1-bulk-seq-data/lentivrial_GOF/analysis_bact_corrected/res1_mouse_kid.csv")


subres_f <- res1_mouse_kid %>% filter(log2FoldChange > 2 & padj < 0.05)
subres_f <- na.omit(subres_f)

library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
enriched <- enrichr(c(unique(subres_f$gene)), dbs)
plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")



enriched[[3]]$log_adj_p <- -log(enriched[[3]]$Adjusted.P.value)


enriched[[3]]$sig <- "NO"
enriched[[3]]$sig[enriched[[3]]$Adjusted.P.value > 0.05] <- "no"
enriched[[3]]$sig[enriched[[3]]$Adjusted.P.value < 0.05] <- "yes"



enriched[[3]] %>%  ggplot(aes(x= Odds.Ratio, y= log_adj_p, label=Term)) + 
                  geom_point() +theme_ArchR()+
    geom_point(data = enriched[[3]] %>% filter(Term == "inflammatory response (GO:0006954)" ), color = "red")
colnames(enriched[[3]])

  

```

```{r}




library(EnhancedVolcano)

genes=c("CSF1R","SEMA7A","CXCL8","TNFAIP6","GPR68","IL20RB","LY96","CXCL3","PTGS2","IL1A","HYAL1","IL1B","CCRL2","SPP1","AOX1","PTX3","SCG2","CCL26")

EnhancedVolcano(res1_mouse_kid,
    lab = res1_mouse_kid$gene,
    x = 'log2FoldChange',
    y = 'padj',
        selectLab = genes,
    title = 'title',
    subtitle = "GEX, red=p_adj<0.05 & FC > 2",
    pCutoff = 0.05,
    FCcutoff = 2,
    pointSize = 2.0,
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
    )+xlim(-10, 10)

```





```{r}
Mm_Matrisome_Masterlist_Naba_et_al_2012 <- read_excel("Mm_Matrisome_Masterlist_Naba et al_2012.xlsx")

subres_f
inhibitor_DEGs_f
DEGs_human_RA

Mm_Matrisome_Masterlist_Naba_et_al_2012_f <- Mm_Matrisome_Masterlist_Naba_et_al_2012[Mm_Matrisome_Masterlist_Naba_et_al_2012$`Gene Symbol` %in% DEGs_human_RA$gene,]


table(Mm_Matrisome_Masterlist_Naba_et_al_2012_f$`Matrisome Category`) %>% as.data.frame() %>% 
ggplot(aes(x = Freq, y = reorder(Var1, Freq), fill = Var1)) +
  geom_col()  +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "off")+theme_ArchR()+NoLegend()



```

```{r}
subres_f <- all_GEX_human_kidney %>% filter(log2FoldChange > 2 & padj < 0.05)

kids_f_med_fibros_overlap <- all_markers_medfibros[all_markers_medfibros$gene %in% subres_f$gene_name,]
kids_f_med_fibros_overlap <- kids_f_med_fibros_overlap %>% filter(p_val_adj < 0.05)
all_markers_medfibros_f <- all_markers_medfibros %>% filter(p_val_adj<0.05)
df <- table(kids_f_med_fibros_overlap$cluster) %>% as.data.frame()
table(all_markers_medfibros_f$cluster) %>% as.data.frame()
df$total <- c(488,281,146,319,283,240,333,409,237,298,219,290,258)


df$pct <- df$Freq/df$total*100

ggplot(df, aes(x=Freq, y=pct)) + 
  geom_point()+theme_ArchR()+xlim(0,17)+ geom_text_repel(aes(label = df$Var1),
                    size = 3.5,  ylim = c(0, 5.4)) +
  geom_hline(yintercept=4, linetype='dotted', col = 'red', size=0.5)+geom_vline(xintercept = 10, linetype="dotted", 
                color = "red", size=0.5)+
    geom_point(data = df %>% filter(Freq > 10 & pct > 4), color = "red")
```



```{r}

stim_obj$sample_tissue_condition<-paste(stim_obj$DonorID, stim_obj$Tissue, stim_obj$Condition, sep=".")
Idents(stim_obj)<-'sample_tissue_condition'
levels(stim_obj)

cts_fibs<-AggregateExpression(stim_obj, group.by = c("sample_tissue_condition"), assays = "RNA", slot = "counts", return.seurat = F)

cts_fibs<-cts_fibs$RNA
cts_fibs<-as.data.frame(cts_fibs)
meta_data=colnames(cts_fibs)
meta_data<-as.data.frame(meta_data)
library(splitstackshape)
meta_data$to_split<-meta_data$meta_data
meta_data<-cSplit(meta_data, splitCols = "to_split", sep=".")
colnames(meta_data)<-c("all","sample", "tissue", "condition")
meta_data$all<-as.factor(meta_data$all)
meta_data$sample<-as.factor(meta_data$sample)
meta_data$tissue<-as.factor(meta_data$tissue)
meta_data$condition<-as.factor(meta_data$condition)
meta_data$comparison <- paste(meta_data$tissue, meta_data$condition, sep="_")
meta_data$comparison <- as.factor(meta_data$comparison)


dds <- DESeqDataSetFromMatrix(countData = cts_fibs,
                                  colData = meta_data,
                                  design = ~1)


dds <- scran::computeSumFactors(dds)
print(dds)
print(quantile(rowSums(counts(dds))))

#mingenecount <- quantile(rowSums(counts(dds)), 0.5)
mingenecount <- 200
maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount & rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

levels= c(     "Lung_ECs"    ,     "Lung_Tcells" ,     
 "Synovium_ECs" ,    "Synovium_Tcells", "Lung_Control", "Synovium_Control" )

dds@colData[['comparison']] <- factor(dds@colData[['comparison']],
                                     levels = levels)

design(dds) <- formula(~ comparison)
print(design(dds))
dds <- DESeq(dds, test = "Wald")


targetvar <- "comparison"

comps1 <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))
      head(comps1)
      
      ress <- apply(comps1, 1, function(cp) {
        print(cp)
        res <- data.frame(results(dds, contrast=c(targetvar, cp[1], cp[2])))
        res[["gene"]] <- rownames(res)
        res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
        res
      })
      
      res1 <- Reduce(rbind, ress)

res1 %>% 
      filter(padj < 0.01) %>%
      mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
      arrange(desc(abs(score))) -> subres

res1_sig <- res1 %>% 
      filter(padj < 0.05)

res2_sig<-res1_sig %>% 
  select(log2FoldChange, gene, comparison) %>%  
  pivot_wider(names_from = comparison, values_from = log2FoldChange) %>% 
  as.data.frame() 

res3_sig <- res2_sig[,c(1,2,3,14,16)]

res3_sig$Lung_Control_vs_Lung_ECs <- res3_sig$Lung_Control_vs_Lung_ECs*-1
res3_sig$Lung_Control_vs_Lung_Tcells <- res3_sig$Lung_Control_vs_Lung_Tcells *-1
res3_sig$Synovium_Control_vs_Synovium_ECs <- res3_sig$Synovium_Control_vs_Synovium_ECs *-1
res3_sig$Synovium_ECs_vs_Synovium_Tcells <- res3_sig$Synovium_ECs_vs_Synovium_Tcells *-1

res2<-res1 %>% 
  select(log2FoldChange, gene, comparison) %>%  
  pivot_wider(names_from = comparison, values_from = log2FoldChange) %>% 
  as.data.frame() 

colnames(res2)

res3 <- res2[,c(1,2,3,7,14,16)]

res3$Lung_ECs_vs_Lung_Tcells <- res3$Lung_ECs_vs_Lung_Tcells*-1
res3$Lung_Control_vs_Lung_ECs <- res3$Lung_Control_vs_Lung_ECs*-1
res3$Lung_Control_vs_Lung_Tcells <- res3$Lung_Control_vs_Lung_Tcells *-1
res3$Synovium_Control_vs_Synovium_ECs <- res3$Synovium_Control_vs_Synovium_ECs *-1
res3$Synovium_ECs_vs_Synovium_Tcells <- res3$Synovium_ECs_vs_Synovium_Tcells *-1


DEGs_human_RA <- subres_human_RA %>% filter(comparison == "EV_vs_R1C" & log2FoldChange < -2 &padj < 0.05)


res_f <- res3_sig[res3_sig$gene %in% DEGs_human_RA$gene_name,]
rownames(res_f) <- res_f$gene

res_f <- res3[res3$gene %in% DEGs_human_RA$gene_name,]
rownames(res_f) <- res_f$gene

genes_to_highlihgt <- res_f %>% filter(Synovium_ECs_vs_Synovium_Tcells >2 & Lung_ECs_vs_Lung_Tcells> 2) %>% rownames()

rownames(res3) <- res3$gene
p3 <- ggplot(data=res3, aes(Synovium_ECs_vs_Synovium_Tcells,Lung_ECs_vs_Lung_Tcells)) + geom_point(alpha=0.6, color="grey", size=0.1) + ggtitle("Lung") +theme_minimal()+
  geom_hline(yintercept=0, linetype='dotted', col = 'red', size=0.5)+geom_vline(xintercept = 0, linetype="dotted", 
                color = "red", size=0.3) +
    geom_point(data = res_f, color = "darkred",size=0.5,  max.overlaps = 20)+labs(y= "Lung (logFC)") 

p3 <- LabelPoints(plot = p3, points = c("CXCL11", "CCL1", "IL4I1", "TNFSF10", "IRF7"), repel = TRUE, xnudge = 0.5,
  ynudge = 1.2)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme_ArchR()+ theme(axis.title.x = element_blank())

p3

```

```

## Custom upset plot

res2<-res1 %>% 
  filter(padj < 0.05) %>% 
  select(log2FoldChange, gene, comparison) %>%  
  pivot_wider(names_from = comparison, values_from = log2FoldChange) %>% 
  as.data.frame() 

colnames(res2)

res3 <- res2[,c(1,2,3,7,14,16)]

res3$Lung_ECs_vs_Lung_Tcells <- res3$Lung_ECs_vs_Lung_Tcells*-1
res3$Lung_Control_vs_Lung_ECs <- res3$Lung_Control_vs_Lung_ECs*-1
res3$Lung_Control_vs_Lung_Tcells <- res3$Lung_Control_vs_Lung_Tcells *-1
res3$Synovium_Control_vs_Synovium_ECs <- res3$Synovium_Control_vs_Synovium_ECs *-1
res3$Synovium_ECs_vs_Synovium_Tcells <- res3$Synovium_ECs_vs_Synovium_Tcells *-1

DEGs_human_RA_new <- subres_human_RA %>%
  filter(
    comparison == "EV_vs_R1C",
    padj < 0.05,
    abs(log2FoldChange) > log2(2)
  )

res3 <- res3[res3$gene %in% DEGs_human_RA_new$gene_name,]


my_list <- list(
  Sy_ECs_vs_Tcells_UP = res3 %>% filter(Synovium_ECs_vs_Synovium_Tcells > 0) %>% pull(gene), 
  Sy_ECs_vs_Tcells_DOWN = res3 %>% filter(Synovium_ECs_vs_Synovium_Tcells < 0) %>% pull(gene),
  Lu_ECs_vs_Tcells_UP = res3 %>% filter(Lung_ECs_vs_Lung_Tcells > 0) %>% pull(gene),
  Lu_ECs_vs_Tcells_DOWN = res3 %>% filter(Lung_ECs_vs_Lung_Tcells < 0) %>% pull(gene)
  )


comb_mat <- make_comb_mat(my_list)
my_names <- set_name(comb_mat)

my_set_sizes <- set_size(comb_mat) %>% 
  as.data.frame() %>% 
  dplyr::rename(sizes = ".") %>% 
  dplyr::mutate(Set = row.names(.)) 

library(RColorBrewer)

p1 <- my_set_sizes %>% 
  mutate(Set = reorder(Set, sizes)) %>% 
  ggplot(aes(x = Set, y= sizes)) +
  geom_bar(stat = "identity", aes(fill = Set), alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sizes), 
            size = 5, angle = 90, hjust = 0, y = 1) +
  scale_fill_manual(values = brewer.pal(4, "Set1"),  # feel free to use some other colors  
                     limits = my_names) + 
  labs(x = NULL,
       y = "Set size",
       fill = NULL) +
  theme_classic() +
  theme(legend.position = "right",
        text = element_text(size= 14),
        axis.ticks.y = element_blank(),
        axis.text = element_blank()
        ) 


p1 <- my_set_sizes %>% 
  mutate(Set = reorder(Set, sizes)) %>% 
  ggplot(aes(x = Set, y= sizes)) +
  geom_bar(stat = "identity", aes(fill = Set),color="black", alpha = 0.8, width = 0.7)  +
  scale_fill_manual(values = brewer.pal(4, "Set1"),  # feel free to use some other colors  
                     limits = my_names) + 
  labs(x = NULL,
       y = "Set size",
       fill = NULL) +
  theme_classic() +
  theme(legend.position = "right",
        text = element_text(size= 14),
        axis.ticks.y = element_blank(),
        axis.text = element_blank()
        ) 


p1


get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

p2 <- get_legend(p1)

my_overlap_sizes <- comb_size(comb_mat) %>% 
  as.data.frame() %>% 
  dplyr::rename(overlap_sizes = ".") %>% 
  dplyr::mutate(category = row.names(.))

p3 <- my_overlap_sizes %>% 
  mutate(category = reorder(category, -overlap_sizes)) %>% 
  ggplot(aes(x = category, y = overlap_sizes)) +
  geom_bar(stat = "identity", fill = "grey80", color = NA, alpha = 0.8, width = 0.7) +
  geom_text(aes(label = overlap_sizes, y = 0), 
            size = 5, hjust = 0, vjust = 0.5) +
  labs(y = "Intersect sizes",
       x = NULL) +
  theme_classic() +
  theme(text = element_text(size= 14, color = "black"),
        axis.text =element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(hjust = 0),
        ) +
  coord_flip()



p3 <- my_overlap_sizes %>% 
  mutate(category = reorder(category, -overlap_sizes)) %>% 
  ggplot(aes(x = category, y = overlap_sizes)) +
  geom_bar(stat = "identity", fill = "grey", color = "black", alpha = 0.8, width = 0.7) +
  labs(y = "Intersect sizes",
       x = NULL) +
  theme_classic() +
  theme(text = element_text(size= 14, color = "black"),
        axis.text =element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(hjust = 0),
        ) +
  coord_flip()

p3


my_overlap_matrix <- str_split(string = my_overlap_sizes$category, pattern = "", simplify = T) %>% 
  as.data.frame() 

colnames(my_overlap_matrix) <- my_names

my_overlap_matrix_tidy <- my_overlap_matrix %>% 
  cbind(category = my_overlap_sizes$category) %>% 
  pivot_longer(cols = !category, names_to = "Set", values_to = "value") %>% 
  full_join(my_overlap_sizes, by = "category") %>% 
  full_join(my_set_sizes, by = "Set")

p4 <- my_overlap_matrix_tidy %>% 
  mutate(category = reorder(category, -overlap_sizes)) %>%  
  mutate(Set = reorder(Set, sizes)) %>%  
  ggplot(aes(x = Set, y = category))+
  geom_tile(aes(fill = Set, alpha = value), color = "black", size = 0.5) +
  scale_fill_manual(values = brewer.pal(4, "Set1"), # feel free to use other colors 
                    limits = my_names) +
  scale_alpha_manual(values = c(0.8, 0),  # color the grid for 1, don't color for 0. 
                     limits = c("1", "0")) +
  labs(x = "Sets",  
       y = "Overlap") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size= 14),
        panel.grid = element_blank(),
        axis.text = element_blank()
        )

p4

library(RVenn)
my_object <- RVenn::Venn(my_list)

ggvenn(
  my_object, slice = 1:3, 
  thickness = 0.5,
  alpha = 0.5, 
  fill = brewer.pal(4, "Set1")
) +
  theme_void() +
  theme(
    legend.position = "none"
  )

wrap_plots(p1, p2, p4, p3, 
          nrow = 2, 
          ncol = 2,
          heights = c(1, 2), # the more rows in the lower part, the longer it should be
          widths = c(1, 0.8),
          guides = "collect") &
  theme(legend.position = "none")


```

                    







```{r}
subres_f <- all_GEX_human_kidney %>% filter(padj < 0.05 & log2FoldChange > 2)


res_human_kid <- res3_sig[res3_sig$gene %in% subres_f$gene_name,]

p3 <- ggplot(data=res3, aes(Synovium_ECs_vs_Synovium_Tcells,Lung_Control_vs_Lung_Tcells)) + geom_point(alpha=0.6, color="grey", size=0.1) + ggtitle("Lung") +theme_minimal()+
  geom_hline(yintercept=0, linetype='dotted', col = 'red', size=0.5)+geom_vline(xintercept = 0, linetype="dotted", 
                color = "red", size=0.3) +
    geom_point(data = res_human_kid, color = "darkred",size=0.5)+labs(y= "Lung (logFC)")
p3 <- LabelPoints(plot = p3, points = c("CXCL11"), repel = TRUE, xnudge = 0.5,
  ynudge = 1.2)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme_ArchR()+ theme(axis.title.x = element_blank())


p4 <- ggplot(data=res3, aes(Synovium_Control_vs_Synovium_ECs,Lung_Control_vs_Lung_ECs)) + geom_point(alpha=0.6, color="grey", size=0.1) + ggtitle("ECs") +theme_minimal()+
  geom_hline(yintercept=0, linetype='dotted', col = 'red', size=0.5)+geom_vline(xintercept = 0, linetype="dotted", 
                color = "red", size=0.3) +
    geom_point(data = res_human_kid, color = "darkred",size=0.5)+ theme(axis.title.x = element_blank())

p4 <- LabelPoints(plot = p4, points = c("CXCL11"), repel = TRUE, xnudge = 1,
  ynudge = -1)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme_ArchR()+ theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
p4

plot <- plot_grid(p3, NULL, p4, align='vh', scale = 1, rel_widths = c(1, -0.07, 1), nrow=1)

y.grob <- textGrob("Lung (logFC)", 
                   gp=gpar( fontsize=10), rot=90)

x.grob <- textGrob("Synovium (logFC)", 
                   gp=gpar( fontsize=10))

grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))



ggdraw(add_sub(plot, "Synovium (logFC)", vpadding=grid::unit(0,"lines"),y=8, x=0.5, vjust=4.5,size=10))
```




```{r}

stim_obj <- AddModuleScore(stim_obj, features = all_markers_medfibros %>% dplyr::filter(cluster== "SPARC+COL3A1+ C4") %>% head(40) %>% rownames() %>% list(), name="endo_interacting_new")

stim_obj <- AddModuleScore(stim_obj, features = all_markers_medfibros %>% dplyr::filter(cluster== "CXCL10+CCL19+ C11") %>% head(40) %>% rownames() %>% list(), name="Tcell_interacting_new")


stim_obj$Condition
Idents(stim_obj)<-'tissue_condition'
Idents(stim_obj)<-'Condition'

DotPlot(stim_obj, features=c("RUNX1", "CXCL11", "IGF1"))


Idents(stim_obj)<-'tissue_condition'
DotPlot(stim_obj, features = c("RUNX1", "MMP14", "IGF1", "endo_interacting_new1", "Tcell_interacting_new1"), idents = c("Lung_Control", "Lung_Tcells", "Lung_ECs"))+RotatedAxis()

DotPlot(stim_obj, features = c("RUNX1", "MMP14", "IGF1", "endo_interacting_new1", "Tcell_interacting_new1"), idents = c("Synovium_Control", "Synovium_Tcells", "Synovium_ECs"))+RotatedAxis()

```


