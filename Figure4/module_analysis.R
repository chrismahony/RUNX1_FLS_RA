
```{r}
rm(list=ls()[! ls() %in% c("stia2021_rna_only","row_cl")])

gc()

#just Runx1 modules, across all clusters
```


```{r}

stia2021_rna_only$cluster_sample_condition <- paste(stia2021_rna_only$Key.marker.genes, stia2021_rna_only$sample_id, stia2021_rna_only$condition)
Idents(stia2021_rna_only) <- 'cluster_sample_condition'

stia2021_rna_only <- AddModuleScore(stia2021_rna_only, features=list(row_cl %>% filter(gene_cluster == "K3") %>% rownames()), name="K3_mod")

stia2021_rna_only <- stia2021_rna_only %>% ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA()

dotplot <- DotPlot(stia2021_rna_only, features= "K3_mod1")

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
dotplot$id_1_1 <- gsub("Chodl", "SL_Chodl", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Clu", "SL_Clu", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Col23a1", "SL_Col23a1", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("F13a1", "LL_Col22a1", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Fmo2", "SL_Fmo2", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Pi16", "SL_Pi16", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Runx2", "SL_Bglap", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Serpina3c", "SL_C3", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Sfrp1", "SL_Cfb", dotplot$id_1_1)

unique(stia2021_rna_only$cluster.name)
dotplot %>% filter(!id_1_1 %in% c("SL_Bglap", "SL_Clu")) %>% 
ggplot(aes(x=id_4, y=avg.exp.scaled)) + 
  geom_point(shape = 21)+
  geom_smooth(alpha = .6, fill="#ffb09c", color="black")+RotatedAxis()+theme(
            axis.line = element_line(),panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
        geom_vline(xintercept = c(0,3,7,15,22,28), linetype = 2, color = 'grey')+ 
        guides(color = FALSE, fill = FALSE) + 
        scale_x_continuous(breaks = c(0,3,7,15,22,28))+theme_ArchR()+facet_wrap(~factor(id_1_1, c("SL_Ccl2", "SL_Col8a1", "LL_Col22a1", "SL_Cfb","SL_Chodl", "SL_Col23a1", "SL_Pi16", "SL_C3", "SL_Ccl11", "SL_Fmo2")), nrow=2)



```


```{r}
Idents(stia2021_rna_only) <- 'cluster_sample_condition'


for (i in 1:length(unique(row_cl$gene_cluster))){
stia2021_rna_only <- AddModuleScore(stia2021_rna_only, features=list(row_cl %>% filter(gene_cluster == unique(row_cl$gene_cluster)[[i]]) %>% rownames()), name=paste0(unique(row_cl$gene_cluster)[[i]],"_mod"))
}



for (i in 1:length(unique(row_cl$gene_cluster))){
dotplot <- DotPlot(stia2021_rna_only, features= paste0(unique(row_cl$gene_cluster)[[i]],"_mod1"))

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
dotplot$id_1_1 <- gsub("Chodl", "SL_Chodl", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Clu", "SL_Clu", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Col23a1", "SL_Col23a1", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("F13a1", "LL_Col22a1", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Fmo2", "SL_Fmo2", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Pi16", "SL_Pi16", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Runx2", "SL_Bglap", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Serpina3c", "SL_C3", dotplot$id_1_1)
dotplot$id_1_1 <- gsub("Sfrp1", "SL_Cfb", dotplot$id_1_1)

print(dotplot %>% filter(!id_1_1 %in% c("SL_Bglap", "SL_Clu")) %>% 
ggplot(aes(x=id_4, y=avg.exp.scaled)) + 
  geom_point(shape = 21)+
  geom_smooth(alpha = .6, fill="#ffb09c", color="black")+RotatedAxis()+theme(
            axis.line = element_line(),panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
        geom_vline(xintercept = c(0,3,7,15,22,28), linetype = 2, color = 'grey')+ 
        guides(color = FALSE, fill = FALSE) + 
        scale_x_continuous(breaks = c(0,3,7,15,22,28))+theme_ArchR()+facet_wrap(~factor(id_1_1, c("SL_Ccl2", "SL_Col8a1", "LL_Col22a1", "SL_Cfb","SL_Chodl", "SL_Col23a1", "SL_Pi16", "SL_C3", "SL_Ccl11", "SL_Fmo2")), nrow=2)+ggtitle(paste0(unique(row_cl$gene_cluster)[[i]],"_mod")))
}

```



```{r}

cols <- ArchR::paletteDiscrete(row_cl$gene_cluster) %>% as.data.frame()


table(row_cl$gene_cluster) %>% as.data.frame() %>% ggplot(aes(y=Freq, x=Var1, fill=Var1))  + 
        geom_bar(stat = 'identity', position = position_dodge(), aes(group = Var1),width = 1) + 
               theme(
            axis.text.x = element_text(angle = 45, hjust=1),
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
            # strip.text = element_blank()
        ) + 
        guides(color = 'none', fill = 'none') + 
        labs(y = '# cells')+
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = c(0, 0)))+theme_ArchR()+scale_fill_manual(values =cols$.)+
  theme(strip.background = element_rect(fill="white", size=1, color="white"))+RotatedAxis()+ coord_flip()
```

```{r}
library(RColorBrewer)


pallete2 <- c(brewer.pal(8,"Set2"),brewer.pal(9,"Set1"))
n <- 30

row_cl


module_genes <- list()
subres_sub <- list()
for(i in 1:length(unique(row_cl$gene_cluster))){
  module_genes[[i]] <- row_cl %>% filter(gene_cluster == unique(row_cl$gene_cluster)[[i]]) %>% rownames()
  subres_sub[[i]] <-subres[subres$gene %in% module_genes[[i]],] 
    top_clust <- sort(as.matrix(subres_sub[[i]])[,"padj"]) %>% as.data.frame()
    top_clust$. <- -log(as.double(top_clust$.))
    top_clust$. <- sort(top_clust$., decreasing=T)
    top_clust$gene <- top_clust %>% rownames()
    top_clust <- top_clust[30:1,]
    print(ggplot(data=top_clust, aes(x = reorder(gene, .), y = .), color=pallete2[[i]]) +
  geom_bar(stat="identity", fill=pallete2[[i]])+ coord_flip()+theme_ArchR())}


    
  


```

```{r}


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


#PERFORM PCA ON DEgenes AND CACLULATE %VARIANCE EXPLAINED
pca = prcomp(t(scale_sub_vsd[subres $gene,]))
percent <- pca$sdev^2/sum(pca$sdev^2)*100
labs <- sapply(seq_along(percent), function(i) {paste("PC ", i, " (", round(percent[i], 1), "%)", sep="")})

#PLOT VARIANCE EXPLAINED USING BARPLOTS AND PIE CHART
par(mar=c(6,4,4,2),mfrow=c(1,3))
barplot(percent,names.arg=colnames(pca$x),las=2,xaxs="i",ylim=c(0,round(max(percent)*1.2,0)),ylab="% variance explained",main = "Variance explained by\neach Principal Component",col=colorRampPalette(c("grey85","red"))(max(round(percent,0))*1.5)[round(percent,0)+1],xpd=F,border=NA)
abline(h=c(0,5),lty=c(1,2),lwd=c(3,1))
pie(percent,labels = colnames(pca$x),col = colorRampPalette(c("grey85","Red"))(max(round(percent,0))*1.5)[round(percent,0)+1],clockwise = T,border = "white",main = "Variance explained by\neach Principal Component")

#PLOT PRINCIPAL COMPONENTS
meta <- dds@colData %>% as.data.frame()
meta$day <- c(rep(0,3), rep(2,3), rep(15,3), rep(22,3), rep(28,3), rep(7,3))

pca$x <- scale(pca$x, center=T,scale = T)
plot(pca$x[,1], pca$x[,2], xlab=labs[1], ylab=labs[2], las=1, xlim=c(-1.9,2.2), ylim=c(-2.3,2.5), cex.main=0.8,main="",type="n")
points(pca$x[,1], pca$x[,2], pch=20, cex=4, col=paste0(colorRampPalette(c("grey60",2,3))(length(meta$sample)),"95")[factor(meta$sample)])
sx <- spline(smooth.spline(1:nrow(pca$x), pca$x[,1],spar = 0.4,tol=2))
sy <- spline(smooth.spline(1:nrow(pca$x), pca$x[,2],spar = 0.4,tol=2))
lines(sx[[2]], sy[[2]], col = 2, lwd = 2)
text(pca$x[,1], pca$x[,2], meta$day, pch=20, cex=1, col="black",tck=-.05)

```
```{r}
par(mar=c(2,2,2,1),mfrow=c(2,3))
for (i in 1:6){
  plot(meta$day,pca$x[,i],bg="grey",pch=21,ylim=c(-3,3),main=labs[i],xlab="",ylab="",xaxt="n",las=1)
  myline <- smooth.spline(meta$day,pca$x[,i], spar=0.4)$y
  lines(unique(meta$day),myline,col="blue",lwd=3)
  abline(h=c(mean(pca$x[meta$day==0,i])),col="red")
}; mypar()




plot(meta$day,pca$x[,1],bg="grey",pch=21,ylim=c(-3,3),main=labs[i],xlab="",ylab="",xaxt="n",las=1)


```
```{r}
https://www.czarnewski.com/uc_classification/docs/analysis_mouse_colitis.html

cluster_contribution <- c()
groups <- row_cl$gene_cluster %>% unique() %>% length()
module_genes <- list()

for(i in 1:max(groups)){
  for(j in 1:ncol(pca$x)){
        module_genes[[i]] <- row_cl %>% filter(gene_cluster == unique(row_cl$gene_cluster)[[i]]) %>% rownames()
    cluster_contribution <- c(cluster_contribution, sum(abs(leading)[subres$gene[groups==i],j]*percent[j] ))}}
        
cluster_contribution_matrix <- matrix(cluster_contribution, nrow=ncol(pca$x), dimnames = list(c(paste0("PC",ifelse(1:ncol(pca$x)<=9,"0",""),1:ncol(pca$x))),c(paste0("Cluster",ifelse(1:max(groups)<=9,"0",""),1:max(groups)))))
```


```{r}

top_PCs <- 5

#Creating visualization plot Sankey Diagram
a <- sort(colSums(cluster_contribution_matrix[1:top_PCs,]))
a2 <- sort(cumsum(a/sum(a)),decreasing = T)

b <- sort(rowSums(cluster_contribution_matrix[1:top_PCs,]))
b2 <- sort(cumsum(b/sum(b)),decreasing = T)

nodes = data.frame(ID = c(names(b2), names(a2)), stringsAsFactors = FALSE)
nodes$x = c(rep(2,length(b)),rep(1,length(a)))
nodes$y = c((b2+c(b2[2:length(b2)],0))/2 , (a2+c(a2[2:length(a2)],0))/2)
rownames(nodes) = c(names(b2), names(a2))

edges <- data.frame(N1 = rownames(cluster_contribution_matrix[1:top_PCs,]),
                    N2 = sort(rep(colnames(cluster_contribution_matrix[1:top_PCs,]),top_PCs)),
                    Value = c(cluster_contribution_matrix[1:top_PCs,]) /sum(c(cluster_contribution_matrix[1:top_PCs,])))
edges <- edges[order(edges[,3]),]

palette = c(paste0(brewer.pal(9, "Set1"), "60"), paste0(brewer.pal(8, "Set2"), "60") )
palette = paste0(colorRampPalette(c("grey85","red"))(101),"60")
styles <- lapply(nodes$y[1:top_PCs], function(n) {
  list(col = palette[(n*10)^2+1],lty = 0, textcol = "black",srt=0)
})
styles <- c( styles, lapply(as.numeric(sub("Cluster","",names(a2))), function(n) {
  list(col = pallete2[n],lty = 0, textcol = "black",srt=0)
}))
names(styles) = nodes$ID


rp <- list(nodes = nodes, edges = edges, styles = styles)
class(rp) <- c(class(rp), "riverplot")
library(rafalib)
mypar()
plot(rp, plot_area = 0.95, yscale=0.95,line=1)


rp$edges


```
```{r}


```

