

```{r}
options(bitmapType='cairo')

setwd("/rds/projects/c/croftap-celldive01/detectionresults/matricies")

results <- dir("./", pattern = "*txt", 
    full.names = TRUE)


data = list()
for (i in 1:length(results)) {
    data[[i]] <- read.delim(results[[i]])
      
}
names(data) <- sub('.//', '', results)

library(Seurat)

library(dplyr)

library(sctransform)

for (i in 1:length(data)) {
data[[i]] <- data[[i]] %>% select(contains('.mean'))
data[[i]] <- data[[i]] %>% select(contains('cell'))
rownames(data[[i]])<-seq_along(data[[i]][,1])
data[[i]]<-as.data.frame(t(data[[i]]))
data[[i]]<-CreateSeuratObject(data[[i]])
data[[i]]<- subset(data[[i]], subset = nCount_RNA > 0)
data[[i]]<-SCTransform(data[[i]])
}
  

BX020_PI_EarlyRA.txt<-data[["BX020_PI_EarlyRA.txt"]]
BX028_PI_Resolver.txt<-data[["BX028_PI_Resolver.txt"]]
BX031_PI_EarlyRA.txt<-data[["BX031_PI_EarlyRA.txt"]]
BX054_PI_Resolver.txt<-data[["BX054_PI_Resolver.txt"]]
BX086_Lymphoid_EstablishedRA.txt<-data[["BX086_Lymphoid_EstablishedRA.txt"]]
BX115_Diffuse_EarlyRA.txt<-data[["BX115_Diffuse_EarlyRA.txt"]]
BX127_Diffuse_EstablishedRA.txt<-data[["BX127_Diffuse_EstablishedRA.txt"]]
BX202_Lymphoid_Resolvers.txt<-data[["BX202_Lymphoid_Resolvers.txt"]]
BX240_Lymphoid_EstablishedRA.txt<-data[["BX240_Lymphoid_EstablishedRA.txt"]]
JRP115_OA.txt<-data[["JRP115_OA.txt"]]
JRP117_OA.txt<-data[["JRP117_OA.txt"]]
JRP127_OA.txt<-data[["JRP127_OA.txt"]]



all<-merge(x=BX054_PI_Resolver.txt, y=c(BX086_Lymphoid_EstablishedRA.txt, BX115_Diffuse_EarlyRA.txt, BX127_Diffuse_EstablishedRA.txt, BX202_Lymphoid_Resolvers.txt, BX240_Lymphoid_EstablishedRA.txt, BX031_PI_EarlyRA.txt, BX028_PI_Resolver.txt, BX020_PI_EarlyRA.txt, JRP115_OA.txt,JRP117_OA.txt, JRP127_OA.txt))


all<-FindVariableFeatures(all, assay = 'RNA')
all<-ScaleData(all, assay = 'RNA')
all<-RunPCA(all, assay = 'RNA')
all<-FindNeighbors(all, dims = 1:20)
all<-FindClusters(all, resolution = c(0.5), graph.name="RNA_snn")
all<-FindClusters(all, resolution = c(0.3), graph.name="RNA_snn")
all<-FindClusters(all, resolution = c(0.7), graph.name="RNA_snn")
all<-FindClusters(all, resolution = c(0.1), graph.name="RNA_snn")
all<-FindClusters(all, resolution = c(0.2), graph.name="RNA_snn")

Idents(all)<-'RNA_snn_res.0.7'
DotPlot(all, features =rownames(all)) +RotatedAxis()
table(all$RNA_snn_res.0.2)
```


```{r}
levels(all)

current.sample.ids<-c( "0" , "1",  "10", "11", "12","13", "14", "15", "16", "17", "18", "19", "2", "20", "21", "3" , "4" , "5" , "6" , "7" , "8",  "9" )
new.sample.ids<-c("mac" , "fib",  "fib", "mac", "fib","LL", "fib", "CD31_COL4A1", "fib", "Bcell", "Tcell", "pericytes", "mac", "fib", "fib", "fib" , "Tcell" , "fib" , "fib" , "mac" , "fib",  "fib" )

all$named<-all@meta.data[["RNA_snn_res.0.2"]]
                   
all@meta.data[["named"]] <- plyr::mapvalues(x = all@meta.data[["named"]], from = current.sample.ids, to = new.sample.ids)

Idents(all) <- 'named'
DotPlot(all, features =rownames(all)) +RotatedAxis()



data_x_y = list()
for (i in 1:length(results)) {
    data_x_y[[i]] <- read.delim(results[[i]])
    rownames(data_x_y[[i]])<-seq_along(data_x_y[[i]][,1])
    data_x_y[[i]] <- data_x_y[[i]] %>% select(c('Centroid.X.µm', 'Centroid.Y.µm'))
      }
names(data_x_y) <- sub('.//', '', results)


cellcodes <- as.data.frame(all@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cellcodes) <- "barcodes"
rownames(cellcodes) <- cellcodes$barcodes

library(splitstackshape)
cellcodes<-cSplit(cellcodes, splitCols = "barcodes", sep="_")

samples<-as.data.frame(c("BX054_PI_Resolver.txt", "BX086_Lymphoid_EstablishedRA.txt", "BX115_Diffuse_EarlyRA.txt", "BX127_Diffuse_EstablishedRA.txt", "BX202_Lymphoid_Resolvers.txt", "BX240_Lymphoid_EstablishedRA.txt", "BX031_PI_EarlyRA.txt", "BX028_PI_Resolver.txt", "BX020_PI_EarlyRA.txt", "JRP115_OA.txt", "JRP117_OA.txt", "JRP127_OA.txt"))
samples$library_id<-c("1" , "2" , "3" , "4" , "5" , "6" , "7" , "8",  "9",  "10", "11", "12")
colnames(samples)<-c("samples", "library_id")

cellcodes$samples <- as.vector(samples$samples[cellcodes$barcodes_2])
cellcodes<-as.data.frame(cellcodes)
cellcodes<-select(cellcodes, c(samples))
rownames(cellcodes)=colnames(all)
all<-AddMetaData(all, cellcodes)


Idents(all)<-'samples'
levels(all)

current.sample.ids<-c("BX054_PI_Resolver.txt", "BX086_Lymphoid_EstablishedRA.txt", "BX115_Diffuse_EarlyRA.txt", "BX127_Diffuse_EstablishedRA.txt", "BX202_Lymphoid_Resolvers.txt", "BX240_Lymphoid_EstablishedRA.txt", "BX031_PI_EarlyRA.txt", "BX028_PI_Resolver.txt", "BX020_PI_EarlyRA.txt", "JRP115_OA.txt", "JRP117_OA.txt", "JRP127_OA.txt")
new.sample.ids<-c("Res" , "EstablishedRA" , "Early_RA" , "EstablishedRA" , "Res" , "EstablishedRA" , "Early_RA" , "Res",  "Early_RA",  "JRP", "JRP", "JRP")

all$condition<-all$samples

all@meta.data[["condition"]] <- plyr::mapvalues(x = all@meta.data[["condition"]], from = current.sample.ids, to = new.sample.ids)

all$niche_condition<-paste(all$named, all$condition, sep = "_")

Idents(all)<-'niche_condition'
levels(all)

DotPlot(all, features = rownames(all), idents = c("fib_EstablishedRA", "fib_Early_RA"))+RotatedAxis()


```
```{r}
#best section#
BX240_Lymphoid_EstablishedRA.txt_x_y<-data_x_y[["BX240_Lymphoid_EstablishedRA.txt"]]

Idents(all)<-'samples'
BX240<-subset(all, idents = 'BX240_Lymphoid_EstablishedRA.txt')

BX240_meta<-BX240@meta.data
BX240_Lymphoid_EstablishedRA.txt_x_y$named<-BX240_meta$RNA_snn_res.0.7

cols <- ArchR::paletteDiscrete(BX240@meta.data[, "RNA_snn_res.0.7"])
ggplot(BX240_Lymphoid_EstablishedRA.txt_x_y, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols))


BX240_Lymphoid_EstablishedRA.txt_x_y_zoom<-BX240_Lymphoid_EstablishedRA.txt_x_y[BX240_Lymphoid_EstablishedRA.txt_x_y$Centroid.X.µm < 2500,]
BX240_Lymphoid_EstablishedRA.txt_x_y_zoom<-BX240_Lymphoid_EstablishedRA.txt_x_y_zoom[BX240_Lymphoid_EstablishedRA.txt_x_y_zoom$Centroid.Y.µm > 3500,]

ggplot(BX240_Lymphoid_EstablishedRA.txt_x_y_zoom, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=1.3)+theme_classic()+scale_color_manual(values = c(cols))
```


```{r}

#"BX054_PI_Resolver.txt", "BX086_Lymphoid_EstablishedRA.txt", #"BX115_Diffuse_EarlyRA.txt", "BX127_Diffuse_EstablishedRA.txt", #"BX202_Lymphoid_Resolvers.txt", "BX240_Lymphoid_EstablishedRA.txt", #"BX031_PI_EarlyRA.txt", "BX028_PI_Resolver.txt", "BX020_PI_EarlyRA.txt", #"JRP115_OA.txt"

obj="BX031_PI_EarlyRA.txt"

obj_x_y<-data_x_y[[obj]]

Idents(all)<-'samples'
s_obj<-subset(all, idents = obj)

s_obj_meta<-s_obj@meta.data
obj_x_y$named<-s_obj_meta$RNA_snn_res.0.5

cols <- ArchR::paletteDiscrete(s_obj@meta.data[, "RNA_snn_res.0.5"])
ggplot(obj_x_y, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=1)+theme_classic()+scale_color_manual(values = c(cols))



```
```{r}
Idents(all)<-'named'
vasc<-subset(all, idents = 'vasc')
VlnPlot(vasc, features = "Cell..CD31.mean", pt.size = 0)
VlnPlot(vasc, features = "Cell..COL4A1.mean", pt.size = 0)


vasc_true<-subset(all, idents = 'vasc', subset = Cell..CD31.mean > 7 )
vasc_f<-subset(all, idents = 'vasc', subset = Cell..CD31.mean <= 7 )
ncol(vasc_true)
ncol(vasc_f)
ncol(vasc)

vasc_true<-vasc_true@meta.data
vasc_true<-subset(vasc_true, select=c("named"))
vasc_true$named<-'vasc_true'

vasc_f<-vasc_f@meta.data
vasc_f<-subset(vasc_f, select=c("named"))
vasc_f$named<-'COL4A1'

vasc_meta<-rbind(vasc_true, vasc_f)

all_meta<-all@meta.data
all_meta<-subset(all_meta, select=c("orig.ident", "named"))
all_meta<-all_meta[!rownames(all_meta) %in% rownames(vasc_meta),]
all_meta$orig.ident<-NULL
all_meta<-rbind(all_meta, vasc_meta)

nrow(all_meta)
ncol(all)

colnames(all_meta)<-"named_new"
all<-AddMetaData(all, all_meta)


```


```{r}
#best section#
BX240_Lymphoid_EstablishedRA.txt_x_y<-data_x_y[["BX240_Lymphoid_EstablishedRA.txt"]]

Idents(all)<-'samples'
BX240<-subset(all, idents = 'BX240_Lymphoid_EstablishedRA.txt')

BX240_meta<-BX240@meta.data
BX240_Lymphoid_EstablishedRA.txt_x_y$named<-BX240_meta$named_new

cols <- ArchR::paletteDiscrete(BX240@meta.data[, "named_new"])
cols2

ggplot(BX240_Lymphoid_EstablishedRA.txt_x_y, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = c(cols))
```
```{r}
Idents(all) <- "named_new"
DotPlot(all, features = rownames(all))+RotatedAxis(
)



```


```{r}
#"BX054_PI_Resolver.txt", "BX086_Lymphoid_EstablishedRA.txt", #"BX115_Diffuse_EarlyRA.txt", "BX127_Diffuse_EstablishedRA.txt", #"BX202_Lymphoid_Resolvers.txt", "BX240_Lymphoid_EstablishedRA.txt", #"BX031_PI_EarlyRA.txt", "BX028_PI_Resolver.txt", "BX020_PI_EarlyRA.txt", #"JRP115_OA.txt"

obj="BX240_Lymphoid_EstablishedRA.txt"

obj_x_y<-data_x_y[[obj]]

Idents(all)<-'samples'
s_obj<-subset(all, idents = obj)

s_obj_meta<-s_obj@meta.data
obj_x_y$named<-s_obj_meta$named_new

cols <- ArchR::paletteDiscrete(s_obj@meta.data[, "named_new"])
cols2=c("green", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red", "blue")
ggplot(obj_x_y, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=1)+theme_classic()+scale_color_manual(values = c(cols))


#zoom
obj_x_y_zoom<-obj_x_y[obj_x_y$Centroid.X.µm > 1200,]
obj_x_y_zoom<-obj_x_y_zoom[obj_x_y_zoom$Centroid.Y.µm < 5900,]

obj_x_y_zoom<-obj_x_y_zoom[obj_x_y_zoom$Centroid.X.µm < 2000,]
obj_x_y_zoom<-obj_x_y_zoom[obj_x_y_zoom$Centroid.Y.µm > 4400,]

cols3=c("green", "green", "grey", "grey", "grey", "grey","grey", "red")

ggplot(obj_x_y_zoom, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=1.3)+theme_classic()+scale_color_manual(values = c(cols3))

```


```{r}
all$nichenew_condition<-paste(all$named_new, all$condition, sep = "_")

Idents(all)<-'nichenew_condition'
levels(all)

rownames(all)

DotPlot(all, features = c( "Cell..RUNX1.mean", "Cell..COL6A1.mean", "Cell..SPARC.mean"), idents = c("Fib_EstablishedRA", "Fib_Early_RA", "Fib_Res"))+RotatedAxis()+coord_flip()

DotPlot(all, features = "Cell..RUNX1.mean", idents = c("COL4A1_EstablishedRA", "COL4A1_Early_RA", "COL4A1_Res", "Fib_EstablishedRA", "Fib_Early_RA", "Fib_Res"))+RotatedAxis()+coord_flip()

DotPlot(all, features = "Cell..DKK3.mean", idents = c("COL4A1_EstablishedRA", "COL4A1_Early_RA", "COL4A1_Res", "Fib_EstablishedRA", "Fib_Early_RA", "Fib_Res"))+RotatedAxis()+coord_flip()

DotPlot(all, features = rownames(all), idents = c("COL4A1_EstablishedRA", "COL4A1_Early_RA"))+RotatedAxis()


Idents(all)<-'named_new'

DotPlot(all, features = rownames(all), idents = c("Fib", "COL4A1"))+RotatedAxis()



Idents(all)<-'nichenew_condition'

DotPlot(all, features = c("Cell..RUNX1.mean"), idents = c("Fib_JRP"                  ,   "lymphatics_JRP"      ,       
 "Mac_and_fibLL_JRP"       ,    "LL_JRP"            ,          "macs_JRP"  ,                  "COL4A1_JRP"     ,            
 "vasc_true_JRP"           ,    "Tcell_JRP"         ,          "Bcell_JRP"    ,               "MertK_NOtch_JRP"  ) )+RotatedAxis()




```
```{r}


#distance calculations

library(spatstat)

BX240_Lymphoid_EstablishedRA.txt_x_y_f<-subset(BX240_Lymphoid_EstablishedRA.txt_x_y, select=c(Centroid.X.µm, Centroid.Y.µm,named))

colnames(BX240_Lymphoid_EstablishedRA.txt_x_y_f)<-c("x", "y", "class")

xlim = range(BX240_Lymphoid_EstablishedRA.txt_x_y_f$x)
ylim = range(BX240_Lymphoid_EstablishedRA.txt_x_y_f$y)

ln = with(BX240_Lymphoid_EstablishedRA.txt_x_y_f, 
  ppp(x = x, y = y, marks = class, xrange = xlim, yrange = ylim)
)

ln

#gln = Gest(ln)
#gln
library(RColorBrewer)
#plot(gln, lty = 1, col = brewer.pal(4, "Set1"), main = "")


ANN <- apply(nndist(ln, k=1:100),2,FUN=mean, by=marks(MertK_NOtch))
plot(ANN ~ eval(1:100), type="b", main=NULL, las=1)

fib<-BX240_Lymphoid_EstablishedRA.txt_x_y_f[BX240_Lymphoid_EstablishedRA.txt_x_y_f$class=="Fib",]

xlim = range(fib$x)
ylim = range(fib$y)

ln_fib = with(fib, 
  ppp(x = x, y = y, marks = class, xrange = xlim, yrange = ylim)
)


plot(ln %mark% nndist(fib), markscale=1)


d<-nndist(ln)

nnda<-nndist(ln, by=marks(ln))

dist<-aggregate(nnda, by=list(from=marks(ln)), mean)

```
```{r}
BX020_PI_EarlyRA.txt_x_y<-data_x_y[["BX020_PI_EarlyRA.txt"]]
BX028_PI_Resolver.txt_x_y<-data_x_y[["BX028_PI_Resolver.txt"]]
BX031_PI_EarlyRA.txt_x_y<-data_x_y[["BX031_PI_EarlyRA.txt"]]
BX054_PI_Resolver.txt_x_y<-data_x_y[["BX054_PI_Resolver.txt"]]
BX086_Lymphoid_EstablishedRA.txt_x_y<-data_x_y[["BX086_Lymphoid_EstablishedRA.txt"]]
BX115_Diffuse_EarlyRA.txt_x_y<-data_x_y[["BX115_Diffuse_EarlyRA.txt"]]
BX127_Diffuse_EstablishedRA.txt_x_y<-data_x_y[["BX127_Diffuse_EstablishedRA.txt"]]
BX202_Lymphoid_Resolvers.txt_x_y<-data_x_y[["BX202_Lymphoid_Resolvers.txt"]]
BX240_Lymphoid_EstablishedRA.txt_x_y<-data_x_y[["BX240_Lymphoid_EstablishedRA.txt"]]
JRP115_OA.txt_x_y<-data_x_y[["JRP115_OA.txt"]]


BX028_PI_Resolver.txt_x_y<-BX028_PI_Resolver.txt_x_y[rownames(BX028_PI_Resolver.txt_x_y) %in% colnames(BX028_PI_Resolver.txt),]


all_x_y<-rbind(BX054_PI_Resolver.txt_x_y, BX086_Lymphoid_EstablishedRA.txt_x_y, BX115_Diffuse_EarlyRA.txt_x_y, BX127_Diffuse_EstablishedRA.txt_x_y, BX202_Lymphoid_Resolvers.txt_x_y, BX240_Lymphoid_EstablishedRA.txt_x_y, BX031_PI_EarlyRA.txt_x_y, BX028_PI_Resolver.txt_x_y, BX020_PI_EarlyRA.txt_x_y, JRP115_OA.txt_x_y)

all_count<-as.data.frame(t(as.matrix(all@assays[["RNA"]]@counts)))
head(all_count)

meta_data_all<-all@meta.data

rownames(all_x_y)=rownames(meta_data_all)
all_x_y$celltype<-meta_data_all$named_new
all_x_y$condition<-meta_data_all$condition
all_x_y$RUNX1<-all_count$Cell..RUNX1.mean
all_x_y$sample<-meta_data_all$samples


all_x_y_fibs<-all_x_y[all_x_y$celltype=="Fib",]
all_x_y_vasc<-all_x_y[all_x_y$celltype=="vasc_true",]
all_x_y_vasc_fibs<-rbind(all_x_y_fibs, all_x_y_vasc)

xlim = range(all_x_y_vasc_fibs$Centroid.X.µm)
ylim = range(all_x_y_vasc_fibs$Centroid.Y.µm)

ln_all = with(all_x_y_vasc_fibs, 
  ppp(x = Centroid.X.µm, y = Centroid.Y.µm, marks = celltype, xrange = xlim, yrange = ylim)
)


d<-nndist(ln_all)

nnda<-nndist(ln_all, by=marks(ln_all))
nnda<-as.data.frame(nnda)
nnda$RUNX1<-all_x_y_vasc_fibs$RUNX1
nnda$condition<-all_x_y_vasc_fibs$condition
nnda$samples<-all_x_y_vasc_fibs$sample
nnda$celltype<-all_x_y_vasc_fibs$celltype

nnda_BX020_PI_EarlyRA.txt<-nnda[nnda$samples=='BX020_PI_EarlyRA.txt',]
nnda_BX028_PI_Resolver.txt<-nnda[nnda$samples=='BX028_PI_Resolver.txt',]
nnda_BX031_PI_EarlyRA.txt<-nnda[nnda$samples=='BX031_PI_EarlyRA.txt',]
nnda_BX054_PI_Resolver.txt<-nnda[nnda$samples=='BX054_PI_Resolver.txt',]
nnda_BX086_Lymphoid_EstablishedRA.txt<-nnda[nnda$samples=='BX086_Lymphoid_EstablishedRA.txt',]
nnda_BX115_Diffuse_EarlyRA.txt<-nnda[nnda$samples=='BX115_Diffuse_EarlyRA.txt',]
nnda_BX127_Diffuse_EstablishedRA.txt<-nnda[nnda$samples=='BX127_Diffuse_EstablishedRA.txt',]
nnda_BX202_Lymphoid_Resolvers.txt<-nnda[nnda$samples=='BX202_Lymphoid_Resolvers.txt',]

nnda_BX240_Lymphoid_EstablishedRA.txt<-nnda[nnda$samples=='BX240_Lymphoid_EstablishedRA.txt',]
nnda_JRP115_OA.txt<-nnda[nnda$samples=='JRP115_OA.txt',]


all_list<-list(nnda_BX020_PI_EarlyRA.txt, nnda_BX028_PI_Resolver.txt, nnda_BX031_PI_EarlyRA.txt, nnda_BX054_PI_Resolver.txt, nnda_BX086_Lymphoid_EstablishedRA.txt,nnda_BX115_Diffuse_EarlyRA.txt, nnda_BX127_Diffuse_EstablishedRA.txt, nnda_BX202_Lymphoid_Resolvers.txt, nnda_BX240_Lymphoid_EstablishedRA.txt, nnda_JRP115_OA.txt )

for (i in 1:length(all_list)) {
    all_list[[i]]$Fib <- all_list[[i]]$Fib/max(all_list[[i]]$Fib)
        }


names(all_list)=c("nnda_BX020_PI_EarlyRA.txt", "nnda_BX028_PI_Resolver.txt", "nnda_BX031_PI_EarlyRA.txt", "nnda_BX054_PI_Resolver.txt", "nnda_BX086_Lymphoid_EstablishedRA.txt","nnda_BX115_Diffuse_EarlyRA.txt", "nnda_BX127_Diffuse_EstablishedRA.txt", "nnda_BX202_Lymphoid_Resolvers.txt", "nnda_BX240_Lymphoid_EstablishedRA.txt", "nnda_JRP115_OA.txt" )


all_data_norm<-rbind(all_list)


library(plyr)
df <- ldply(all_list, data.frame)

df_fibs<-df[df$celltype=='Fib',]
df_fibs_nojrp<-df_fibs[!df_fibs$condition=='JRP',]


ggplot(df_fibs_nojrp, aes(Fib, RUNX1, color=condition)) + geom_smooth(alpha=.1) + labs(x="Distance (towards lumen)", y="Cell Type Signal") + theme_light(base_size = 16)+theme_classic()




nnda

nnda$Fib<- nnda$Fib/max(nnda$Fib)

ggplot(df, aes(Fib, RUNX1), color=celltype) + geom_smooth(alpha=.1) + labs(x="Distance (towards lumen)", y="Cell Type Signal") + theme_light(base_size = 16)+theme_classic()

dist<-aggregate(nnda, by=list(from=marks(ln_all)), mean)

#need to pull each sample and normlaize for total distance!
```
```{r}

all_x_y_vasc<-all_x_y[all_x_y$celltype=='vasc_true',]

library(sp)
for ( spot in 1:nrow(all_x_y)){
  
  dists <- spDistsN1(as.matrix(all_x_y_vasc[, 2:1]), pt=as.numeric(all_x_y[spot, c( "Centroid.X.µm", "Centroid.Y.µm")]))
  all_x_y$dist1[spot] <- min(dists)
  all_x_y$dist.x1[spot] <- all_x_y_vasc[which(dists == min(dists))[1], 1]
  all_x_y$dist.y1[spot] <- all_x_y_vasc[which(dists == min(dists))[1], 2]
    

}


```
```{r}


Idents(all)<-'named_new'
levels(all)<-c( "COL4A1", "Fib"   ,   "Mac_and_fibLL", "macs" ,         "lymphatics" , "LL"  ,          "Tcell"  ,       "MertK_NOtch"  , "Bcell" , "vasc_true"         )
dotplot<-DotPlot(all, features='Cell..RUNX1.mean', idents = c( "COL4A1", "Fib" ))
dotplot<-dotplot$data

ggplot(dotplot, aes(x=avg.exp, y=id, fill=id))+
geom_bar(stat="identity", color="black")+
scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme_minimal()


median.stat <- function(x){
    out <- quantile(x, probs = c(0.5))
    names(out) <- c("ymed")
    return(out) 
}


VlnPlot(all, features='Cell..RUNX1.mean', idents = c( "COL4A1", "Fib" ), pt.size = 0) +stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black") 

VlnPlot(all, features='Cell..DKK3.mean', idents = c( "COL4A1", "Fib" ), pt.size = 0) +stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black") 
```
```{r}

Idents(all)<-'nichenew_condition'
levels(all)<-c( "Fib_Early_RA"             ,"Fib_EstablishedRA"       ,
                "Fib_Res"                  ,
                "COL4A1_Early_RA"          ,
                "COL4A1_EstablishedRA"      ,
                "COL4A1_Res"                  ,
  
  "macs_Res"   ,                 "lymphatics_Res",             
  "LL_Res"                    ,    
  "MertK_NOtch_Res"            , "Mac_and_fibLL_Res"        ,  
  "Tcell_Res"                ,  
  "Bcell_Res"             ,      "vasc_true_Res"            ,  
     "LL_EstablishedRA"          , 
 "macs_EstablishedRA"       ,    
 "Tcell_EstablishedRA"       ,  "MertK_NOtch_EstablishedRA" , 
 "Mac_and_fibLL_EstablishedRA", "lymphatics_EstablishedRA"  , 
 "Bcell_EstablishedRA"      ,   "vasc_true_EstablishedRA"   , 
 "macs_Early_RA"             , 
    "LL_Early_RA"               , 
 "vasc_true_Early_RA"       ,   "Tcell_Early_RA"            , 
 "Mac_and_fibLL_Early_RA"   ,   "MertK_NOtch_Early_RA"      , 
 "Bcell_Early_RA"           ,   "lymphatics_Early_RA"       , 
 "Fib_JRP"                  ,   "lymphatics_JRP"            , 
 "Mac_and_fibLL_JRP"        ,   "LL_JRP"                    , 
 "macs_JRP"                 ,   "COL4A1_JRP"                , 
 "vasc_true_JRP"            ,   "Tcell_JRP"                  ,
 "Bcell_JRP"                ,   "MertK_NOtch_JRP")

rownames(all)

DotPlot(all, features = c( "Cell..RUNX1.mean", "Cell..COL6A1.mean", "Cell..SPARC.mean"), idents = c("Fib_EstablishedRA", "Fib_Early_RA", "Fib_Res"))+RotatedAxis()+coord_flip()

VlnPlot(all, features = "Cell..RUNX1.mean", idents = c("COL4A1_EstablishedRA", "COL4A1_Early_RA", "COL4A1_Res"), pt.size = 0, cols = c("#90D5E4", "#90D5E4", "#90D5E4"))+stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black")

VlnPlot(all, features = "Cell..RUNX1.mean", idents = c( "Fib_EstablishedRA", "Fib_Early_RA", "Fib_Res"), pt.size = 0, cols = c("#D51F26", "#D51F26", "#D51F26"))+stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black")



VlnPlot(all, features = "Cell..DKK3.mean", idents = c("COL4A1_EstablishedRA", "COL4A1_Early_RA", "COL4A1_Res"), pt.size = 0, cols = c("#90D5E4", "#90D5E4", "#90D5E4"))+stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black")

VlnPlot(all, features = "Cell..DKK3.mean", idents = c( "Fib_EstablishedRA", "Fib_Early_RA", "Fib_Res"), pt.size = 0, cols = c("#D51F26", "#D51F26", "#D51F26"))+stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black")


```


```{r}
table(all_x_y$condition)
table(all_x_y$celltype)

all_x_y_fibs<-all_x_y[all_x_y$celltype=="Fib",]
all_x_y_vasc<-all_x_y[all_x_y$celltype=="vasc_true",]
all_vasc_fibs<-rbind(all_x_y_fibs, all_x_y_vasc)

all_x_y_early<-all_vasc_fibs[all_vasc_fibs$condition=="Early_RA",]
all_x_y_est<-all_vasc_fibs[all_vasc_fibs$condition=="EstablishedRA",]


xlim = range(all_x_y_early$Centroid.X.µm)
ylim = range(all_x_y_early$Centroid.Y.µm)

ln_all_early = with(all_x_y_early, 
  ppp(x = Centroid.X.µm, y = Centroid.Y.µm, marks = celltype, xrange = xlim, yrange = ylim)
)


ANN <- apply(nndist(ln_all_early, k=1:100),2,FUN=mean, by=marks(vasc_true))
plot(ANN ~ eval(1:100), type="b", main=NULL, las=1)

d<-nndist(ln_all_early)

nnda<-nndist(ln_all_early, by=marks(ln_all))

```

```{r}

library(dbscan)

all_split<-SplitObject(all, split.by = "samples")

BX240_Lymphoid_EstablishedRA.txt_x_y_nn<-BX240_Lymphoid_EstablishedRA.txt_x_y
BX240_Lymphoid_EstablishedRA.txt_x_y_nn$named<-NULL

nn <- frNN(x= BX240_Lymphoid_EstablishedRA.txt_x_y_nn, eps = 50)

dat<-BX240_Lymphoid_EstablishedRA.txt_x_y_nn[nn$id$`7722`, ]

BX240_Lymphoid_EstablishedRA.txt_x_y_nn['7722',]


BX240_Lymphoid_EstablishedRA.txt_x_y_nn %>%
        filter(rownames(BX240_Lymphoid_EstablishedRA.txt_x_y_nn) %in% rownames(dat))%>%
        ggplot(aes(x=Centroid.X.µm, y = Centroid.Y.µm)) +
        geom_point() +
        ggforce::geom_circle(aes(x0 = 5212.2 , y0 = 1555, r = 50)) +
        geom_point(data = as.data.frame(BX240_Lymphoid_EstablishedRA.txt_x_y_nn['7722', , drop=FALSE]), aes(x=Centroid.X.µm, y=Centroid.Y.µm), color = "red", size = 3) +
        coord_fixed()

```
```{r}
x<- purrr::map(nn$id, ~all_split[["JRP115_OA.txt"]]@meta.data[["named_new"]][.x] %>% table())

nn_matrix<- do.call(rbind,x)

head(nn_matrix)

nn_obj<- CreateSeuratObject(counts = t(nn_matrix),  min.features = 1)
nn_obj<- SCTransform(nn_obj, vst.flavor = "v2")
nn_obj <- RunPCA(nn_obj, npcs = 30, features = rownames(nn_obj))
nn_obj <- FindNeighbors(nn_obj, reduction = "pca", dims = 1:5)
nn_obj <- FindClusters(nn_obj, resolution = 0.1)



old_meta<- all_split[["BX240_Lymphoid_EstablishedRA.txt"]]@meta.data %>% 
        tibble::rownames_to_column(var= "cell_id")

nn_meta<- nn_obj@meta.data %>%
        tibble::rownames_to_column(var= "cell_id") %>%
        select(cell_id, SCT_snn_res.0.001)

table(nn_meta$SCT_snn_res.0.1)

## note, we filtered out some cells for the neighborhood analysis

new_meta<- old_meta

new_meta$SCT_snn_res.0.1<-nn_meta$SCT_snn_res.0.1

new_meta<- as.data.frame(new_meta)
rownames(new_meta)<- old_meta$cell_id

all_split[["BX240_Lymphoid_EstablishedRA.txt"]]@meta.data<- new_meta

BX240_Lymphoid_EstablishedRA.txt_x_y$nn_group<-new_meta$SCT_snn_res.0.001



table(BX240_Lymphoid_EstablishedRA.txt_x_y$nn_group)
```


```{r}
Idents(all)<-'RNA_snn_res.0.3'
DotPlot(all, features =rownames(all)) +RotatedAxis()
```


```{r}

Idents(all)<-'RNA_snn_res.0.3'

current.sample.ids<-levels(all)


new.sample.ids <- c( "DKK_fibs" , "mac_fibLL",  "COL4A1fibs", "COL6A1fibs", "vascfibs","macs", "SL_fibs", "lymphatics","COL6A1fibs", "sparc_fibs", "CLIC5_fibs", "Bcells", "Tcells", "MertK_notch", "CD34fibs", "PRG4_fibs", "Bcells" , "Fap_fibs" , "Tcells" ,"CD34fibs", "Tcell_interacting_fibs" )


all$highres_named <-all$RNA_snn_res.0.3


all@meta.data[["highres_named"]] <- plyr::mapvalues(x = all@meta.data[["highres_named"]], from = current.sample.ids, to = new.sample.ids)
```


```{r}
Idents(all)<-'highres_named'
vascfibs<-subset(all, idents = 'vascfibs')
VlnPlot(vascfibs, features = "Cell..CD31.mean", pt.size = 0)
VlnPlot(vascfibs, features = "Cell..COL4A1.mean", pt.size = 0)


vasc_true2<-subset(all, idents = 'vascfibs', subset = Cell..CD31.mean > 7 )
vasc_f2<-subset(all, idents = 'vascfibs', subset = Cell..CD31.mean <= 7 )


vasc_true2<-vasc_true2@meta.data
vasc_true2<-subset(vasc_true2, select=c("highres_named"))
vasc_true2$highres_named<-'vasc_true'

vasc_f2<-vasc_f2@meta.data
vasc_f2<-subset(vasc_f2, select=c("highres_named"))
vasc_f2$highres_named<-'COL4A1'

vasc_meta2<-rbind(vasc_true2, vasc_f2)

all_meta<-all@meta.data
all_meta<-subset(all_meta, select=c("orig.ident", "highres_named"))
all_meta<-all_meta[!rownames(all_meta) %in% rownames(vasc_meta2),]
all_meta$orig.ident<-NULL
vasc_meta2$named <- NULL
all_meta<-rbind(all_meta, vasc_meta2)

nrow(all_meta)
ncol(all)

colnames(all_meta)<-"highres_named"
all<-AddMetaData(all, all_meta)

Idents(all) <- 'highres_named'
DotPlot(all, features =rownames(all)) +RotatedAxis()


Idents(all) <- 'RNA_snn_res.0.3'
DotPlot(all, features =rownames(all)) +RotatedAxis()
```


```{r}
Idents(all)<-'highres_named'
vascfibs<-subset(all, idents = 'mac_fibLL')
VlnPlot(vascfibs, features = "Cell..CD68.mean", pt.size = 0)
VlnPlot(vascfibs, features = "Cell..PDPN.mean", pt.size = 0)


vasc_true2<-subset(all, idents = 'mac_fibLL', subset = Cell..CD68.mean > 7.5 )
vasc_f2<-subset(all, idents = 'mac_fibLL', subset = Cell..CD68.mean <= 7.5 )


vasc_true2<-vasc_true2@meta.data
vasc_true2<-subset(vasc_true2, select=c("highres_named"))
vasc_true2$highres_named<-'LLmacs'

vasc_f2<-vasc_f2@meta.data
vasc_f2<-subset(vasc_f2, select=c("highres_named"))
vasc_f2$highres_named<-'PDPN_fibs'

vasc_meta2<-rbind(vasc_true2, vasc_f2)

all_meta<-all@meta.data
all_meta<-subset(all_meta, select=c("orig.ident", "highres_named"))
all_meta<-all_meta[!rownames(all_meta) %in% rownames(vasc_meta2),]
all_meta$orig.ident<-NULL
vasc_meta2$named <- NULL
all_meta<-rbind(all_meta, vasc_meta2)

nrow(all_meta)
ncol(all)

colnames(all_meta)<-"highres_named"
all<-AddMetaData(all, all_meta)

Idents(all) <- 'highres_named'
DotPlot(all, features =rownames(all)) +RotatedAxis()
```






```{r}
Idents(all) <- 'highres_named'

levels(all)

current.sample.ids <- c( "DKK_fibs"    ,          "COL4A1fibs"   ,          "COL6A1fibs"     ,        "macs" ,         "SL_fibs"     ,           "lymphatics"          ,   "sparc_fibs"     ,       "CLIC5_fibs"    ,         "Bcells"    ,            
"Tcells"        ,         "MertK_notch"     ,       "CD34fibs"     ,          "PRG4_fibs"   ,           "Fap_fibs",              
"Tcell_interacting_fibs" ,"vasc_true"      ,        "COL4A1"  ,  "LLmacs",  "PDPN_fibs" )


new.sample.ids <- c( "SL_fibs1"    ,             "pericytes"   ,          "SL_fibs2"     ,        "LLmacs" ,         "SL_fibs2"     ,           "lymphatics"          ,   "SL_fibs2"     ,       "LL_fibs"    ,         "Bcells"    ,            
"Tcells"        ,         "pericytes"     ,       "SL_fibs2"     ,          "LL_fibs"   ,           "SL_fibs1",              
"SL_fibs2" ,"vasc_true"      ,        "pericytes" ,  "LLmacs",  "LL_fibs" )


all$highres_grouped <-all$highres_named


all@meta.data[["highres_grouped"]] <- plyr::mapvalues(x = all@meta.data[["highres_grouped"]], from = current.sample.ids, to = new.sample.ids)

BX240_Lymphoid_EstablishedRA.txt_x_y<-data_x_y[["BX240_Lymphoid_EstablishedRA.txt"]]

Idents(all)<-'samples'
BX240<-subset(all, idents = 'BX240_Lymphoid_EstablishedRA.txt')

BX240_meta<-BX240@meta.data
BX240_Lymphoid_EstablishedRA.txt_x_y$named<-BX240_meta$highres_grouped

cols <- ArchR::paletteDiscrete(BX240@meta.data[, "highres_grouped"])
cols2

ggplot(BX240_Lymphoid_EstablishedRA.txt_x_y, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = c(cols))+xlim(2000,4000)+ylim(3500,6000)

#data_x_y_safe <- data_x_y
#data_x_y <- data_x_y_safe
samples <- list()
Idents(all)<-'samples'
subseting <- names(data_x_y)
meta_list <- list()
cols <- ArchR::paletteDiscrete(BX240@meta.data[, "highres_grouped"])
ggplots <- list()
for (i in 1:length(data_x_y)){
  samples[[i]]<-subset(all, idents = subseting[[i]])
  meta_list[[i]] <- samples[[i]]@meta.data}

#data_x_y[[2]] <- data_x_y[[2]][-c(23131), ]

for (i in 1:length(data_x_y)){
  #rownames(data_x_y[[2]]) <- paste(rownames(data_x_y[[2]]), "_8", sep="")
  #data_x_y[[2]] <- data_x_y[[2]][rownames(data_x_y[[2]]) %in% rownames(meta_list[[2]]),]
  #data_x_y[[i]] <- data_x_y[[i]][rownames(data_x_y[[i]]) %in% rownames(meta_list[[i]]),]
   # meta_list[[i]] <- meta_list[[i]][rownames(meta_list[[i]]) %in% rownames(data_x_y[[i]])]
data_x_y[[i]]$named<-meta_list[[i]]$highres_grouped
ggplots[[i]] <-   ggplot(data_x_y[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = c(cols))
print(ggplots[[i]])}


cols <- as.data.frame(cols)

Idents(BX240) <- 'highres_grouped'
levels(BX240) <- c("SL_fibs1" ,"LLmacs"  , "pericytes" , "SL_fibs2",   "Tcells",  "vasc_true",    "Bcells"    ,  "LL_fibs",  "lymphatics")
BX240$highres_grouped <- BX240@active.ident

cols <- ArchR::paletteDiscrete(BX240@meta.data[, "highres_grouped"])

ggplot(data_x_y[[8]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = cols)+xlim(4400,5600)+ylim(3700,5200)
```


```{r}
ggplot(data_x_y[[9]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.3)+theme_classic()+scale_color_manual(values = cols)+xlim(2700,3800)+ylim(3600,5200)

ggplot(data_x_y[[9]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.3)+theme_classic()+scale_color_manual(values = cols_fibs)+xlim(2700,3800)+ylim(3600,5200)


ggplot(data_x_y[[9]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = cols)+xlim(450,2000)+ylim(3800,5600)


ggplot(data_x_y[[8]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = cols)+xlim(2800,6000)+ylim(3700,6600)



ggplot(data_x_y[[1]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = cols)+xlim(2500,4000)+ylim(3600,5200)


cols_fibs <- c("#D51F26" ,  "#208A42" , "#C06CAB" ,"grey" , "grey" , "grey",  "grey" , "grey" , "#FEE500" )


ggplot(data_x_y[[8]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = cols)+xlim(4400,5600)+ylim(3700,5200)

ggplot(data_x_y[[8]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = cols_fibs)+xlim(4400,5600)+ylim(3700,5200)


ggplot(data_x_y[[9]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.5)+theme_classic()+scale_color_manual(values = cols_fibs)+xlim(2500,4000)+ylim(3600,5200)

```


```{r}
Idents(all) <- 'highres_grouped'
dotplot<-DotPlot(all, features =c("Cell..CD34.mean", "Cell..CD90.mean", "Cell..CLIC5.mean", "Cell..COL4A1.mean", "Cell..COL6A1.mean", "Cell..DKK3.mean", "Cell..FAP.mean", "Cell..NOTCH3.mean", "Cell..PDGFRA.mean", "Cell..PDPN.mean", "Cell..PRG4.mean", "Cell..SPARC.mean"), idents=levels(all)[c(1,2,3,6)]) +RotatedAxis()




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

Heatmap(dotplot)



Idents(all) <- 'RNA_snn_res.0.3'
DotPlot(all, features =rownames(all)) +RotatedAxis()


Idents(stia2021_rna) <- 'cluster.name'
genes <- c("CD34", "THY1", "CLIC5", "COL4A1", "COL6A1", "DKK3", "FAP", "NOTCH3", "PDGFRA", "PDPN", "PRG4", "RUNX1", "SPARC", "TC")
genes <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(genes), perl=TRUE)
dotplot<-DotPlot(stia2021_rna, features =genes, idents=levels(stia2021_rna)[-c(2,5)]) +RotatedAxis()


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

Heatmap(dotplot)

```
```{r}
Idents(all) <- 'highres_grouped'
dotplot<-DotPlot(all, features =c("Cell..PDGFRA.mean", "Cell..COL4A1.mean", "Cell..COL6A1.mean", "Cell..TREM2.mean", "Cell..LYVE1.mean", "Cell..PRG4.mean", "Cell..CD206.mean", "Cell..CD4.mean", "Cell..CD31.mean")) +RotatedAxis()

rownames(all)


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

Heatmap(dotplot, cluster_rows = F, cluster_columns = F)

```


```{r}

levels(all)

VlnPlot(all, features = "Cell..RUNX1.mean", idents=levels(all)[c(1,3)], pt.size = 0)


```






```{r}


all$niche_grouped_condition <- paste(all$highres_grouped, all$condition, sep="_")
Idents(all) <- 'niche_grouped_condition'
levels(all)
VlnPlot(all, features = "Cell..RUNX1.mean", idents=levels(all)[c(5,10,21)], pt.size = 0)+stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black") 

VlnPlot(all, features = "Cell..RUNX1.mean", idents=levels(all)[c(3,11,24)], pt.size = 0)+stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black") 
VlnPlot(all, features = "Cell..RUNX1.mean", idents=levels(all)[c(4,12,21)], pt.size = 0)+stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black") 


VlnPlot(all, features='Cell..RUNX1.mean', idents = c( "COL4A1", "Fib" ), pt.size = 0) +stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black") 

Idents(all) <- 'highres_grouped'
levels(all)
VlnPlot(all, features = "Cell..RUNX1.mean", idents=levels(all)[c(1,3)], pt.size = 0)

```

```{r}

#run proximity analysis and look at expression, hopefully show sl_fibs2 are close to vessels

library(devtools)
library(spatula)
library(furrr)
 
all_x_y <- rbindlist(data_x_y)
 
#need to add a column to all_x_y with the ident for each cell. in this case mine is called 'named'
 
coloc_res_coarse<- coloc_all_types(
        index_type = unique(all_x_y$named),
        coords = all_x_y[, c("Centroid.X.µm", "Centroid.Y.µm")],
        y = all_x_y$named,
        compartments = NULL,
        max_dist = 40,
        nperm = 1000,
        parallel = TRUE
    )
 
 save.image("/rds/projects/c/croftap-celldive01/detectionresults/matricies/NEWannotation.RData")
library(data.table)
library(splitstackshape)
 
library(ggpubr)
library(rstatix)
library(DescTools)
 
 
 
plt_df<-coloc_res_coarse %>%
    subset(pval < 0.05) %>%
    dplyr::select(index_type, type, zscore) %>%
    spread(type, zscore, fill = 0) %>%
    column_to_rownames('index_type') %>%
    as.matrix
 
plt_df %>%
    Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = T, cluster_rows = T)
 
pltnew <- plt_df+t(plt_df)

pltnew %>%
    Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = T, cluster_rows = T)

```
```{r}

prop_table <- table(all$samples, all$highres_grouped) %>% as.data.frame() %>% 
    pivot_wider(names_from = Var2, values_from = Freq) %>% 
  as.data.frame() 

rownames(prop_table) <- prop_table$Var1
prop_table$Var1 <- NULL

prop_table <- prop_table/rowSums(prop_table)


Idents(all) <- 'samples'

dotplot <- DotPlot(all, features = "Cell..RUNX1.mean")

dotplot <- dotplot$data

table(all$samples) %>% as.data.frame()

dotplot$Var1 <- dotplot$id
prop_table$Var1 <- rownames(prop_table)


combined_df <- merge(prop_table, dotplot, by = "Var1")


ggplot(combined_df, aes(x=Tcells, y=avg.exp.scaled)) + 
  geom_point()+
  geom_smooth(method=lm, se=T)




ggscatter(combined_df[-c(2,4,8,10),], x = "Tcells", y = "avg.exp.scaled",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
          )+
  stat_cor(method = "pearson")+theme_ArchR()


```
```{r}


df_new <- combined_df[-c(2,4,8,10),]

df_new$condition <- c("early", "early", "est", "early", "est", "est")

ggplot(df_new, aes(x=condition, y=Tcells))+
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Cluster1")+theme_classic()+RotatedAxis()
cluster1
res.aov_cluster1 <- aov(cluster1` ~ condition_2, data = pt)
summary(res.aov_cluster1)
stats_cluster1 <- TukeyHSD(res.aov_cluster1)
stats_cluster1


```
```{r}

all$condition %>% unique()
Idents(all) <- 'highres_grouped'
levels(all)

all$condition_grouped <- paste(all$condition, all$highres_grouped, sep="_")

Idents(all) <- 'condition_grouped'

to_plot <- levels(all)[grep("SL_fibs2",levels(all))]

VlnPlot(all, features="Cell..RUNX1.mean", idents=to_plot[c(2,3)], pt.size = 0)+stat_summary(fun.y = median.stat, geom='point', size = 1, colour = "black") 


markers_df <- FindMarkers(all, ident.1 = "EstablishedRA_SL_fibs2", ident.2 = "Early_RA_SL_fibs2", recorrect_umi=FALSE)



Idents(all) <- 'condition'
all_f2 <- subset(all, idents=levels(all)[c(2,3)])
test <- sc_utils(all_f2)


library(scProportionTest)
prop.test_1 <- permutation_test(test, cluster_identity = "highres_grouped", sample_1="Early_RA", sample_2="EstablishedRA", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test_1, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = T)#+
  #geom_point(size=5,color = c("grey", "grey", "grey", "red", "red"))
##




```

