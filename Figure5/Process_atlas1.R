library(ArchR)#needs r4.1
library(Nebulosa)#needs r4.1
library(scMEGA)#needs to be r4.1
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
#library(EnsDb.Mmusculus.v75)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(cowplot)
library(patchwork) #latest version is required!
library(TFBSTools)
library(JASPAR2020)
library(gsfisher)
library(harmony)
library(clustree)
library(plyr)
library(dplyr)
library(EnhancedVolcano)
library(monocle)
library(SoupX)
library(DoubletFinder)

setwd("/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool1_count")
pool1<-list.dirs(path = "/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool1_count", recursive = FALSE)
pool2<-list.dirs(path = "/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool2_count", recursive = FALSE)
pool3<-list.dirs(path = "/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool3_count", recursive = FALSE)
pool4<-list.dirs(path = "/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool4_count", recursive = FALSE)
pool5<-list.dirs(path = "/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool5_count", recursive = FALSE)

reseq<-list.dirs(path = "/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_second_seq_round_count", recursive = FALSE)
reseq <- reseq[reseq != c("/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_second_seq_round_count/old")]

controls_old<-c("/rds/projects/c/croftap-sitia-cite-seq-tc/Kemble_CITE_seq_OLD/CD45N/CD45N", "/rds/projects/c/croftap-sitia-cite-seq-tc/Kemble_CITE_seq_OLD/CD45P/CD45P")
  


dirs<-c(pool1, pool2, pool3, pool4, pool5, reseq, controls_old)

samples<-dirs


samples <- gsub('/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool1_count/', '', samples)
samples <- gsub('/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool2_count/', '', samples)
samples <- gsub('/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool3_count/', '', samples)
samples <- gsub('/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool4_count/', '', samples)
samples <- gsub('/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_count/Pool5_count/', '', samples)
samples <- gsub('/rds/projects/c/croftap-sitia-cite-seq-tc/STIA_CITEseq_second_seq_round_count/', '', samples)
samples <- gsub('/rds/projects/c/croftap-sitia-cite-seq-tc/Kemble_CITE_seq_OLD/CD45N/', '', samples)
samples <- gsub('/rds/projects/c/croftap-sitia-cite-seq-tc/Kemble_CITE_seq_OLD/CD45P/', '', samples)
samples <- gsub('/', '', samples)

samples_STIA<-samples


dirs_STIA <-paste0( dirs, "/outs/") 

db_pct_STIA<-as.numeric(c("0.024", "0.01", "0.01", "0.017", "0.008", "0.008", "0.008", "0.026", "0.008", "0.016", "0.014", "0.01", "0.014", "0.023", "0.014", "0.08", "0.02", "0.016", "0.009", "0.008", "0.08", "0.061", "0.1", "0.1", "0.048", "0.1", "0.061", "0.035", "0.1", "0.039", "0.039", "0.039", "0.054", "0.046", "0.01", "0.009"))





#cd45neg
setwd("/rds/projects/c/croftap-aia-seq-data/mono_arthritis_count (first and second seq)")
dirs_Day4 <-dir("./", pattern = "Day4")
dirs_Day2neg3 <-dir("./", pattern = "Day2_neg3")
dirs_Day2pos3 <-dir("./", pattern = "Day2_pos3")
dirs_Day7neg <-dir("./", pattern = "Day7_neg")

setwd("/rds/projects/c/croftap-aia-seq-data/D7_count_reseq_combined")
dirs_Day7pos <-dir("./", pattern="Day7")

setwd("/rds/projects/c/croftap-aia-seq-data/ThirdSeq")
dirs_Day0_2 <-dir("./", pattern = "Day0")
dirs_Day0_2 <- dirs_Day0_2[dirs_Day0_2!= c("Day0neg", "Day0pos")]
dirs_Day0_2 <- dirs_Day0_2[dirs_Day0_2!= c("Day0neg", "Day0pos")]


#B5
dirs_b5 <-dir("/rds/projects/c/croftap-mono01/AIA_count")

#B6
setwd("/rds/projects/c/croftap-mono01/AIA_count2")
dirs_b6 <-dir("/rds/projects/c/croftap-mono01/AIA_count2")


setwd("/rds/projects/c/croftap-aia-seq-data/ALL_data/all_pipline2")

all_dirs<-c(dirs_b5   ,    dirs_b6    ,   dirs_Day0_2 ,  dirs_Day2neg3, dirs_Day2pos3 ,dirs_Day4    , dirs_Day7neg , dirs_Day7pos)

samples<- all_dirs
samples_AIA<- all_dirs


dirs_b5 <-paste0("/rds/projects/c/croftap-mono01/AIA_count/", dirs_b5, "/outs/") 
dirs_b6 <-paste0("/rds/projects/c/croftap-mono01/AIA_count2/", dirs_b6, "/outs/")
dirs_Day0_2 <-paste0("/rds/projects/c/croftap-aia-seq-data/ThirdSeq/", dirs_Day0_2, "/outs/")
dirs_Day7pos <-paste0("/rds/projects/c/croftap-aia-seq-data/D7_count_reseq_combined/", dirs_Day7pos, "/outs/")
dirs_Day4 <-paste0("/rds/projects/c/croftap-aia-seq-data/mono_arthritis_count (first and second seq)/", dirs_Day4, "/outs/") 

dirs_Day2neg3 <-paste0("/rds/projects/c/croftap-aia-seq-data/mono_arthritis_count (first and second seq)/", dirs_Day2neg3, "/outs/") 
dirs_Day2pos3 <-paste0("/rds/projects/c/croftap-aia-seq-data/mono_arthritis_count (first and second seq)/", dirs_Day2pos3, "/outs/") 
dirs_Day7neg <-paste0("/rds/projects/c/croftap-aia-seq-data/mono_arthritis_count (first and second seq)/", dirs_Day7neg, "/outs/") 

dirs_AIA<-c(dirs_b5   ,    dirs_b6    ,   dirs_Day0_2 ,  dirs_Day2neg3, dirs_Day2pos3 ,dirs_Day4    , dirs_Day7neg , dirs_Day7pos)


db_pct_AIA<-as.numeric(c("0.069", "0.08", "0.076", "0.08", "0.076", "0.058", "0.069", "0.073", "0.04", "0.05", "0.069", "0.077", "0.085", "0.004", "0.058", "0.031", "0.012", "0.08", "0.01", "0.032", "0.039", "0.042", "0.026", "0.023", "0.023","0.05", "0.061", "0.12", "0.11", "0.16"))

CD45pos <-dir("/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45_pos", pattern = "pos")
CD45pos <- c("Con1_CD45pos"     ,        "Con3_CD45pos"  ,           "Infla1_CD45pos")

CD45pos_control4<-dir("/rds/projects/c/croftap-ktr-cia-sc-01/Third_seq", pattern = "CIA_Control_pos_CR5")

CD45pos_ctrl2<-dir("/rds/projects/c/croftap-ktr-cia-sc-01/SecondSeq_Pos_and_neg/COUNT", pattern = "Con2_CD45pos")

CD45pos_infla2<-dir("/rds/projects/c/croftap-ktr-cia-sc-01/SecondSeq_Pos_and_neg/COUNT", pattern = "Infla2_CD45pos")

CD45pos_infla3<-dir("/rds/projects/c/croftap-ktr-cia-sc-01/SecondSeq_Pos_and_neg/COUNT", pattern = "Infla3_CD45pos")

CD45pos_infla4<-dir("/rds/projects/c/croftap-ktr-cia-sc-01/SecondSeq_Pos_and_neg/COUNT", pattern = "Infla4_CD45pos")

all_dirs<-c(CD45pos   ,    CD45pos_control4    ,   CD45pos_ctrl2 ,  CD45pos_infla2, CD45pos_infla3 ,CD45pos_infla4)

samples<- all_dirs
samples_CIA<- all_dirs


CD45pos <-paste0("/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45_pos/", CD45pos, "/outs/") 

CD45pos_control4 <-paste0("/rds/projects/c/croftap-ktr-cia-sc-01/Third_seq/", CD45pos_control4, "/outs/") 
CD45pos_ctrl2<-paste0("/rds/projects/c/croftap-ktr-cia-sc-01/SecondSeq_Pos_and_neg/COUNT/", CD45pos_ctrl2, "/outs/")

CD45pos_infla2<-paste0("/rds/projects/c/croftap-ktr-cia-sc-01/SecondSeq_Pos_and_neg/COUNT/", CD45pos_infla2, "/outs/")

CD45pos_infla3<-paste0("/rds/projects/c/croftap-ktr-cia-sc-01/SecondSeq_Pos_and_neg/COUNT/", CD45pos_infla3, "/outs/")

CD45pos_infla4<-paste0("/rds/projects/c/croftap-ktr-cia-sc-01/SecondSeq_Pos_and_neg/COUNT/", CD45pos_infla4, "/outs/")

all_dirs<-c(CD45pos   ,    CD45pos_control4    ,   CD45pos_ctrl2 ,  CD45pos_infla2, CD45pos_infla3 ,CD45pos_infla4)


CD45neg <-c("Con1_CD45neg", "Con2_CD45neg", "Con3_CD45neg", "Infla1_CD45neg", "Infla2_CD45neg")

CD45neg2<-c("Con4_CD45neg", "Infla3_CD45neg")

all_dirs_neg<-c(CD45neg, CD45neg2)

samples_neg<- all_dirs_neg

CD45neg <-paste0("/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45neg/", CD45neg, "/outs/") 
CD45neg2 <-paste0("/rds/projects/c/croftap-ktr-cia-sc-01/SecondSeq_Pos_and_neg/COUNT/", CD45neg2, "/outs/") 

all_dirs_neg<-c(CD45neg, CD45neg2)


samples_CIA<-c(samples, samples_neg)

dirs_CIA<-c(all_dirs, all_dirs_neg)


db_pct_CIA<-as.numeric(c("0.016", "0.038", "0.012", "0.031", "0.12", "0.039", "0.09", "0.02", "0.007", "0.03" , "0.03", "0.02", "0.02", "0.024", "0.032" ))

all_samples<-c(samples_STIA, samples_AIA, samples_CIA)
all_dirs<-c(dirs_STIA, dirs_AIA, dirs_CIA)
db_pct<-c(db_pct_STIA, db_pct_AIA, db_pct_CIA)

rm(list=ls()[! ls() %in% c("all_samples","all_dirs", "db_pct")])


save.image("/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/pipeline_all_samples/all_analysis.RData")

#to configure doublet finder

dirs=all_dirs
samples=all_samples

#create seurat list and QC
data.10x = list()
for (i in 1:length(dirs)) {
data.10x[[i]] = load10X(dataDir =dirs[[i]]);
data.10x[[i]] = autoEstCont(data.10x[[i]], tfidfMin=0.5, forceAccept=TRUE);
data.10x[[i]] = adjustCounts(data.10x[[i]], roundToInt=T );
data.10x[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=3, min.features=0, project=samples[i]);
data.10x[[i]][["percent.mt"]] = PercentageFeatureSet(object=data.10x[[i]], pattern = "^mt-");
data.10x[[i]]<-subset(x = data.10x[[i]], subset = nFeature_RNA > 500 & nFeature_RNA <7000 & percent.mt < 10);
Idents(data.10x[[i]])<-'orig.ident';
data.10x[[i]]<-subset(x = data.10x[[i]], downsample = 3000);
data.10x[[i]] =NormalizeData(object = data.10x[[i]]);
data.10x[[i]] =ScaleData(object = data.10x[[i]]);
data.10x[[i]] =FindVariableFeatures(object = data.10x[[i]]);
data.10x[[i]] =RunPCA(object = data.10x[[i]], verbose = FALSE);
data.10x[[i]] =RunUMAP(object = data.10x[[i]], dims=1:30)
}

#calclate optimal pk
sweep.stats.list <- list()
optimal.pk.list<-list()
for (i in 1:length(data.10x)) {
  seu_temp <- data.10x[[i]]
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
  1:30)
  sweep.stats <- summarizeSweep(sweep.res.list)
  bcmvn <- find.pK(sweep.stats)
  sweep.stats.list[[i]] <- sweep.stats
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  optimal.pk.list[[i]] <- optimal.pk
}

#claisfy doublets

optimal.pk.list<-as.numeric(optimal.pk.list)
cleaned.list<-list()
for (i in 1:length(data.10x)) {
  seu_temp <- data.10x[[i]]
  nExp_poi <-  db_pct[[i]]*nrow(seu_temp@meta.data)
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
  seu_temp@commands$RunUMAP.RNA.pca$dims, pN = 0.25, pK = optimal.pk.list[i], nExp = nExp_poi, reuse.pANN = FALSE)
  cleaned.list[[i]] <- seu_temp
}

names(cleaned.list)<-samples


#remove doublets and re process data
for (i in 1:length(cleaned.list)) {
DF.name = colnames(cleaned.list[[i]]@meta.data)[grepl("DF.classification", colnames(cleaned.list[[i]]@meta.data))];
cleaned.list[[i]] = cleaned.list[[i]][, cleaned.list[[i]]@meta.data[, DF.name] == "Singlet"];
cleaned.list[[i]] =ScaleData(object = cleaned.list[[i]]);
cleaned.list[[i]] =FindVariableFeatures(object = cleaned.list[[i]]);
cleaned.list[[i]] =RunPCA(object = cleaned.list[[i]], verbose = FALSE)
}


anchors <- FindIntegrationAnchors(object.list = cleaned.list, reduction = "rpca",   dims = 1:30, anchor.features = 2000)
aggr <- IntegrateData(anchorset = anchors, dims = 1:30)
VlnPlot(aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
aggr <- FindVariableFeatures(aggr)
aggr <- ScaleData(aggr, verbose = FALSE)
aggr <- RunPCA(aggr, verbose = FALSE)
aggr <- RunUMAP(aggr, dims = 1:30)

aggr <- FindNeighbors(aggr, dims = 1:20)
aggr <- FindClusters(aggr, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3), graph.name = 'integrated_snn')
Idents(aggr)<-'integrated_snn_res.0.01'  
res0.01markers<-FindAllMarkers(aggr, only.pos = T)


Idents(aggr)<-'integrated_snn_res.0.05'
res0.05markers<-FindAllMarkers(aggr, only.pos = T)


Idents(aggr)<-'integrated_snn_res.0.1'
res0.1markers<-FindAllMarkers(aggr, only.pos = T)


Idents(aggr)<-'integrated_snn_res.0.2'
res0.2markers<-FindAllMarkers(aggr, only.pos = T)


save.image("/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/pipeline_all_samples/all_analysis.RData")
