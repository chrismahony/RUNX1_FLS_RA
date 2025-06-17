library(Signac)
library(Seurat)
library(sporkforlife)



data.10x = list()
dirs <- c("/rds/projects/m/mahonyc-cesar-data/inhibitor_organoids/count/DMSO", "/rds/projects/m/mahonyc-cesar-data/inhibitor_organoids/count/cbfbi")
dirs <- paste(dirs, "/outs/filtered_feature_bc_matrix/", sep="")

sample_names <- c("DMSO", "cbfbi")

process_scrna_data(dirs, sample_names, target_n_clusters = 5,
                               resolution_range = seq(0.05, 0.3, by = 0.5),
                               min_nFeature_RNA = 500, max_nFeature_RNA = 7000,
                               max_percent_mt = 10, n_dims=50)



save.image("/rds/projects/m/mahonyc-cesar-data/inhibitor_organoids/analysis/analysis.RData")





