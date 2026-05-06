qc_counts_fine_seur_clin_ctap_anno_miles_81924 <- readRDS("/rds/projects/c/croftap-mapjagdata/AMP2_new/qc_counts_fine_seur_clin_ctap_anno_miles_81924.rds")

# Get the Assay5 object
assay5_obj <- qc_counts_fine_seur_clin_ctap_anno_miles_81924[["RNA"]]

# Extract gene and cell names
gene_names <- dimnames(assay5_obj@features)[[1]]
cell_names <- dimnames(assay5_obj@cells)[[1]]

# Pull the counts matrix
counts_matrix <- assay5_obj@layers[["counts"]]

# Set dimnames manually
dimnames(counts_matrix) <- list(gene_names, cell_names)

meta <- qc_counts_fine_seur_clin_ctap_anno_miles_81924@meta.data %>% as.data.frame()

amp2_new <- CreateSeuratObject(counts = counts_matrix, meta.data = meta)

rm(counts_matrix, meta, qc_counts_fine_seur_clin_ctap_anno_miles_81924)

amp2_new <- amp2_new %>% NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:30)




load("/rds/projects/c/croftap-mapjagdata/AMP2_new/amp2_new.RData")
rm(assay5_obj)
rm(cell_names, gene_names)


DimPlot(amp2_new, raster=F, group.by = "cell_type")

all_cells_reference <- readRDS("/rds/projects/c/croftap-celldive01/amp2/processed_output_04-11-2023/all_cells_reference.rds")
umap <- all_cells_reference$umap %>% as.data.frame()
rownames(umap) <- all_cells_reference[["meta_data"]][["cell"]]
colnames(umap) <- c("UMAP_1", "UMAP_2")
umap <- umap %>% as.matrix()
umap %>% nrow()
amp2_new[['umap']] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", global = T, assay = "RNA")

DimPlot(amp2_new, raster=F, group.by = "cell_type")
rm(list=ls()[! ls() %in% c("amp2_new")])
gc()

DimPlot(amp2_new, raster=F, group.by = "CTAP")



DimPlot(amp2_new, raster=F, group.by = "cell_type")

DotPlot(amp2_new, features = c("RUNX1", "RUNX2", "RUNX3"), group.by="cell_type")
