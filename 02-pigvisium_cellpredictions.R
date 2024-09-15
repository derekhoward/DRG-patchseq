library(here)
library(Seurat)
library(readr)
library(dplyr)
library(patchwork)
library(ggplot2)
source(here("utils.R"))

options(matrixStats.useNames.NA = 'deprecated')

results_path <- here("results/pig_spatial_integration")
dir.create(results_path)

pig_ref <- readRDS(here("data/pig_snRNAseq.RDS"))

# Load in pig spatial DRG data
spatial_drg <- readRDS(here('./data/price_spatial_transcriptomics/Pig_visium_rds_file/pig_clustering3_res1_REGRESS_seq123_04262023_test05302023.rds'))
spatial_drg[["percent.mt"]] <- PercentageFeatureSet(spatial_drg, features = c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB"))
spatial_pig_qc <- VlnPlot(spatial_drg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(spatial_drg@meta.data)

table(spatial_drg@meta.data$id)

# How many genes per dataset?
DefaultAssay(pig_ref)
n_distinct(rownames(pig_ref))
DefaultAssay(spatial_drg)
n_distinct(rownames(spatial_drg))

genes_intersection <- Reduce(intersect, list(rownames(pig_ref), rownames(spatial_drg)))
n_distinct(genes_intersection)

# object generation
subset_genes_pig_ref <- GetAssayData(object = pig_ref, assay = "RNA", slot = "counts")[genes_intersection, ]

subset_genes_spatial <- GetAssayData(object = spatial_drg, assay = "RNA", slot = "counts")[genes_intersection, ]

pig_atlas <- CreateSeuratObject(subset_genes_pig_ref, meta.data = pig_ref@meta.data)
spatial <- CreateSeuratObject(subset_genes_spatial, meta.data = spatial_drg@meta.data)
spatial$labels <- spatial_drg$seurat_clusters

pig_atlas$batch <- pig_atlas$orig.ident
pig_atlas$dataset <- "atlas"
spatial$batch <- spatial$id
spatial$dataset <- "spatial"

integrate_with_reference <- function(batch_data, reference_data) {
  merged_object <- merge(reference_data, y = batch_data)
  DefaultAssay(merged_object)
  data.list <- SplitObject(merged_object, split.by = "batch")
  
  # SCT normalization of batches
  for (i in names(data.list)) {
    data.list[[i]] <- SCTransform(data.list[[i]], vst.flavor = "v2", verbose = TRUE, vars.to.regress="percent.mt", return.only.var.genes = FALSE)
  }
  drg.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 6000)
  drg.list <- PrepSCTIntegration(object.list = data.list, anchor.features = drg.features)
  drg.anchors <- FindIntegrationAnchors(object.list = drg.list, normalization.method = "SCT", anchor.features = drg.features)
  # the following only works in k.weight is decreased to min number of cells in smallest sample
  drg.integrated <- IntegrateData(anchorset = drg.anchors, normalization.method = "SCT", k.weight = 75) %>% 
    RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:30)
  
  drg.integrated
}

process_animal_data <- function(animal_id, full_query_data, reference_data, max_cells_per_batch = 250) {
  min_cells_per_batch <- 75
  reference_labels_df <- as_tibble(reference_data$labels, rownames = "cell_id") %>% rename(celltype = value)
  # Subset data for the specific animal
  animal_data <- subset(x = full_query_data, subset = id == animal_id)
  
  # Total number of cells and number of batches needed
  total_cells <- ncol(animal_data)
  optimal_batch_size <- max(min_cells_per_batch, ceiling(total_cells / ceiling(total_cells / max_cells_per_batch)))
  # num_batches <- ceiling(total_cells / max_cells_per_batch)
  num_batches <- ceiling(total_cells / optimal_batch_size)
  
  # List to store processed batches
  processed_batches <- list()
  
  for (batch_idx in 1:num_batches) {
    print(strrep('*', times = 100))
    print(str_glue("     PROCESSING BATCH ", batch_idx, ' of ', num_batches))
    print(strrep('*', times = 100))
    # Calculate start and end indices for the batch
    start_idx <- (batch_idx - 1) * optimal_batch_size + 1
    end_idx <- min(batch_idx * optimal_batch_size, total_cells)
    
    # Subset batch based on indices
    batch_cells <- colnames(animal_data)[start_idx:end_idx]
    batch_data <- subset(x = animal_data, cells = batch_cells)
    
    # Integration with reference data
    integrated_data <- integrate_with_reference(batch_data, reference_data) # Define this function
    
    # Cell-type prediction
    batch_predictions <- predict_celltype_by_neighbours(integrated_data, cell_ids_to_predict = colnames(batch_data), labels_df = reference_labels_df)
    
    # Store integrated and annotated data
    processed_batches[[batch_idx]] <- batch_predictions
  }
  
  # Combine processed batches
  combined_data <- do.call(rbind, processed_batches)
  combined_data$animal_id <- animal_id
  return(combined_data)
}

library(foreach)
library(doParallel)

n_cores <- detectCores() - 1 # Reserve one core for system processes
registerDoParallel(cores = n_cores)


animal_ids <- c('pig1', 'pig2', 'pig3', 'pig4', 'pig5', 'pig6', 'pig7')
# Parallel processing
processed_animals <- foreach(animal_id = animal_ids, .combine = 'rbind', .packages = c("Seurat")) %dopar% {
  process_animal_data(animal_id = animal_id, full_query_data = spatial, reference_data = pig_atlas, max_cells_per_batch = 250)
}

# Stop the parallel cluster
stopImplicitCluster()

write_csv(processed_animals, file = here(results_path, 'iterative_spatialpig_preds.csv'))
