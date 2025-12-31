library(here)
library(dplyr)
library(Seurat)
library(SummarizedExperiment)

pig_ref <- readRDS(here("data/pig_snRNAseq.RDS"))

drg.integrated <- readRDS(here("results/SCT-integration/drg_integrated.RDS"))
patchseq_subset <- subset(drg.integrated, subset = dataset == 'patch')
source(here("utils.R"))
Patchseqobject <- prep_patchseq_obj_all_batches(convert_genes = "pig")

spatial_drg <- readRDS(here('./data/price_spatial_transcriptomics/Pig_visium_rds_file/pig_clustering3_res1_REGRESS_seq123_04262023_test05302023.rds'))
spatial_drg[["percent.mt"]] <- PercentageFeatureSet(spatial_drg, features = c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB"))
pigspatial_preds <- read.csv(here('./results/pig_spatial_integration/iterative_spatialpig_preds.csv'))
rownames(pigspatial_preds) <- pigspatial_preds$cell_id
pigspatial_preds <- pigspatial_preds %>% dplyr::select(labels.p)
spatial_drg <- AddMetaData(object = spatial_drg, metadata = pigspatial_preds)

human_spatial_drg <- readRDS(here('./data/price_spatial_transcriptomics/human.drg.combined.rds'))
human_spatial_drg[["percent.mt"]] <- PercentageFeatureSet(human_spatial_drg, pattern = "^MT-")
human_spatial_drg$orig.labels <- Idents(human_spatial_drg)
humanspatial_preds <- read.csv(here("results/human_spatial_integration/iterative_humanspatial_preds.csv"))
rownames(humanspatial_preds) <- humanspatial_preds$cell_id
humanspatial_preds <- humanspatial_preds %>% dplyr::select(labels.p)
human_spatial_drg <- AddMetaData(object = human_spatial_drg, metadata = humanspatial_preds)


convertSeuratToSummarizedExperiment <- function(seurat_object) {
  # Extract the relevant data from the Seurat object
  counts <- as.matrix(seurat_object@assays$RNA@counts) # or use the appropriate assay
  meta_data <- seurat_object@meta.data
  row_data <- rownames(seurat_object)
  
  # Create a SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = counts),
                             colData = meta_data,
                             rowData = row_data)
  
  # Return the SummarizedExperiment object
  return(se)
}

output_path <- here("data/summarized_experiments")
dir.create(output_path)
DefaultAssay(patchseq_subset) <- "RNA"
patchseq_subset_SCE <- as.SingleCellExperiment(patchseq_subset)
saveRDS(patchseq_subset_SCE, file = here(output_path, 'patchseq_subset_SCE.RDS'))
patchseq_subset_se <- convertSeuratToSummarizedExperiment(patchseq_subset)
saveRDS(patchseq_subset_se, file = here(output_path, 'patchseq_subset_se.RDS'))

Patchseqobject$labels.p <- patchseq_subset$labels.p
patchseqobject_SCE <- as.SingleCellExperiment(Patchseqobject)
saveRDS(patchseqobject_SCE, file = here(output_path, 'patchseq_object_SCE.RDS'))
patchseqobject_se <- convertSeuratToSummarizedExperiment(Patchseqobject)
saveRDS(patchseqobject_se, file = here(output_path, 'patchseq_object_se.RDS'))


pig_ref_SCE <- as.SingleCellExperiment(pig_ref)
saveRDS(pig_ref_SCE, file = here(output_path, 'pig_ref_SCE.RDS'))
pig_ref_se <- convertSeuratToSummarizedExperiment(pig_ref)
saveRDS(pig_ref_se, file = here(output_path, 'pig_ref_se.RDS'))

pig_visium_SCE <- as.SingleCellExperiment(spatial_drg)
saveRDS(pig_visium_SCE, file = here(output_path, 'pig_visium_SCE.RDS'))

pig_visium_se <- convertSeuratToSummarizedExperiment(spatial_drg)
saveRDS(pig_visium_se, file = here(output_path, 'pig_visium_se.RDS'))

human_visium_SCE <- as.SingleCellExperiment(human_spatial_drg)
saveRDS(human_visium_SCE, file = here(output_path, 'human_visium_SCE.RDS'))

human_visium_se <- convertSeuratToSummarizedExperiment(human_spatial_drg)
saveRDS(human_visium_se, file = here(output_path, 'human_visium_se.RDS'))
