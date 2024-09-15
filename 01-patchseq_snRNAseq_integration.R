library(here)
library(Seurat)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(purrr)

# ======================== Setup and Configuration ========================

source(here("utils.R"))
set.seed(1)

theme_set(theme_cowplot())

results_path <- here("results/SCT-integration")
dir.create(results_path, recursive = TRUE)
figs_path <- here(results_path, "figures")
dir.create(figs_path)

# Define color palettes and other constants
dataset_colors <- c("atlas" = "darkgrey", "patch" = "firebrick")
mt_genes <- c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB")

# ======================== Data Loading and Preprocessing ========================

# Load reference data
pig_ref <- readRDS(here("data/pig_snRNAseq.RDS"))

# Calculate mitochondrial percentage
pig_ref[["percent.mt"]] <- PercentageFeatureSet(pig_ref, features = mt_genes)

# Prepare patch-seq data
Patchseqobject <- prep_patchseq_obj_all_batches(convert_genes = "pig")
Patchseqobject[["percent.mt"]] <- PercentageFeatureSet(Patchseqobject, features = mt_genes)

# ======================== Quality Control ========================

# Visualize QC metrics for reference data
Idents(pig_ref) <- pig_ref$orig.ident
ref_qc_plot <- VlnPlot(pig_ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), cols = reference_celltype_colors)
ggsave(here(figs_path, "01-QC_snrnaseq_batch.png"), ref_qc_plot, width = 9, height = 6, dpi = 300)

# Visualize QC metrics for patch-seq data
Patchseqobject$pseq_batch <- factor(Idents(Patchseqobject), levels = c("p1", "p2", "p3", "p4", "p5", "p6"))
Idents(Patchseqobject) <- Patchseqobject$pseq_batch
patch_qc_plot <- VlnPlot(Patchseqobject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 1)
ggsave(here(figs_path, "01-QC_patchseq_batch.png"), patch_qc_plot, width = 18, height = 8, dpi = 300)

# ======================== Integration ========================

# Find common genes
genes_intersection <- intersect(rownames(Patchseqobject), rownames(pig_ref))

# Create Seurat objects with common genes
pig_atlas <- CreateSeuratObject(GetAssayData(pig_ref, assay = "RNA", slot = "counts")[genes_intersection, ],
                                meta.data = pig_ref@meta.data)
patchseq <- CreateSeuratObject(GetAssayData(Patchseqobject, assay = "RNA", slot = "counts")[genes_intersection, ],
                               meta.data = Patchseqobject@meta.data)

# Add metadata
pig_atlas$batch <- pig_atlas$orig.ident
pig_atlas$dataset <- "atlas"
patchseq$batch <- "p"
patchseq$dataset <- "patch"
patchseq$labels <- "patchseq"

# Merge objects
merged_object <- merge(pig_atlas, y = patchseq)

# Perform integration
merged_object <- SCTransform(merged_object, vars.to.regress = "percent.mt", verbose = FALSE)
merged_object <- RunPCA(merged_object, verbose = FALSE)
merged_object <- RunUMAP(merged_object, dims = 1:30, verbose = FALSE)

# Visualize merged data before integration
p_merged <- DimPlot(merged_object, group.by = "dataset", cols = dataset_colors, label = FALSE) +
  ggtitle("Before Integration")
ggsave(here(figs_path, "01-UMAP_merged_before_integration.png"), p_merged, width = 10, height = 8, dpi = 300)
ggsave(here(figs_path, "01-UMAP_merged_before_integration.svg"), p_merged, width = 10, height = 8, dpi = 300)

# Perform integration
data.list <- SplitObject(merged_object, split.by = "batch")
data.list <- data.list[c(1, 2, 3, 4, 5, "p")]

options(matrixStats.useNames.NA = 'deprecated')

# SCT normalization of batches
data.list <- lapply(data.list, function(x) {
  SCTransform(x, vst.flavor = "v2", verbose = TRUE, vars.to.regress="percent.mt")
})

# Integration
integration_features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 6000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = integration_features)
integration_anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = integration_features)
integrated_data <- IntegrateData(anchorset = integration_anchors, normalization.method = "SCT", k.weight = 75)

# ======================== Post-integration Analysis ========================

# Dimension reduction
integrated_data <- RunPCA(integrated_data, verbose = FALSE)
integrated_data <- RunUMAP(integrated_data, dims = 1:30)

# Find neighbors
integrated_data <- FindNeighbors(integrated_data, return.neighbor = TRUE, reduction = "pca", dims = 1:30)

# Identify neighbors for patch-seq cells
pseq_neighbours <- map(colnames(patchseq), ~TopNeighbors(integrated_data@neighbors$integrated.nn, cell = ., n = 20))
pseq_neighbours <- as_tibble(t(as.data.frame(pseq_neighbours))) %>%
  pivot_longer(cols = V2:V20, names_to = "neighbor") %>%
  rename(pseq = V1)

# Assign cell types based on neighbors
labels_df <- as_tibble(pig_atlas$labels, rownames = "cell_id") %>% rename(celltype = value)
pseq_neighbours <- left_join(pseq_neighbours, labels_df, by = c("value" = "cell_id"))

results <- pseq_neighbours %>%
  drop_na() %>%
  group_by(pseq) %>%
  count(celltype) %>%
  mutate(type_fraction = n / sum(n), proportion_neighbours_atlas_cells = sum(n)/100) %>%
  arrange(pseq, desc(n)) %>%
  slice(1) %>%
  mutate(confidence = case_when(
    type_fraction > 0.7 ~ "strong",
    type_fraction >= 0.4 ~ "mid",
    TRUE ~ "low"
  ))

# Save predictions
predictions_meta <- results %>%
  select(cell_id = pseq, labels.p = celltype, type_fraction, confidence, proportion_neighbours_atlas_cells) %>%
  as.data.frame()

write_csv(predictions_meta, here(results_path, 'predictions_meta.csv'))

# Add predictions to integrated data
rownames(predictions_meta) <- predictions_meta$cell_id
predictions_meta <- predictions_meta %>%
  select(labels.p, type_fraction, confidence, proportion_neighbours_atlas_cells)

pig_meta <- pig_atlas@meta.data %>%
  select(labels.p = labels)
pig_meta$type_fraction <- 0
pig_meta$confidence <- "strong"
pig_meta$proportion_neighbours_atlas_cells <- 0

new_meta <- rbind(pig_meta, predictions_meta)
integrated_data <- AddMetaData(integrated_data, metadata = new_meta)

# Save integrated data
saveRDS(integrated_data, file = here(results_path, "drg_integrated.RDS"))

# ======================== Visualization ========================

# Generate UMAP plots
umap_plots <- list(
  DimPlot(integrated_data, group.by = c("labels.p"), cols = reference_celltype_colors),
  DimPlot(integrated_data, group.by = c("labels"), cols = reference_celltype_colors, label = TRUE),
  DimPlot(integrated_data, group.by = c("dataset"), cols = dataset_colors),
  DimPlot(integrated_data, group.by = c("batch"))
)

annotation_umaps <- wrap_plots(umap_plots)
ggsave(here(figs_path, "02-post_integration_umaps.png"), annotation_umaps, width = 20, height = 16)

# Generate supplementary UMAP plots
forSuppLabels <- integrated_data@meta.data %>%
  dplyr::select(dataset, labels.p, Animal, pseq_batch, age, weight) %>% 
  mutate(Animal = as.character(Animal), batch=as.character(pseq_batch), age=as.character(age), weight=as.character(weight)) %>%
  mutate(supp_label = if_else(dataset == "atlas", "atlas", labels.p)) %>% 
  mutate(supp_animal = if_else(dataset == "atlas", "atlas", coalesce(Animal, "NA"))) %>% 
  mutate(supp_batch = if_else(dataset == "atlas", "atlas", pseq_batch)) %>% 
  mutate(supp_age = if_else(dataset == "atlas", "atlas", age)) %>% #NA_integer_
  mutate(supp_weight = if_else(dataset == "atlas", "atlas", weight)) #NA_integer_

integrated_data <- AddMetaData(integrated_data, metadata = forSuppLabels)
supplabels_reference_celltype_colors <- c(
  "C-TAC1-KCNQ5" = "#F8766D", "C-TAC1-LRP1B" = "#E88526", "C-LTMR" = "#D39200", "Aβ-HTMR" = "#B79F00",
  "C-OSMR-SST" = "#93AA00", "Aδ-COOL" = "#5EB300", "C-TAC1-TRPA1" = "#00BA38", "Aβ-LTMR-PALM" = "#00BF74",
  "Aδ-TRPV1" = "#00C19F", "C-COLD" = "#00BFC4", "Aβ-PROPRIO" = "#00B9E3", "Aβ-LTMR-FGF14" = "#00ADFA",
  "C-OSMR-GFRA1_2" = "#619CFF", "Aδ-LTMR" = "#DB72FB", "Aβ-LTMR-ALDH1A1" = "#F564E3",
  "Aδ-HTMR" = "#FF61C3", "patchseq" = "#c00000", "atlas" = "#E0E0E0"
)

# Function to generate color vectors
generate_colors <- function(unique_values) {
  # Remove NA and "atlas" from unique values for color generation
  unique_non_special <- unique_values[!is.na(unique_values) & unique_values != "atlas"]
  
  # Number of unique non-special levels
  n <- length(unique_non_special)
  
  # Generate color vector for non-"atlas"/non-NA levels
  set.seed(123)  # for reproducibility, you can remove this line
  colors <- rainbow(n)  # you can replace `rainbow` with any other color palette function
  
  # Create named color vector
  names(colors) <- unique_non_special
  
  # Add colors for "atlas" and NA
  colors <- c("atlas" = "#E0E0E0", "NA" = "#E0E0E0", colors)
  
  return(colors)
}

# Use the function for all your new columns
supp_animal_colors <- generate_colors(unique(integrated_data@meta.data$supp_animal))
supp_batch_colors <- generate_colors(unique(integrated_data@meta.data$supp_batch))
supp_age_colors <- generate_colors(unique(integrated_data@meta.data$supp_age))
supp_weight_colors <- generate_colors(unique(integrated_data@meta.data$supp_weight))

supp_plots <- list(
  DimPlot(integrated_data, group.by = c("supp_label"), cols = supplabels_reference_celltype_colors),
  DimPlot(integrated_data, group.by = c("supp_animal"), cols = generate_colors(unique(integrated_data$supp_animal))),
  DimPlot(integrated_data, group.by = c("supp_batch"), cols = generate_colors(unique(integrated_data$supp_batch))),
  DimPlot(integrated_data, group.by = c("supp_age"), cols = generate_colors(unique(integrated_data$supp_age))),
  DimPlot(integrated_data, group.by = c("supp_weight"), cols = generate_colors(unique(integrated_data$supp_weight)))
)

supp_umaps <- wrap_plots(supp_plots)
ggsave(here(figs_path, "03-Supp-QC-umaps.png"), supp_umaps, width = 25, height = 15)

# Generate QC violin plots
qc_violin_plots <- VlnPlot(integrated_data, 
                           features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                           group.by = 'labels.p', 
                           ncol = 3, 
                           pt.size = 1, 
                           assay = 'SCT')

ggsave(here(figs_path, "04-QC_integrated_celltypes.png"), qc_violin_plots, width = 18, height = 6, dpi = 300)

# Generate snRNA-seq batch UMAP
snRNAseq_batches <- c("1", "2", "3", "4", "5")
cells_highlight <- lapply(snRNAseq_batches, function(batch) {
  rownames(integrated_data@meta.data[integrated_data@meta.data$orig.ident == batch, ])
})
names(cells_highlight) <- snRNAseq_batches

snRNAseq_QC_UMAP <- DimPlot(
  integrated_data,
  reduction = "umap",
  group.by = "orig.ident",
  cells.highlight = cells_highlight,
  cols.highlight = RColorBrewer::brewer.pal(n = 5, name = 'Set1')
) + ggtitle('snRNAseq batches')

ggsave(here(figs_path, "05-Supp-QC-snRNAseq_batch.png"), snRNAseq_QC_UMAP, width = 8, height = 7)

# Combine all UMAP plots
all_umaps <- wrap_plots(c(umap_plots, supp_plots, list(snRNAseq_QC_UMAP)), ncol = 3) + 
  plot_annotation(tag_levels = 'A')
ggsave(here(figs_path, "06-All-QC-umaps.svg"), all_umaps, width = 30, height = 30)
ggsave(here(figs_path, "06-All-QC-umaps.png"), all_umaps, width = 30, height = 30)

