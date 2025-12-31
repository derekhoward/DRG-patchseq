library(here)
library(readr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
source(here("utils.R"))

results_path <- here("results/figures/fig2")
dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------

load_atlas <- function(atlas_path) {
  atlas <- readRDS(atlas_path)
  VariableFeatures(atlas@assays$RNA) <- VariableFeatures(atlas@assays$integrated)
  DefaultAssay(atlas) <- "RNA"
  Idents(atlas) <- 'Atlas_annotation'
  atlas <- ScaleData(atlas)
  atlas <- RunPCA(atlas, reduction.name = "pca.rna")
  return(atlas)
}

#' Convert gene symbols from human (HGNC) to mouse (MGI) using ortholog mapping
convert_gene_symbols <- function(counts, orthologs) {
  gene_map <- setNames(orthologs$mgi_symbol, orthologs$hgnc_symbol)
  convert_genes <- function(genes) {
    ifelse(genes %in% names(gene_map), gene_map[genes], genes)
  }
  new_rownames <- convert_genes(rownames(counts))
  counts_df <- as.data.frame(as.matrix(counts))
  counts_df$gene <- new_rownames
  counts_df <- aggregate(. ~ gene, data = counts_df, sum)
  rownames(counts_df) <- counts_df$gene
  counts_df$gene <- NULL
  return(as.matrix(counts_df))
}

preprocess_seurat <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj)
  return(seurat_obj)
}

label_transfer <- function(reference, query) {
  transfer.anchors <- FindTransferAnchors(
    reference = reference, 
    query = query, 
    reference.assay = "RNA", 
    query.assay = "RNA", 
    reduction = "cca"
  )
  predictions <- TransferData(
    anchorset = transfer.anchors, 
    refdata = reference$Atlas_annotation, 
    dims = 1:30, 
    weight.reduction = "cca"
  )
  query <- AddMetaData(query, metadata = predictions)
  reference <- RunUMAP(reference, dims = 1:30, reduction = "pca", return.model = TRUE)
  query <- MapQuery(
    anchorset = transfer.anchors, 
    reference = reference, 
    query = query, 
    refdata = list(celltype = "Atlas_annotation"), 
    reference.reduction = "pca", 
    reduction.model = "umap"
  )
  return(query)
}

#' Prepare dataset for label transfer (convert genes and preprocess)
prep_dataset <- function(seurat_obj, orthologs) {
  counts <- GetAssayData(seurat_obj, slot = "counts")
  counts_converted <- convert_gene_symbols(counts, orthologs)
  seurat_obj_converted <- CreateSeuratObject(counts = counts_converted)
  seurat_obj_converted@meta.data <- seurat_obj@meta.data
  seurat_obj_converted <- preprocess_seurat(seurat_obj_converted)
  return(seurat_obj_converted)
}

# ------------------------------------------------------------------------------
# Plotting Functions
# ------------------------------------------------------------------------------

#' Plot reference and query UMAPs side by side
plot_dimplot_results <- function(reference, query, img_path) {
  plot1 <- DimPlot(
    reference, 
    reduction = "umap", 
    group.by = "Atlas_annotation", 
    label = TRUE, 
    label.size = 3, 
    repel = TRUE
  ) + NoLegend() + ggtitle("Reference annotations")
  
  plot2 <- DimPlot(
    query, 
    reduction = "ref.umap", 
    group.by = "predicted.celltype", 
    label = TRUE, 
    label.size = 3, 
    repel = TRUE
  ) + NoLegend() + ggtitle("Query transferred labels")
  
  compound_plot <- plot1 + plot2
  ggsave(img_path, compound_plot, height = 9, width = 12)
  return(compound_plot)
}

#' Plot confusion matrix comparing original vs predicted cell types
plot_confusion_matrix <- function(query, img_path) {
  celltype_order <- c(
    "Pvalb", "Ntrk3high+Ntrk2", "Ntrk3high+S100a16", "Ntrk3low+Ntrk2", 
    "Calca+Bmpr1b", "Calca+Dcn", "Rxfp1", "Trpm8", "Sst", 
    "Mrgpra3+Mrgprb4", "Mrgprd", "Th", "Calca+Smr2", "Calca+Sstr2", 
    "Calca+Oprk1", "Calca+Adra2a"
  )
  
  confusion_matrix <- table(query$labels, query$predicted.celltype)
  confusion_matrix_prop <- prop.table(confusion_matrix, margin = 1)
  confusion_matrix_df <- as.data.frame(as.table(confusion_matrix_prop))
  names(confusion_matrix_df) <- c("Original", "Predicted", "Proportion")
  
  # Ensure all cell types are present
  all_celltypes <- union(celltype_order, unique(confusion_matrix_df$Predicted))
  celltype_order <- union(celltype_order, setdiff(all_celltypes, celltype_order))
  confusion_matrix_df$Predicted <- factor(confusion_matrix_df$Predicted, levels = celltype_order)
  
  conf_mat <- ggplot(confusion_matrix_df, aes(x = Predicted, y = Original, fill = Proportion)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", Proportion)), color = "black", size = 3) +
    scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
    labs(x = "Predicted cell types", y = "Original cell types") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(limits = celltype_order)
  
  ggsave(img_path, conf_mat, height = 12, width = 14)
  return(conf_mat)
}

#' Create custom UMAP with species-colored points and pig overlay
create_species_umap <- function(plotData, custom_colors) {
  ggplot() +
    geom_point(
      data = plotData %>% filter(Species_group == 'Mouse'), 
      aes(x = UMAP_1, y = UMAP_2, color = Species_group), 
      alpha = 0.4, size = 1.2
    ) +
    geom_point(
      data = plotData %>% filter(Species_group == 'Other'), 
      aes(x = UMAP_1, y = UMAP_2, color = Species_group), 
      alpha = 0.4, size = 1.2
    ) +
    geom_point(
      data = plotData %>% filter(Species_group == 'Primate'), 
      aes(x = UMAP_1, y = UMAP_2, color = Species_group), 
      alpha = 0.4, size = 1.2
    ) +
    geom_point(
      data = plotData %>% filter(Species_group == 'Human'), 
      aes(x = UMAP_1, y = UMAP_2, color = Species_group), 
      alpha = 0.4, size = 1.2
    ) +
    geom_point(
      data = plotData %>% filter(Species_group == 'Pig'), 
      aes(x = UMAP_1, y = UMAP_2, fill = Species_group),
      color = "black", alpha = 0.6, size = 1.6, shape = 21, stroke = 0.3
    ) +
    scale_color_manual(values = custom_colors, name = NULL) +
    scale_fill_manual(values = custom_colors, name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme_cowplot() +
    labs(x = 'UMAP 1', y = 'UMAP 2') +
    theme(legend.position = "right") +
    ggtitle("Projected pig snRNAseq")
}

# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------

atlas <- load_atlas(here('data/DRG_neurons_complete.Rds'))
pig_snRNA <- readRDS(here("data/pig_snRNAseq.RDS"))
pig_snRNA_converted <- readRDS(here("data/pig_snRNAseq_converted_human_genes_anchored_to_humans_atlas_correct.Rds"))
patchseq <- prep_patchseq_obj_all_batches(convert_genes = "pig")
human_mouse_orthologs <- read_csv(here('data/processed/mouse_to_human_orthologs.csv'))

# Add cell type labels to patchseq from previous integration
patchseq_metadata <- read.csv(
  here('results/SCT-integration/predictions_meta.csv'), 
  row.names = 'cell_id'
) %>% 
  select('labels.p') %>% 
  rename(labels = labels.p)
patchseq <- AddMetaData(patchseq, metadata = patchseq_metadata)

# ------------------------------------------------------------------------------
# Process pig snRNA-seq Data
# ------------------------------------------------------------------------------

pig_snRNA_prepped <- prep_dataset(pig_snRNA, human_mouse_orthologs)
pig_snRNA_labeled <- label_transfer(atlas, pig_snRNA_prepped)

pig_snRNA_labeled@meta.data %>% 
  dplyr::select(labels, predicted.celltype, predicted.celltype.score) %>% 
  write.csv(here(results_path, 'snRNAseq_predicted_celltypes.csv'))

# ------------------------------------------------------------------------------
# Prep Data for Visualization
# ------------------------------------------------------------------------------

# Add species grouping to atlas
atlas$Species_group <- case_when(
  atlas$Species == "Mouse" ~ "Mouse",
  atlas$Species == "Human" ~ "Human",
  atlas$Species %in% c("C. Macaque", "R. Macaque") ~ "Primate",
  TRUE ~ "Other"
)
atlas$Species_group <- factor(
  atlas$Species_group, 
  levels = c("Mouse", "Human", "Primate", "Other")
)

# Species color palette
custom_colors <- c(
  "Mouse" = "grey70",
  "Human" = "#E69F00",
  "Primate" = "#56B4E9",
  "Pig" = "#FF6B6B",
  "Other" = "#009E73"
)

# Prepare pig UMAP coordinates from pre-processed data (pig_snRNA_converted)
pig_umap_coords <- Embeddings(pig_snRNA_converted, reduction = "ref.umap") %>%
  as_tibble(rownames = "cell_id")

pig_meta <- pig_snRNA_converted@meta.data %>% 
  as_tibble(rownames = 'cell_id') %>% 
  dplyr::select(cell_id, Atlas_annotation = predicted.celltype)

pig_plotData <- left_join(pig_umap_coords, pig_meta, by = "cell_id") %>% 
  dplyr::select(cell_id, UMAP_1 = refUMAP_1, UMAP_2 = refUMAP_2, Atlas_annotation) %>% 
  dplyr::mutate(Species = 'Pig', Species_group = 'Pig', dataset = 'snRNA')

# Prepare atlas UMAP coordinates
atlas_umap_coords <- Embeddings(atlas, reduction = "umap") %>%
  as_tibble(rownames = "cell_id")

atlas_meta <- atlas@meta.data %>% 
  as_tibble(rownames = 'cell_id') %>% 
  dplyr::select(cell_id, Atlas_annotation, Species, Species_group) %>% 
  dplyr::mutate(dataset = 'atlas')

atlas_plotData <- left_join(atlas_umap_coords, atlas_meta, by = "cell_id")

# Combine for plotting
plotData <- bind_rows(pig_plotData, atlas_plotData)

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

# Reference atlas UMAP
p1 <- DimPlot(
  atlas, 
  reduction = "umap", 
  group.by = "Atlas_annotation", 
  label = FALSE
) + 
  theme_cowplot() + 
  ggtitle("Reference atlas celltypes")

# Projected pig snRNA-seq UMAP
custom_UMAP <- create_species_umap(plotData, custom_colors)

# Save individual and combined plots
ggsave(
  here(results_path, 'umap-half_species_projection.png'), 
  custom_UMAP, height = 7, width = 7, dpi = 300
)
ggsave(
  here(results_path, 'umap-half_species_projection.svg'), 
  custom_UMAP, height = 7, width = 7
)

compound_plot <- p1 + custom_UMAP
ggsave(
  here(results_path, 'umap-reference_snRNA_projected.png'), 
  compound_plot, height = 7, width = 15
)
ggsave(
  here(results_path, 'umap-reference_snRNA_projected.svg'), 
  compound_plot, height = 7, width = 15
)

# Confusion matrix
snrna_confmat <- plot_confusion_matrix(
  pig_snRNA_labeled, 
  here(results_path, 'conf_mat-snRNA_mapquery_predictions.png')
)

# Combined figure with confusion matrix
final <- compound_plot / snrna_confmat + plot_annotation(tag_levels = 'A')
ggsave(
  here(results_path, 'fig2_w_confmat.png'), 
  final, height = 12, width = 12
)

# ------------------------------------------------------------------------------
# Process Patch-seq Data
# ------------------------------------------------------------------------------
patchseq_prepped <- prep_dataset(patchseq, human_mouse_orthologs)
patchseq_labeled <- label_transfer(atlas, patchseq_prepped)

# Save predictions
patchseq_labeled@meta.data %>% 
  dplyr::select(labels, predicted.celltype, predicted.celltype.score) %>% 
  write.csv(here(results_path, 'patchseq_predicted_celltypes.csv'))

# Generate plots
plot_dimplot_results(
  atlas, 
  patchseq_labeled, 
  here(results_path, 'umap-patchseq_mapquery_predictions.png')
)
