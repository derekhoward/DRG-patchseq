library(here)
library(Seurat)
library(readr)
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
source(here("utils.R"))

options(matrixStats.useNames.NA = 'deprecated')

results_path <- here("results/pig_spatial_integration")
dir.create(results_path)

# full integration of pig snRNAseq with pig visium  data
pig_ref <- readRDS(here("data/pig_snRNAseq.RDS"))

spatial_drg <- readRDS(here('./data/price_spatial_transcriptomics/Pig_visium_rds_file/pig_clustering3_res1_REGRESS_seq123_04262023_test05302023.rds'))
spatial_drg[["percent.mt"]] <- PercentageFeatureSet(spatial_drg, features = c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB"))

DefaultAssay(pig_ref)
n_distinct(rownames(pig_ref))
DefaultAssay(spatial_drg)
n_distinct(rownames(spatial_drg))

genes_intersection <- intersect(rownames(pig_ref), rownames(spatial_drg))
n_distinct(genes_intersection)

# subset_genes_pig_ref <- pig_ref@assays$RNA[genes_intersection, ]
subset_genes_pig_ref <- GetAssayData(object = pig_ref, assay = "RNA", slot = "counts")[genes_intersection, ]

# Warning message:
# The `slot` argument of `GetAssayData()` is deprecated as of        
# SeuratObject 5.0.0.      
# subset_genes_spatial <- spatial_drg@assays$RNA[genes_intersection, ]
subset_genes_spatial <- GetAssayData(object = spatial_drg, assay = "RNA", slot = "counts")[genes_intersection, ]


pig_atlas <- CreateSeuratObject(subset_genes_pig_ref, meta.data = pig_ref@meta.data)
spatial <- CreateSeuratObject(subset_genes_spatial, meta.data = spatial_drg@meta.data)

pig_atlas$batch <- pig_atlas$orig.ident
pig_atlas$dataset <- "snRNAseq"

spatial$batch <- spatial$id
spatial$dataset <- "Visium"

merged_object <- merge(pig_atlas, y = spatial)

DefaultAssay(merged_object)
data.list <- SplitObject(merged_object, split.by = "batch")
# drop batch #2 (from snRNAseq)
data.list <- data.list[c(1, 2, 3, 4, 5, "pig1", "pig2", "pig3", "pig4", "pig5", "pig6", "pig7")]

# SCT normalization of batches
for (i in names(data.list)) {
  data.list[[i]] <- SCTransform(data.list[[i]], vst.flavor = "v2", verbose = TRUE, vars.to.regress="percent.mt", return.only.var.genes = FALSE)
}

# integration
drg.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 6000)
drg.list <- PrepSCTIntegration(object.list = data.list, anchor.features = drg.features)
drg.anchors <- FindIntegrationAnchors(object.list = drg.list, normalization.method = "SCT", anchor.features = drg.features)
# the following only works in k.weight is decreased to min number of cells in smallest sample
drg.integrated <- IntegrateData(anchorset = drg.anchors, normalization.method = "SCT", k.weight = 75)

# read in label preds for spatialdrg
pigspatial_preds <- read.csv(here('./results/pig_spatial_integration/iterative_spatialpig_preds.csv'), row.names = 'cell_id')
pigspatial_preds <- pigspatial_preds %>% dplyr::select(labels.p)
drg.integrated <- AddMetaData(object = drg.integrated, metada = pigspatial_preds)
drg.integrated$final_labels <- ifelse(drg.integrated$dataset == "snRNAseq",
                                      as.character(drg.integrated$labels),
                                      as.character(drg.integrated$labels.p))
table(drg.integrated$dataset, drg.integrated$final_labels)


drg.integrated <- RunPCA(object = drg.integrated, verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)

drg.integrated$final_labels <- factor(drg.integrated$final_labels, levels = reference_celltype_levels)

saveRDS(drg.integrated, file = here(results_path, "drg_snRNAseq_visium_reintegrated.RDS"))

drg.integrated <- readRDS(here(results_path, "drg_snRNAseq_visium_reintegrated.RDS"))
spatial_subset <- subset(drg.integrated, dataset == 'Visium')

DefaultAssay(spatial_subset) <- "RNA"

# QC Plots
## violin plots of batches, showing n_genes, n_reads, percent.mt
## done by batch, then also by celltype assignment

qc_batch_spatial <- VlnPlot(spatial_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by ='batch', ncol = 3, pt.size = 1, assay = 'RNA') & theme(axis.title.x = element_blank())
qc_celltype_spatial <- VlnPlot(spatial_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by ='final_labels', ncol = 3, pt.size = 1, assay = 'RNA', cols = reference_celltype_colors) & theme(axis.title.x = element_blank())


qc_spatial_nFeature <- VlnPlot(spatial_subset, features = "nFeature_RNA", group.by ='final_labels', pt.size = 1, assay = 'SCT', cols = reference_celltype_colors) + NoLegend()+ 
  ylab("Number of Features") & theme(axis.title.x = element_blank())
qc_spatial_nCount <- VlnPlot(spatial_subset, features = "nCount_RNA", group.by ='final_labels', pt.size = 1, assay = 'SCT', log = T, cols = reference_celltype_colors) + NoLegend()+ 
  ylab("Total RNA Counts (log scale)") & theme(axis.title.x = element_blank())
qc_spatial_percent_mt <- VlnPlot(spatial_subset, features = "percent.mt", group.by ='final_labels', pt.size = 1, assay = 'SCT', cols = reference_celltype_colors) + NoLegend() + 
  ylab("Percent Mitochondrial Content") & theme(axis.title.x = element_blank())
QC_pig_spatial <- (qc_spatial_nFeature | qc_spatial_nCount | qc_spatial_percent_mt) + plot_annotation(tag_levels = 'A')
ggsave(here("results/pig_spatial_integration/", "02-QC_pig_spatial.svg"), QC_pig_spatial, width = 14, height = 6, dpi = 300)
ggsave(here("results/pig_spatial_integration/", "02-QC_pig_spatial.png"), QC_pig_spatial, width = 14, height = 6, dpi = 300)


dataset_colors <- c("snRNAseq" = "firebrick", "Visium" = "steelblue")


p2 <- DimPlot(drg.integrated, group.by = c("dataset"), cols = dataset_colors, pt.size = 1)+ ggtitle('Integrated datasets')

# create a dimplot highlighting the different pig/batches
UMAP_animals <- DimPlot(
  drg.integrated,
  reduction = "umap",
  group.by = "id",
  pt.size = 1,
  cells.highlight = list(
    'pig1' = drg.integrated@meta.data %>%
      filter(id == "pig1") %>%
      row.names(.),
    'pig2' = drg.integrated@meta.data %>%
      filter(id == "pig2") %>%
      row.names(.),
    'pig3' = drg.integrated@meta.data %>%
      filter(id == "pig3") %>%
      row.names(.),
    'pig4' = drg.integrated@meta.data %>%
      filter(id == "pig4") %>%
      row.names(.),
    'pig5' = drg.integrated@meta.data %>%
      filter(id == "pig5") %>%
      row.names(.),
    'pig6' = drg.integrated@meta.data %>%
      filter(id == "pig6") %>%
      row.names(.),
    'pig7' = drg.integrated@meta.data %>%
      filter(id == "pig7") %>%
      row.names(.)
  ),
  cols.highlight = RColorBrewer::brewer.pal(n = 7, name = 'Set1'),sizes.highlight = 1
) + ggtitle('Visium data: animals')

forSuppLabels <- drg.integrated@meta.data %>%
  dplyr::select(dataset, final_labels) %>% 
  mutate(supp_label = if_else(dataset == "snRNAseq", "snRNAseq", final_labels))

drg.integrated <- AddMetaData(drg.integrated, metadata = forSuppLabels)

cells.highlight <- setNames(lapply(reference_celltype_levels, function(celltype) {
  drg.integrated@meta.data %>%
    dplyr::filter(supp_label == celltype) %>%
    row.names()
}), reference_celltype_levels)

supplabels_reference_celltype_colors <- c(
  "C-TAC1-KCNQ5" = "#F8766D", "C-TAC1-LRP1B" = "#E88526", "C-LTMR" = "#D39200", "Aβ-HTMR" = "#B79F00",
  "C-OSMR-IL31RA" = "#93AA00", "Aδ-COOL" = "#5EB300", "C-TAC1-TRPA1" = "#00BA38", "Aβ-LTMR-PALM" = "#00BF74",
  "Aδ-TRPV1" = "#00C19F", "C-COLD" = "#00BFC4", "Aβ-PROPRIO" = "#00B9E3", "Aβ-LTMR-FGF14" = "#00ADFA",
  "C-OSMR-GFRA1_2" = "#619CFF", "Aδ-LTMR" = "#DB72FB", "Aβ-LTMR-ALDH1A1" = "#F564E3",
  "Aδ-HTMR" = "#FF61C3", "patchseq" = "#c00000", "snRNAseq" = "#E0E0E0"
)


UMAP_celltypes <- DimPlot(
  drg.integrated,
  reduction = "umap",
  group.by = "supp_label", 
  cells.highlight = cells.highlight, pt.size = 1,
  cols.highlight = supplabels_reference_celltype_colors,
  sizes.highlight = 1) + ggtitle('Visium data: mapped cells') #RColorBrewer::brewer.pal(n = 17, name = 'Set1')

qc_violins <- qc_batch_spatial + qc_celltype_spatial
QC_umaps <- p2+UMAP_celltypes+UMAP_animals #p1+p2+p3+UMAP_animals

qc_violins / QC_umaps

visium_QC_plots <- (qc_batch_spatial / qc_celltype_spatial / QC_umaps) + plot_layout(heights = c(1, 1, 2))+ plot_annotation(tag_levels = 'A')
ggsave(here("results/pig_spatial_integration/", "02-Supp-QC-visium.png"), visium_QC_plots, width = 15, height = 20)
ggsave(here("results/pig_spatial_integration/", "02-Supp-QC-visium.svg"), visium_QC_plots, width = 15, height = 20)
