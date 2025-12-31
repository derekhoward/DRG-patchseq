library(here)
library(Seurat)
library(edgeR)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)
source(here("utils.R"))

theme_set(theme_cowplot(font_size = 26))

results_path <- here("results/SCT-integration")
figs_path <- here("results/figures", "fig1")
dir.create(figs_path, recursive = TRUE)

##########################################################################
# Load data
##########################################################################
drg.integrated <- readRDS(here(results_path, "drg_integrated.RDS"))
drg.integrated$labels <- factor(drg.integrated$labels, levels = reference_celltype_levels)

pig_ref <- readRDS(here("data/pig_snRNAseq.RDS"))
pig_ref$labels <- factor(pig_ref$labels, levels = reference_celltype_levels)

spatial_drg <- readRDS(here('./data/price_spatial_transcriptomics/Pig_visium_rds_file/pig_clustering3_res1_REGRESS_seq123_04262023_test05302023.rds'))
spatial_drg[["percent.mt"]] <- PercentageFeatureSet(spatial_drg, features = c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB"))
pigspatial_preds <- read.csv(here('results/pig_spatial_integration/iterative_spatialpig_preds.csv'))
rownames(pigspatial_preds) <- pigspatial_preds$cell_id
pigspatial_preds <- pigspatial_preds %>% dplyr::select(labels.p)
spatial_drg <- AddMetaData(object = spatial_drg, metadata = pigspatial_preds)
spatial_drg$labels.p <- factor(spatial_drg$labels.p, levels = reference_celltype_levels)

##########################################################################
# Fig 1A: UMAP with cluster labels
##########################################################################
integrated_umap_coords <- Embeddings(object = drg.integrated, reduction = "umap")
integrated_umap_coords <- as_tibble(integrated_umap_coords, rownames = "cell_id")

snrnaseq_ids <- rownames(subset(drg.integrated@meta.data, subset = dataset == "atlas"))

reference_coords <- integrated_umap_coords %>% filter(cell_id %in% snrnaseq_ids)
patchseq_coords <- integrated_umap_coords %>% filter(!cell_id %in% snrnaseq_ids)

reference_coords$dataset <- "reference"
patchseq_coords$dataset <- "patchseq"

reference_plotdata <- as_tibble(subset(drg.integrated@meta.data, subset = dataset == "atlas"), rownames = "cell_id") %>%
  dplyr::select(cell_id, celltype = labels)
reference_plotdata <- left_join(reference_coords, reference_plotdata)

patchseq_plotdata <- as_tibble(subset(drg.integrated@meta.data, subset = dataset == "patch"), rownames = "cell_id") %>%
  dplyr::select(cell_id, celltype = "labels.p")
patchseq_plotdata <- left_join(patchseq_coords, patchseq_plotdata)

plotData <- bind_rows(reference_plotdata, patchseq_plotdata)
plotData <- plotData %>% rename(`Cell type` = celltype)

cluster_centers <- plotData %>%
  group_by(`Cell type`) %>%
  summarise(center_UMAP_1 = mean(UMAP_1), center_UMAP_2 = mean(UMAP_2))

label_colors <- reference_celltype_colors[match(cluster_centers$`Cell type`, names(reference_celltype_colors))]

custom_UMAP2 <- ggplot() +
geom_point(data = plotData %>% filter(dataset == "reference"),
             aes(x = UMAP_1, y = UMAP_2, color = `Cell type`), alpha = 0.4, size = 1.5) +
  geom_point(data = plotData %>% filter(dataset == "patchseq"),
             aes(x = UMAP_1, y = UMAP_2, fill = `Cell type`), shape = 21, stroke = 0.5, size = 3, alpha = 1) +
  scale_fill_manual(values = reference_celltype_colors) +
  scale_color_manual(values = reference_celltype_colors) +
  theme(legend.position = "none") +
  geom_label_repel(
    data = cluster_centers,
    aes(x = center_UMAP_1, y = center_UMAP_2, label = `Cell type`),
    min.segment.length = Inf,
    fill = label_colors,
    colour = "black",
    size = 6,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = "grey50",
    inherit.aes = FALSE,
    seed = 42
  )

ggsave(here(figs_path, "1A-annotated_UMAP.png"), custom_UMAP2, height = 10, width = 10)

##########################################################################
# Prepare ephys metadata for subsequent plots
##########################################################################
obj.list <- SplitObject(drg.integrated, split.by = "dataset")
int_patchseq <- obj.list[[2]]

meta_table <- FetchData(object = int_patchseq, vars = c(
  "Animal", "age", "weight", "labels.p", "confidence",
  "days.in.dish", "diameter", "cslow", "RMP", "pA_start",
  "pA_threshold", "Rheobase", "input_res", "Input_res_new",
  "AP_threshhold", "AP_max", "Apamp", "Min_AHP", "Loc_APHW",
  "APD", "max_slope_raw", "max_slope_smth", "slope_oversh",
  "slope_oversh_0.5", "subthresh_slope", "subthresh_slope_0.5",
  "TTP", "Shoulder", "AHP_new", "AHP_Tau", "RMP_beforeAP",
  "X2nd_Infliction", "Spontan_active", "X2Hz", "X4Hz", "X10Hz",
  "X20Hz", "X50Hz", "X100Hz", "HS_Threshold", "AP_HS_2TH",
  "Sinusscore_old", "Sinusscore_new", "HS_Quotient", "batch",
  "Delta.I.Initial", "Delta.ind.time", "Delta.Rmem"
))

feats_to_select <- meta_table %>% select_if(is.numeric) %>% names()

tidy_meta_table <- meta_table %>%
  pivot_longer(cols = all_of(feats_to_select), names_to = "efeature", values_to = "measure")

tidy_meta_table$labels.p <- factor(tidy_meta_table$labels.p, levels = reference_celltype_levels)

##########################################################################
# Fig 1B: Double barplot (snRNAseq + Patchseq counts)
##########################################################################
theme_set(theme_cowplot(font_size = 28))

patchseq_subset <- subset(drg.integrated, subset = dataset == 'patch')
nucseq_subset <- subset(drg.integrated, subset = dataset == 'atlas')

patchseq_count <- as.data.frame(table(patchseq_subset$labels.p))
names(patchseq_count) <- c("labels", "Patchseq")

nucseq_count <- as.data.frame(table(nucseq_subset$labels.p))
names(nucseq_count) <- c("labels", "snRNAseq")

combined_count <- left_join(patchseq_count, nucseq_count)

long_combined_count <- combined_count %>%
  pivot_longer(cols = c(Patchseq, snRNAseq), names_to = "dataset", values_to = "count")

long_combined_count$labels <- factor(long_combined_count$labels, levels = rev(reference_celltype_levels))

patchseq_data <- subset(long_combined_count, dataset == "Patchseq")
snRNAseq_data <- subset(long_combined_count, dataset == "snRNAseq")

patchseq_plot <- ggplot(patchseq_data, aes(x = labels, y = count, fill = labels)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = count), vjust = 0.5, hjust = 0, position = position_dodge(width = 0)) +
  scale_fill_manual(values = reference_celltype_colors) +
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  NoLegend() +
  ggtitle('Patch-seq')

snRNAseq_plot <- ggplot(snRNAseq_data, aes(x = labels, y = count, fill = labels)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = count), vjust = 0.5, hjust = 1, position = position_dodge(width = 0)) +
  scale_fill_manual(values = reference_celltype_colors) +
  scale_y_reverse() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.length.x = unit(0.2, "cm"),
        plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  ggtitle('snRNAseq')

counts_double_barplot <- snRNAseq_plot + patchseq_plot + plot_layout(ncol = 2)

ggsave(here(figs_path, "1B-snrna_patchseq_counts_barcharts.png"), counts_double_barplot, height = 8, width = 12)

##########################################################################
# Fig 1C: Patchseq cell counts (single barplot)
##########################################################################
labels_count <- as.data.frame(table(patchseq_subset$labels.p))
names(labels_count) <- c("labels.p", "count")
labels_count$labels.p <- factor(labels_count$labels.p, levels = reference_celltype_levels)

patchseq_cellcounts <- ggplot(labels_count, aes(x = labels.p, y = count, fill = labels.p)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = reference_celltype_colors) +
  scale_y_continuous(breaks = seq(0, max(labels_count$count), by = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "Number of cells", title = "") +
  NoLegend()

ggsave(here(figs_path, "1C-patchseq_cell_counts.png"), patchseq_cellcounts, height = 8, width = 8, bg = "white")

##########################################################################
# Fig 1D: AP Duration
##########################################################################
theme_set(theme_cowplot(font_size = 26))

common_theme <- function() {
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_blank()
  )
}

AP_duration <- tidy_meta_table %>%
  filter(efeature == "APD") %>%
  ggplot(aes(x = labels.p, y = measure, color = labels.p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15) +
  scale_colour_manual(values = reference_celltype_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank()) +
  labs(x = "", y = "AP duration (ms)") +
  guides(color = "none") +
  ggtitle("Action potential duration")

ggsave(here(figs_path, "1D-ap_duration.png"), AP_duration, height = 8, width = 8, bg = "white")

##########################################################################
# Fig 1E: AP Duration + SST Expression Comparison
##########################################################################
# Prepare patchseq object with CPM
patchseq <- prep_patchseq_obj_all_batches(convert_genes = "pig")
metadata <- read.csv(here('results/SCT-integration/predictions_meta.csv'), row.names = 'cell_id')
metadata <- metadata %>% select('labels.p') %>% rename(labels = labels.p)
patchseq <- AddMetaData(patchseq, metadata = metadata)
patchseq[["percent.mt"]] <- PercentageFeatureSet(patchseq, features = c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB"))

counts_matrix <- GetAssayData(patchseq, assay = "RNA", layer = "counts")
cpm_matrix <- edgeR::cpm(counts_matrix)
patchseq[["CPM"]] <- CreateAssayObject(data = cpm_matrix)
patchseq$labels <- factor(patchseq$labels, levels = reference_celltype_levels)

# C-fiber subsets for SST comparison
cfibers <- c("C-TAC1-LRP1B", "C-TAC1-TRPA1", "C-TAC1-KCNQ5", "C-COLD", "C-LTMR", "C-OSMR-GFRA1_2", "C-OSMR-SST")
pig_snRNA_subset <- subset(pig_ref, subset = labels %in% cfibers)
patchseq_subset <- subset(patchseq, subset = labels %in% cfibers)
spatial_drg_subset <- subset(spatial_drg, subset = labels.p %in% cfibers)

pig_snRNA_subset$labels <- factor(pig_snRNA_subset$labels, levels = cfibers)
patchseq_subset$labels <- factor(patchseq_subset$labels, levels = cfibers)

# Gene expression comparison function
plot_gene_expression_comparison <- function(
    gene_symbol,
    snrna_data,
    patchseq_data,
    visium_data,
    celltype_levels,
    celltype_colors,
    output_dir,
    snrna_group_by_col = "labels",
    patchseq_group_by_col = "labels",
    visium_group_by_col = "labels.p",
    snrna_assay = "RNA",
    visium_assay = "RNA",
    snrna_slot = "data",
    patchseq_slot_cpm = "data",
    visium_slot = "data",
    log_transform_snrna = TRUE,
    log_transform_patchseq = TRUE,
    log_transform_visium = TRUE,
    pt_size_snrna = 0.25,
    pt_size_patchseq = 0.25,
    pt_size_visium = 0.25,
    save_plot = TRUE,
    plot_filename_suffix = "_expr_comparison.png",
    plot_width = 12,
    plot_height = 7,
    xaxis_label_size = 10,
    verbose = TRUE
) {
  
  common_theme_local <- function() {
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = xaxis_label_size),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_blank(),
      axis.title.x = element_blank()
    )
  }
  
  snrna_data_local <- snrna_data
  patchseq_data_local <- patchseq_data
  visium_data_local <- visium_data
  
  if (!("CPM" %in% Assays(patchseq_data_local))) {
    if (verbose) message("CPM assay not found in patchseq_data. Calculating from RNA assay counts.")
    counts_matrix <- GetAssayData(patchseq_data_local, assay = "RNA", layer = "counts")
    cpm_matrix <- edgeR::cpm(counts_matrix)
    patchseq_data_local[["CPM"]] <- CreateAssayObject(data = cpm_matrix)
  }
  
  if (snrna_group_by_col %in% colnames(snrna_data_local@meta.data)) {
    snrna_data_local@meta.data[[snrna_group_by_col]] <- factor(snrna_data_local@meta.data[[snrna_group_by_col]], levels = intersect(celltype_levels, unique(snrna_data_local@meta.data[[snrna_group_by_col]])))
  }
  
  if (patchseq_group_by_col %in% colnames(patchseq_data_local@meta.data)) {
    patchseq_data_local@meta.data[[patchseq_group_by_col]] <- factor(patchseq_data_local@meta.data[[patchseq_group_by_col]], levels = intersect(celltype_levels, unique(patchseq_data_local@meta.data[[patchseq_group_by_col]])))
  }
  
  if (visium_group_by_col %in% colnames(visium_data_local@meta.data)) {
    visium_data_local@meta.data[[visium_group_by_col]] <- factor(visium_data_local@meta.data[[visium_group_by_col]], levels = intersect(celltype_levels, unique(visium_data_local@meta.data[[visium_group_by_col]])))
  }
  
  gene_present_snrna <- gene_symbol %in% rownames(snrna_data_local[[snrna_assay]])
  gene_present_patchseq <- gene_symbol %in% rownames(patchseq_data_local[["CPM"]])
  gene_present_visium <- gene_symbol %in% rownames(visium_data_local[[visium_assay]])
  
  if (gene_present_snrna) {
    ylab_snrna <- bquote(.(gene_symbol) ~ "expr. ("*log[2]*" norm.)")
    plot_snrna <- VlnPlot(snrna_data_local, features = gene_symbol, group.by = snrna_group_by_col,
                          assay = snrna_assay, slot = snrna_slot, cols = celltype_colors,
                          log = log_transform_snrna, pt.size = pt_size_snrna, raster = FALSE) +
      ggtitle("snRNAseq") + NoLegend() + ylab(ylab_snrna) + common_theme_local()
  } else {
    plot_snrna <- ggplot() + annotate("text", x=0.5, y=0.5, label=paste(gene_symbol, "\nnot in snRNAseq")) + theme_void() + ggtitle("snRNAseq") + common_theme_local()
  }
  
  if (gene_present_patchseq) {
    ylab_patchseq <- bquote(.(gene_symbol) ~ "expr. ("*log[2]*" CPM)")
    plot_patchseq <- VlnPlot(patchseq_data_local, features = gene_symbol, group.by = patchseq_group_by_col,
                             assay = "CPM", slot = patchseq_slot_cpm, cols = celltype_colors,
                             log = log_transform_patchseq, pt.size = pt_size_patchseq, raster = FALSE) +
      ggtitle("Patch-seq") + NoLegend() + ylab(ylab_patchseq) + common_theme_local()
  } else {
    plot_patchseq <- ggplot() + annotate("text", x=0.5, y=0.5, label=paste(gene_symbol, "\nnot in Patch-seq CPM")) + theme_void() + ggtitle("Patch-seq") + common_theme_local()
  }
  
  if (gene_present_visium) {
    ylab_visium <- bquote(.(gene_symbol) ~ "expr. ("*log[2]*" norm.)")
    plot_visium <- VlnPlot(visium_data_local, features = gene_symbol, group.by = visium_group_by_col,
                           assay = visium_assay, slot = visium_slot, cols = celltype_colors,
                           log = log_transform_visium, pt.size = pt_size_visium, raster = FALSE) +
      ggtitle("Pig visium") + NoLegend() + ylab(ylab_visium) + common_theme_local()
  } else {
    plot_visium <- ggplot() + annotate("text", x=0.5, y=0.5, label=paste(gene_symbol, "\nnot in Visium")) + theme_void() + ggtitle("Pig visium") + common_theme_local()
  }
  
  combined_plot <- (plot_snrna | plot_patchseq | plot_visium)
  
  if (save_plot) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    plot_filename <- paste0(gsub("[^a-zA-Z0-9_]", "_", gene_symbol), plot_filename_suffix)
    save_path <- file.path(output_dir, plot_filename)
    ggsave(filename = save_path, plot = combined_plot, height = plot_height, width = plot_width, bg = "white")
    if (verbose) message(paste("Plot saved to:", save_path))
  }
  
  return(combined_plot)
}

# SST expression comparison (C-fiber subset)
subset_sst_exp_plot <- plot_gene_expression_comparison(
  gene_symbol = "SST",
  snrna_data = pig_snRNA_subset,
  patchseq_data = patchseq_subset,
  visium_data = spatial_drg_subset,
  celltype_levels = reference_celltype_levels,
  celltype_colors = reference_celltype_colors,
  output_dir = figs_path,
  save_plot = FALSE,
  xaxis_label_size = 12
)

fig1e <- (AP_duration | subset_sst_exp_plot) + plot_layout(widths = c(1, 2))
ggsave(here(figs_path, "1E-AP_duration-SST_expr.png"), fig1e, height = 8, width = 14, bg = "white")

##########################################################################
# Fig 1F: Expression Dotplot
##########################################################################
Idents(pig_ref) <- pig_ref$labels

# Marker gene lists
cfiber_genelist <- c('TRPM8', 'KIT', 'TAC1', 'KCNQ5', 'CALCB', 'ADCYAP', 'LRP1B', 'TRPA1', 'CDH9','GFRA2', 'GFRA1', 'LPAR3', 'SYNPR', 'TMC3', 'IL31RA', 'OSMR', 'HRH1', 'JAK1', 'NPPB')
adelta_genelist <- c('PCDH7', 'KCNQ3', 'NTRK2', 'B4GALT6', 'DCC', 'CREB5', 'CACNA1E', 'UNC5D', 'ZNF521', 'TRPV1', 'RORA', 'FOXP2', 'ANTXR3')
abeta_genelist <- c('PVALB', 'SPP1', 'ETV1', 'ALDH1A1', 'SLIT2', 'PALM', 'CPNE4', 'PLXNA2', 'S100A100')
scn_channels_list <- c('PIEZO2', 'SCN1A', 'SCN2A', 'SCN8A', 'SCN9A', 'SCN10A', 'SCN11A')

allmarkerslist <- c(abeta_genelist, adelta_genelist, cfiber_genelist, scn_channels_list)

pig_ref_for_flipped <- pig_ref
desired_order <- rev(reference_celltype_levels)
pig_ref_for_flipped$labels <- factor(pig_ref_for_flipped$labels, levels = desired_order)
Idents(pig_ref_for_flipped) <- pig_ref_for_flipped$labels

p2a_flipped <- DotPlot(pig_ref_for_flipped, assay = "SCT", features = allmarkerslist) +
  scale_colour_viridis_c(option = "magma", begin = 0, end = 0.9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab('') + ylab('')

ggsave(here(figs_path, "1F-dotplot.png"), p2a_flipped, height = 8, width = 16, bg = "white")

##########################################################################
# Supp: All electrophysiology features
##########################################################################
theme_set(theme_cowplot(font_size = 14))

supp_feats <- c(
  "diameter", "cslow", "RMP", "Rheobase", "Input_res_new", "AP_threshhold",
  "AP_max", "Apamp", "Min_AHP", "APD", "max_slope_smth", "slope_oversh_0.5",
  "subthresh_slope_0.5", "TTP", "AHP_new", "X2Hz", "X4Hz", "X10Hz", "X20Hz",
  "X50Hz", "X100Hz", "HS_Quotient", "Sinusscore_old", "Spontan_active"
)

meta_table_supp <- FetchData(object = int_patchseq, vars = c('labels.p', 'confidence', 'batch', supp_feats))

tidy_meta_table_supp <- meta_table_supp %>%
  pivot_longer(cols = all_of(supp_feats), names_to = "efeature", values_to = "measure")

tidy_meta_table_supp$labels.p <- factor(tidy_meta_table_supp$labels.p, levels = reference_celltype_levels)
tidy_meta_table_supp$batch <- as.factor(tidy_meta_table_supp$batch)
tidy_meta_table_supp$efeature <- as.factor(tidy_meta_table_supp$efeature)

name_mapping <- c(
  "diameter" = "Diameter",
  "cslow" = "Cell capacitance",
  "RMP" = "Resting membrane potential",
  "Input_res_new" = "Input resistance",
  "AP_threshhold" = "Action potential threshold",
  "AP_max" = "AP maximum",
  "Apamp" = "AP amplitude",
  "Min_AHP" = "Afterhyperpolarization minimum",
  "APD" = "Action potential duration",
  "max_slope_smth" = "Maximum slope",
  "slope_oversh_0.5" = "Maximum overshoot slope",
  "subthresh_slope_0.5" = "Maximum subthreshold slope",
  "TTP" = "Time to peak",
  "AHP_new" = "Afterhyperpolarization Tau",
  "X2Hz" = "Follow frequency 2Hz",
  "X4Hz" = "Follow frequency 4Hz",
  "X10Hz" = "Follow frequency 10Hz",
  "X20Hz" = "Follow frequency 20Hz",
  "X50Hz" = "Follow frequency 50Hz",
  "X100Hz" = "Follow frequency 100Hz",
  "HS_Quotient" = "Half sinus threshold",
  "Sinusscore_old" = "Sinus score",
  "Spontan_active" = "Spontaneously active"
)

tidy_meta_table_supp$efeature <- recode(tidy_meta_table_supp$efeature, !!!name_mapping)

efeat_order <- c(
  "Afterhyperpolarization Tau", "AP maximum", "Action potential threshold",
  "AP amplitude", "Action potential duration", "Cell capacitance", "Diameter",
  "Half sinus threshold", "Input resistance", "Maximum slope",
  "Afterhyperpolarization minimum", "Rheobase", "Resting membrane potential",
  "Sinus score", "Maximum overshoot slope", "Spontaneously active",
  "Maximum subthreshold slope", "Time to peak", "Follow frequency 2Hz",
  "Follow frequency 4Hz", "Follow frequency 10Hz", "Follow frequency 20Hz",
  "Follow frequency 50Hz", "Follow frequency 100Hz"
)
tidy_meta_table_supp$efeature <- factor(tidy_meta_table_supp$efeature, levels = efeat_order)

tidy_meta_table_supp <- tidy_meta_table_supp %>%
  mutate(y_labels = case_when(
    efeature == "Diameter" ~ "µm",
    efeature == "Cell capacitance" ~ "pF",
    efeature == "Resting membrane potential" ~ "mV",
    efeature == "Rheobase" ~ "pA",
    efeature == "Input resistance" ~ "GΩ",
    efeature == "Action potential threshold" ~ "mV",
    efeature == "AP maximum" ~ "mV",
    efeature == "AP amplitude" ~ "mV",
    efeature == "Afterhyperpolarization minimum" ~ "mV",
    efeature == "Action potential duration" ~ "ms",
    efeature == "Maximum slope" ~ "V/s",
    efeature == "Maximum overshoot slope" ~ "V/s",
    efeature == "Maximum subthreshold slope" ~ "V/s",
    efeature == "Time to peak" ~ "ms",
    efeature == "Afterhyperpolarization Tau" ~ "ms",
    efeature == "Follow frequency 2Hz" ~ "#APs",
    efeature == "Follow frequency 4Hz" ~ "#APs",
    efeature == "Follow frequency 10Hz" ~ "#APs",
    efeature == "Follow frequency 20Hz" ~ "#APs",
    efeature == "Follow frequency 50Hz" ~ "#APs",
    efeature == "Follow frequency 100Hz" ~ "#APs",
    efeature == "Half sinus threshold" ~ "Half sinus threshold",
    efeature == "Sinus score" ~ "Sinus score",
    efeature == "Spontaneously active" ~ "Spontaneously active",
    TRUE ~ "NA"
  ))

plot_feature <- function(data, feature_name) {
  if (feature_name == 'Spontaneously active') {
    spontaneous_proportions <- data %>%
      filter(efeature == 'Spontaneously active') %>%
      group_by(labels.p) %>%
      summarise(
        Active = sum(measure == 1),
        Not_Active = sum(measure == 0),
        Total = n()
      ) %>%
      mutate(Proportion_Active = Active / Total * 100) %>%
      dplyr::select(labels.p, Proportion_Active) %>%
      pivot_longer(cols = starts_with("Proportion"), names_to = "Activity", values_to = "Proportion")

    plot <- ggplot(spontaneous_proportions, aes(x = labels.p, y = Proportion, fill = labels.p)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = reference_celltype_colors) +
      labs(x = "Cell Type", y = "Proportion (%)", title = "Proportion of Spontaneously Active Cells") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
      ylim(0, 100)

    return(plot)
  } else {
    plot <- ggplot(data %>% filter(efeature == feature_name), aes(x = labels.p, y = measure, fill = labels.p)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
      scale_fill_manual(values = reference_celltype_colors) +
      scale_color_manual(values = reference_celltype_colors) +
      labs(x = "", y = data$y_labels[data$efeature == feature_name][1], title = feature_name) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none")

    return(plot)
  }
}

features <- unique(tidy_meta_table_supp$efeature)
plots_list <- list()

for (feature_name in features) {
  plots_list[[feature_name]] <- plot_feature(tidy_meta_table_supp, feature_name)
}

combined_efeats_plot <- wrap_plots(plots_list, ncol = 5)

ggsave(here(figs_path, "Supp-efeats.png"), combined_efeats_plot, height = 24, width = 24, bg = "white")

##########################################################################
# Supp: TRPV1 Expression Comparison
##########################################################################
trpv1_exp_plot <- plot_gene_expression_comparison(
  gene_symbol = "TRPV1",
  snrna_data = pig_ref,
  patchseq_data = patchseq,
  visium_data = spatial_drg,
  celltype_levels = reference_celltype_levels,
  celltype_colors = reference_celltype_colors,
  output_dir = figs_path,
  plot_filename_suffix = "_expr_comparison.png",
  verbose = TRUE
)

ggsave(here(figs_path, "Supp-TRPV1_expr.png"), trpv1_exp_plot, height = 7, width = 12, bg = "white")

##########################################################################
# Supp: SCN11A Expression Comparison
##########################################################################
scn11a_exp_plot <- plot_gene_expression_comparison(
  gene_symbol = "SCN11A",
  snrna_data = pig_ref,
  patchseq_data = patchseq,
  visium_data = spatial_drg,
  celltype_levels = reference_celltype_levels,
  celltype_colors = reference_celltype_colors,
  output_dir = figs_path,
  plot_filename_suffix = "_expr_comparison.png",
  verbose = TRUE
)

ggsave(here(figs_path, "Supp-SCN11A_expr.png"), scn11a_exp_plot, height = 7, width = 12, bg = "white")

##########################################################################
# Supp: Injury and Culturing Marker Violins
##########################################################################
theme_set(theme_cowplot(font_size = 14))

injury_theme <- function() {
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank()
  )
}

# ATF1
snrna_atf1_vln <- VlnPlot(pig_ref, features = 'ATF1', group.by = 'labels', assay = 'RNA', slot = 'data', cols = reference_celltype_colors, log = T, pt.size = 0.2, raster = FALSE) +
  ggtitle("snRNAseq") + NoLegend() + ylab(expression("ATF1 expr. ("*log[2]*" norm.)")) + injury_theme()

cpm_patchseq_atf1_vln <- VlnPlot(patchseq, features = 'ATF1', group.by = 'labels', assay = 'CPM', slot = 'data', cols = reference_celltype_colors, log = T, raster = F, pt.size = 0.2) +
  ggtitle("Patch-seq") + NoLegend() + ylab(expression("ATF1 expr. ("*log[2]*" CPM)")) + injury_theme()

atf1_dataset_comparison <- (snrna_atf1_vln | cpm_patchseq_atf1_vln) & theme(axis.title.x = element_blank())

# SOX11
snrna_sox11_vln <- VlnPlot(pig_ref, features = 'SOX11', group.by = 'labels', assay = 'RNA', slot = 'data', cols = reference_celltype_colors, log = T, pt.size = 0.2, raster = FALSE) +
  ggtitle("snRNAseq") + NoLegend() + ylab(expression("SOX11 expr. ("*log[2]*" norm.)")) + injury_theme()

cpm_patchseq_sox11_vln <- VlnPlot(patchseq, features = 'SOX11', group.by = 'labels', assay = 'CPM', slot = 'data', cols = reference_celltype_colors, log = T, raster = F, pt.size = 0.2) +
  ggtitle("Patch-seq") + NoLegend() + ylab(expression("SOX11 expr. ("*log[2]*" CPM)")) + injury_theme()

sox11_dataset_comparison <- (snrna_sox11_vln | cpm_patchseq_sox11_vln) & theme(axis.title.x = element_blank())

# CLIC1
snrna_clic1_vln <- VlnPlot(pig_ref, features = 'ENSSSCG00000039071', group.by = 'labels', assay = 'RNA', slot = 'data', cols = reference_celltype_colors, log = T, pt.size = 0.2, raster = FALSE) +
  ggtitle("snRNAseq") + NoLegend() + ylab(expression("CLIC1 expr. ("*log[2]*" norm.)")) + injury_theme()

cpm_patchseq_clic1_vln <- VlnPlot(patchseq, features = 'ENSSSCG00000039071', group.by = 'labels', assay = 'CPM', slot = 'data', cols = reference_celltype_colors, log = T, raster = F, pt.size = 0.2) +
  ggtitle("Patch-seq") + NoLegend() + ylab(expression("CLIC1 expr. ("*log[2]*" CPM)")) + injury_theme()

clic1_dataset_comparison <- (snrna_clic1_vln | cpm_patchseq_clic1_vln) & theme(axis.title.x = element_blank())

# CLIC4
snrna_clic4_vln <- VlnPlot(pig_ref, features = 'CLIC4', group.by = 'labels', assay = 'RNA', slot = 'data', cols = reference_celltype_colors, log = T, pt.size = 0.2, raster = FALSE) +
  ggtitle("snRNAseq") + NoLegend() + ylab(expression("CLIC4 expr. ("*log[2]*" norm.)")) + injury_theme()

cpm_patchseq_clic4_vln <- VlnPlot(patchseq, features = 'CLIC4', group.by = 'labels', assay = 'CPM', slot = 'data', cols = reference_celltype_colors, log = T, raster = F, pt.size = 0.2) +
  ggtitle("Patch-seq") + NoLegend() + ylab(expression("CLIC4 expr. ("*log[2]*" CPM)")) + injury_theme()

clic4_dataset_comparison <- (snrna_clic4_vln | cpm_patchseq_clic4_vln) & theme(axis.title.x = element_blank())

# Combined injury/culturing violins
inj_cult_vlns <- atf1_dataset_comparison / sox11_dataset_comparison / clic1_dataset_comparison / clic4_dataset_comparison + plot_annotation(tag_levels = 'A')
ggsave(here(figs_path, "Supp-injury_culturing_violins.png"), inj_cult_vlns, height = 20, width = 10, bg = "white")

##########################################################################
# Generate and save reference markers
##########################################################################
DefaultAssay(drg.integrated) <- "SCT"

if (!file.exists(here(figs_path, "subclass_reference_markers.csv"))) {
  Idents(object = drg.integrated) <- "labels.p"
  drg.integrated <- PrepSCTFindMarkers(drg.integrated)
  subclass.reference.markers <- FindAllMarkers(drg.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write_csv(subclass.reference.markers, here(figs_path, "subclass_reference_markers.csv"))
  cat("Reference markers saved to:", here(figs_path, "subclass_reference_markers.csv"), "\n")
} else {
  cat("Loading pre-computed reference markers from:", here(figs_path, "subclass_reference_markers.csv"), "\n")
  subclass.reference.markers <- read_csv(here(figs_path, "subclass_reference_markers.csv"))
}

subclass.reference.markers$cluster <- factor(subclass.reference.markers$cluster, levels = reference_celltype_levels)

top_3_subclass_markers <- subclass.reference.markers %>%
  filter(str_starts(gene, pattern = 'ENSSSCG', negate = TRUE)) %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC) %>%
  arrange(factor(cluster, levels = reference_celltype_levels))

write_csv(top_3_subclass_markers, here(figs_path, "top3_subclass_reference_markers.csv"))

