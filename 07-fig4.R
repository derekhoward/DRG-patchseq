library(here)
library(Seurat)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(viridis)
library(cluster)
library(factoextra)
library(impute)
library(tidymodels)
library(recipes)
library(ggrepel)
library(limma)
library(edgeR)
library(ggVennDiagram)

theme_set(theme_cowplot())
# Source utility functions
source(here("utils.R"))

base_font_size <- 24

# Set up results path
results_path <- here("results/figures/fig4")
dir.create(results_path, recursive = TRUE)

# Load and prepare data
patchseq <- prep_patchseq_obj_all_batches(convert_genes = "pig")
patchseq_meta <- read.csv(here('results/SCT-integration/predictions_meta.csv'), row.names = 'cell_id') %>% 
  dplyr::select(labels.p)
patchseq_meta$labels.p[patchseq_meta$labels.p == "C-OSMR-IL31RA"] <- "C-OSMR-SST"
patchseq <- AddMetaData(patchseq, metadata = patchseq_meta)

# Prepare ephys data
ephys_table <- FetchData(object = patchseq, vars = c(
  "labels.p", "APD", "X50Hz", "Sinusscore_old", "Delta.ind.time"))
rownames(ephys_table) <- Cells(patchseq)

# Impute missing values
ephys_table$id <- rownames(ephys_table)
recipe_obj <- recipe(~ ., data = ephys_table[, -1]) %>%
  step_impute_knn(all_predictors(), neighbors = 5)
prep_obj <- prep(recipe_obj)
ephys_data_imputed <- bake(prep_obj, new_data = NULL)
rownames(ephys_data_imputed) <- ephys_data_imputed$id
ephys_data_imputed$id <- NULL

# Normalize and cluster data
ephys_data_scaled <- scale(ephys_data_imputed)
distance_matrix <- dist(ephys_data_scaled, method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")
clusters <- cutree(hc, k = 5)
ephys_table$Cluster <- as.factor(clusters)

# Perform PCA
pca_result <- prcomp(ephys_data_scaled, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$PC1 <- pca_scores$PC1 * -1

# Prepare for plotting
cell_types <- levels(factor(ephys_table$labels.p))
custom_colors <- viridis(16, option = "turbo")
osmr_index <- which(cell_types == "C-OSMR-SST")
tac1_index <- which(cell_types == "C-TAC1-LRP1B")
custom_colors[c(osmr_index, tac1_index)] <- custom_colors[c(tac1_index, osmr_index)]
custom_colors[osmr_index] <- "black"

# Create PCA plot
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = ephys_table$labels.p)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = reference_celltype_colors, name = "Cell-type") +
  labs(title = "PCA of selected electrophysiology features",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(here(results_path, 'ephys_pca.png'), pca_plot)

# Create PCA biplot
loadings <- as.data.frame(pca_result$rotation)
loadings$feature <- rownames(loadings)
loadings$PC1 <- loadings$PC1 * -1
scaling_factor <- 5
loadings$PC1 <- loadings$PC1 * scaling_factor
loadings$PC2 <- loadings$PC2 * scaling_factor

# Rename the features
loadings$feature <- dplyr::case_when(
  loadings$feature == "X50Hz" ~ "50Hz",
  loadings$feature == "Sinusscore_old" ~ "Sinus score",
  loadings$feature == "Delta.ind.time" ~ "Rel. ADS",
  TRUE ~ loadings$feature  # Keep any other features unchanged
)

pca_biplot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = ephys_table$labels.p), 
             size = 3, alpha = 0.7) +  # Reduced point size
  scale_color_manual(values = reference_celltype_colors, name = "Cell-type") +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.1, "cm")), color = "black", alpha = 0.5) +  # Smaller arrows
  geom_text_repel(data = loadings, aes(x = PC1, y = PC2, label = feature),
                  color = "black", fontface = "bold", size = base_font_size/5, max.overlaps = 5) +  # Limit label overlaps
  labs(#title = "PCA Biplot of CMi-related\nelectrophysiology features",
    x = paste0("CMi ephys PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("CMi ephys PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_cowplot(font_size = base_font_size) +
  theme(legend.position = "none",
        plot.title = element_text(size = base_font_size * 1.1),
        axis.title = element_text(size = base_font_size * 1))

ggsave(here(results_path, '4A-ephys_pca_biplot.png'), pca_biplot)

# Add PC1 to metadata
PC1_meta <- data.frame(PC1 = pca_result$x[, 1] * -1)
rownames(PC1_meta) <- Cells(patchseq)
patchseq <- AddMetaData(patchseq, metadata = PC1_meta)
write.csv(patchseq@meta.data %>% 
            dplyr::select(labels.p, PC1), 
          here(results_path, 'PC1meta.csv'))

# Prepare data for differential expression analysis
prepare_data <- function(seurat_obj) {
  counts <- GetAssayData(seurat_obj, slot = "counts")
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge)
  return(dge)
}

run_de_analysis <- function(dge, design, coef_name) {
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = coef_name, number = Inf)
  return(results)
}

identify_degs <- function(results, fc_threshold = 1, p_threshold = 0.05) {
  degs <- rownames(results)[abs(results$logFC) > fc_threshold & results$adj.P.Val < p_threshold]
  return(degs)
}

ic_list <- readLines(here("data/processed/human_IC_genes.csv"))

# Define the specific genes you want to label 
genes_of_interest <- c('OSMR', 'IL31RA', 'SCN11A', 'TMEM45B')

create_volcano_plot <- function(results, fc_threshold = 1, p_threshold = 0.05) {
  results$diffexpressed <- "NO"
  results$diffexpressed[results$adj.P.Val < p_threshold & results$logFC > fc_threshold ] <- "YES"
  
  # Create a new column for coloring and labeling
  results$label_color <- "none"
  results$label_color[results$diffexpressed == "YES" & rownames(results) %in% ic_list] <- "black"
  
  # Create label column
  results$label <- ""
  results$label[results$label_color != "none"] <- rownames(results)[results$label_color != "none"]
  
  ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
    # Plot all points
    geom_point(color = "grey", alpha = 0.5) +
    # Plot significant points
    geom_point(data = subset(results, diffexpressed == "YES"), 
               color = "red", alpha = 0.5) +
    # Plot ion channel genes
    geom_point(data = subset(results, label_color == "black"), 
               color = "black", alpha = 1) +
    # # Label ion channel genes
    geom_text_repel(data = subset(results, label_color == "black", size=5),
                    aes(label = label),
                    color = "black", fontface = "bold", max.overlaps = Inf) +
    # Label genes of interest
    geom_text_repel(data = subset(results, label_color == "blue"),
                    aes(label = label),
                    color = "blue", fontface = "bold", max.overlaps = Inf) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), col = "blue", linetype = "dashed") +
    geom_hline(yintercept = -log10(p_threshold), col = "blue", linetype = "dashed") +
    labs(#title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value") +
    theme_cowplot() +
    theme(legend.position = "none",
          plot.title = element_text(size = base_font_size * 1.1),
          axis.title = element_text(size = base_font_size * 1))
}


# Prepare data for DE analysis
dge <- prepare_data(patchseq)

# Analysis 1: PC1-associated genes
design_pc1 <- model.matrix(~ PC1, data = patchseq@meta.data)
results_pc1 <- run_de_analysis(dge, design_pc1, "PC1")
write.csv(results_pc1, file = here(results_path, "PC1_DE_results.csv"))


volcano_plot_pc1 <- create_volcano_plot(results_pc1, fc_threshold = 0.5) +
  theme(legend.position = "none")
ggsave(here(results_path, "4B-volcano_plot_PC1.png"), volcano_plot_pc1, width = 12, height = 10)

degs_to_visualize <- c("TRPC6", "SCN11A", "SCN2A", "KCNK7", "CACNA1B", "SCN9A",
                       "KCNH8", "KCNK2", "KCNB1", "KCNT2", "SCN10A", "KCNA2", 
                       "KCNH7", "TRPV2", "KCNG3", "MCOLN2", "KCNU1", "KCNA3",
                       "KCNN3", "KCNS2", "TRPV1")

###################################################################################
# Compare DEGs with snRNAseq markers
###################################################################################

snrna_markers <- read.csv(here("results/figures/fig1v3/subclass_reference_markers.csv"))

# Filter for C-OSMR-SST markers
cosmr_markers <- snrna_markers %>%
  filter(cluster == 'C-OSMR-SST')

# Identify overlapping genes
degs_pc1 <- identify_degs(results_pc1, fc_threshold = 0.5, p_threshold = 0.05)
overlapping_genes <- intersect(degs_pc1, cosmr_markers$gene)


# Print results
cat("Number of DEGs associated with PC1:", length(degs_pc1), "\n")
cat("Number of C-OSMR-SST markers from snRNAseq:", nrow(cosmr_markers), "\n")
cat("Number of overlapping genes:", length(overlapping_genes), "\n")

# Save DEG results
write.csv(results_pc1, file = here(results_path, "PC1_DE_results.csv"))
write.table(degs_pc1, file = here(results_path, 'degs_pc1.txt'), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(overlapping_genes, file = here(results_path, 'overlapping_genes.txt'), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.csv(cosmr_markers, file = here(results_path, "OSMR_IL31RA_markers_snRNAseq.csv"))

# Load integrated data for gene overlap check
drg.integrated <- readRDS(here("results/SCT-integration", "drg_integrated.RDS"))
DefaultAssay(drg.integrated) <- 'RNA'
genes_overlap <- intersect(rownames(drg.integrated), rownames(patchseq))
message("Number of overlapping genes: ", length(genes_overlap))

###################################################################################
# Create Venn diagram
###################################################################################

# Create a list of gene sets
gene_sets <- list(
  "C-OSMR-SST markers" = cosmr_markers$gene,
  "CMi-associated genes" = degs_pc1
)

# Create Venn diagram with proper aspect ratio
venn_plot <- ggVennDiagram(gene_sets, label_alpha = 0,
                           label = "count",
                           category.names = c("C-OSMR-SST\nmarkers","CMi-associated\ngenes"),
                           label_size = 6,
                           set_size = 6) +
  scale_fill_gradient(low = "white", high = "lightblue") +
  theme_cowplot(font_size = base_font_size) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = base_font_size * 1.1),
        aspect.ratio = 1,
        plot.margin = margin(20, 5, 20, 5, "pt")) +
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  scale_y_continuous(expand = expansion(mult = 0.3))

# Perform hypergeometric test
universe_size <- 20159
set1_size <- length(degs_pc1)
set2_size <- nrow(cosmr_markers)
overlap_size <- length(overlapping_genes)

hyper_test <- phyper(overlap_size - 1, set1_size, universe_size - set1_size, set2_size, lower.tail = FALSE)

format_scientific <- function(x) {
  parts <- strsplit(format(x, scientific = TRUE), "e")[[1]]
  base <- as.numeric(parts[1])
  exponent <- as.numeric(parts[2])
  sprintf("%.2f × 10⁻%d", base, abs(exponent))
}

venn_plot <- venn_plot +
  labs(caption = paste("Hypergeometric test p-value:", format_scientific(hyper_test))) +
  theme(plot.caption = element_text(size = 12, margin = margin(t = 20)))

# Save the Venn diagram with increased size
ggsave(here(results_path, "4C-venn_diagram_PC1_COSMR_markers.png"), venn_plot, width = 12, height = 10, dpi = 300)



plot_gene_expression_vs_pc1 <- function(seurat_obj, gene, reference_celltype_colors) {
  expression_data <- GetAssayData(seurat_obj, slot = "data")[gene, ]
  plot_data <- data.frame(
    PC1 = seurat_obj$PC1,
    Expression = log1p(expression_data),
    Celltype = seurat_obj$labels.p
  )
  
  # Calculate mean values for each cell type
  celltype_means <- plot_data %>%
    group_by(Celltype) %>%
    summarize(
      mean_PC1 = mean(PC1),
      mean_Expression = mean(Expression)
    )
  
  correlation <- cor.test(plot_data$PC1, plot_data$Expression)
  r_value <- round(correlation$estimate, 3)
  p_value <- correlation$p.value
  
  ggplot(plot_data, aes(x = Expression, y = PC1, color = Celltype)) +
    geom_point(alpha = 0.5, size = 3) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    geom_point(data = celltype_means, aes(x = mean_Expression, y = mean_PC1, fill = Celltype), 
               size = 5, shape = 23, color = "black") +
    scale_color_manual(values = reference_celltype_colors) +
    scale_fill_manual(values = reference_celltype_colors) +
    labs(x = paste0(gene, " log(Expression + 1)"),
         y = "CMi score") +
    theme_cowplot(font_size = base_font_size) +
    theme(
      legend.position = "right",
      legend.justification = "center",
      legend.box.just = "right",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0, "pt"),
      legend.spacing.y = unit(0.1, "cm"),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
      axis.title = element_text(size = base_font_size * 1),
      legend.title = element_text(size = base_font_size * 0.6),
      legend.text = element_text(size = base_font_size * 0.6),
      legend.key.size = unit(0.8, "lines")  # Reduce the size of legend keys
    ) +
    guides(
      color = guide_legend(
        title = "Celltype",
        title.position = "top",
        title.hjust = 0.5,
        ncol = 1,
        byrow = TRUE,
        override.aes = list(alpha = 1, size = 3)  # Make legend points more visible
      ),
      fill = "none"  # Hide the fill legend since it's redundant
    ) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("R = ", r_value, "; p = ", format.pval(p_value, digits = 2)),
             size = base_font_size/4)
}

gene_plot <- plot_gene_expression_vs_pc1(patchseq, "SCN11A", reference_celltype_colors) +
  theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing = unit(0.1, "cm"),
    plot.margin = margin(5.5, 15.5, 5.5, 5.5)
  )

ggsave(here(results_path, "4D-gene_expression_vs_CMi_PC1_SCN11A.png"),
       gene_plot,
       width = 8,
       height = 6,
       dpi = 300)

# Create final composite figure using top_row layout
top_row <- (
  pca_biplot + 
    volcano_plot_pc1 + 
    plot_layout(widths = c(1.1, 1.1))
) | (
  venn_plot + 
    plot_layout(width = 0.7)
)

final_f4 <- top_row / (
  gene_plot + 
    gene_plot + 
    plot_layout(guides = 'collect')
) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A',
                  theme = theme(plot.tag = element_text(size = 16, face = "bold")))

# Save final figure
ggsave(
  here(results_path, "fig4.png"),
  final_f4,
  width = 18,
  height = 16,
  bg = "white",
  dpi = 300,
  limitsize = FALSE
)

###################################################################################
# Generate individual gene expression vs CMi plots
###################################################################################

cmi_exp_corr_results_path <- here(results_path, "cmi_geneExp_corrs")
dir.create(cmi_exp_corr_results_path, recursive = TRUE)

purrr::walk(degs_to_visualize, function(current_gene) {
  # Generate the plot for the current gene
  gene_plot <- plot_gene_expression_vs_pc1(patchseq, current_gene, reference_celltype_colors) +
    theme(
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing = unit(0.1, "cm"),
      plot.margin = margin(5.5, 15.5, 5.5, 5.5)
    )
  
  # Define the filename using file.path for robustness
  plot_filename <- file.path(cmi_exp_corr_results_path, paste0(current_gene, "_cmi_vs_expression.png"))
  
  ggsave(
    filename = plot_filename,
    plot = gene_plot,
    width = 8,  # Adjust width as needed
    height = 6, # Adjust height as needed
    dpi = 300   # Adjust DPI as needed
  )
  
  message("Saved plot for gene: ", current_gene, " to ", plot_filename)
})

message("All plots have been generated and saved to: ", cmi_exp_corr_results_path)

###################################################################################
# Create Faceted Figure of All CMI vs. Expression Plots (Alphabetical)
###################################################################################

message("Starting generation of alphabetical faceted CMI vs. Expression plot...")

# 1. Define and validate genes to plot
genes_to_plot <- intersect(degs_to_visualize, rownames(patchseq))
if (length(genes_to_plot) == 0) {
  stop("None of the genes in 'degs_to_visualize' were found in the patchseq object.")
}
genes_to_plot_alphabetical <- sort(genes_to_plot)

# 2. Prepare combined data frame for plotting
plot_meta <- patchseq@meta.data %>%
  dplyr::select(PC1, labels.p)
plot_meta$cell_id <- rownames(plot_meta)

# Get expression data and reshape to long format
expression_matrix <- GetAssayData(patchseq, slot = "data")[genes_to_plot_alphabetical, , drop = FALSE]
expression_long <- t(expression_matrix) %>%
  as.data.frame() %>%
  log1p()
expression_long$cell_id <- rownames(expression_long)

expression_long_tidy <- expression_long %>%
  pivot_longer(
    cols = -cell_id,
    names_to = "Gene",
    values_to = "Expression"
  )

# Join metadata and expression data
combined_plot_data <- plot_meta %>%
  left_join(expression_long_tidy, by = "cell_id") %>%
  rename(Celltype = labels.p) %>%
  mutate(Gene = factor(Gene, levels = genes_to_plot_alphabetical))

# 3. Calculate summary statistics for plotting
combined_celltype_means <- combined_plot_data %>%
  group_by(Celltype, Gene) %>%
  summarize(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    mean_Expression = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  )

correlation_stats <- combined_plot_data %>%
  group_by(Gene) %>%
  summarize(
    correlation = list(cor.test(PC1, Expression)),
    .groups = "drop"
  ) %>%
  mutate(
    r_value = round(map_dbl(correlation, "estimate"), 3),
    p_value = map_dbl(correlation, "p.value"),
    label = paste0("R = ", r_value, "\np = ", format.pval(p_value, digits = 2))
  )

# 4. Create the faceted plot
faceted_gene_plot <- ggplot(combined_plot_data, aes(x = Expression, y = PC1, color = Celltype)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  geom_point(data = combined_celltype_means,
             aes(x = mean_Expression, y = mean_PC1, fill = Celltype),
             size = 2.5, shape = 23, color = "black") +
  scale_color_manual(values = reference_celltype_colors) +
  scale_fill_manual(values = reference_celltype_colors) +
  geom_text(data = correlation_stats,
            aes(x = -Inf, y = Inf, label = label),
            color = "black",
            hjust = -0.1, vjust = 1.1,
            size = base_font_size / 6,
            inherit.aes = FALSE) +
  labs(x = "Gene Expression [log(counts + 1)]",
       y = "CMi score") +
  theme_cowplot(font_size = base_font_size * 0.9) +
  facet_wrap(~ Gene, scales = "free_x", ncol = 6) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.box.just = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0.1, "cm"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
    axis.title = element_text(size = base_font_size * 0.9),
    legend.title = element_text(size = base_font_size * 0.6),
    legend.text = element_text(size = base_font_size * 0.6),
    legend.key.size = unit(0.8, "lines"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = base_font_size * 0.7, face = "bold")
  ) +
  guides(
    color = guide_legend(
      title = "Celltype",
      title.position = "top",
      title.hjust = 0.5,
      ncol = 1,
      byrow = TRUE,
      override.aes = list(alpha = 1, size = 3)
    ),
    fill = "none"
  )

# 5. Save the faceted plot
num_rows <- ceiling(length(genes_to_plot_alphabetical) / 6)
plot_height <- num_rows * 4
plot_width <- 22

ggsave(
  filename = here(results_path, "cmi_vs_expression_FACETED_alphabetical.png"),
  plot = faceted_gene_plot,
  width = plot_width,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)

message(paste0("Alphabetized faceted plot saved to: ", here(results_path, "cmi_vs_expression_FACETED_alphabetical.png")))
