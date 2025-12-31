library(here)
library(ggplot2)
library(MetaNeighbor)
library(kableExtra)
library(SummarizedExperiment)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)

output_path <- here("results/figures/fig2_metaneighbor/")
dir.create(output_path, recursive = TRUE)

pig_ref <- readRDS(here('./data/summarized_experiments/pig_ref_se.RDS'))
atlas <- readRDS(here('./data/summarized_experiments/drg_complete_SE.RDS'))

my_data <- list(
  pig_snRNA = pig_ref,
  integrated_atlas = atlas
)

colData(my_data$integrated_atlas)$labels <- colData(my_data$integrated_atlas)$Atlas_annotation

lapply(my_data, function(x) head(rownames(x), 3))
lapply(my_data, function(x) colnames(colData(x)))
lapply(my_data, function(x) names(assays(x)))

# Merge data
fused_data <- mergeSCE(my_data)
dim(fused_data)
head(colData(fused_data))
table(fused_data$labels, fused_data$study_id)

# =============================================================================
# MetaNeighbor analysis
# =============================================================================
global_hvgs <- variableGenes(dat = fused_data, exp_labels = fused_data$study_id)
length(global_hvgs)

aurocs <- MetaNeighborUS(
  var_genes = global_hvgs,
  dat = fused_data,
  study_id = fused_data$study_id,
  cell_type = fused_data$labels,
  fast_version = TRUE
)

# Save heatmaps
metaneighborheatmap <- ggPlotHeatmap(aurocs)
ggsave(
  filename = here(output_path, 'ggplot_metaneighborheatmap_small_origlabels.png'),
  plot = metaneighborheatmap,
  width = 12,
  height = 12,
  bg = "white"
)

png(here(output_path, "small_metaneighborHeatmap_origlabels.png"), width = 14, height = 13, units = "in", res = 300)
plotHeatmap(aurocs, cex = 0.5)
dev.off()

# =============================================================================
# Top hits
# =============================================================================
top_hits_by_study <- topHitsByStudy(aurocs, threshold = 0.7, n_digits = 2, collapse_duplicates = TRUE)
print(top_hits_by_study, n = 17)

uniquehits <- top_hits_by_study %>%
  filter(str_detect(`Study_ID|Celltype_1`, 'pig')) %>%
  kbl() %>%
  kable_styling()

top_hits_by_study %>%
  filter(str_detect(`Study_ID|Celltype_1`, 'pig')) %>%
  write_csv(here(output_path, 'top_hits_by_study.csv'))

top_hits_by_study %>%
  filter(str_detect(`Study_ID|Celltype_1`, 'atlas'))

topHits(aurocs, dat = fused_data, study_id = fused_data$study_id,
        cell_type = fused_data$labels, threshold = 0.9) %>%
  filter(str_detect(`Study_ID|Celltype_1`, 'pig'))

# =============================================================================
# Hierarchical splitting
# =============================================================================
full_labels <- makeClusterName(fused_data$study_id, fused_data$labels)

# Level 1 split
level1_split <- splitClusters(aurocs, k = 2)
first_split <- level1_split[[2]]

subdata <- fused_data[, full_labels %in% first_split]
dim(subdata)

var_genes <- variableGenes(dat = subdata, exp_labels = subdata$study_id)

sub_aurocs <- MetaNeighborUS(
  var_genes = var_genes,
  dat = subdata,
  fast_version = TRUE,
  study_id = subdata$study_id,
  cell_type = subdata$labels
)
plotHeatmap(sub_aurocs, cex = 0.7)

# Level 2 split
level2_split <- splitClusters(sub_aurocs, k = 3)
my_split <- level2_split[[3]]
subdata <- fused_data[, full_labels %in% my_split]
var_genes <- variableGenes(dat = subdata, exp_labels = subdata$study_id)
length(var_genes)

# =============================================================================
# Best hits analysis
# =============================================================================
best_hits <- MetaNeighborUS(
  var_genes = global_hvgs,
  dat = fused_data,
  study_id = fused_data$study_id,
  cell_type = fused_data$labels,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

write_csv(best_hits %>% as.data.frame() %>% tibble::rownames_to_column("rowname"), 
          here(output_path, 'best_hits_matrix.csv'))

tidy_best_hits_df <- best_hits %>%
  as_tibble(rownames = "from") %>%
  pivot_longer(cols = -from, names_to = "to", values_to = "similarity") %>%
  filter(!is.na(similarity)) %>%
  arrange(desc(similarity))

write_csv(tidy_best_hits_df, here(output_path, 'best_hits_df.csv'))

png(here(output_path, "best_hits_Heatmap.png"), width = 7, height = 6.5, units = "in", res = 1200)
plotHeatmap(best_hits, cex = 0.5)
dev.off()

# =============================================================================
# Cluster graph visualization
# =============================================================================
cluster_graph <- makeClusterGraph(best_hits, low_threshold = 0.3)

png(here(output_path, "cluster_graph_lowthr3_sizef45.png"), width = 8, height = 8, units = "in", res = 300)
plotClusterGraph(cluster_graph, fused_data$study_id, fused_data$labels, size_factor = 3.5)
dev.off()