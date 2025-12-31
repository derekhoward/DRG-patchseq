library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(ggsankey)
library(scales)

data <- read_csv(here("results/figures/fig2/snRNAseq_predicted_celltypes.csv"))


# Define the reference cell type levels and colors
reference_celltype_levels <- c("Aβ-PROPRIO", "Aβ-LTMR-ALDH1A1", "Aβ-LTMR-SLIT2", "Aβ-LTMR-PALM", "Aβ-HTMR", 
                               "Aδ-LTMR", "Aδ-CACNA1E", "Aδ-TRPV1", "Aδ-COOL", "C-COLD", 
                               "C-TAC1-KCNQ5", "C-TAC1-LRP1B", "C-TAC1-TRPA1", "C-LTMR", 
                               "C-OSMR-GFRA1_2", "C-OSMR-SST")

reference_celltype_colors <- c(
  "C-TAC1-KCNQ5" = "#F8766D", "C-TAC1-LRP1B" = "#E88526", "C-LTMR" = "#D39200", "Aβ-HTMR" = "#B79F00",
  "C-OSMR-SST" = "#93AA00", "Aδ-COOL" = "#5EB300", "C-TAC1-TRPA1" = "#00BA38", "Aβ-LTMR-PALM" = "#00BF74",
  "Aδ-TRPV1" = "#00C19F", "C-COLD" = "#00BFC4", "Aβ-PROPRIO" = "#00B9E3", "Aβ-LTMR-SLIT2" = "#00ADFA",
  "C-OSMR-GFRA1_2" = "#619CFF", "Aδ-LTMR" = "#DB72FB", "Aβ-LTMR-ALDH1A1" = "#F564E3",
  "Aδ-CACNA1E" = "#FF61C3", "patchseq" = "#c00000"
)

cross_species_order <- rev(c("Pvalb", "Ntrk3high+Ntrk2", "Ntrk3high+S100a16", "Ntrk3low+Ntrk2", 
                             "Calca+Bmpr1b", "Calca+Dcn", "Rxfp1", "Trpm8", "Calca+Oprk1", "Calca+Adra2a", 
                             "Calca+Smr2", "Calca+Sstr2", "Mrgpra3+Mrgprb4", "Mrgpra3+Trpv1", "Th", "Mrgprd", "Sst", 
                             "Atf3"))

# Generate colors in alphabetical order
seurat_colors <- setNames(hue_pal()(length(cross_species_order)), sort(cross_species_order))

sankey_data <- data %>%
  select(labels, predicted.celltype) %>%
  mutate('Cross-species taxonomy labels (projected)' = factor(predicted.celltype, levels = cross_species_order),
         'Pig snRNAseq taxonomy labels (original)' = factor(labels, levels = rev(reference_celltype_levels))) %>%
  make_long('Cross-species taxonomy labels (projected)', 'Pig snRNAseq taxonomy labels (original)')

# Manually reorder the factors after make_long
sankey_data <- sankey_data %>%
  mutate(node = factor(node, levels = c(cross_species_order, rev(reference_celltype_levels))),
         x = factor(x, levels = c("Cross-species taxonomy labels (projected)", "Pig snRNAseq taxonomy labels (original)")))

# Combine both color palettes
combined_colors <- c(seurat_colors, reference_celltype_colors)

sankey_plot <- ggplot(sankey_data, aes(x = x, 
                                       next_x = next_x, 
                                       node = node, 
                                       next_node = next_node,
                                       fill = factor(node),
                                       label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black") +
  geom_sankey_label(size = 5, color = "black", fill = "white") +
  theme_sankey(base_size = 20) +
  labs(title = NULL,
       subtitle = NULL,
       x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_fill_manual(values = combined_colors)

ggsave(here("results/figures/fig2/sankey_plot.png"), sankey_plot, width = 12, height = 10, dpi = 300)