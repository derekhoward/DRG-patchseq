library(here)
library(Seurat)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggpubr)
source(here("utils.R"))

theme_set(theme_cowplot(font_size = 20))
figs_path <- here('results/figures/fig3')
dir.create(figs_path, recursive = TRUE)

# Color palette for microneurography plots
microneurography_colors <- list(
  'c_fibres' = '#FF6347',
  'nocis' = '#93AA00'
)

# Cell types of interest
celltypes_of_interest <- c(
  "C-COLD", "C-TAC1-KCNQ5", "C-TAC1-LRP1B", "C-TAC1-TRPA1",
  "C-LTMR", "C-OSMR-GFRA1_2", "C-OSMR-SST"
)

drg.integrated <- readRDS(here("results/SCT-integration/drg_integrated.RDS"))
obj.list <- SplitObject(drg.integrated, split.by = "dataset")
reference <- obj.list[[1]]
int_patchseq <- obj.list[[2]]

# Extract metadata
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

feats_to_select <- meta_table %>%
  select_if(is.numeric) %>%
  names()

tidy_meta_table <- meta_table %>%
  tibble::rownames_to_column(var = 'cell_id') %>%
  pivot_longer(cols = all_of(feats_to_select), names_to = "efeature", values_to = "measure")

# Filter for cell types of interest
filtered_for_celltypes <- tidy_meta_table %>%
  filter(labels.p %in% celltypes_of_interest) %>%
  droplevels()

filtered_for_celltypes$labels.p <- factor(
  filtered_for_celltypes$labels.p,
  levels = celltypes_of_interest
)

# =============================================================================
# Figure 3E: Pruriceptor Electrophysiology Features
# =============================================================================

# Sinus score plot
my_comparisons <- list(
  c("C-OSMR-SST", "C-OSMR-GFRA1_2"),
  c("C-OSMR-SST", "C-TAC1-TRPA1"),
  c("C-OSMR-SST", "C-TAC1-LRP1B"),
  c("C-OSMR-SST", "C-TAC1-KCNQ5")
)
# Sinus score plot
sinusscore <- filtered_for_celltypes %>%
  filter(efeature == "Sinusscore_old") %>%
  ggplot(aes(x = labels.p, y = measure)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = labels.p), size = 2, alpha = 0.75, position = position_jitter(width = 0.1)) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif') +
  scale_colour_manual(values = reference_celltype_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank()) +
  labs(x = "", y = "Sinus score (AP difference)") +
  guides(color = FALSE) +
  ggtitle("Sinus score")

# 50 Hz following frequency plot
my_comparisons <- list(
  c("C-OSMR-SST", "C-LTMR"),
  c("C-OSMR-SST", "C-TAC1-TRPA1"),
  c("C-OSMR-SST", "C-TAC1-LRP1B"),
  c("C-OSMR-SST", "C-TAC1-KCNQ5")
)
follow_freq <- filtered_for_celltypes %>%
  filter(efeature == "X50Hz") %>%
  ggplot(aes(x = labels.p, y = measure)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = labels.p), size = 2, alpha = 0.75, position = position_jitter(width = 0.1)) +
  scale_colour_manual(values = reference_celltype_colors) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank()) +
  labs(x = "", y = "AP responses (Count)") +
  guides(color = FALSE) +
  ggtitle("50 Hz following frequency")

# Peak delay plot
my_comparisons <- list(
  c("C-OSMR-SST", "C-TAC1-TRPA1"),
  c("C-OSMR-SST", "C-TAC1-LRP1B"),
  c("C-OSMR-SST", "C-TAC1-KCNQ5")
)
peak_delay <- filtered_for_celltypes %>%
  filter(efeature == "Delta.ind.time") %>%
  filter((measure > -1) & (measure < 2)) %>%
  ggplot(aes(x = labels.p, y = measure)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = labels.p), size = 2, alpha = 0.75, position = position_jitter(width = 0.1)) +
  scale_colour_manual(values = reference_celltype_colors) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank()) +
  labs(x = "", y = "Time (ms)") +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  guides(color = FALSE) +
  ggtitle("Activity dependent\ntime to peak delay")

# Save Figure 3E
fig3e <- sinusscore + follow_freq
ggsave(here(figs_path, "3e-pruri_ephys_feats.png"), fig3e, height = 8, width = 12, bg = "white")

# ANOVA statistics


anova_sinus <- aov(formula = measure ~ labels.p, data = filtered_for_celltypes %>% filter(efeature == 'Sinusscore_old'))
summary(anova_sinus)
TukeyHSD(anova_sinus)

anova_50hz <- aov(formula = measure ~ labels.p, data = filtered_for_celltypes %>% filter(efeature == 'X50Hz'))
summary(anova_50hz)
TukeyHSD(anova_50hz)

# =============================================================================
# Figure 3E (with ANOVA labels version)
# =============================================================================

my_comparisons <- list(
  c("C-OSMR-SST", "C-OSMR-GFRA1_2"),
  c("C-OSMR-SST", "C-TAC1-TRPA1"),
  c("C-OSMR-SST", "C-TAC1-LRP1B"),
  c("C-OSMR-SST", "C-TAC1-KCNQ5")
)
sinusscore <- filtered_for_celltypes %>%
  filter(efeature == "Sinusscore_old") %>%
  ggplot(aes(x = labels.p, y = measure)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = labels.p), size = 2, alpha = 0.75, position = position_jitter(width = 0.1)) +
  stat_compare_means(method = "anova", label.y = 25) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif') +
  scale_colour_manual(values = reference_celltype_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank()) +
  labs(x = "", y = "Sinus score") +
  guides(color = FALSE) +
  ggtitle("Sinus score")

my_comparisons <- list(
  c("C-OSMR-SST", "C-TAC1-TRPA1"),
  c("C-OSMR-SST", "C-TAC1-LRP1B"),
  c("C-OSMR-SST", "C-TAC1-KCNQ5"),
  c("C-OSMR-SST", "C-LTMR")
)
follow_freq <- filtered_for_celltypes %>%
  filter(efeature == "X50Hz") %>%
  ggplot(aes(x = labels.p, y = measure)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = labels.p), size = 2, alpha = 0.75, position = position_jitter(width = 0.1)) +
  scale_colour_manual(values = reference_celltype_colors) +
  stat_compare_means(method = "anova", label.y = 31) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank()) +
  labs(x = "", y = "AP responses") +
  guides(color = FALSE) +
  ggtitle("50 Hz following frequency")

my_comparisons <- list(
  c("C-OSMR-SST", "C-TAC1-TRPA1"),
  c("C-OSMR-SST", "C-TAC1-LRP1B"),
  c("C-OSMR-SST", "C-TAC1-KCNQ5")
)
peak_delay <- filtered_for_celltypes %>%
  filter(efeature == "Delta.ind.time") %>%
  filter((measure > -1) & (measure < 2)) %>%
  ggplot(aes(x = labels.p, y = measure)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = labels.p), size = 2, alpha = 0.75, position = position_jitter(width = 0.1)) +
  scale_colour_manual(values = reference_celltype_colors) +
  stat_compare_means(method = "anova", label.y = 3.1) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_blank()) +
  labs(x = "", y = "Time (ms)") +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  guides(color = FALSE) +
  ggtitle("Activity dependent\ntime to peak delay")

fig3e <- sinusscore + follow_freq + peak_delay
ggsave(here(figs_path, "3e-pruri_ephys_feats_ANOVA.png"), fig3e, height = 8, width = 18, bg = "white")

# =============================================================================
# Activity-Dependent Slowing (ADS) Analysis
# =============================================================================

# -----------------------------------------------------------------------------
# Prepare Data (excluding outliers)
# -----------------------------------------------------------------------------
drg.integrated$sample_ids <- Cells(drg.integrated)
to_exclude <- c('16287', '16730', '16740', '16757')
sample_ids_to_keep <- setdiff(drg.integrated$sample_ids, to_exclude)
drg.integrated <- subset(drg.integrated, subset = sample_ids %in% sample_ids_to_keep)

obj.list <- SplitObject(drg.integrated, split.by = "dataset")

reference <- obj.list[[1]]
int_patchseq <- obj.list[[2]]

all_meta <- names(int_patchseq@meta.data)
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

# -----------------------------------------------------------------------------
# Create ADS Feature Vectors
# -----------------------------------------------------------------------------
ads_rmp_vector <- sapply(1:75, function(i) paste0("ADS_RMP_Pulse", i))
ads_slope_vector <- sapply(1:75, function(i) paste0("ADS_Slope_Pulse", i))
ads_peak_vector <- sapply(1:75, function(i) paste0("ADS_Peak_Pulse", i))
ads_ttp_vector <- sapply(1:75, function(i) paste0("ADS_TTP_Pulse", i))
ads_RT_vector <- sapply(1:75, function(i) paste0("ADS_RT_ABS", i))

# Extract ADS metadata tables
rmp_meta_table <- FetchData(object = int_patchseq, vars = c("labels.p", ads_rmp_vector)) %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = all_of(ads_rmp_vector), names_to = "efeature", values_to = "measure")

slope_meta_table <- FetchData(object = int_patchseq, vars = c("labels.p", ads_slope_vector)) %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = all_of(ads_slope_vector), names_to = "efeature", values_to = "measure")

peak_meta_table <- FetchData(object = int_patchseq, vars = c("labels.p", ads_peak_vector)) %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = all_of(ads_peak_vector), names_to = "efeature", values_to = "measure")

ttp_meta_table <- FetchData(object = int_patchseq, vars = c("labels.p", ads_ttp_vector)) %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = all_of(ads_ttp_vector), names_to = "efeature", values_to = "measure")

RT_meta_table <- FetchData(object = int_patchseq, vars = c("labels.p", ads_RT_vector)) %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = all_of(ads_RT_vector), names_to = "efeature", values_to = "measure")

# =============================================================================
# ADS RMP Plot
# =============================================================================

adsrmp <- rmp_meta_table %>%
  filter(labels.p %in% celltypes_of_interest) %>%
  filter(!is.na(measure)) %>%
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_RMP_Pulse"))) %>%
  mutate(measure = if_else(f_idx == 1, 0, measure)) %>%
  mutate(efeature = factor(efeature, levels = paste0("ADS_RMP_Pulse", 1:75))) %>%
  group_by(rowname) %>%
  mutate(Celltype = ifelse(labels.p == "C-OSMR-SST", "C-OSMR-SST", "C fibres")) %>%
  group_by(Celltype, efeature) %>%
  summarise(
    mean = mean(measure),
    se = sqrt(sum((measure - mean(measure))^2 / (length(measure) - 1))) / sqrt(length(measure))
  ) %>%
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_RMP_Pulse"))) %>%
  arrange(f_idx) %>%
  ggplot(aes(x = f_idx, y = mean, color = Celltype)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, linewidth = 1.5) +
  scale_color_manual(
    values = c("C-OSMR-SST" = microneurography_colors$nocis, "C fibres" = microneurography_colors$c_fibres),
    labels = c("C fibres" = "C-fibres, n=101", "C-OSMR-SST" = "C-OSMR-SST, n=34")
  ) +
  scale_x_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 15), labels = seq(0, 75, by = 15)) +
  labs(x = "Pulse (repetition at 2Hz)", y = "RMP decrease (mV)") +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.2, 0.97)
  )

# ADS RMP Statistics
rmp_table_for_stats <- rmp_meta_table %>%
  filter(labels.p %in% celltypes_of_interest) %>%
  filter(!is.na(measure)) %>%
  mutate(efeature = factor(efeature, levels = paste0("ADS_RMP_Pulse", 1:75))) %>%
  group_by(rowname) %>%
  mutate(relative_slowing = (measure / first(measure)) * 100) %>%
  mutate(Celltype = ifelse(labels.p == "C-OSMR-SST", "C-OSMR-SST", "C fibres")) %>%
  mutate(Celltype = factor(Celltype, levels = c("C fibres", "C-OSMR-SST"))) %>% 
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_RMP_Pulse"))) %>% 
  dplyr::select(rowname, Celltype, f_idx, measure, relative_slowing)

RMPdata_final_f_idx <- filter(rmp_table_for_stats, f_idx == max(rmp_table_for_stats$f_idx))

RMPdata_final_f_idx %>%
  group_by(Celltype) %>%
  summarise(
    mean_value = mean(measure, na.rm = TRUE),
    sem = sd(measure, na.rm = TRUE) / sqrt(n())
  )

RMPdata_final_f_idx %>%
  group_by(Celltype) %>%
  summarise(
    mean_value = mean(relative_slowing, na.rm = TRUE),
    sem = sd(relative_slowing, na.rm = TRUE) / sqrt(n())
  )

wilcox.test(relative_slowing ~ Celltype, data = RMPdata_final_f_idx, exact = FALSE)
wilcox.test(measure ~ Celltype, data = RMPdata_final_f_idx, exact = FALSE)

# =============================================================================
# ADS Slope Plot
# =============================================================================

adsslope <- slope_meta_table %>%
  filter(labels.p %in% celltypes_of_interest) %>%
  filter(!is.na(measure)) %>%
  filter(!rowname %in% c("16283")) %>%
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_Slope_Pulse"))) %>%
  mutate(efeature = factor(efeature, levels = paste0("ADS_Slope_Pulse", 1:75))) %>%
  mutate(measure = measure * 100) %>%
  group_by(rowname) %>%
  mutate(Celltype = ifelse(labels.p == "C-OSMR-SST", "C-OSMR-SST", "C fibres")) %>%
  group_by(Celltype, efeature) %>%
  summarise(
    mean = mean(measure),
    se = sqrt(sum((measure - mean(measure))^2 / (length(measure) - 1))) / sqrt(length(measure))
  ) %>%
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_Slope_Pulse"))) %>%
  arrange(f_idx) %>%
  ggplot(aes(x = f_idx, y = mean, color = Celltype)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, linewidth = 1.5) +
  scale_color_manual(values = c("C-OSMR-SST" = microneurography_colors$nocis, "C fibres" = microneurography_colors$c_fibres)) +
  scale_x_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 15), labels = seq(0, 75, by = 15)) +
  labs(x = "Pulse (repetition at 2Hz)", y = "Relative slope decrease (% of baseline)") +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.6, 0.97)
  )

# ADS Slope Statistics
slope_table_for_stats <- slope_meta_table %>%
  filter(labels.p %in% celltypes_of_interest) %>%
  filter(!is.na(measure)) %>%
  filter(!rowname %in% c("16283")) %>% 
  mutate(efeature = factor(efeature, levels = paste0("ADS_Slope_Pulse", 1:75))) %>%
  group_by(rowname) %>%
  mutate(relative_slowing = (measure / first(measure)) * 100) %>%
  mutate(Celltype = ifelse(labels.p == "C-OSMR-SST", "C-OSMR-SST", "C fibres")) %>%
  mutate(Celltype = factor(Celltype, levels = c("C fibres", "C-OSMR-SST"))) %>% 
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_Slope_Pulse"))) %>% 
  dplyr::select(rowname, Celltype, f_idx, measure, relative_slowing)

slopedata_final_f_idx <- filter(slope_table_for_stats, f_idx == max(slope_table_for_stats$f_idx))

slopedata_final_f_idx %>%
  group_by(Celltype) %>%
  summarise(
    mean_value = mean(measure, na.rm = TRUE),
    sem = sd(measure, na.rm = TRUE) / sqrt(n())
  )

slopedata_final_f_idx %>%
  group_by(Celltype) %>%
  summarise(
    mean_value = mean(relative_slowing, na.rm = TRUE),
    sem = sd(relative_slowing, na.rm = TRUE) / sqrt(n())
  )

wilcox.test(relative_slowing ~ Celltype, data = slopedata_final_f_idx, exact = FALSE)
wilcox.test(measure ~ Celltype, data = slopedata_final_f_idx, exact = FALSE)

# =============================================================================
# ADS Peak Plot
# =============================================================================

adspeak <- peak_meta_table %>%
  filter(labels.p %in% celltypes_of_interest) %>%
  filter(!is.na(measure)) %>%
  filter(!rowname %in% c("16283")) %>%
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_Peak_Pulse"))) %>%
  mutate(efeature = factor(efeature, levels = paste0("ADS_Peak_Pulse", 1:75))) %>%
  group_by(rowname) %>%
  mutate(Celltype = ifelse(labels.p == "C-OSMR-SST", "C-OSMR-SST", "C fibres")) %>%
  group_by(Celltype, efeature) %>%
  summarise(
    mean = mean(measure),
    se = sqrt(sum((measure - mean(measure))^2 / (length(measure) - 1))) / sqrt(length(measure))
  ) %>%
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_Peak_Pulse"))) %>%
  arrange(f_idx) %>%
  ggplot(aes(x = f_idx, y = mean, color = Celltype)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, linewidth = 1.5) +
  scale_color_manual(values = c("C-OSMR-SST" = microneurography_colors$nocis, "C fibres" = microneurography_colors$c_fibres)) +
  scale_x_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 15), labels = seq(0, 75, by = 15)) +
  labs(x = "Pulse (repetition at 2Hz)", y = "Action potential peak difference (mV)") +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.1, 0.97)
  )

# ADS Peak Statistics
peak_table_for_stats <- peak_meta_table %>%
  filter(labels.p %in% celltypes_of_interest) %>%
  filter(!is.na(measure)) %>%
  filter(!rowname %in% c("16283")) %>%
  mutate(efeature = factor(efeature, levels = paste0("ADS_Peak_Pulse", 1:75))) %>%
  group_by(rowname) %>%
  mutate(relative_slowing = (measure / first(measure)) * 100) %>%
  mutate(Celltype = ifelse(labels.p == "C-OSMR-SST", "C-OSMR-SST", "C fibres")) %>%
  mutate(Celltype = factor(Celltype, levels = c("C fibres", "C-OSMR-SST"))) %>% 
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_Peak_Pulse"))) %>% 
  dplyr::select(rowname, Celltype, f_idx, measure, relative_slowing)

peakdata_final_f_idx <- filter(peak_table_for_stats, f_idx == max(peak_table_for_stats$f_idx))

peakdata_final_f_idx %>%
  group_by(Celltype) %>%
  summarise(
    mean_value = mean(measure, na.rm = TRUE),
    sem = sd(measure, na.rm = TRUE) / sqrt(n())
  )

wilcox.test(relative_slowing ~ Celltype, data = peakdata_final_f_idx, exact = FALSE)
wilcox.test(measure ~ Celltype, data = peakdata_final_f_idx, exact = FALSE)

# =============================================================================
# ADS RT (Response Time) Plot
# =============================================================================

adsRT <- RT_meta_table %>%
  filter(labels.p %in% celltypes_of_interest) %>%
  filter(!is.na(measure)) %>%
  mutate(efeature = factor(efeature, levels = paste0("ADS_RT_ABS", 1:75))) %>%
  group_by(rowname) %>%
  mutate(relative_slowing = (measure / first(measure)) * 100) %>%
  mutate(Celltype = ifelse(labels.p == "C-OSMR-SST", "C-OSMR-SST", "C fibres")) %>%
  mutate(Celltype = factor(Celltype, levels = c("C-OSMR-SST", "C fibres"))) %>%
  group_by(Celltype, efeature) %>%
  summarise(
    mean = mean(relative_slowing),
    se = sqrt(sum((relative_slowing - mean(relative_slowing))^2 / (length(relative_slowing) - 1))) / sqrt(length(relative_slowing))
  ) %>%
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_RT_ABS"))) %>%
  arrange(f_idx) %>%
  ggplot(aes(x = f_idx, y = mean, color = Celltype)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, linewidth = 1) +
  scale_color_manual(
    values = c("C fibres" = microneurography_colors$c_fibres, "C-OSMR-SST" = microneurography_colors$nocis),
    labels = c("C fibres" = "Other C-fibres, n=101", "C-OSMR-SST" = "C-OSMR-SST, n=34")
  ) +
  scale_y_continuous(limits = c(100, 114)) +
  scale_x_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 15), labels = seq(0, 75, by = 15)) +
  labs(x = "Pulse (repetition at 2Hz)", y = "Relative slowing (% baseline)") +
  theme(legend.title = element_blank(), legend.position = c(0.03, 0.97)) +
  ggtitle('In vitro patch-clamp (pig)')

# ADS RT Statistics
RT_table_for_stats <- RT_meta_table %>%
  filter(labels.p %in% celltypes_of_interest) %>%
  filter(!is.na(measure)) %>%
  mutate(efeature = factor(efeature, levels = paste0("ADS_RT_ABS", 1:75))) %>%
  group_by(rowname) %>%
  mutate(relative_slowing = (measure / first(measure)) * 100) %>%
  mutate(Celltype = ifelse(labels.p == "C-OSMR-SST", "C-OSMR-SST", "C fibres")) %>%
  mutate(Celltype = factor(Celltype, levels = c("C fibres", "C-OSMR-SST"))) %>% 
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_RT_ABS"))) %>% 
  dplyr::select(rowname, Celltype, f_idx, measure, relative_slowing)

final_f_idx_data <- filter(RT_table_for_stats, f_idx == max(RT_table_for_stats$f_idx))

final_f_idx_data %>%
  group_by(Celltype) %>%
  summarise(
    mean_value = mean(measure, na.rm = TRUE),
    sem = sd(measure, na.rm = TRUE) / sqrt(n())
  )

final_f_idx_data %>%
  group_by(Celltype) %>%
  summarise(
    mean_value = mean(relative_slowing, na.rm = TRUE),
    sem = sd(relative_slowing, na.rm = TRUE) / sqrt(n())
  )

mw_ADS_test_result <- wilcox.test(relative_slowing ~ Celltype, data = final_f_idx_data, exact = FALSE)
wilcox.test(measure ~ Celltype, data = final_f_idx_data, exact = FALSE)

# =============================================================================
# Microneurography Data
# =============================================================================

# Pig microneurography
pig_microneurography <- read_csv(here('./data/microneurography/pig_microneurography.csv'))

pig_microneurography_adsRT <- pig_microneurography %>%
  tidyr::pivot_longer(!Pulse, names_to = 'fibre', values_to = 'relative_slowing') %>%
  arrange(fibre) %>%
  mutate(Celltype = case_when(
    str_starts(string = fibre, pattern = 'CMi') ~ 'CMi, n=8',
    TRUE ~ "CM, n=26"
  )) %>%
  mutate(Celltype = factor(Celltype, levels = c("CMi, n=8", "CM, n=26"))) %>%
  group_by(Celltype, Pulse) %>%
  summarise(
    mean = mean(relative_slowing),
    se = sqrt(sum((relative_slowing - mean(relative_slowing))^2 / (length(relative_slowing) - 1))) / sqrt(length(relative_slowing))
  ) %>%
  arrange(Pulse) %>%
  ggplot(aes(x = Pulse, y = mean, color = Celltype)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, linewidth = 1) +
  scale_color_manual(values = c("CMi, n=8" = microneurography_colors$nocis, "CM, n=26" = microneurography_colors$c_fibres)) +
  scale_x_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 15), labels = seq(0, 75, by = 15)) +
  labs(x = "Pulse", y = "Relative slowing (% baseline)") +
  theme(legend.title = element_blank(), legend.position = c(0.03, 0.97)) +
  ggtitle('In vivo microneurography (pig)')

# Human microneurography
human_microneurography <- read_csv(here('./data/microneurography/human_microneurography.csv'))

human_microneurography_adsRT <- human_microneurography %>%
  tidyr::pivot_longer(!Pulse, names_to = 'fibre', values_to = 'relative_slowing') %>%
  arrange(fibre) %>%
  mutate(Celltype = case_when(
    str_starts(string = fibre, pattern = 'CMi') ~ 'CMi, n=22',
    TRUE ~ 'CM, n=69'
  )) %>%
  mutate(Celltype = factor(Celltype, levels = c("CMi, n=22", "CM, n=69"))) %>%
  group_by(Celltype, Pulse) %>%
  summarise(
    mean = mean(relative_slowing, na.rm = TRUE),
    se = sqrt(sum((relative_slowing - mean(relative_slowing, na.rm = TRUE))^2, na.rm = TRUE) /
                (sum(!is.na(relative_slowing)) - 1)) /
      sqrt(sum(!is.na(relative_slowing)))
  ) %>%
  arrange(Pulse) %>%
  ggplot(aes(x = Pulse, y = mean, color = Celltype)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, linewidth = 1) +
  scale_color_manual(values = c("CMi, n=22" = microneurography_colors$nocis, "CM, n=69" = microneurography_colors$c_fibres)) +
  scale_x_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 15), labels = seq(0, 75, by = 15)) +
  labs(x = "Pulse", y = "Relative slowing (% baseline)") +
  theme(legend.title = element_blank(), legend.position = c(0.03, 0.97)) +
  ggtitle('In vivo microneurography (human)')

# Outlier check for human microneurography
human_microneurography %>%
  tidyr::pivot_longer(!Pulse, names_to = 'fibre', values_to = 'relative_slowing') %>%
  arrange(fibre) %>%
  mutate(Celltype = case_when(
    str_starts(string = fibre, pattern = 'CMi') ~ 'Silent, n=22',
    TRUE ~ 'Polymodal, n=69'
  )) %>%
  filter(Pulse == 2) %>%
  arrange(Celltype, relative_slowing)

# Compound plot
compound_RT <- human_microneurography_adsRT + pig_microneurography_adsRT + adsRT

# =============================================================================
# Save ADS Plots
# =============================================================================

ggsave(filename = here(figs_path, "3D-ads_rt_c-osmr_patchseq.png"), plot = adsRT, height = 6, width = 8)
ggsave(filename = here(figs_path, "ads_RMP_c-osmr_patchseq.png"), plot = adsrmp, height = 6, width = 8)
ggsave(filename = here(figs_path, "ads_slope_c-osmr_patchseq.png"), plot = adsslope, height = 6, width = 8)
ggsave(filename = here(figs_path, "ads_peak_c-osmr_patchseq.png"), plot = adspeak, height = 6, width = 8)

# Supplementary figure with 3 ADS plots
suppA <- adsrmp + adspeak + adsslope + plot_annotation(tag_levels = 'A')
ggsave(filename = here(figs_path, "supp-ADS_3plots.png"), plot = suppA, height = 6, width = 22)

# Microneurography plots
ggsave(filename = here(figs_path, "3B-human_microneurography_ADS_RT.png"), plot = human_microneurography_adsRT, height = 6, width = 8)
ggsave(filename = here(figs_path, "3E-pig_microneurography_adsRT.png"), plot = pig_microneurography_adsRT, height = 6, width = 8)
ggsave(filename = here(figs_path, "3BCD-compound_ADS_RT.png"), plot = compound_RT, height = 7, width = 18)

# =============================================================================
# Cell Heterogeneity Analysis Functions
# =============================================================================

#' Generate ADS RT heterogeneity plot for a specific cell type
#' @param celltype Character string specifying the cell type
#' @return ggplot object or NULL if no data available
ads_RT_cell_heterogeneity <- function(celltype = "C-OSMR-SST") {
  # Filter and prepare the data
  data_filtered <- RT_meta_table %>%
    filter(!is.na(measure), labels.p == celltype) %>%
    mutate(f_idx = as.numeric(str_remove(efeature, "ADS_RT_ABS"))) %>%
    arrange(rowname, f_idx) %>%
    group_by(rowname) %>%
    mutate(relative_slowing = (measure / first(measure)) * 100) %>%
    ungroup()
  
  # Check if the data is empty
  if (nrow(data_filtered) == 0) {
    warning(paste("No data available for cell type:", celltype))
    return(NULL)  # Or return an alternative output
  }
  
  # Proceed with plotting
  data_filtered <- data_filtered %>%
    mutate(efeature = factor(efeature, levels = paste0("ADS_RT_ABS", 1:75)))
  
  fig <- ggplot(data_filtered, aes(x = efeature, y = relative_slowing, color = labels.p, group = rowname)) +
    geom_line(alpha = .5) +
    facet_wrap(~rowname, scales = "free") +
    scale_x_discrete(
      breaks = paste0("ADS_RT_ABS", seq(0, 70, by = 10)),
      labels = seq(0, 70, by = 10)
    )
  
  return(fig)
}

#' Generate ADS RMP heterogeneity plot for a specific cell type
#' @param celltype Character string specifying the cell type
#' @return ggplot object or NULL if no data available
ads_RMP_cell_heterogeneity <- function(celltype = "C-OSMR-SST") {
  data_filtered <- rmp_meta_table %>%
    filter(!is.na(measure), labels.p == celltype)
  
  if (nrow(data_filtered) == 0) {
    warning(paste("No data available for cell type:", celltype))
    return(NULL)
  }
  
  fig <- data_filtered %>%
    mutate(f_idx = as.numeric(str_remove(efeature, "ADS_RMP_Pulse")),
           efeature = factor(efeature, levels = paste0("ADS_RMP_Pulse", 1:75))) %>%
    ggplot(aes(x = efeature, y = measure, color = labels.p, group = rowname)) +
    geom_line(alpha = .5) +
    facet_wrap(~rowname, scales = "free") +
    scale_x_discrete(breaks = paste0("ADS_RMP_Pulse", seq(0, 70, by = 10)),
                     labels = seq(0, 70, by = 10)) + 
    labs(y = "ADS RMP")
  
  return(fig)
}

#' Generate ADS Peak heterogeneity plot for a specific cell type
#' @param celltype Character string specifying the cell type
#' @return ggplot object or NULL if no data available
ads_peak_cell_heterogeneity <- function(celltype = "C-OSMR-SST") {
  data_filtered <- peak_meta_table %>%
    filter(!is.na(measure), labels.p == celltype)
  
  if (nrow(data_filtered) == 0) {
    warning(paste("No data available for cell type:", celltype))
    return(NULL)
  }
  
  fig <- data_filtered %>%
    mutate(f_idx = as.numeric(str_remove(efeature, "ADS_Peak_Pulse")),
           efeature = factor(efeature, levels = paste0("ADS_Peak_Pulse", 1:75))) %>%
    ggplot(aes(x = efeature, y = measure, color = labels.p, group = rowname)) +
    geom_line(alpha = .5) +
    facet_wrap(~rowname, scales = "free") +
    scale_x_discrete(breaks = paste0("ADS_Peak_Pulse", seq(0, 70, by = 10)),
                     labels = seq(0, 70, by = 10))
  
  return(fig)
}

#' Generate ADS Slope heterogeneity plot for a specific cell type
#' @param celltype Character string specifying the cell type
#' @return ggplot object or NULL if no data available
ads_slope_cell_heterogeneity <- function(celltype = "C-OSMR-SST") {
  data_filtered <- slope_meta_table %>%
    filter(!is.na(measure), labels.p == celltype)
  
  if (nrow(data_filtered) == 0) {
    warning(paste("No data available for cell type:", celltype))
    return(NULL)
  }
  
  fig <- data_filtered %>%
    mutate(f_idx = as.numeric(str_remove(efeature, "ADS_Slope_Pulse")),
           efeature = factor(efeature, levels = paste0("ADS_Slope_Pulse", 1:75))) %>%
    ggplot(aes(x = efeature, y = measure, color = labels.p, group = rowname)) +
    geom_line(alpha = .5) +
    facet_wrap(~rowname, scales = "free") +
    scale_x_discrete(breaks = paste0("ADS_Slope_Pulse", seq(0, 70, by = 10)),
                     labels = seq(0, 70, by = 10))
  
  return(fig)
}

# =============================================================================
# Generate Heterogeneity Plots for All Cell Types
# =============================================================================

# Create output directories
ads_RT_path <- here(figs_path, "ads_plots", "ads_RT")
ads_RMP_path <- here(figs_path, "ads_plots", "ads_RMP")
ads_peak_path <- here(figs_path, "ads_plots", "ads_peak")
ads_slope_path <- here(figs_path, "ads_plots", "ads_slope")

dir.create(ads_RT_path, recursive = TRUE)
dir.create(ads_RMP_path, recursive = TRUE)
dir.create(ads_peak_path, recursive = TRUE)
dir.create(ads_slope_path, recursive = TRUE)

# Generate heterogeneity plots for each cell type
for (celltype in unique(RT_meta_table$labels.p)) {
  ads_RT_plot <- ads_RT_cell_heterogeneity(celltype = celltype)
  ggsave(
    filename = here(ads_RT_path, str_glue("ads_RT_{celltype}_heterogeneity.png")),
    plot = ads_RT_plot, height = 12, width = 16
  )
  
  ads_RMP_plot <- ads_RMP_cell_heterogeneity(celltype = celltype)
  ggsave(
    filename = here(ads_RMP_path, str_glue("ads_RMP_{celltype}_heterogeneity.png")),
    plot = ads_RMP_plot, height = 12, width = 16
  )
  
  ads_peak_plot <- ads_peak_cell_heterogeneity(celltype = celltype)
  ggsave(
    filename = here(ads_peak_path, str_glue("ads_peak_{celltype}_heterogeneity.png")),
    plot = ads_peak_plot, height = 12, width = 16
  )
  
  ads_slope_plot <- ads_slope_cell_heterogeneity(celltype = celltype)
  ggsave(
    filename = here(ads_slope_path, str_glue("ads_slope_{celltype}_heterogeneity.png")),
    plot = ads_slope_plot, height = 12, width = 16
  )
}

# =============================================================================
# ADS RT by Cell Type Summary Plots
# =============================================================================

ads_RT_by_celltype <- RT_meta_table %>%
  filter(measure < 250) %>%
  filter(!is.na(measure)) %>%
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_RT_ABS"))) %>%
  arrange(rowname, f_idx) %>%
  group_by(rowname) %>%
  mutate(relative_slowing = (measure / first(measure)) * 100) %>%
  ungroup() %>%
  mutate(efeature = factor(efeature, levels = paste0("ADS_RT_ABS", 1:75))) %>%
  group_by(labels.p, efeature) %>%
  summarise(mean = mean(relative_slowing), se = sqrt(sum((relative_slowing - mean(relative_slowing))^2 / (length(relative_slowing) - 1))) / sqrt(length(relative_slowing))) %>%
  ggplot(aes(x = efeature, y = mean, color = labels.p, group = labels.p)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
  facet_wrap(~labels.p) +
  scale_x_discrete(
    breaks = paste0("ADS_RT_ABS", seq(0, 70, by = 10)),
    labels = seq(0, 70, by = 10)
  ) +
  scale_color_manual(values = reference_celltype_colors)

ggsave(filename = here(figs_path, "ads_RT_by_celltype.png"), plot = ads_RT_by_celltype, height = 12, width = 18)

# Version 2: Free y-axis scales
ads_RT_by_celltypev2 <- RT_meta_table %>%
  filter(measure < 250) %>%
  filter(!is.na(measure)) %>%
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_RT_ABS"))) %>%
  arrange(rowname, f_idx) %>%
  group_by(rowname) %>%
  mutate(relative_slowing = (measure / first(measure)) * 100) %>%
  ungroup() %>%
  mutate(efeature = factor(efeature, levels = paste0("ADS_RT_ABS", 1:75))) %>%
  group_by(labels.p, efeature) %>%
  summarise(mean = mean(relative_slowing), se = sqrt(sum((relative_slowing - mean(relative_slowing))^2 / (length(relative_slowing) - 1))) / sqrt(length(relative_slowing)), .groups = "drop") %>%
  ggplot(aes(x = efeature, y = mean, color = labels.p, group = labels.p)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  facet_wrap(~labels.p, scales = "free_y") + # This line has been updated
  scale_x_discrete(
    breaks = paste0("ADS_RT_ABS", seq(0, 70, by = 10)),
    labels = seq(0, 70, by = 10)
  ) +
  scale_color_manual(values = reference_celltype_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = here(figs_path, "ads_RT_by_celltypev2.png"), plot = ads_RT_by_celltypev2, height = 12, width = 18)

# Version 3: Fixed y-axis range
ads_RT_by_celltypev3 <- RT_meta_table %>%
  filter(measure < 250) %>%
  filter(!is.na(measure)) %>%
  mutate(f_idx = as.numeric(str_remove(efeature, "ADS_RT_ABS"))) %>%
  arrange(rowname, f_idx) %>%
  group_by(rowname) %>%
  mutate(relative_slowing = (measure / first(measure)) * 100) %>%
  ungroup() %>%
  mutate(efeature = factor(efeature, levels = paste0("ADS_RT_ABS", 1:75))) %>%
  group_by(labels.p, efeature) %>%
  summarise(mean = mean(relative_slowing), se = sqrt(sum((relative_slowing - mean(relative_slowing))^2 / (length(relative_slowing) - 1))) / sqrt(length(relative_slowing)), .groups = "drop") %>%
  ggplot(aes(x = efeature, y = mean, color = labels.p, group = labels.p)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  facet_wrap(~labels.p) + # No need for scales="free_y"
  scale_x_discrete(
    breaks = paste0("ADS_RT_ABS", seq(0, 70, by = 10)),
    labels = seq(0, 70, by = 10)
  ) +
  scale_color_manual(values = reference_celltype_colors) +
  scale_y_continuous(limits = c(80, 140)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = here(figs_path, "ads_RT_by_celltypev3.png"), plot = ads_RT_by_celltypev3, height = 12, width = 18)
