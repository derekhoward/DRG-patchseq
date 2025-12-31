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

results_path <- here("results/human_spatial_integration")
dir.create(results_path)

# full integration of pig snRNAseq with human visium  data
pig_ref <- readRDS(here("data/pig_snRNAseq.RDS"))

human_spatial_drg <- readRDS(here('./data/price_spatial_transcriptomics/human.drg.combined.rds'))
human_spatial_drg[["percent.mt"]] <- PercentageFeatureSet(human_spatial_drg, pattern = "^MT-")
human_spatial_QC <- VlnPlot(human_spatial_drg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
human_spatial_dimplot <- DimPlot(human_spatial_drg, label = T)
human_spatial_drg$sp.labels <- Idents(human_spatial_drg)

genes_intersection <- intersect(rownames(pig_ref), rownames(human_spatial_drg))
n_distinct(genes_intersection)

## object generation
# subset_genes_pig_ref <- pig_ref@assays$RNA[genes_intersection, ]
subset_genes_pig_ref <- GetAssayData(object = pig_ref, assay = "RNA", slot = "counts")[genes_intersection, ]

# subset_genes_patchseq <- Patchseqobject@assays$RNA[genes_intersection, ]

# subset_genes_spatial <- human_spatial_drg@assays$RNA[genes_intersection, ]
subset_genes_spatial <- GetAssayData(object = human_spatial_drg, assay = "RNA", slot = "counts")[genes_intersection, ]

pig_atlas <- CreateSeuratObject(subset_genes_pig_ref, meta.data = pig_ref@meta.data)
# patchseq <- CreateSeuratObject(subset_genes_patchseq, meta.data = Patchseqobject@meta.data)
# patchseq$labels <- "patchseq"
spatial <- CreateSeuratObject(subset_genes_spatial, meta.data = human_spatial_drg@meta.data)
glimpse(human_spatial_drg@meta.data)

pig_atlas$batch <- pig_atlas$orig.ident
pig_atlas$dataset <- "snRNAseq"
# patchseq$batch <- "p"
# patchseq$dataset <- "patch"
spatial$batch <- spatial$id
spatial$dataset <- "H_Visium"
# merged_object <- merge(pig_atlas, y = c(patchseq, spatial))
merged_object <- merge(pig_atlas, y = spatial)

# remove old objects so integration can happen smoothly with less RAM
rm(pig_ref, pig_atlas, spatial_drg, spatial)#, patchseq, Patchseqobject)

DefaultAssay(merged_object)
data.list <- SplitObject(merged_object, split.by = "batch")

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
humanspatial_preds <- read.csv(here('./results/human_spatial_integration/iterative_humanspatial_preds.csv'))
humanspatial_preds %>% head
rownames(humanspatial_preds) <- humanspatial_preds$cell_id
humanspatial_preds <- humanspatial_preds %>% dplyr::select(labels.p)
drg.integrated <- AddMetaData(object = drg.integrated, metadata = humanspatial_preds)

drg.integrated$final_labels <- ifelse(drg.integrated$dataset == "snRNAseq",
                                      as.character(drg.integrated$labels),
                                      as.character(drg.integrated$labels.p))
table(drg.integrated$final_labels, drg.integrated$dataset)
drg.integrated <- RunPCA(object = drg.integrated, verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)

saveRDS(drg.integrated, file = here(results_path, "drg_snRNAseq_humanvisium_reintegrated.RDS"))
