library(Seurat)
library(here)
library(dplyr)
library(SummarizedExperiment)
library(biomaRt)

# Function to convert mouse gene symbols to human gene symbols using biomaRt
convert_mouse_to_human <- function(gene_list) {
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genes <- getLDS(attributes = c("mgi_symbol"), 
                  filters = "mgi_symbol", 
                  values = gene_list, 
                  mart = mouse, 
                  attributesL = c("hgnc_symbol"), 
                  martL = human, 
                  uniqueRows = TRUE)
  
  return(genes)
}

# Function to remove duplicate gene symbols
remove_duplicates <- function(counts, genes) {
  genes <- genes[complete.cases(genes), ]
  unique_genes <- genes[!duplicated(genes$HGNC.symbol), ]
  matched_genes <- unique_genes[unique_genes$MGI.symbol %in% rownames(counts), ]
  unique_counts <- counts[matched_genes$MGI.symbol, , drop = FALSE]
  rownames(unique_counts) <- matched_genes$HGNC.symbol
  return(unique_counts)
}

# Read in the data
drg <- readRDS(here('data/DRG_neurons_complete.Rds'))
DefaultAssay(drg) <- "RNA"

# Extract the raw counts data for all cells
raw_counts <- GetAssayData(drg, layer = "counts")

# Convert mouse gene symbols to human gene symbols
mouse_genes <- rownames(raw_counts)
converted_genes <- convert_mouse_to_human(mouse_genes)

# Filter the raw counts data to only include the converted genes and remove duplicates
converted_counts <- remove_duplicates(raw_counts, converted_genes)

# Create a new Seurat object with the converted counts
drg_converted <- CreateSeuratObject(counts = converted_counts)

# Add the original metadata back to the new Seurat object
drg_converted@meta.data <- drg@meta.data[colnames(converted_counts),]

# Function to convert Seurat object to SummarizedExperiment
convertSeuratToSummarizedExperiment <- function(seurat_object) {
  counts <- as.matrix(GetAssayData(seurat_object, slot = "counts"))
  meta_data <- seurat_object@meta.data
  colnames(counts) <- rownames(meta_data)
  row_data <- DataFrame(gene = rownames(counts))
  rownames(row_data) <- rownames(counts)
  
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = meta_data,
    rowData = row_data
  )
  
  return(se)
}

# Convert the entire Seurat object to SummarizedExperiment
drg_se <- convertSeuratToSummarizedExperiment(drg_converted)

# Save the SummarizedExperiment object
dir.create(here('./data/summarized_experiments'), recursive = TRUE, showWarnings = FALSE)
saveRDS(drg_se, file = here('./data/summarized_experiments/drg_complete_SE.RDS'))
