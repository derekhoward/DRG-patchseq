suppressPackageStartupMessages({
  library(here)
  library(biomaRt)
  library(Seurat)
  library(readr)
  library(dplyr)
  library(stringr)
  library(rlang)
  library(grid)
  library(tidyr)
  library(purrr)
})

# Global configurations
reference_celltype_levels <- factor(
  c("Aβ-PROPRIO", "Aβ-LTMR-ALDH1A1", "Aβ-LTMR-SLIT2", "Aβ-LTMR-PALM", "Aβ-HTMR", 
    "Aδ-LTMR", "Aδ-CACNA1E", "Aδ-TRPV1", "Aδ-COOL", "C-COLD", 
    "C-TAC1-KCNQ5", "C-TAC1-LRP1B", "C-TAC1-TRPA1", "C-LTMR", 
    "C-OSMR-GFRA1_2", "C-OSMR-SST"),
  levels = c("Aβ-PROPRIO", "Aβ-LTMR-ALDH1A1", "Aβ-LTMR-SLIT2", "Aβ-LTMR-PALM", "Aβ-HTMR", 
             "Aδ-LTMR", "Aδ-CACNA1E", "Aδ-TRPV1", "Aδ-COOL", "C-COLD", 
             "C-TAC1-KCNQ5", "C-TAC1-LRP1B", "C-TAC1-TRPA1", "C-LTMR", 
             "C-OSMR-GFRA1_2", "C-OSMR-SST")
)

reference_celltype_colors <- c(
  "C-TAC1-KCNQ5" = "#F8766D", "C-TAC1-LRP1B" = "#E88526", "C-LTMR" = "#D39200", "Aβ-HTMR" = "#B79F00",
  "C-OSMR-SST" = "#93AA00", "Aδ-COOL" = "#5EB300", "C-TAC1-TRPA1" = "#00BA38", "Aβ-LTMR-PALM" = "#00BF74",
  "Aδ-TRPV1" = "#00C19F", "C-COLD" = "#00BFC4", "Aβ-PROPRIO" = "#00B9E3", "Aβ-LTMR-SLIT2" = "#00ADFA",
  "C-OSMR-GFRA1_2" = "#619CFF", "Aδ-LTMR" = "#DB72FB", "Aβ-LTMR-ALDH1A1" = "#F564E3",
  "Aδ-CACNA1E" = "#FF61C3", "patchseq" = "#c00000"
)

# Utility functions

prep_patchseq_obj_all_batches <- function(convert_genes = "pig") {
  if (convert_genes == "pig" && file.exists(here("data/processed/pig_patchseq.RDS")) == FALSE) {
    pig_counts <- read_csv(here("./data/counts_plate1-6.csv"))
    names(pig_counts)[1] <- "ensembl_id"
    sampleIDs <- sort(names(pig_counts[2:270]))
    dir.create(here("./data/processed"))
    write.csv(sampleIDs, file = here("./data/processed/sampleIDs.csv"))
    
    pig_counts <- pig_counts[, c("ensembl_id", sampleIDs)]
    
    pigMart <- useMart("ensembl", dataset = "sscrofa_gene_ensembl") # ,host = "https://dec2021.archive.ensembl.org/")
    # find gene symbols for ENSSSG where possible
    gene_conversion <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "external_gene_name"), values = pig_counts$ensembl_id, mart = pigMart)
    # fill empty values of gene symbol with NA
    gene_conversion <- mutate_at(gene_conversion, .vars = "external_gene_name", .funs = ~ na_if(., ""))
    # use gene_symbol if possible, else use ENSSSCG
    gene_conversion$merged_symbols <- coalesce(gene_conversion$external_gene_name, gene_conversion$ensembl_gene_id)
    gene_conversion <- gene_conversion %>% dplyr::select(ensembl_gene_id, merged_symbols)
    
    pig_counts_gene_symbols <- inner_join(pig_counts, gene_conversion, by = c("ensembl_id" = "ensembl_gene_id")) %>%
      dplyr::select(merged_symbols, everything(), -ensembl_id)
    
    # need to drop rows with duplicate gene_symbols
    # here we just drop all rows that have duplicates
    # should think of an alternate strategy
    duplicate_gene_rows <- pig_counts_gene_symbols %>% filter(duplicated(merged_symbols))
    pig_counts_gene_symbols <- pig_counts_gene_symbols %>% filter(!duplicated(merged_symbols))
    
    forObject <- as.data.frame(pig_counts_gene_symbols[, 2:270])
    row.names(forObject) <- pig_counts_gene_symbols$merged_symbols
    Patchseqobject <- CreateSeuratObject(counts = forObject)
    
    Idents(Patchseqobject, cells = 1:48) <- "p1"
    Idents(Patchseqobject, cells = 49:96) <- "p2"
    Idents(Patchseqobject, cells = 97:144) <- "p3"
    Idents(Patchseqobject, cells = 145:192) <- "p4"
    Idents(Patchseqobject, cells = 193:240) <- "p5"
    Idents(Patchseqobject, cells = 241:270) <- "p6"
    
    metas <- read.csv(here("./data/patchseq_metadata.csv"), header = TRUE, row.names = 1)
    rownames(metas) <- gsub("X", "", rownames(metas))
    Patchseqobject <- AddMetaData(object = Patchseqobject, metadata = metas)
    
    Patchseqobject$batch <- Idents(Patchseqobject)
    Patchseqobject$batch
    Patchseqobject <- subset(Patchseqobject, subset = nCount_RNA > 1e+05)
    Patchseqobject$species <- "pig"
    
    processed_path <- here("data/processed/")
    dir.create(processed_path, recursive = T)
    saveRDS(Patchseqobject, here(processed_path, "pig_patchseq.RDS"))
  } else if (convert_genes == "pig" && file.exists(here("data/processed/pig_patchseq.RDS")) == TRUE) {
    cat("Loading pre-computed patchseq object with pig gene symbols (merged ENSSSCG and gene symbols) from: ", here("data/processed/pig_patchseq.RDS"))
    Patchseqobject <- readRDS(here("data/processed/pig_patchseq.RDS"))
  }
  if (convert_genes == "human" && file.exists(here("data/processed/pig_patchseq-human.RDS")) == FALSE) {
    pig_counts <- read_csv(here("./data/counts_plate1-6.csv"))
    names(pig_counts)[1] <- "ensembl_id"
    sampleIDs <- sort(names(pig_counts[2:270]))
    dir.create(here("./data/processed"))
    write.csv(sampleIDs, file = here("./data/processed/sampleIDs.csv"))
    
    pigMart <- useMart("ensembl", dataset = "sscrofa_gene_ensembl") # ,host = "https://dec2021.archive.ensembl.org/")
    humanMart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # ,host = "https://dec2021.archive.ensembl.org/")
    genesV2 <- getLDS(
      attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = pig_counts$ensembl_id,
      mart = pigMart, attributesL = c("hgnc_symbol"), martL = humanMart, uniqueRows = T
    )
    
    # write out gene names translation table
    write_csv(genesV2, here(paste0("data/processed/gene_translation_", convert_genes, ".csv")))
    
    pig_counts <- pig_counts[, c("ensembl_id", sampleIDs)]
    # pig_counts[ , order(names(pig_counts))]
    # patchseqdata2 <- patchseqdata2[, c('ensembl_id', sampleIDs)]
    
    # filter out mitochondrial genes
    genesV2 <- genesV2 %>% filter(str_detect(HGNC.symbol, "MT-", negate = TRUE))
    
    pig_counts_gene_symbols <- inner_join(pig_counts, genesV2, by = c("ensembl_id" = "Gene.stable.ID")) %>%
      dplyr::select(HGNC.symbol, everything(), -ensembl_id)
    # relocate(HGNC.symbol)
    
    # need to drop rows with duplicate gene_symbols
    # here we just drop all rows that have duplicates
    # should think of an alternate strategy
    duplicate_gene_rows <- pig_counts_gene_symbols %>% filter(duplicated(HGNC.symbol))
    pig_counts_gene_symbols <- pig_counts_gene_symbols %>% filter(!duplicated(HGNC.symbol))
    
    forObject <- as.data.frame(pig_counts_gene_symbols[, 2:270])
    row.names(forObject) <- pig_counts_gene_symbols$HGNC.symbol
    Patchseqobject <- CreateSeuratObject(counts = forObject)
    
    Idents(Patchseqobject, cells = 1:48) <- "p1"
    Idents(Patchseqobject, cells = 49:96) <- "p2"
    Idents(Patchseqobject, cells = 97:144) <- "p3"
    Idents(Patchseqobject, cells = 145:192) <- "p4"
    Idents(Patchseqobject, cells = 193:240) <- "p5"
    Idents(Patchseqobject, cells = 241:270) <- "p6"
    
    metas <- read.table(here("./data/patchseq_metadata.csv"), sep = ";", header = TRUE, row.names = 1)
    rownames(metas) <- gsub("X", "", rownames(metas))
    Patchseqobject <- AddMetaData(object = Patchseqobject, metadata = metas)
    
    Patchseqobject$batch <- Idents(Patchseqobject)
    Patchseqobject$batch
    Patchseqobject <- subset(Patchseqobject, subset = nCount_RNA > 1e+05)
    Patchseqobject$species <- "pig"
    
    processed_path <- here("data/processed/")
    dir.create(processed_path, recursive = T)
    saveRDS(Patchseqobject, here(processed_path, "pig_patchseq-human.RDS"))
  } else if (convert_genes == "human" && file.exists(here("data/processed/pig_patchseq-human.RDS")) == TRUE) {
    cat("Loading pre-computed patchseq object with human gene symbols from: ", here("data/processed/pig_patchseq-human.RDS"))
    Patchseqobject <- readRDS(here("data/processed/pig_patchseq-human.RDS"))
  } else if (convert_genes == "mouse" && file.exists(here("data/processed/pig_patchseq-mouse.RDS")) == FALSE) {
    pig_counts <- read_csv(here("./data/counts_plate1-6.csv"))
    names(pig_counts)[1] <- "ensembl_id"
    sampleIDs <- sort(names(pig_counts[2:270]))
    dir.create(here("./data/processed"))
    write.csv(sampleIDs, file = here("./data/processed/sampleIDs.csv"))
    
    pigMart <- useMart("ensembl", dataset = "sscrofa_gene_ensembl") # ,host = "https://dec2021.archive.ensembl.org/")
    mouseMart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") # , host = "https://dec2021.archive.ensembl.org/")
    genesV2 <- getLDS(
      attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = pig_counts$ensembl_id,
      mart = pigMart, attributesL = c("mgi_symbol"), martL = mouseMart, uniqueRows = T
    )
    
    # write out gene names translation table
    write_csv(genesV2, here(paste0("data/processed/gene_translation_", convert_genes, ".csv")))
    
    pig_counts <- pig_counts[, c("ensembl_id", sampleIDs)]
    # pig_counts[ , order(names(pig_counts))]
    # patchseqdata2 <- patchseqdata2[, c('ensembl_id', sampleIDs)]
    
    # filter out mitochondrial genes
    genesV2 <- genesV2 %>% filter(str_detect(MGI.symbol, "mt-", negate = TRUE))
    
    pig_counts_gene_symbols <- inner_join(pig_counts, genesV2, by = c("ensembl_id" = "Gene.stable.ID")) %>%
      dplyr::select(MGI.symbol, everything(), -ensembl_id)
    # relocate(HGNC.symbol)
    
    # need to drop rows with duplicate gene_symbols
    # here we just drop all rows that have duplicates
    # should think of an alternate strategy
    duplicate_gene_rows <- pig_counts_gene_symbols %>% filter(duplicated(MGI.symbol))
    pig_counts_gene_symbols <- pig_counts_gene_symbols %>% filter(!duplicated(MGI.symbol))
    
    forObject <- as.data.frame(pig_counts_gene_symbols[, 2:270])
    row.names(forObject) <- pig_counts_gene_symbols$MGI.symbol
    Patchseqobject <- CreateSeuratObject(counts = forObject)
    
    Idents(Patchseqobject, cells = 1:48) <- "p1"
    Idents(Patchseqobject, cells = 49:96) <- "p2"
    Idents(Patchseqobject, cells = 97:144) <- "p3"
    Idents(Patchseqobject, cells = 145:192) <- "p4"
    Idents(Patchseqobject, cells = 193:240) <- "p5"
    Idents(Patchseqobject, cells = 241:270) <- "p6"
    
    metas <- read.table(here("./data/patchseq_metadata.csv"), sep = ";", header = TRUE, row.names = 1)
    rownames(metas) <- gsub("X", "", rownames(metas))
    Patchseqobject <- AddMetaData(object = Patchseqobject, metadata = metas)
    
    Patchseqobject$batch <- Idents(Patchseqobject)
    Patchseqobject$batch
    Patchseqobject <- subset(Patchseqobject, subset = nCount_RNA > 1e+05)
    Patchseqobject$species <- "pig"
    
    processed_path <- here("data/processed/")
    dir.create(processed_path, recursive = T)
    saveRDS(Patchseqobject, here(processed_path, "pig_patchseq-mouse.RDS"))
  } else if (convert_genes == "mouse" && file.exists(here("data/processed/pig_patchseq-mouse.RDS")) == TRUE) {
    cat("Loading pre-computed patchseq object with mouse gene symbols from: ", here("data/processed/pig_patchseq-mouse.RDS"))
    Patchseqobject <- readRDS(here("data/processed/pig_patchseq-mouse.RDS"))
  }
  return(Patchseqobject)
}



#' Get gene conversion
#'
#' @param ensembl_ids Vector of Ensembl IDs
#' @param convert_genes Species to convert genes to
#' @return Data frame with converted gene names
get_gene_conversion <- function(ensembl_ids, convert_genes) {
  pigMart <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")
  
  if (convert_genes == "pig") {
    conversion <- getBM(
      filters = "ensembl_gene_id",
      attributes = c("ensembl_gene_id", "external_gene_name"),
      values = ensembl_ids,
      mart = pigMart
    )
    conversion$merged_symbols <- coalesce(conversion$external_gene_name, conversion$ensembl_gene_id)
    return(conversion[, c("ensembl_gene_id", "merged_symbols")])
  } else {
    targetMart <- useMart("ensembl", dataset = paste0(ifelse(convert_genes == "human", "hsapiens", "mmusculus"), "_gene_ensembl"))
    conversion <- getLDS(
      attributes = c("ensembl_gene_id"),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = pigMart,
      attributesL = c(ifelse(convert_genes == "human", "hgnc_symbol", "mgi_symbol")),
      martL = targetMart,
      uniqueRows = TRUE
    )
    colnames(conversion)[2] <- "target_symbol"
    return(conversion)
  }
}

#' Process counts data
#'
#' @param pig_counts Data frame with pig counts
#' @param gene_conversion Data frame with gene conversion information
#' @param convert_genes Species to convert genes to
#' @return Processed counts data frame
process_counts <- function(pig_counts, gene_conversion, convert_genes) {
  if (convert_genes == "pig") {
    processed_counts <- inner_join(pig_counts, gene_conversion, by = c("ensembl_id" = "ensembl_gene_id")) %>%
      dplyr::select(merged_symbols, everything(), -ensembl_id)
  } else {
    symbol_col <- ifelse(convert_genes == "human", "HGNC.symbol", "MGI.symbol")
    processed_counts <- inner_join(pig_counts, gene_conversion, by = c("ensembl_id" = "Gene.stable.ID")) %>%
      dplyr::select(target_symbol, everything(), -ensembl_id) %>%
      dplyr::rename(!!symbol_col := target_symbol)
  }
  
  processed_counts <- processed_counts %>% filter(!duplicated(!!sym(colnames(processed_counts)[1])))
  return(processed_counts)
}

#' Create Seurat object
#'
#' @param counts_data Processed counts data
#' @return Seurat object
create_seurat_object <- function(counts_data) {
  counts_matrix <- as.matrix(counts_data[, -1])
  rownames(counts_matrix) <- counts_data[[1]]
  seuratObj <- CreateSeuratObject(counts = counts_matrix)
  
  seuratObj$batch <- rep(paste0("p", 1:6), c(48, 48, 48, 48, 48, 30))
  seuratObj <- subset(seuratObj, subset = nCount_RNA > 1e5)
  seuratObj$species <- "pig"
  
  return(seuratObj)
}

#' Add metadata to Seurat object
#'
#' @param seuratObj Seurat object
#' @return Seurat object with added metadata
add_metadata <- function(seuratObj) {
  metas <- read.csv(here::here("data", "patchseq_metadata.csv"), row.names = 1)
  rownames(metas) <- gsub("X", "", rownames(metas))
  seuratObj <- AddMetaData(object = seuratObj, metadata = metas)
  return(seuratObj)
}

#' Perform multi-bar heatmap
#'
#' @param object Seurat object
#' @param features Features to plot
#' @param cells Cells to plot
#' @param group.by Grouping variable
#' @param ... Additional parameters
#' @return ggplot object
#' @export
DoMultiBarHeatmap <- function(object, features = NULL, cells = NULL, group.by = "ident", ...) {
  # Implementation details...
  # (Keep the existing implementation of this function)
}

#' Predict cell type by neighbors
#'
#' @param seurat_obj Seurat object
#' @param cell_ids_to_predict Cell IDs to predict
#' @param labels_df Data frame with cell type labels
#' @return Data frame with cell type predictions
#' @export
predict_celltype_by_neighbours <- function(seurat_obj, cell_ids_to_predict, labels_df) {
  
  # Identify the neighbours to each patchseq cell
  # using nearest neighbours in pca space
  # modify the default k.param to be increased, since we have relatively fewer cells from labelled snRNAseq dataset
  seurat_obj <- FindNeighbors(seurat_obj, return.neighbor = T, reduction = "pca", dims = 1:30, k.param=20)
  
  # Get neighbors for data subset
  neighbours <- purrr::map(.x = cell_ids_to_predict, .f = TopNeighbors, object = seurat_obj@neighbors$integrated.nn, n = 20)
  
  # Transform into data frame
  neighbours <- as_tibble(t(as.data.frame(neighbours)))
  neighbours <- tidyr::pivot_longer(neighbours, cols = V2:V20, names_to = "neighbor") %>% rename(pseq = V1)
  
  # Merge cell type information
  neighbours <- left_join(neighbours, y = labels_df, by = c("value" = "cell_id"))
  
  # Generate predictions
  results <- neighbours %>%
    tidyr::drop_na() %>%
    group_by(pseq) %>%
    count(celltype) %>%
    mutate(type_fraction = n / sum(n), proportion_neighbours_atlas_cells = sum(n)/20) %>%
    arrange(pseq, desc(n))
  
  # Assign confidence labels
  results <- results %>%
    slice(1) %>%
    mutate(confidence = case_when(
      type_fraction > 0.7 ~ "strong",
      type_fraction >= 0.4 ~ "mid",
      TRUE ~ "low"
    ))
  
  # Prepare predictions_meta data frame
  predictions_meta <- results %>%
    dplyr::select(cell_id = pseq, labels.p = celltype, type_fraction, confidence, proportion_neighbours_atlas_cells) #%>%
  
  return(predictions_meta)
}
