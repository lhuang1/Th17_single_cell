#########################################################################
######## Differential Expression: spleen EAE cluster 0 vs. 1 ############
#########################################################################

## Identify differentially expressed genes comparing cluster 0 and 1 in SPL (EAE, GFP all)

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
# setwd("~/Google Drive/research/Th17_SC")

source("1_code/utils.R")
source("1_code/differential_expr/utils.R")

library(Seurat)
library(dplyr)
library(tibble)
library(edgeR)
set.seed(1)

# configure output directories
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2020-07-22"
  prep_date <- "2020-03-24"
  clustering_date <- "2020-03-25"
} else {
  today <- cargs[1]
}

dir_out <- paste0("2_pipeline/differential_expression/SPL_EAE_0_1/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
so_all <- readRDS(paste0("2_pipeline/preprocessing/so_processed_dominant_TCR_", prep_date, ".rds"))
meta <- read.table(paste0("2_pipeline/clustering/EAE_GFPall_SPL/", clustering_date, "/FILES/meta_data.txt"), sep = " ")
so <- so_all[,rownames(meta)[meta$seurat_clusters %in% c(0, 1)]]
so@meta.data <- droplevels(meta[colnames(so),])
rm(so_all)
DefaultAssay(object = so) <- "RNA"

########Differential Expression#########
# prefilter genes
# gene_prefilter <- 0.1 # a gene should be detected by this proportion in at least one of the clusters
# prop_mat <- sapply(unique(so$seurat_clusters), function(x){
#   Matrix::rowMeans(so@assays$RNA@counts[,so$seurat_clusters == x] > 0) 
# }) 
# keep <- (prop_mat > gene_prefilter) %>% Matrix::rowSums(.) > 0
# so <- subset(so, features = which(keep))

# run edgeR
results <- run_edgeR_spl_0_1(so)
saveRDS(results, file = paste0(dir_out, "results.rds"))

# ####### process results into signature ######
results <- readRDS(file = paste0(dir_out, "results.rds"))
cxcr6 <- results$table %>% rownames_to_column(var = "gene") %>% filter(logFC < -log2(1.5) & FDR < 0.05) %>% arrange(FDR)
slamf6 <- results$table %>% rownames_to_column(var = "gene") %>% filter(logFC > log2(1.5) & FDR < 0.05) %>% arrange(FDR)

sig <- list(
  'Cxcr6_sc-plus' = cxcr6$gene,
  'Slamf6_sc-plus' = slamf6$gene
)

saveRDS(sig, file = paste0(dir_out, "sc_signature.rds"))
