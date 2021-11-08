################################################################
######## Differential Expression: cluster markers   ############
################################################################

## Identify differentially expressed genes comparing each cluster to all other clusters in UT, GFP all, intra tissue clustering

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")

source("1_code/utils.R")
source("1_code/differential_expr/utils.R")

library(Seurat)
library(dplyr)
library(edgeR)
set.seed(1)

# configure output directories
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2021-05-18"
  clustering_date <- "2020-03-25"
} else {
  today <- cargs[1]
  clustering_date <- cargs[2]
}

dir_out <- paste0("2_pipeline/differential_expression/Cluster_UT_GFPall_intra/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}


so_all <- readRDS("for_paper/results/preprocessing/so_processed_dominant_TCR_2020-03-24.rds")
DefaultAssay(so_all) <- "RNA"
so_all[["integrated"]] <- NULL
so_all[["SCT"]] <- NULL

tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL")
for (tissue in tissue_vec) {
  cat(tissue, "\n")

  ### load data
  if (tissue == "SPL") {
    fname <- paste0("2_pipeline/clustering/UT_GFPall_intra/", clustering_date, "/FILES/SPL/meta_data_switched.txt")
    meta <- read.table(fname, sep = " ", row.names = 1, header = T)
  } else {
    fname <- paste0("2_pipeline/clustering/UT_GFPall_intra/", clustering_date, "/FILES/", tissue, "/meta_data.csv")
    meta <- read.table(fname, sep = ",", row.names = 1, header = T)
  }
  
  so <- so_all[,rownames(meta)]
  so@meta.data <- meta

  ########Differential Expression#########
  # prefilter genes
  gene_prefilter <- 0.1 # a gene should be detected by this proportion in at least one of the clusters
  prop_mat <- sapply(unique(so$seurat_clusters), function(x){
    Matrix::rowMeans(so@assays$RNA@counts[,so$seurat_clusters == x] > 0) 
  }) 
  keep <- Matrix::rowSums(prop_mat > gene_prefilter) > 0
  so <- subset(so, features = which(keep))
  
  # run edgeR
  so$batch <- factor(so$batch)
  so$seurat_clusters <- factor(so$seurat_clusters)
  # Create dummy variables
  contrasts(so$batch) = contr.sum(nlevels(so$batch))
  # Create a DGEList data object.
  dgeFull <- DGEList(so[['RNA']]@counts)
  # Estimate the normalization factors
  dgeFull <- calcNormFactors(dgeFull, method="TMM")
  # compute detection rate [here I use scaled number of detected genes]
  so$cdr <- scale(so$nFeature_RNA)
  
  
  results <- list()
  for (cluster_to_test in sort(unique(so$seurat_clusters))) {
    so$in_cluster <- factor(so$seurat_clusters == cluster_to_test, levels = c(FALSE, TRUE))
    # create design matrix
    design <- model.matrix(~ in_cluster + batch + cdr, data = so@meta.data)
    # Estimate dispersion
    dgeFull <- estimateDisp(dgeFull, design = design)
    # Perform QLF tests
    fit <- glmQLFit(dgeFull, design = design)
    qlf <- glmQLFTest(fit, coef = 2)
    results[[paste0(tissue, "_", cluster_to_test)]] <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)$table
  }
  saveRDS(results, file = paste0(dir_out, tissue, "_results.rds"))
}
