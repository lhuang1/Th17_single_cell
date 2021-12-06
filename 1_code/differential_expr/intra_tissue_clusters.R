################################################################
######## Differential Expression: cluster markers   ############
################################################################

## Identify differentially expressed genes comparing each cluster to all other clusters in UT/EAE, GFP all, intra tissue clustering
## edgeR/QLF with cellular detection rate as covariate

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")

source("1_code/utils.R")
source("1_code/differential_expr/utils.R")

library(Seurat)
library(dplyr)
library(edgeR)
set.seed(1)

# configuration
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- Sys.Date()
  trt_gfp <- "EAE_GFPall" # or "UT_GFPall"; treatment and GFP combination. 
  prep_date <- "2020-03-24"
  clustering_date <- "2020-03-25"
} else {
  today <- cargs[1]
  trt_gfp <- cargs[2]
  prep_date <- cargs[3]
  clustering_date <- cargs[4]
}
dir_out <- paste0("2_pipeline/differential_expression/Cluster_", trt_gfp, "_intra/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
so_all <- readRDS(paste0("2_pipeline/preprocessing/so_processed_dominant_TCR_", prep_date, ".rds"))

tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL")
if (trt_gfp == "EAE_GFPall") {
  tissue_vec <- c(tissue_vec, "CNS", "DLN")
}
for (tissue in tissue_vec) {
  meta <- read.table(paste0("2_pipeline/clustering/", trt_gfp, "_intra/", clustering_date, "/FILES/", tissue, "/meta_data.csv"),
                     sep = ",", row.names = 1, header = T)
  so <- so_all[,rownames(meta)]
  so@meta.data <- meta
  rm(so_all)
  DefaultAssay(object = so) <- "RNA"
  
  ########Differential Expression#########
  # prefilter genes
  gene_prefilter <- 0.1 # a gene should be detected by this proportion in at least one of the tissues
  prop_mat <- sapply(unique(so$seurat_clusters), function(x){
    Matrix::rowMeans(so@assays$RNA@counts[,so$seurat_clusters == x] > 0) 
  }) 
  keep <- (prop_mat > gene_prefilter) %>% Matrix::rowSums(.) > 0
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
    results[[paste0("C_", cluster_to_test)]] <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)$table
  }
  saveRDS(results, file = paste0(dir_out, tissue, "_results.rds"))
}
