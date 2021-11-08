################################################################
######## Differential Expression: tissue vs. spleen ############
################################################################

## Identify differentially expressed genes comparing each tissue to spleen (UT, GFPall)

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
# setwd("~/Google Drive/research/Th17_SC")

source("1_code/utils.R")
source("1_code/differential_expr/utils.R")

library(Seurat)
library(dplyr)
library(edgeR)
set.seed(1)

# configure output directories
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2020-05-18_downsample"
  trt_gfp <- "EAE_UT_GFPall"
} else {
  today <- cargs[1]
  trt_gfp <- cargs[2]
}

dir_out <- paste0("2_pipeline/differential_expression/", trt_gfp, "_vs_SPL/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
so_all <- readRDS("2_pipeline/preprocessing/so_processed_dominant_TCR_2020-03-24.rds")
### fix batch coding #####
so_all$batch[so_all$batch == "b8" & so_all$treatment == "UT"] <- "b9"
so_all$batch[so_all$batch == "b7" & so_all$treatment == "UT"] <- "b8"

if (trt_gfp == "UT_GFPpos") {
  keep_cells_sub <- colnames(so_all)[so_all$treatment == "UT" & so_all$batch %in% paste0("b", 1:5) & so_all$GFP_positive]
} else if (trt_gfp == "UT_GFPall") {
  keep_cells_sub <- colnames(so_all)[so_all$treatment == "UT" & so_all$batch %in% paste0("b", 1:5)]
} else if (trt_gfp == "EAE_UT_GFPall") {
  keep_cells_sub <- colnames(so_all)
  so_all$tissue_treatment <- paste(so_all$tissue, so_all$treatment, sep = "_")
}
so <- so_all[,keep_cells_sub]
so@meta.data <- droplevels(so@meta.data)
rm(so_all)
DefaultAssay(object = so) <- "RNA"

########Differential Expression#########
## prefilter genes
gene_prefilter <- 0.1 # a gene should be detected by this proportion in at least one of the tissues
if (trt_gfp == "UT_GFPpos") {
  prop_mat <- sapply(unique(so$tissue), function(x){
    Matrix::rowMeans(so@assays$RNA@counts[,so$tissue == x] > 0) 
  }) 
  keep <- (prop_mat > gene_prefilter) %>% Matrix::rowSums(.) > 0
  so <- subset(so, features = which(keep))
} else if (trt_gfp == "UT_GFPall") {
  prop_mat <- sapply(unique(so$tissue), function(x){
    Matrix::rowMeans(so@assays$RNA@counts[,so$tissue == x] > 0) 
  }) 
  keep <- (prop_mat > gene_prefilter) %>% Matrix::rowSums(.) > 0
  ## make sure genes used in GFP pos are also included so can compare results
  prop_mat_gfppos <- sapply(unique(so$tissue), function(x){
    Matrix::rowMeans(so@assays$RNA@counts[,so$tissue == x & so$GFP_positive] > 0) 
  }) 
  keep_gfppos <- (prop_mat_gfppos > gene_prefilter) %>% Matrix::rowSums(.) > 0
  so <- subset(so, features = which(keep | keep_gfppos))
} else if (trt_gfp == "EAE_UT_GFPall") {
  prop_mat <- sapply(unique(so$tissue_treatment), function(x){
    Matrix::rowMeans(so@assays$RNA@counts[,so$tissue_treatment == x] > 0) 
  }) 
  keep <- (prop_mat > gene_prefilter) %>% Matrix::rowSums(.) > 0
  so <- subset(so, features = which(keep))
}


## downsample??
downsample <- TRUE
if (downsample) {
  n_downsample_cell <- ncol(so)/10
  if (trt_gfp == "EAE_UT_GFPall" & !is.na(n_downsample_cell)) {
    so <- so[,sample(1:ncol(so), round(n_downsample_cell, 0), replace = FALSE)]
    so@meta.data <- droplevels(so@meta.data)
  }
}



# run edgeR
if (trt_gfp == "UT_GFPpos" | trt_gfp == "UT_GFPall") {
  results <- run_edgeR_vs_spleen(so)
} else if (trt_gfp == "EAE_UT_GFPall") {
  results <- run_edgeR_vs_spleen_tissue_trt(so)
}
saveRDS(results, file = paste0(dir_out, "results.rds"))

####### process DE results ###########
# (1) obtain and plot (stain, violin) genes for each tissue; 
# (2) compute sample fold changes; 
# results <- readRDS(file = paste0(dir_out, "results.rds"))
# process_de_vs_spleen_results(so, results, prop_mat, dir_out, thresh_fdr = 0.05, thresh_fc = 1.5)
